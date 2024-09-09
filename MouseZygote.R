###---------------------------------------------------------------------------------------------------------.
###--------------------------------------- Mouse Zygote Analysis -------------------------------------------------------
###---------------------------------------------------------------------------------------------------------.

###---------------------------------------------------------------------------------------------------------.
###--------------------------------------- parameter set-up -------------------------------------------------------
###---------------------------------------------------------------------------------------------------------.

source("libraries_and_functions.R")

# prot to gene names 
convertProts <- read.delim("uniprot_Mus_musculus_ProteinsTOGenes.tab")
convertProts$LeadingGene <- sub(" .*", "", convertProts$Gene.names)
convertProts$uppercase <- toupper(convertProts$LeadingGene)

d2 <- read.delim("ev_updated.txt")
reference_channel = 2

# Load experimental design
design <- read.csv("annotation_parth_cutZygotes.csv")
design.melt <- reshape2::melt(design, id.vars = "Set")
design.melt$set_RI <- paste(design.melt$Set, design.melt$variable)
design.melt$set_RI_celltype <- paste(design.melt$Set, design.melt$variable, design.melt$value)

# obtain only the single cell samples 
ev <- d2 #%>% dplyr::filter(grepl("006|007|008|009", Raw.file))

# load significantly differential proteins: 
sig2c <- read.csv(paste0(WindowsPath, "/Collaborators/Zernicka_Goetz/code/2020_DecPlate01/2-cell-stage-combined-sig-Prots_5percentFDR_2022_02_07.csv"))


#------------------------------------------------------------------------------------.
################## QC #################################
#------------------------------------------------------------------------------------.
parse_row<-grep("|",ev$Leading.razor.protein, fixed=T)

if(length(parse_row)>0){
  split_prot<-str_split(ev$Leading.razor.protein[parse_row], pattern = fixed("|"))
  split_prot2<-unlist(split_prot)[seq(2,3*length(split_prot),3)]
  ev$Leading.razor.protein[parse_row]<-split_prot2
  
}

# Filter out reverse hits, contaminants, and contaminated spectra...
ev<-ev[-which(ev$Reverse=="+"),]
if(length(grep("REV", ev$Leading.razor.protein))>0){ ev<-ev[-grep("REV", ev$Leading.razor.protein),] }
if(length(grep("CON", ev$Leading.razor.protein))>0){ ev<-ev[-grep("CON", ev$Leading.razor.protein),] }
if(length(which(ev$Potential.contaminant=="+"))>0){ ev<-ev[-which(ev$Potential.contaminant=="+"),] }
ev<-ev[!is.na(ev$PIF),]
# create a new column called SeqCharge 
ev$SeqCharge <- paste0(ev$Modified.sequence, ev$Charge)

# filter for PEP and PIF 
ev.1 <- ev %>% dplyr::filter(PEP <= 0.03) %>% dplyr::filter(PIF >=0.8)
ev.1 <- ev %>% dplyr::filter(dart_PEP <= 0.02) %>% dplyr::filter(PIF >=0.8)
# select the most necessary columns
ev.2 <- ev.1 %>% dplyr::select(Raw.file, Leading.razor.protein, SeqCharge, contains("intensity.corrected"))
# change O's to NA's 
ev.2[ev.2 == 0] <- NA 



#------------------------------------------------------------------------------------..
######################## Peptides -> Proteins #################################################
#------------------------------------------------------------------------------------..

ev.3 <- ev.2 %>% dplyr::filter(grepl("eAP008|eAP009", Raw.file))


# normalize to the reference 
RI.ref <- paste0("Reporter.intensity.corrected.", reference_channel)
ev.3[,4:19] <- ev.3[,4:19] / ev.3[,which(colnames(ev.3) == RI.ref)]

# melt 
ev.3.melt <- reshape2::melt(ev.3, id.vars = c("Raw.file", "Leading.razor.protein", "SeqCharge"))
ev.3.melt$RI <- paste0("RI",gsub("(?:[^.]+\\.){3}([^.]+).*", "\\1", ev.3.melt$variable))
ev.3.melt$Raw.file_RI <- paste(ev.3.melt$Raw.file, ev.3.melt$RI)

# now we can match the annotation file with this dataframe
ev.3.melt$celltype <- design.melt$value[match(ev.3.melt$Raw.file_RI, design.melt$set_RI)]
ev.3.melt$Raw.file_RI_celltype <- paste(ev.3.melt$Raw.file_RI, ev.3.melt$celltype)
# filter out the carrier, unused, and control 
if(chooseSets == "parth") { 
  ev.3.melt.quant <- ev.3.melt %>% 
    dplyr::filter(!grepl("carrier|unused|reference|control|sc_x|sc_c", Raw.file_RI_celltype))
}
if(chooseSets == "cuts") { 
  ev.3.melt.quant <- ev.3.melt %>% 
    dplyr::filter(!grepl("carrier|unused|reference|control|sc_x|sc_c", Raw.file_RI_celltype)) %>% 
    dplyr::filter(Raw.file_RI_celltype %in% cells2keep$variable)
  ev.3.melt.quant <- ev.3.melt %>% 
    dplyr::filter(!grepl("carrier|unused|reference|control|sc_x|sc_c", Raw.file_RI_celltype))
}



# take the median of repeat peptide observations 
ev.3.melt.quant.2 <- ev.3.melt.quant %>% dplyr::group_by(Raw.file_RI_celltype, Leading.razor.protein, SeqCharge) %>% 
  dplyr::summarise(medPep = median(value, na.rm = T))
# unmelt so we can do relative quant of the peptide level 
ev.4 <- reshape2::dcast(ev.3.melt.quant.2, Leading.razor.protein + SeqCharge  ~ Raw.file_RI_celltype, 
                        value.var = "medPep")

# normalize for sample loading
ev.6 <- columned(ev.4, start.num = 3, end.num = ncol(ev.4))
# normalize each row by the median 
ev.6$rowMedians <- apply(ev.6[,3:ncol(ev.6)],1, median, na.rm = TRUE)
ev.7a <- rowed(ev.6, start.num = 3, end.num = (ncol(ev.6)-1), static.num = ncol(ev.6))
#remove the column with the rowMedians
ev.7b <- ev.7a[,-ncol(ev.7a)]
# melt the row normalized peptide data 
ev.7b.melt <- reshape2::melt(ev.7b)

ev.7b.melt[ev.7b.melt == 0] <- NA
sc <- ev.7b.melt # this dataframe contains the reference / column / row normalized peptide values for each embryo 
sc$embryo <- sub(".* ", "", sc$variable)
sc$log2Value <- log2(sc$value)




# collapse the peptide data into protein level by the median across peptides for each sample 
sc.collapse <- sc %>% dplyr::group_by(variable, Leading.razor.protein) %>%
  dplyr::summarise(medPeptideLevel = median(log2Value, na.rm = T))

ggplot(data = sc.collapse,  aes(x = variable, y = medPeptideLevel)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# filter out rows with more than % missing data on the protein level - relaxed filter since we want to see differences within embryos
sc.collapse.d <- reshape2::dcast(sc.collapse, Leading.razor.protein ~ variable, value.var = "medPeptideLevel")
sc.collapse.1.missing <- sc.collapse.d

# column normalize the proteins x samples matrix 
sc.collapse.1 <- columned_log2(sc.collapse.1.missing, start.num = 2, end.num = ncol(sc.collapse.1.missing))
sc.collapse.1.m <- reshape2::melt(sc.collapse.1)
sc.collapse.1.m$embryo <- sub(".* ", "", sc.collapse.1.m$variable)
ggplot(data = sc.collapse.1.m,  aes(x = variable, y = value)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# row normalize by controlling for embryo-specific effects 
meanProteinEmbryo <- sc.collapse.1.m %>% dplyr::group_by(embryo, Leading.razor.protein) %>% 
  dplyr::summarise(meanProtein = mean(value, na.rm = T))
meanProteinEmbryo[meanProteinEmbryo == 0] <- NA
# match the protein level data with the mean protein data per embryo 
sc.collapse.2 <- merge(sc.collapse.1.m, meanProteinEmbryo, by = c("embryo", "Leading.razor.protein"))
sc.collapse.2$log2ValueAdj <- sc.collapse.2$value - sc.collapse.2$meanProtein
ggplot(data = sc.collapse.2,  aes(x = variable, y = log2ValueAdj)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# make the dataframe wide, proteins x samples 
sc.collapse.3 <- reshape2::dcast(sc.collapse.2, Leading.razor.protein ~ variable, value.var = "log2ValueAdj")
rownames(sc.collapse.3) <- sc.collapse.3$Leading.razor.protein
sc.collapse.3 <- sc.collapse.3[,-1]

# count number of proteins per column 
numProtein <- as.data.frame(colSums(!is.na(sc.collapse.3), na.rm = T))
colnames(numProtein) <- "n_proteins"
ggplot(data = numProtein, aes(x = rownames(numProtein), y= n_proteins)) + geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 12)) + 
  labs(title = "Number of Proteins per cell") + 
  geom_text(aes(label = n_proteins), vjust = -0.3)

# impute the protein data by knn, used only in plotting purposes
sc.collapse3.impute <- impute.knn(as.matrix(sc.collapse.3), k =3)$data





#------------------------------------------------------------------------------------------------------------.
############### Kmeans clustering + PCA  ########################################
#------------------------------------------------------------------------------------------------------------.
df <- sc.collapse.3 #%>% dplyr::filter(rownames(sc.collapse.3) %in% sig2c$Protein)
# missing data filter - let's observe each protein in x amount of cells
df.NA <- df[rowSums(is.na(df)) <= 0.7*ncol(df),]  
dist.mat <- cor(df.NA, use = "pairwise.complete.obs", method = "spearman")
#dist.mat <- dist(df[,2:ncol(df)], method = "euclidean")
# change the color gradient 
library(circlize)
col_fun = colorRamp2(c(-0.3, 0, 0.3), c("yellow", "white", "purple"))
col_fun(seq(-3, 3))
Heatmap(dist.mat, col = col_fun) #, km=2,  km=2, column_km = 2,
Heatmap(dist.mat) #, km=2,  km=2, column_km = 2,



ht_clust <- Heatmap(dist.mat, column_names_gp = gpar(fontsize = 5))
ht_clust <- draw(ht_clust)
dist.mat.o <- dist.mat[row_order(ht_clust),column_order(ht_clust)]

# dendrogram 
dd <- dist(t(df.NA), method = "euclidean")
hc <- hclust(dd)
hcd <- as.dendrogram(hc)
plot(hc, hang = -1, cex = 0.6)
ggdendrogram(hc, rotate = TRUE)
ddata_x <- dendro_data(hcd)
# make a labels dataframe that will be used to label and color the branches 
hLabs <- label(ddata_x)
hLabs$Embryo <- sub(".* ", "", hLabs$label)
# make a dataframe of which sister is in which cluster 
clust <- cutree(hc, k=2)
clust.df <- data.frame(label = names(clust), cluster = factor(clust))
# add the cluster information back to the labels dataframe
hLabs$cluster <- clust.df$cluster[match(hLabs$label, clust.df$label)]





#assess optimal number of clusters 
library(factoextra)
factoextra::fviz_nbclust(dist.mat, kmeans, method = "silhouette")





# kmeans clustering based on similarity metric 
set.seed(1)
fir <- kmeans(dist.mat, centers = 2, iter.max = 50, nstart = 25) 
fir.cluster <- as.data.frame(fir$cluster)
colnames(fir.cluster) <- "Cluster"

clusplot(dist.mat, fir$cluster)
# count the number of cells per embryo in this clustering
CluK <- as.data.frame(fir$cluster)
CluK$embryo <- sub(".* ", "", rownames(CluK))
colnames(CluK)[1] <- "cluster"
CluK$clusterName <- paste0("Cluster",CluK$cluster)
countCells <- CluK %>% 
  dplyr::group_by(embryo) %>% 
  dplyr::summarise(n_clu1 = sum(cluster == "1"), n_clu2 = sum(cluster == "2"))

# let's do the pca ourselves 
pca.imp.cor <-  dist.mat
sc.pca<-eigen(pca.imp.cor) 
scx<-as.data.frame(sc.pca$vectors)
colnames(scx)<-paste0("PC",1:ncol(scx))
scx$cells<-colnames(pca.imp.cor)
# Percent of variance explained by each principle component
pca_var <- sc.pca$values
percent_var<- pca_var/sum(pca_var)*100
plot(1:length(percent_var), percent_var, xlab="PC", ylab="% of variance explained")
pca.melt <- reshape2::melt(scx)
colnames(pca.melt)<-c("id","pc","value")
pca.display <- reshape2::dcast(pca.melt, id ~ pc, value.var = "value", fill=NA)
pca.display$cluster <- fir.cluster$Cluster[match(pca.display$id, rownames(fir.cluster))]
pca.display$Cluster <- paste0("Cluster",pca.display$cluster)
mycolor2 <- c("#cc0066", "#e2ac00", 
              #"orchid", "darkgoldenrod", 
              "indianred1", "red", "cadetblue", "orange", "khaki")
PCx<-"PC1"
PCy<-"PC2"
pca.display$Sample_Pair <- sub(".* ", "", pca.display$id)
ggscatter(pca.display, x =PCx, y = PCy, color = "Cluster", shape = "Cluster", size = 8, alpha=0.6) +
  xlab(paste0(PCx,"  (", round(percent_var[1],0),"%)")) +
  ylab(paste0(PCy,"  (", round(percent_var[2],0),"%)")) +
  font("ylab",size=30) +
  font("xlab",size=30) +
  font("xy.text", size=20) +
  scale_color_manual(values = mycolor2) +
  theme(legend.background = element_rect(fill = "gray97", color = "darkblue"), 
        legend.text = element_text(size = 25), legend.title = element_text(size = 25), legend.position = "bottom", 
        plot.title = element_text(size = 30)) 


#####------------------------------------------------------------------------------------------------------------.
##### Compare cut zygotes to mouse 2-cell blastomeres -----------------------------------------------------------
#####------------------------------------------------------------------------------------------------------------.
# calculate the Fold changes between sisters
df <- sc.collapse.3 %>% dplyr::filter(rownames(sc.collapse.3) %in% sig2c$Protein)
df.m <- reshape2::melt(as.matrix(df))
df.m$embryo <- sub(".* ", "", df.m$Var2)
df.m$cluster <- paste0("Cluster",fir.cluster$Cluster[match(df.m$Var2, rownames(fir.cluster))])
df.m.sel <- df.m 
FC <- df.m.sel %>% dplyr::group_by(embryo, Var1) %>% 
  dplyr::summarise(Clu1_Clu2 = value[cluster == "Cluster1"] -value[cluster == "Cluster2"])
ggplot(data = FC, aes(x = embryo, y = Clu1_Clu2, color = embryo)) + geom_violin() + geom_point() + 
  theme(legend.position = "none") + 
  labs(title = "Sister Fold changes for each protein")

# compare these fold changes in parthenotes to the fold changes in 2-cell stage 
NormalMouse <- read.csv("2021_12_07_SisterFoldChanges.csv")
# attach age to embryoID 
NormalMouse$embryoID_type <- paste0(NormalMouse$embryoID, ".", NormalMouse$type)
# let's focus on 2-cell stage 
NormalMouse2c <- NormalMouse %>% dplyr::filter(grepl("_2c", type))




# get the significantly differential proteins 
# obtain the proteins that are significantly differential 
FC.shared <- FC %>% dplyr::filter(Var1 %in% sig2c$Protein)
NormalMouse.shared <- NormalMouse2c %>% dplyr::filter(X %in% FC.shared$Var1)




# scatterplot of mean values across cells for each protein
# take the mean across single cells for each protein in the parthenotes data 
FC.shared.mean <- FC.shared %>% dplyr::group_by(Var1) %>% dplyr::summarise(Parthenotes = median(Clu1_Clu2, na.rm = T))
# obtain the FCs in the normal 2cell stage data 
NormalMouse.shared <- NormalMouse2c %>% dplyr::filter(X %in% FC.shared.mean$Var1)
# take the mean across single cells for each protein in the normal 2-cell stage data 
NormalMouse.shared.mean <- NormalMouse.shared %>% dplyr::group_by(X) %>% dplyr::summarise(Normal_2c = median(FC, na.rm = T))
# merge the datasets together 
MergedFC <- merge(FC.shared.mean, NormalMouse.shared.mean, by.x = "Var1", by.y = "X")
MergedFC$GeneName <- convertProts$uppercase[match(MergedFC$Var1, convertProts$Entry)]
# obtain spearman correlation and associated pvalue
testCOR <- cor.test(MergedFC$Parthenotes, MergedFC$Normal_2c, method = "spearman")
testCOR$estimate
testCOR$p.value

# plot a scatterplot 
ggplot(data = MergedFC, aes(x = Parthenotes, y = Normal_2c)) + geom_point(color = "purple3", alpha = 0.7, size = 4) + 
  theme_bw() + 
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 24), 
        plot.title = element_text(size = 20), plot.subtitle = element_text(size = 16)) + 
  labs(title = "All 2-cell stage embryos", 
       x = "Median Protein Fold Change in Cut Zygotes", y = "Median Protein Fold Change in 2-cell stage") + 
  geom_abline(intercept = 0, slope = 1, linetype= "dashed")





# let's get two plots arranged
# one with differences in PC between sisters as a barplot
# the other with distribution of correlations per cut zygote 
pca.diff <- pca.display[,c(1,2, 17, 18)]
pca.diff <- pca.diff %>% dplyr::group_by(Sample_Pair) %>% 
  dplyr::summarise(Difference = PC1[Cluster == "Cluster1"] - PC1[Cluster == "Cluster2"])
pca.diff <- pca.diff[order(-pca.diff$Difference),]
pca.diff$Sample_Pair <- factor(pca.diff$Sample_Pair, levels = pca.diff$Sample_Pair)
plo1 <- ggplot(data = pca.diff, aes(x = Sample_Pair, y = Difference)) + geom_bar(stat = "identity", fill = "purple3") + 
  theme_bw() + 
  theme(axis.text.y = element_text(size = 20, color = "black"), 
        axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_text(size = 24), 
        plot.title = element_text(size = 20), plot.subtitle = element_text(size = 16)) + 
  labs(x = "Zygote ID", y = "PC1\nDifference")
MergedFC.cor.m$Cut_Zygotes <- factor(MergedFC.cor.m$Var1, levels = pca.diff$Sample_Pair)
plo2 <- ggplot(data = MergedFC.cor.m,aes(y = value, x = Cut_Zygotes)) + 
  theme_bw() + 
  geom_beeswarm(size = 2, priority = "descending") + 
  stat_summary(fun = median, na.rm =T, geom = 'point', size = 9, shape = 23, 
               color = "red", fill = "purple3", alpha =0.5) + 
  theme(axis.text = element_text(size = 20, color = "black"), 
        axis.title = element_text(size = 24), 
        plot.title = element_text(size = 20), plot.subtitle = element_text(size = 16)) + 
  labs(x = "Zygote ID", y = "Correlations\nto Alpha - Beta")
ggarrange(plo1, NULL, plo2, ncol = 1, heights = c(1.5, 0.05,3))

ggarrange(plo1, plo2, ncol = 1, heights = c(1.5, 3))

