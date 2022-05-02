#------------------------------------------------------------------------------------.
################## SECTION 0: QC #################################
#------------------------------------------------------------------------------------.

humanGenesConvert <- read.delim("G:/My Drive/MS/Collaborators/Zernicka_Goetz/code/2020_DecPlate01/uniprot_Homo_sapiens_ProteinsTOGenes.tab")
humanGenesConvert$LeadingGene <- sub(" .*", "", humanGenesConvert$Gene.names)
humanGenesConvert$uppercase <- toupper(humanGenesConvert$LeadingGene)
# old human data, run with tmt 
UK <- read.csv("G:/My Drive/MS/Collaborators/Zernicka_Goetz/code/2020_DecPlate01/Exp8_humanUK_ProteinsXsamples.csv")
Clinic <- read.csv("G:/My Drive/MS/Collaborators/Zernicka_Goetz/code/2020_DecPlate01/Exp9_humanClinic_ProteinsXsamples.csv")
# human data from Lisa run with DIA methods (single shot blastomeres), searched with MBR
x <- read.delim("G:/My Drive/MS/Collaborators/Zernicka_Goetz/dat/Experiment13_LF_DIA_human/03_noLIB/LF_DIA_human_report.tsv")
ann <- read_excel("G:/My Drive/MS/Collaborators/Zernicka_Goetz/dat/Experiment13_LF_DIA_human/Exp_13_annotation.xlsx")
# filter on Q.value
x.fil <- x %>% dplyr::filter(Lib.Q.Value <= 0.01 & Lib.PG.Q.Value <= 0.01)
x.fil$SeqCharge <- paste0(x.fil$Modified.Sequence, x.fil$Precursor.Charge)
# count the number of unique protein groups per run
x.count1 <- x.fil %>% dplyr::group_by(Run) %>% dplyr::summarise(n_Genes = length(unique(Genes)))
# plot the number of counts 
ggplot(data = x.count1, aes(x = Run, y = n_Genes)) + geom_bar(stat="identity", fill = "darkolivegreen") + 
  geom_text(aes(label = n_Genes), vjust = -0.3, size = 3.5) + 
  labs(title = "Number of Proteins by GeneNames after Q.value filter") + theme_bw()
# count the number of unique peptides per run and then plot 
x.count2 <- x.fil %>% dplyr::group_by(Run) %>% dplyr::summarise(n_modSeqCharge = length(unique(SeqCharge)))
ggplot(data = x.count2, aes(x = Run, y = n_modSeqCharge)) + geom_bar(stat="identity", fill = "darkolivegreen3") + 
  geom_text(aes(label = n_modSeqCharge), vjust = -0.3, size = 3.5) + 
  labs(title = "Number of Peptides after Q.value filter") + 
  theme_bw()

#------------------------------------------------------------------------------------.
################## SECTION 1: NORMALIZE FROM MS1.Area for DIA  data #################################
#------------------------------------------------------------------------------------.

# now let's use the Ms1.Area column for peptide normalization first 
x.prot <- x.fil %>% dplyr::select(Run, Protein.Group, Genes, SeqCharge, Ms1.Area)
x.prot$Protein <- sub(";.*", "", x.prot$Protein.Group)
# transform df into proteins x cells 
x.prot1 <- reshape2::dcast(x.prot, Protein  + SeqCharge ~ Run, value.var = "Ms1.Area")
# transform zeroes into NA's 
x.prot2 <- x.prot1
x.prot2[x.prot2 == 0] <- NA
# log2 transform the raw intensities
x.prot2[,3:11] <- log2(x.prot2[,3:11])
# filter for missing data
pct.r = 0.5
x.prot2a <- c()
for(k in 1:nrow(x.prot2)){
  
  pct.na<-length(which(is.na(x.prot2[k,]))) / length(x.prot2[k,])
  if(pct.na <= pct.r){ x.prot2a<-rbind(x.prot2a, x.prot2[k,])}
  
}
# column normalize
x.prot2 <- columned_log2(x.prot2, start.num = 3, end.num =11)
# row normalize 
x.prot2$rowMean <- rowMeans(x.prot2[,3:11], na.rm = T)
x.prot2 <- rowed_log2(x.prot2, start.num = 3, end.num = 11, static.num = 12)
x.prot2 <- x.prot2[,-12]
x.prot3 <- reshape2::melt(x.prot2)
# collapse peptide level into protein level matrix by the median 
x.prot4 <- x.prot3 %>% dplyr::group_by(Protein, variable) %>% 
  dplyr::summarise(medPep = median(value, na.rm = T))
ggplot(data = x.prot4,  aes(x = variable, y = medPep)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(title = "Peptide Normalized Protein Matrix")
# median column normalize 
x.prot5 <- reshape2::dcast(data = x.prot4, Protein ~variable, value.var="medPep")
x.prot5 <- columned_log2(x.prot5, start.num = 2, end.num = 10)
x.prot6 <- reshape2::melt(data = x.prot5)
# find the mean value of each protein across cells within each embryo 
x.prot6$embryo <- ann$Embryo[match(x.prot6$variable, ann$Set)]
meanEmbryo <- x.prot6 %>% dplyr::group_by(embryo, Protein) %>% 
  dplyr::summarise(meanProtein = mean(value, na.rm = T))
meanEmbryo[meanEmbryo == 0] <- NA
# match the protein level data with the mean protein data per embryo 
x.prot7 <- merge(x.prot6, meanEmbryo, by = c("embryo", "Protein"))
# row normalize by controlling for embryo-specific effects 
x.prot7$log2ValueAdj <- x.prot7$value - x.prot7$meanProtein
ggplot(data = x.prot7,  aes(x = variable, y = log2ValueAdj)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(title = "Within Embryo Normalized Values")

# make the dataframe wide, proteins x samples 
x.prot8 <- reshape2::dcast(x.prot7, Protein ~ variable, value.var = "log2ValueAdj")
x.prot8 <- x.prot8 %>% dplyr::select(-wAP561)
x.prot8[x.prot8 == 0] <- NA
# add in the Gene Names 
x.prot8$Gene <- humanGenesConvert$uppercase[match(x.prot8$Protein, humanGenesConvert$Entry)]


#------------------------------------------------------------------------------------.
################## SECTION 4: Combine DIA with TMT #################################
#------------------------------------------------------------------------------------.
Clinic$Gene <- humanGenesConvert$uppercase[match(Clinic$X, humanGenesConvert$Entry)]
colnames(Clinic)[1] <- "Protein"
UK$Gene <- humanGenesConvert$uppercase[match(UK$X, humanGenesConvert$Entry)]
colnames(UK)[1] <- "Protein"
human <- plyr::join_all(list(x.prot8, UK, Clinic), by = c("Protein", "Gene"), type = "full")

#------------------------------------------------------------------------------------.
################## SECTION 4a: Global view of ALL human data _ k-means #################################
#------------------------------------------------------------------------------------.
# Calculate similarity metrics for all column pairs (default is Euclidean distance)
allprots <- human %>% dplyr::filter(!is.na(Gene))
allprots$Gene <- ifelse(allprots$Protein =="E9PAV3", "NACA_MUSCLE", allprots$Gene)
rownames(allprots) <- allprots$Gene
allprots <- allprots[,c(-1,-10)]

# OR LOAD from saved file: write.csv(allprots, "2021_09_13_human2cell_protsXcells.csv")
#allprots <- read.csv("2021_09_13_human2cell_protsXcells.csv")
#rownames(allprots) <- allprots$X
#allprots <- allprots[,c(-1)]
#filter for proteins that are present in at least x cells out of 26 
xCELLS <- 13 # 20 -> 224 proteins, 13 -> 369 proteins 
allprotsFIL <- data.frame()
for(i in 1:nrow(allprots)) { 
  if(sum(!is.na(allprots[i,])) >= xCELLS) { 
    allprotsFIL <- rbind(allprotsFIL, allprots[i,])
  }
}

# SELECT WHETHER TO filter on protein groups (like degradation or signaling?)
allprotsFIL <- allprots#[grepl("^PSM|^UB", rownames(allprots)),] # ^PSM|^UB  #^KIF|^RAB|^DYN|^MYO|^VAMP
#dist.mat<-as.matrix( dist(t(allprots), method = "canberra") )
dist.mat<- cor((allprotsFIL), method = "pearson", use = "pairwise.complete.obs")
# kmeans clustering based on similarity metric 
set.seed(1)
fir <- kmeans(dist.mat, centers = 2, iter.max = 50, nstart = 25) 
fir.cluster <- as.data.frame(fir$cluster)
colnames(fir.cluster) <- "Cluster"
# count the number of cells per embryo in this clustering
CluK <- as.data.frame(fir$cluster)
# count how many times an embryo is in each cluster 
CluK$embryo <- sub(".*_", "", rownames(CluK))
CluK$cleavage <- sub("([^_]*)_([^_]*)_([^_]*)", "\\2", rownames(CluK))
colnames(CluK)[1] <- "cluster"
# plot a dendrogram? 
forhclust <- allprotsFIL
forhclustCOLs <- data.frame(colnames(forhclust))
forhclustCOLs$shortName1 <- ifelse(grepl("_",forhclustCOLs$colnames.forhclust.) == TRUE, 
                                   sub(".*_","",forhclustCOLs$colnames.forhclust.), 
                                   ann$Embryo[match(forhclustCOLs$colnames.forhclust., ann$Set)])
forhclustCOLs$shortName2 <- ifelse(grepl("_",forhclustCOLs$shortName1) == TRUE, 
                                   sub(".*_","",forhclustCOLs$shortName1), 
                                   forhclustCOLs$shortName1)
forhclustCOLs$shortName2 <- as.numeric(forhclustCOLs$shortName2)



for(i in 1:nrow(forhclustCOLs)) {
  if((i %% 2) == 0 ) { 
    forhclustCOLs$sister[i] <- "B"
  }
  else{ forhclustCOLs$sister[i] <- "A" 
  }
  
}
forhclustCOLs$shortName3 <- ifelse(grepl("_",forhclustCOLs$colnames.forhclust.) == TRUE, 
                                   paste0("Embryo_",forhclustCOLs$shortName2, "_Sister_", forhclustCOLs$sister),
                                   paste0("Embryo_", forhclustCOLs$shortName2+10, "_Sister_", forhclustCOLs$sister))
forhclustCOLs$Embryo <- ifelse(grepl("_",forhclustCOLs$colnames.forhclust.) == TRUE, 
                               paste0("Embryo_",forhclustCOLs$shortName2),
                               paste0("Embryo_", forhclustCOLs$shortName2+10))

colnames(forhclust) <- forhclustCOLs$shortName3
dd <- dist(t(forhclust), method = "euclidean")
hc <- hclust(dd)
hcd <- as.dendrogram(hc)
plot(hc, hang = -1, cex = 0.6)
ggdendrogram(hc, rotate = TRUE)
ddata_x <- dendro_data(hcd)
# make a labels dataframe that will be used to label and color the branches 
hLabs <- label(ddata_x)
hLabs$Embryo <- forhclustCOLs$Embryo[match(hLabs$label, forhclustCOLs$shortName3)]
# make a dataframe of which sister is in which cluster 
clust <- cutree(hc, k=2)
clust.df <- data.frame(label = names(clust), cluster = factor(clust))
# add the cluster information back to the labels dataframe
hLabs$cluster <- clust.df$cluster[match(hLabs$label, clust.df$label)]
# color dendrogram 
ggplot(segment(ddata_x)) + 
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend))  + 
  # color by embryo: 
  #geom_text(data = hLabs, aes(label = label,x=x,y=0, color = Embryo), angle = 90, hjust = 1, size = 5) + 
  # color by cluster 
  geom_text(data = hLabs, aes(label = label,x=x,y=0, color = cluster), angle = 90, hjust = 1, size = 5) + 
  #geom_text(data = hLabs, aes(label = Embryo,x=x,y=0, color = cluster), angle = 90, hjust = 1, size = 5) + 
  ylim(-42, 65) + 
  theme_classic() + 
  theme(legend.position = "right") + 
  scale_color_manual(labels = c("Alpha", "Beta"),
                     values = c("#cc0066", "#e2ac00"))



#------------------------------------------------------------------------------------.
################## SECTION 4b: Global view -> Sig Proteins #################################
#------------------------------------------------------------------------------------.
# GOAL - build a heatmap plot illustrating the differences between alpha & beta cells 
# similar to what we have done with mouse 
humanPROT <- human[,c(1,10,2:9,11:28)]
humanPROT.m <- reshape2::melt(humanPROT)
humanPROT.m$cluster <- pca.display$clusterPCA[match(humanPROT.m$variable, pca.display$id)]

ProteinsAcrossClusters <- data.frame(matrix(nrow = length(unique(humanPROT.m$Protein)), ncol = 2))
colnames(ProteinsAcrossClusters) <- c("Protein", "p.value")

for(i in 1:nrow(ProteinsAcrossClusters)) { 
  temp.prot <- humanPROT.m %>% dplyr::filter(Protein == unique(humanPROT.m$Protein)[i])
  temp.prot.1 <- temp.prot %>% dplyr::filter(cluster == "Cluster1")
  temp.prot.2 <- temp.prot %>% dplyr::filter(cluster == "Cluster2")
  ProteinsAcrossClusters$Protein[i] <- unique(humanPROT.m$Protein)[i]
  if(sum(!is.na(temp.prot.1$value)) > 4 & sum(!is.na(temp.prot.2$value)) > 4) { 
    #temp.test <- kruskal.test(value ~ cluster, data = temp.prot)
    temp.test <- t.test(temp.prot.1$value, temp.prot.2$value)
    ProteinsAcrossClusters$p.value[i] <- temp.test$p.value
  }
  
  print(i)
}

# filter out proteins that were not tested (have a "NA" in place of pvalue)
ProteinsAcrossClusters.NA <- ProteinsAcrossClusters %>% dplyr::filter(!is.na(p.value))
# correct the pvalues using fdr correction 
ProteinsAcrossClusters.NA$q.value <- p.adjust(ProteinsAcrossClusters.NA$p.value, method = "fdr")
# filter proteins at 5% FDR
ProteinsAcrossClusters.NA.sig <- ProteinsAcrossClusters.NA %>% dplyr::filter(q.value < 0.05)
# obtain the GeneName/Symbol for each Uniprot ID 
ProteinsAcrossClusters.NA.sig$GeneName <- humanGenesConvert$uppercase[match(ProteinsAcrossClusters.NA.sig$Protein,humanGenesConvert$Entry)]
# obtain the Gene Description for each Uniprot ID 
ProteinsAcrossClusters.NA.sig$GeneDescription <- humanGenesConvert$Protein.names[match(ProteinsAcrossClusters.NA.sig$Protein,humanGenesConvert$Entry)]
# filter the protein-level data for the proteins that are significantly differential between the two cell classes
humanPROT.sig <- humanPROT %>% dplyr::filter(Protein %in% ProteinsAcrossClusters.NA.sig$Protein)
write.csv(humanPROT.sig, "2022_04_22_Human_sigProteins.csv")
# filter the MELTED protein-level data for the proteins that are significantly differential between the two cell classes
humanPROT.m.sig <- humanPROT.m %>% dplyr::filter(Protein %in% ProteinsAcrossClusters.NA.sig$Protein)
# Factor the blastomeres so that we can get the two clusters in order 
humanPROT.m.sig$dendro.name <- forhclustCOLs$shortName3[match(humanPROT.m.sig$variable, forhclustCOLs$colnames.forhclust.)]
humanPROT.m.sig$dendro.name <- factor(humanPROT.m.sig$dendro.name, levels = hLabs$label)
# plot the data, unimputed level: 
ggplot(data = humanPROT.m.sig,aes(x = dendro.name, y = Protein, fill = value)) + geom_tile() + 
  scale_fill_gradientn(colors = my_col2) + 
  geom_vline(xintercept = 13.5)

# impute the data for plotting purposes
humanPROT.sig.i <- impute.knn(as.matrix(humanPROT.sig[,3:28]), k =3)$data
# add back the protein names as rownames
rownames(humanPROT.sig.i) <- humanPROT.sig$Protein
# hclust the protein rows
row.order <- hclust(dist(humanPROT.sig.i))$order
humanPROT.sig.i <- humanPROT.sig.i[row.order,]
# melt the imputed dataframe
humanPROT.sig.i.m <- reshape2::melt(humanPROT.sig.i)
# again, factor the blastomeres via dendrogram names 
humanPROT.sig.i.m$dendro.name <- forhclustCOLs$shortName3[match(humanPROT.sig.i.m$Var2, forhclustCOLs$colnames.forhclust.)]
humanPROT.sig.i.m$dendro.name <- factor(humanPROT.sig.i.m$dendro.name, levels = hLabs$label)
# plot the data, imputed level: 
ggplot(data = humanPROT.sig.i.m,aes(x = dendro.name, y = Var1, fill = value)) + geom_tile() + 
  scale_fill_gradientn(colors = my_col2, name = "Log2 Protein \nAbudance") + 
  geom_vline(xintercept = 13.5) + theme_bw()+ 
  theme(#axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
    legend.position = "top", legend.key.height = unit(1, "cm"), legend.key.width = unit(2, "cm"), 
    legend.title =element_text(size = 15), legend.text = element_text(size = 15), 
    axis.text = element_blank(), axis.title = element_blank()) + 
  geom_rect(aes(xmin = 0.5, xmax = 13.5,
                ymin = length(row.order), ymax = length(row.order)+9), 
            fill = "#D60093", alpha = 0.5) +  
  geom_rect(aes(xmin = 13.5, xmax=26.5, 
                ymin = length(row.order), ymax = length(row.order)+9), 
            fill = "#E2AC00", alpha = 0.5) + 
  annotate(geom="text", x = 7, y = length(row.order)+4.5, label = "Cluster 1", color = "white", size = 10)+ 
  annotate(geom="text", x = 20, y = length(row.order)+4.5, label = "Cluster 2", color = "white", size = 10)

# order the unimputed data 
humanPROT.sig.1 <- humanPROT.sig
rownames(humanPROT.sig.1) <- humanPROT.sig.1$Protein
humanPROT.sig.1 <- humanPROT.sig.1[,-1]
humanPROT.sig.1 <- humanPROT.sig.1[row.order,]
humanPROT.sig.1.m <- reshape2::melt(humanPROT.sig.1)
humanPROT.sig.1.m$dendro.name <- forhclustCOLs$shortName3[match(humanPROT.sig.1.m$variable, forhclustCOLs$colnames.forhclust.)]
humanPROT.sig.1.m$dendro.name <- factor(humanPROT.sig.1.m$dendro.name, levels = hLabs$label)
ggplot(data = humanPROT.sig.1.m,aes(x = dendro.name, y = Gene, fill = value)) + geom_tile() + 
  scale_fill_gradientn(colors = my_col2, name = "Log2 Protein \nAbudance") + 
  geom_vline(xintercept = 13.5) + 
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1)) 


# build a legend axis bar
Expmap <- read_excel("G:/My Drive/MS/Collaborators/Zernicka_Goetz/dat/EXP_Map_Human.xlsx", sheet = 1)
humanPROT.sig.i.m$Raw.file <- sub("\\..*", "", humanPROT.sig.i.m$Var2)
humanPROT.sig.i.m$Data.Acq <- Expmap$Data.Acq[match(humanPROT.sig.i.m$Raw.file, Expmap$Raw.file)]
leg.groups <- as.data.frame(unique(humanPROT.sig.i.m[,c(2,4,5,6)]))
leg.groups.mat <- match(levels(humanPROT.sig.i.m$dendro.name), leg.groups$dendro.name)
leg.groups <- leg.groups[leg.groups.mat,]
plotHeat <- ggplot(data = humanPROT.sig.i.m,aes(x = dendro.name, y = Var1, fill = value)) + geom_tile() + 
  scale_fill_gradientn(colors = my_col2, name = "Log2 Protein \nAbudance") + 
  geom_vline(xintercept = 13.5) + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(), 
        legend.position = "right", legend.key.height = unit(2, "cm"), legend.key.width = unit(1, "cm"), 
        legend.title = element_blank(), legend.text = element_text(size = 18)) + 
  labs(x = '', y = '') 
plotHeat
leg <- ggplot(leg.groups,aes(y = 0, x = dendro.name)) + 
  geom_point(aes(color = Data.Acq), shape = 15, size = 7, show.legend = F) + 
  theme_classic() + 
  theme(axis.title = element_blank(), axis.line = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), 
        plot.margin = unit(c(0,0,0,0), "cm")) +  
  scale_color_manual(values = c("#fcad03", "#03b1fc")) #bafc06 #"#56fc03", "#c203fc" #03fcce #03ebfc
leg

plotHeat + annotation_custom(ggplotGrob(leg), 
                             xmin = 0.2, xmax = 26.7, 
                             ymin = -3, ymax = 1)


#------------------------------------------------------------------------------------.
################## SECTION 4b: Global view -> PSEA #################################
#------------------------------------------------------------------------------------.

# download the human GO terms 
GOterms <- read.delim("G:/My Drive/MS/Ribo/Aleks_LungCancer/Code/human.txt", sep = ",", header= F)
GONames <- read.delim("G:/My Drive/MS/Ribo/Aleks_LungCancer/Code/GoNames.txt", sep = ",", header= F)
colnames(GOterms) <- c("Uniprot", "GeneSymbol", "Qualifier", "GO_term","Reference", "Evidence", "With/From","UpLetter", "LongGeneName")
GOterms$annotation <- GONames$V3[match(GOterms$GO_term, GONames$V1)]

allprots.m <- reshape2::melt(as.matrix(allprots))
# label which cluster each cell belongs to
allprots.m$cluster <- pca.display$clusterPCA[match(allprots.m$Var2, pca.display$id)]
colnames(allprots.m)[1] <- "GeneName"
# GSEA
GOinDat <- GOterms %>% dplyr::filter(GeneSymbol %in% allprots.m$GeneName)
GOAcrossClusters <- data.frame(matrix(nrow = length(unique(GOinDat$GO_term)), ncol = 8))
colnames(GOAcrossClusters) <- c("GOterm", "Annotation", "p.value", "Proteins", "NumProteins", "NumProteinsInvolved", 
                                "median_Cluster1", "median_Cluster2")
GOAcrossClusters$GOterm <- unique(GOinDat$GO_term)
i<-2
for(i in 1:nrow(GOAcrossClusters)) { 
  temp.GO <- GOinDat %>% dplyr::filter(GO_term== GOAcrossClusters$GOterm[i])
  temp.GO2 <- GOterms %>% dplyr::filter(GO_term== GOAcrossClusters$GOterm[i])
  temp.prot <- allprots.m %>% dplyr::filter(GeneName %in% temp.GO$GeneSymbol)
  temp1 <- temp.prot %>% dplyr::filter( cluster == "Cluster1") 
  temp2 <- temp.prot %>% dplyr::filter( cluster == "Cluster2") 
  if(sum(!is.na(temp1$value)) >3 & sum(!is.na(temp2$value)) >3 ) { 
    temp.test <- kruskal.test(value ~ cluster, data = temp.prot)
    GOAcrossClusters$p.value[i] <- temp.test$p.value
    GOAcrossClusters$median_Cluster1[i] <- median(temp1$value, na.rm = T)
    GOAcrossClusters$median_Cluster2[i] <- median(temp2$value, na.rm = T)
    temp1med <- temp1 %>% dplyr::group_by(Var2) %>% dplyr::summarise(medLevel = median(value, na.rm = T))
    temp1med.t <- as.data.frame(as.matrix(t(temp1med[,-1])))
    colnames(temp1med.t) <- temp1med$Var2
    temp2med <- temp2 %>% dplyr::group_by(Var2) %>% dplyr::summarise(medLevel = median(value, na.rm = T))
    temp2med.t <- as.data.frame(as.matrix(t(temp2med[,-1])))
    colnames(temp2med.t) <- temp2med$Var2
    tempmed <- cbind(temp1med.t, temp2med.t)
    GOAcrossClusters[i,9:34] <- tempmed
  } 
  
  GOAcrossClusters$GOterm[i] <- unique(temp.GO$GO_term)
  GOAcrossClusters$Annotation[i] <- unique(temp.GO$annotation)
  GOAcrossClusters$Proteins[i] <- paste(unique(temp.prot$GeneName), collapse = " ")
  GOAcrossClusters$NumProteins[i] <- length(unique(temp.prot$GeneName))
  GOAcrossClusters$NumProteinsInvolved[i] <- length(temp.GO2$GeneSymbol)
  print(i)
}

#------------------------------------------------------------------------------------.
################## SECTION 4b: Global view  -> PSEA PLOTTING ver2 #################################
#------------------------------------------------------------------------------------.

# brand new analysis: Feb 7, 2022
GOAcrossClusters$fraction <- GOAcrossClusters$NumProteins / GOAcrossClusters$NumProteinsInvolved
GOAcrossClusters$Clus1_vs_Clus2 <- GOAcrossClusters$median_Cluster1 - GOAcrossClusters$median_Cluster2
# go through to obtain GO terms that are not repeated in terms of proteins used 
# keep only unique protein lists
GOAcrossStages.NA.pre <- GOAcrossClusters %>% dplyr::distinct(Proteins, .keep_all = TRUE)
GOAcrossClusters.NA <- GOAcrossStages.NA.pre %>% dplyr::filter(!is.na(p.value)) %>% dplyr::filter(NumProteins >1)
GOAcrossClusters.NA$q.value <- p.adjust(GOAcrossClusters.NA$p.value, method = "fdr")
GOAcrossClusters.NA.sig <- GOAcrossClusters.NA %>% dplyr::filter(q.value < 0.05)
GO2plot <- GOAcrossClusters.NA.sig %>% 
  #CHOOSE: degradation or transport: ubiquit|proteas|neddyl|catab OR endocy|vesicl|clathrin
  dplyr::filter(grepl("endocy|vesicl|clathrin", Annotation) & !grepl("regulation", Annotation)) %>% 
  dplyr::filter(fraction > 0.15 & NumProteins > 2)
GO2plot.m <- reshape2::melt((GO2plot %>% dplyr::select(Annotation, contains("wAP"))))
GO2plot.m$cluster <- pca.display$cluster[match(GO2plot.m$variable, pca.display$id)]
GO2plot.m$group <- GOAcrossClusters.sig$group[match(GO2plot.m$Annotation, GOAcrossClusters.sig$Annotation)]
GO2plot.mPLOT <- GO2plot.m  
ggplot(data = GO2plot.mPLOT, aes(y = Annotation, x = value, fill = cluster)) + geom_boxplot() + 
  theme_bw() + 
  coord_cartesian(xlim = c(-0.3, 0.3)) + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18),
        legend.text = element_text(size = 20), legend.position = "bottom", 
        legend.title = element_text(size = 20)) + 
  scale_fill_manual(labels = c("Alpha", "Beta"),
                    #values = c("orchid", "darkgoldenrod")) + 
                    values = c("#cc0066", "#e2ac00"))+
  labs(x = "GO term abundance, log[2]", fill = "Cell Classification")
# find the proteins associated with GO term and plot that 
Go.plot1 <- c()
i<-2
for(i in 1:length(unique(GO2plot$Annotation))) {
  temp.GO <- GOinDat %>% dplyr::filter(GO_term == GO2plot$GOterm[i])
  temp.gene <- temp.GO$GeneSymbol
  temp.prot <- allprots.m %>% dplyr::filter(GeneName %in% temp.gene)
  temp.prot$GO_Annotation <- GO2plot$Annotation[i]
  temp.prot$GO_ID <- GO2plot$GOterm[i]
  if(i == 1) { Go.plot1 <- temp.prot}
  if(i != 1) { Go.plot1 <- rbind(Go.plot1, temp.prot)}
}
ggplot(data = Go.plot1, aes(y = GO_Annotation, x = value, fill = cluster)) + geom_boxplot() + 
  theme_bw() + 
  coord_cartesian(xlim = c(-0.5, 0.5)) + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18),
        legend.text = element_text(size = 20), legend.position = "bottom", 
        legend.title = element_text(size = 20)) + 
  scale_fill_manual(labels = c("Alpha", "Beta"),
                    #values = c("orchid", "darkgoldenrod")) + 
                    values = c("#cc0066", "#e2ac00"))+
  labs(x = "GO term abundance, log[2]", fill = "Cell Classification")


# rather, we should find the fold changes between the proteins for each sister-pair
# after obtaining GO terms with high FC
GOAcrossClusters.NA <- GOAcrossStages.NA.pre %>% dplyr::filter(!is.na(p.value)) %>% 
  #dplyr::filter(fraction > 0.5)
  dplyr::filter(NumProteins >=2)
GOAcrossClusters.NA$q.value <- p.adjust(GOAcrossClusters.NA$p.value, method = "fdr")
GOAcrossClusters.NA.sig <- GOAcrossClusters.NA %>% dplyr::filter(q.value < 0.05)
GOAcrossClusters.NA.sig2 <- GOAcrossClusters.NA.sig %>% dplyr::filter(abs(Clus1_vs_Clus2) >= 0.2) %>% 
  dplyr::filter(grepl("endocy|vesicl|clathrin|ubiquit|proteas|neddyl|catab", Annotation) & !grepl("regulation|protein K", Annotation))
# melt the dataframe with the median GO levels in each blastomere 
GOAcrossClusters.NA.sig3 <- reshape2::melt(GOAcrossClusters.NA.sig2[,c(1,2,9:34)])
# add the sister/cluster information 
GOAcrossClusters.NA.sig3$EmbryoID <- forhclustCOLs$Embryo[match(GOAcrossClusters.NA.sig3$variable,forhclustCOLs$colnames.forhclust.)]
GOAcrossClusters.NA.sig3$cluster <- fir.cluster$Cluster[match(GOAcrossClusters.NA.sig3$variable, row.names(fir.cluster))]
# obtain the FCs between sister blastomeres 
GOAcrossClusters.NA.sig4 <- GOAcrossClusters.NA.sig3 %>% dplyr::group_by(GOterm, Annotation, EmbryoID) %>% 
  dplyr::summarise(FC = value[cluster == "2"] - value[cluster == "1"])
GOAcrossClusters.NA.sig4$Annotation_trim <- str_trunc(GOAcrossClusters.NA.sig4$Annotation, 30, "right")
GOAcrossClusters.NA.sig5 <- GOAcrossClusters.NA.sig4 %>% dplyr::group_by(Annotation_trim) %>% 
  dplyr::summarise(medianFC = median(FC, na.rm = T))
GOAcrossClusters.NA.sig5 <- GOAcrossClusters.NA.sig5[order(GOAcrossClusters.NA.sig5$medianFC),]
GOAcrossClusters.NA.sig4$Annotation_trim <- factor(GOAcrossClusters.NA.sig4$Annotation_trim, 
                                                   levels = GOAcrossClusters.NA.sig5$Annotation_trim)
GOAcrossClusters.NA.sig5 <- GOAcrossClusters.NA.sig4 %>% dplyr::filter(!(grepl("exonucleo", Annotation_trim)))
# plot distributions of FCs per GO term 
ggplot(data = GOAcrossClusters.NA.sig5, aes(y = Annotation_trim, x = FC)) + 
  geom_boxplot(aes(fill = Annotation_trim)) + 
  theme_bw() + 
  coord_cartesian(xlim = c(-1.2, 1.2)) + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  theme(axis.text = element_text(size = 24), axis.title = element_blank(),
        legend.text = element_text(size = 20), legend.position = "none", 
        legend.title = element_text(size = 20)) + 
  #scale_fill_brewer(palette = "RdYlGn") +
  scale_fill_brewer(palette = "Set2") + 
  labs(x = "Alpha / Beta,GO term abundance, log[2]")



#------------------------------------------------------------------------------------.
################## SECTION 4c: Compare to mouse (proteins) - VERSION 3   #################################
#------------------------------------------------------------------------------------.
# load in mouse uniprot --> GeneNames 
convertProts <- read.delim("G:/My Drive/MS/Collaborators/Zernicka_Goetz/code/2020_DecPlate01/uniprot_Mus_musculus_ProteinsTOGenes.tab")
convertProts$LeadingGene <- sub(" .*", "", convertProts$Gene.names)
convertProts$uppercase <- toupper(convertProts$LeadingGene)

# load in all the proteins at the mouse 2-cell stage 
x1.mouse <- read.csv("G:/My Drive/MS/Collaborators/Zernicka_Goetz/code/2020_DecPlate01/Exp4_early2cell_ProteinsXSamples_WithinEmbryo.csv")
x2.mouse <- read.csv("G:/My Drive/MS/Collaborators/Zernicka_Goetz/code/2020_DecPlate01/Exp2_late2cell_ProteinsXSamples_WithinEmbryo.csv")
df.mouse <- plyr::join_all(list(x1.mouse, x2.mouse), by = "X", type = "full")
df.mouse$Gene <- convertProts$uppercase[match(df.mouse$X, convertProts$Entry)]
df.mouse$Gene <- ifelse(df.mouse$X == "P70670", "NACA.m", df.mouse$Gene)
#P70670

# intersect the lists
# Option 2 : All intersecting proteins between datasets, regardless of significant 
#E9PAV3
human$Gene <- ifelse(human$Protein == "E9PAV3", "NACA.m", human$Gene)
human1 <- human %>% dplyr::filter(!is.na(Gene)) %>% dplyr::filter(Gene != "")
mouse <- df.mouse %>% dplyr::filter(!is.na(Gene)) %>% dplyr::filter(Gene != "")
human1 <- human1 %>% dplyr::filter(Gene %in% mouse$Gene)
mouse <- mouse %>% dplyr::filter(Gene %in% human1$Gene)
merged2plot <- merge(human1, mouse, by = "Gene")
merged2plot.A <- merged2plot #%>% dplyr::filter(Gene %in% df.mouse3$GeneName)
noccur <- data.frame(table(merged2plot.A$Gene))

# prepare matrix for pairwise correlations between cells 
merged2plot.1 <- merged2plot.A
colnames(merged2plot.1)[3:28] <- paste0(colnames(merged2plot.1)[3:28], "_h")
colnames(merged2plot.1)[30:101] <- paste0(colnames(merged2plot.1)[30:101], "_m")
merged2plot.2 <- merged2plot.1[,c(1,3:28,30:101)]
merCOR <- cor(merged2plot.2[,2:99], method = "pearson", use = "pairwise.complete.obs")
#Heatmap(merCOR)
# this correlation heatmap does not make any sense 
# perhaps let's have mouse on x axis and human on y axis
merCOR.m <- reshape2::melt(merCOR)
merCOR.m <- merCOR.m %>% dplyr::filter(grepl("_m", Var1)) %>% dplyr::filter(grepl("_h", Var2))
merCOR.m.d <- reshape2::dcast(merCOR.m, Var1 ~Var2, value.var = "value")
rownames(merCOR.m.d) <- merCOR.m.d$Var1
merCOR.m.d <- merCOR.m.d[,-1]
row.order <- hclust(dist(merCOR.m.d))$order
merCOR.m.d <- merCOR.m.d[row.order,]
merCOR.m.d.m <- reshape2::melt(as.matrix(merCOR.m.d))
# Factor the blastomeres so that we can get the two clusters in order 
merCOR.m.d.m$human <- sub("_h.*", "", merCOR.m.d.m$Var2)
merCOR.m.d.m$dendro.name.huma <- forhclustCOLs$shortName3[match(merCOR.m.d.m$human, forhclustCOLs$colnames.forhclust.)]
merCOR.m.d.m$dendro.name.huma <- factor(merCOR.m.d.m$dendro.name.huma, levels = hLabs$label)

ggplot(data = merCOR.m.d.m, aes(y = dendro.name.huma, x = Var1, fill = value)) + geom_tile() + 
  scale_fill_gradientn(colors = my_col3) + 
  geom_hline(yintercept = 13.5) + geom_vline(xintercept = ((72/2)+0.5)) + 
  theme(#axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text = element_blank(),
    legend.position = "right", legend.key.height = unit(2, "cm"), legend.key.width = unit(1, "cm"), 
    legend.title = element_blank(), legend.text = element_text(size = 18), 
    axis.title = element_text(size = 26)) + 
  labs(x = "Mouse Blastomeres", y = "Human Blastomeres")

# we should actually use the clusters we use in the entire paper probably: 
mouseClusDF <- read.csv("G:/My Drive/MS/Collaborators/Zernicka_Goetz/code/2020_DecPlate01/2022_01_11_fircluster.csv")
mouseClusDF <- mouseClusDF[order(mouseClusDF$Cluster),]
merCOR.m <- reshape2::melt(merCOR)
merCOR.m <- merCOR.m %>% dplyr::filter(grepl("_m", Var1)) %>% dplyr::filter(grepl("_h", Var2))
# Factor the HUMAN blastomeres so that we can get the two clusters in order 
merCOR.m$human <- sub("_h.*", "", merCOR.m$Var2)
merCOR.m$dendro.name.huma <- forhclustCOLs$shortName3[match(merCOR.m$human, forhclustCOLs$colnames.forhclust.)]
merCOR.m$dendro.name.huma <- factor(merCOR.m$dendro.name.huma, levels = hLabs$label)
# Factor the MOUSE blastomeres so that we can get the two clusters in order 
merCOR.m$mouse <- sub("_m.*", "", merCOR.m$Var1)
mouseClusDF.in <- mouseClusDF %>% dplyr::filter(X %in% merCOR.m$mouse)
merCOR.m$mouse <- factor(merCOR.m$mouse, levels = rev(mouseClusDF.in$X))
my_col4 = c("purple", "white", "darkolivegreen4")
ggplot(data = merCOR.m, aes(x = dendro.name.huma, y = mouse, fill = value)) + geom_tile() + 
  scale_fill_gradientn(colors = my_col4) + 
  geom_vline(xintercept = 13.5) + geom_hline(yintercept = ((72/2)+0.5)) + 
  theme(#axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text = element_blank(),
    axis.title = element_blank(),
    #axis.title = element_text(size = 26),
    #axis.text.x = element_blank(),
    #axis.text.y = element_text(),
    legend.position = "right", legend.key.height = unit(2, "cm"), legend.key.width = unit(1, "cm"), 
    legend.title = element_blank(), legend.text = element_text(size = 18)) + 
  labs(y = "Mouse Blastomeres", x = "Human Blastomeres")


#------------------------------------------------------------------------------------.
################## SECTION 4d: Compare to mouse (GO terms) version 2 #################################
#------------------------------------------------------------------------------------.
HumanGO <- GOAcrossClusters.NA.sig
HumanGO.a <- HumanGO #%>% dplyr::filter(NumProteins > 4)
# read in significantly differential GO terms from mouse 
MouseGO <- read.csv("G:/My Drive/MS/Collaborators/Zernicka_Goetz/code/2020_DecPlate01/2021_09_14_mouse2cell_5%FDRsigGOterms.csv")
MouseGOlevels <- read.csv("G:/My Drive/MS/Collaborators/Zernicka_Goetz/code/2020_DecPlate01/2021_09_14_mouse2cell_5%FDRsigGOterms_MedianValues.csv")
colnames(MouseGOlevels)[2:3] <- c("GOterm", "Annotation")
MouseGO1 <- merge(MouseGO, MouseGOlevels, by = c("GOterm", "Annotation"))
MouseGO1.a <- MouseGO1 #%>% dplyr::filter(NumProteins > 4)
interGO <- intersect(HumanGO.a$GOterm, MouseGO1.a$GOterm)
interGO2 <- interGO


HumanGO.fil <- HumanGO %>% dplyr::filter(GOterm %in% interGO2) %>% 
  dplyr::filter(!(grepl("sperm|zona", Annotation)))
HumanGO.fil2 <- HumanGO.fil %>% dplyr::select(GOterm, Annotation, contains("AP"))
MouseGO.fil <- MouseGO1 %>% dplyr::filter(GOterm %in% interGO2) %>% 
  dplyr::filter(!(grepl("sperm|zona", Annotation)))
MouseGO.fil2 <- MouseGO.fil %>% dplyr::select(GOterm, Annotation, contains("AP"))


MergeGO <- merge(HumanGO.fil2, MouseGO.fil2, by = c("GOterm"))
MergeGO2 <- MergeGO[,c(1,2,29,3:28,30:ncol(MergeGO))]
MergeGO3 <- MergeGO2 %>% dplyr::filter(Annotation.x != Annotation.y) # this looks ok, can continue 
MergeGO4 <- MergeGO2
rownames(MergeGO4) <- MergeGO2$Annotation.x
MergeGO4 <- MergeGO4[,-c(1,2,3)]
Heatmap(as.matrix(MergeGO4))

# download the cluster information for mouse 
mouseClusDF <- read.csv("G:/My Drive/MS/Collaborators/Zernicka_Goetz/code/2020_DecPlate01/2022_01_11_fircluster.csv")
# order the clusters 
mouseClusDF <- mouseClusDF[order(mouseClusDF$Cluster),]
# get cluster 1 and cluster 2 samples: 
mouseClusDF.1 <- mouseClusDF %>% dplyr::filter(Cluster == 1)
mouseClusDF.2 <- mouseClusDF %>% dplyr::filter(Cluster == 2)
# obtain just the mouse samples in the merged GO dataframe
MergeGO4.mouse <- MergeGO4 %>% dplyr::select(contains("Cluster"))
# obtain blastomeres in separate cluster and hierarchial cluster the columns
MergeGO4.mouse.selecting <- MergeGO4.mouse
colnames(MergeGO4.mouse.selecting) <- sub("^.*?\\.","", colnames(MergeGO4.mouse.selecting))
MergeGO4.mouse.1 <- MergeGO4.mouse.selecting %>% dplyr::select(any_of(mouseClusDF.1$X))
MergeGO4.mouse.2 <- MergeGO4.mouse.selecting %>% dplyr::select(any_of(mouseClusDF.2$X))
column.order.1 <- hclust(dist(t(MergeGO4.mouse.1), method = "euclidean"))$order
column.order.2 <- hclust(dist(t(MergeGO4.mouse.2), method = "euclidean"))$order
MergeGO4.mouse.1 <- MergeGO4.mouse.1[,column.order.1]
MergeGO4.mouse.2 <- MergeGO4.mouse.2[,column.order.2]
MergeGO4.mouse.3 <- cbind(MergeGO4.mouse.2, MergeGO4.mouse.1)
# obtain a rational row order 
row.orderGO <- hclust(dist((MergeGO4.mouse.3), method = "euclidean"))$order
# implement the new row order
MergeGO4.mouse.3 <- MergeGO4.mouse.3[row.orderGO,]
MergeGO4.mouse.m <- reshape2::melt(as.matrix(MergeGO4.mouse.3))

# obtain the 2-cell stage blastomeres that are in the clustering pattern dataframe 
mouseClusDF.in <- mouseClusDF %>% dplyr::filter(X %in% MergeGO4.mouse.m$Var2)
# match the cluster information 
MergeGO4.mouse.m$cluster <- mouseClusDF.in$Cluster[match(MergeGO4.mouse.m$Var2, mouseClusDF.in$X)]
# obtain the median value in each cluster
mouseMedianGo <- MergeGO4.mouse.m %>% dplyr::group_by(cluster, Var1) %>% 
  dplyr::summarise(medianGO = median(value, na.rm = T))
# obtain the median fold change
mouseMedianGo.1 <- mouseMedianGo %>% dplyr::group_by(Var1) %>% 
  dplyr::summarise(medianGOfc = medianGO[cluster == "1"] - medianGO[cluster == "2"])
mouseMedianGo.1$absmedianGOfc <- abs(mouseMedianGo.1$medianGOfc)
mouseMedianGo.1 <- mouseMedianGo.1[order(mouseMedianGo.1$absmedianGOfc),]
mouseMedianGo.2 <- mouseMedianGo.1[1:50,]
MergeGO4.mouse.m.2 <- MergeGO4.mouse.m %>% dplyr::filter(Var1 %in% mouseMedianGo.2$Var1)
#MergeGO4.mouse.m.2$blastomere <- factor(MergeGO4.mouse.m.2$blastomere, levels = rev(mouseClusDF.in$X))
numLimit <- 0.5
MergeGO4.mouse.m.2$value[MergeGO4.mouse.m.2$value >=numLimit] <- numLimit
MergeGO4.mouse.m.2$value[MergeGO4.mouse.m.2$value <=-numLimit] <- -numLimit
MergeGO4.mouse.m.2.d <- reshape2::dcast(data =MergeGO4.mouse.m.2,Var1 ~Var2,value.var = "value")
# plot 
ggplot(data = MergeGO4.mouse.m.2,
       aes(x = Var2, y = Var1, fill = (value))) + 
  geom_tile() +
  scale_fill_gradientn(colors = my_col2) +
  geom_vline(xintercept = 36.5) + 
  theme(legend.position = "bottom",legend.key.height = unit(0.75, "cm"), legend.key.width = unit(1.25, "cm"),
        legend.text = element_text(size = 12), legend.title = element_blank(),
        axis.text.x =element_blank()) + 
  labs(title = paste("Data\nShared GO terms"), fill = "Median GO",
       x = paste("Blastomeres"), y = "GO terms")




plotTOBE2 <- MergeGO4 %>% dplyr::select(!contains("Cluster")) %>% 
  dplyr::filter(row.names(MergeGO4) %in% mouseMedianGo.2$Var1)
row.orderGO2 <- match( MergeGO4.mouse.m.2.d$Var1, rownames(plotTOBE2))
plotTOBE2.o <- plotTOBE2[row.orderGO2,]
plotTOBE2.m <- reshape2::melt(as.matrix(plotTOBE2.o))
# Factor the HUMAN blastomeres so that we can get the two clusters in order 
plotTOBE2.m$dendro.name.huma <- forhclustCOLs$shortName3[match(plotTOBE2.m$Var2, forhclustCOLs$colnames.forhclust.)]
plotTOBE2.m$dendro.name.huma <- factor(plotTOBE2.m$dendro.name.huma, levels = hLabs$label)


plotTOBE2.m$value[plotTOBE2.m$value > 1] <- 1
plotTOBE2.m$value[plotTOBE2.m$value < -1] <- -1
ggplot(data = plotTOBE2.m,
       aes(y = Var1, x = dendro.name.huma, fill = (value))) + 
  geom_tile() +
  scale_fill_gradientn(colors = my_col2) + 
  geom_vline(xintercept = 13.5) + #theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(legend.position = "bottom",legend.key.height = unit(0.75, "cm"), legend.key.width = unit(1.25, "cm"),
        legend.text = element_text(size = 12), legend.title = element_blank(),
        axis.text.x =element_blank()) +
  labs(title = paste("Data\nShared GO terms"), 
       #subtitle = "Using Canberra for clustering", 
       x = paste( "Blastomeres"), y = "GO terms")
