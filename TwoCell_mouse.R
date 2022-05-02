###---------------------------------------------------------------------------------------------------------.
###--------------------------------------- 2-cell embryo analysis -------------------------------------------------------
###---------------------------------------------------------------------------------------------------------.

###---------------------------------------------------------------------------------------------------------.
###--------------------------------------- parameter set-up -------------------------------------------------------
###---------------------------------------------------------------------------------------------------------.

# compare the late and early 2-cell stage samples 
WindowsPath <- "G:/My Drive/MS"
source("libraries_and_functions.R")


# let's try to combine the data post protein normalization 
protEarly <- read.csv("Exp4_early2cell_ProteinsXSamples_WithinEmbryo.csv")
protLate <- read.csv("Exp2_late2cell_ProteinsXSamples_WithinEmbryo.csv")
interProt <- intersect(rownames(protEarly), rownames(protLate))

mergeSC.byX<- merge(protEarly, protLate, by="X")
mergeSC.all<- merge(protEarly, protLate, by="X", all = T)

# protein to gene names 
prot2gene <- read.delim("uniprot_Mus_musculus_ProteinsTOGenes.tab")
convertProts <- read.delim("uniprot_Mus_musculus_ProteinsTOGenes.tab")
convertProts$LeadingGene <- sub(" .*", "", convertProts$Gene.names)
convertProts$uppercase <- toupper(convertProts$LeadingGene)

# load in mouse GO terms 
GOterms <- fread("gene_association.mgi.gz")
colnames(GOterms) <- c("DatabaseDes", "MGIMarker", "MouseMarkerSymbol", "NOTdesign", "GOtermID", "MGIRefAccesssionID",
                       "GOEviCode", "InferredFrom", 
                       "Ontology", "MouseMarkerName", "MouseMarkerSynonyms", "MouseMarkerType", "Taxon", "ModDate", 
                       "AssignedBy", "AnnotationExt", "GeneProduct")
# make all letters uppercase
GOterms$GeneSymbol <- toupper(GOterms$MouseMarkerSymbol)
GOannotations <- fread("go_terms.mgi.txt", header = F)
GOterms$annotation <- GOannotations$V3[match(GOterms$GOtermID, GOannotations$V2)]

# set the color bar for heatmap 
my_col2<-c("blue",rgb(0,0,1,0.5),"white",rgb(1,0,0,0.5),"red")


###---------------------------------------------------------------------------------------------------------.
###------------------- 2-cell stage cells - Form Clusters using kmeans ----------------------------
###---------------------------------------------------------------------------------------------------------.

df <- mergeSC.all
colnames(df)[1] <- "Protein"
dist.mat <- cor(df[,2:ncol(df)], use = "pairwise.complete.obs", method = "spearman") 
# kmeans clustering based on similarity metric 
set.seed(1)
fir <- kmeans(dist.mat, centers = 2, iter.max = 50, nstart = 25) 
fir.cluster <- as.data.frame(fir$cluster)
colnames(fir.cluster) <- "Cluster"
fir.cluster$Cluster <- paste0("Cluster", fir.cluster$Cluster)
Clusters <- reshape2::melt(df)
Clusters$cluster <- fir.cluster$Cluster[match(Clusters$variable, rownames(fir.cluster))]
colnames(Clusters)[1] <- "X"

###---------------------------------------------------------------------------------------------------------.
###------------------- 2-cell stage cells - Find Differential Proteins between clusters ----------------------------
###---------------------------------------------------------------------------------------------------------.

ProteinsAcrossClusters <- data.frame(matrix(nrow = length(unique(Clusters$X)), ncol = 2))
colnames(ProteinsAcrossClusters) <- c("Protein", "p.value")

for(i in 1:nrow(ProteinsAcrossClusters)) { 
  temp.prot <- Clusters %>% dplyr::filter(X == unique(Clusters$X)[i])
  temp.test <- kruskal.test(value ~ cluster, data = temp.prot)
  ProteinsAcrossClusters$Protein[i] <- unique(Clusters$X)[i]
  ProteinsAcrossClusters$p.value[i] <- temp.test$p.value
  print(i)
}

ProteinsAcrossClusters$q.value <- p.adjust(ProteinsAcrossClusters$p.value, method = "fdr")
ProteinsAcrossClusters.sig <- ProteinsAcrossClusters %>% dplyr::filter(q.value <= 0.05)

Clusters.sig <- Clusters %>% dplyr::filter(X %in% ProteinsAcrossClusters.sig$Protein)
Clusters.sig$GeneName <- convertProts$uppercase[match(Clusters.sig$X, convertProts$Entry)]

Clusters.sig.d <- reshape2::dcast(Clusters.sig, X ~ variable, value.var = "value")
rownames(Clusters.sig.d) <- Clusters.sig.d$X
Clusters.sig.d <- Clusters.sig.d[,-1]
Clusters.sig.d.i <- impute.knn(as.matrix(Clusters.sig.d), k =3)$data

# save this list 
#Clusters.sig.d1 <- reshape2::dcast(Clusters.sig, X + GeneName ~ variable, value.var = "value")
#write.csv(Clusters.sig.d1, "2-cell-stage-combined-sig-Prots_5percentFDR_2021_09_16_mergeALL.csv", row.names= F)



# analysis for obtaining list of proteins that are increased in abundance in alpha or in beta
ProteinsAcrossClusters <- data.frame(matrix(nrow = length(unique(Clusters$X)), ncol = 4))
colnames(ProteinsAcrossClusters) <- c("Protein", "p.value", "medLevel_Cluster1", "medLevel_Cluster2")

for(i in 1:nrow(ProteinsAcrossClusters)) { 
  temp.prot <- Clusters %>% dplyr::filter(X == unique(Clusters$X)[i])
  temp.test <- kruskal.test(value ~ cluster, data = temp.prot)
  ProteinsAcrossClusters$Protein[i] <- unique(Clusters$X)[i]
  ProteinsAcrossClusters$p.value[i] <- temp.test$p.value
  temp.1 <- temp.prot %>% dplyr::filter(cluster == "Cluster1")
  ProteinsAcrossClusters$medLevel_Cluster1[i] <- median(temp.1$value, na.rm= T)
  temp.2 <- temp.prot %>% dplyr::filter(cluster == "Cluster2")
  ProteinsAcrossClusters$medLevel_Cluster2[i] <- median(temp.2$value, na.rm= T)
  print(i)
}

ProteinsAcrossClusters$q.value <- p.adjust(ProteinsAcrossClusters$p.value, method = "fdr")
ProteinsAcrossClusters.sig <- ProteinsAcrossClusters %>% dplyr::filter(q.value <= 0.05)
ProteinsAcrossClusters.sig$Protein_cut <- sub("-.*", "", ProteinsAcrossClusters.sig$Protein)
ProteinsAcrossClusters.sig$Gene <- convertProts$uppercase[match(ProteinsAcrossClusters.sig$Protein_cut, 
                                                                convertProts$Entry)]
ProteinsAcrossClusters.sig$Gene[ProteinsAcrossClusters.sig$Protein_cut == "Q4VA45"] <- "ZFTA"
Proteins.Cluster1 <- ProteinsAcrossClusters.sig %>% dplyr::filter(medLevel_Cluster1 > medLevel_Cluster2)
Proteins.Cluster1$GreaterIN <- "Cluster1_Alpha"
Proteins.Cluster2 <- ProteinsAcrossClusters.sig %>% dplyr::filter(medLevel_Cluster1 < medLevel_Cluster2)
Proteins.Cluster2$GreaterIN <- "Cluster2_Beta"
tbsProteinsAcrossClusters <- rbind(Proteins.Cluster1, Proteins.Cluster2)
#write.csv(tbsProteinsAcrossClusters, "2021_11_03_2cellstage_SignificantProteins.csv")


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# FINAL FIGURE (MAYBE) V2 plot early and late stage separately Significant proteins
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
highnum = 1
lownum = -1
Clusters.sig.d.1 <- Clusters.sig.d.i
Clusters.sig.d.1[Clusters.sig.d.1 > highnum] <- highnum
Clusters.sig.d.1[Clusters.sig.d.1 < lownum] <- lownum
Clusters.sig.d.1.m <- reshape2::melt(Clusters.sig.d.1)
Clusters.sig.d.1.m$cluster <- Clusters$cluster[match(Clusters.sig.d.1.m$Var2, Clusters$variable)]
Clusters.sig.d.1.m$Stages <-  ifelse(grepl("406|407|408|409|410|411",Clusters.sig.d.1.m$Var2) == TRUE, "early", "late")

# obtain just the early stage embryos 
Clusters.sig.d.1.m.early <- Clusters.sig.d.1.m %>% dplyr::filter(Stages == "early")
meanProt_inCluster <- Clusters.sig.d.1.m.early %>% dplyr::group_by(cluster, Var1) %>% 
  dplyr::summarise(meanProt = mean(value, na.rm = T))
FCProt_Cluster <- meanProt_inCluster %>% dplyr::group_by(Var1) %>% 
  dplyr::summarise(FC = meanProt[cluster == "Cluster1"] - meanProt[cluster == "Cluster2"])
FCProt_Cluster <- FCProt_Cluster[order(FCProt_Cluster$FC),]
FCProt_Cluster$absFC <- abs(FCProt_Cluster$FC)
# also make a column that contains the clustering based on the protein FC 
FCProt_Cluster$proteincluster <- ifelse(FCProt_Cluster$FC > 0, "ProtCluster1", "ProtCluster2")
# now you can use this number to delinate the two clusters of proteins within a cluster of single cells 
Clusters.sig.d.1.m.early$proteincluster <- FCProt_Cluster$proteincluster[match(Clusters.sig.d.1.m.early$Var1, FCProt_Cluster$Var1)]
meanCell_inCluster <- Clusters.sig.d.1.m.early %>% dplyr::group_by(proteincluster, Var2) %>% dplyr::summarise(meanProt = mean(value, na.rm = T))
FCCell_Cluster <- meanCell_inCluster %>% dplyr::group_by(Var2) %>% 
  dplyr::summarise(FC = meanProt[proteincluster == "ProtCluster1"] - meanProt[proteincluster == "ProtCluster2"])
FCCell_Cluster <- FCCell_Cluster[order(FCCell_Cluster$FC),]
col.order.early <- match(as.character(FCCell_Cluster$Var2),colnames(Clusters.sig.d.1))
#row.order <- hclust(dist(Clusters.sig.d.1))$order
row.order <- match(as.character(FCProt_Cluster$Var1),rownames(Clusters.sig.d.1))
# select just the early embryos
Clusters.sig.d.2 <- Clusters.sig.d.1[row.order,col.order.early]
Clusters.sig.d.2.m <- reshape2::melt(as.matrix(Clusters.sig.d.2))

ggplot(data = Clusters.sig.d.2.m,
       aes(y = Var1, x = Var2, fill = (value))) + 
  geom_tile() +
  scale_fill_gradientn(colors = my_col2) + 
  theme_classic() + 
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(), 
        #axis.text.x = element_text(angle = 45),
        title = element_text(size = 16), legend.position = "bottom") + 
  #ggtitle("2-cell stage Differential Proteins") + 
  labs(y = "Protein", x = "Individual Blastomeres", caption = "Early Stages, imputed data", fill = "Log2(Protein Abundance)") + 
  geom_vline(xintercept = (length(col.order.early)/2)+0.5) +
  geom_rect(aes(xmin = 0.5, xmax = (length(col.order.early)/2)+0.5, 
                ymin = length(row.order)+3, ymax = length(row.order)+8), 
            fill = "orchid", alpha = 0.5) +  ##00BFC4
  geom_rect(aes(xmin = (length(col.order.early)/2)+0.5, xmax=length(col.order.early)+0.5, 
                ymin = length(row.order)+3, ymax = length(row.order)+8), 
            fill = "darkgoldenrod", alpha = 0.5) #"#F8766D"
#scale_fill_manual(values = c("orchid", "darkgoldenrod")) + 



# obtain the late stage embryos 
Clusters.sig.d.1.m.late <- Clusters.sig.d.1.m %>% dplyr::filter(Stages == "late")

meanProt_inCluster <- Clusters.sig.d.1.m.late %>% dplyr::group_by(cluster, Var1) %>% 
  dplyr::summarise(meanProt = mean(value, na.rm = T))
FCProt_Cluster <- meanProt_inCluster %>% dplyr::group_by(Var1) %>% 
  dplyr::summarise(FC = meanProt[cluster == "Cluster1"] - meanProt[cluster == "Cluster2"])
FCProt_Cluster <- FCProt_Cluster[order(FCProt_Cluster$FC),]
FCProt_Cluster$absFC <- abs(FCProt_Cluster$FC)
# what is the location of the smallest fold change
minFCProt <- which(FCProt_Cluster$absFC == min(FCProt_Cluster$absFC))
# also make a column that contains the clustering based on the protein FC 
FCProt_Cluster$proteincluster <- ifelse(FCProt_Cluster$FC > 0, "ProtCluster1", "ProtCluster2")
# now you can use this number to delinate the two clusters of proteins within a cluster of single cells 
Clusters.sig.d.1.m.late$proteincluster <- FCProt_Cluster$proteincluster[match(Clusters.sig.d.1.m.late$Var1, FCProt_Cluster$Var1)]
meanCell_inCluster <- Clusters.sig.d.1.m.late %>% dplyr::group_by(proteincluster, Var2) %>%
  dplyr::summarise(meanProt = mean(value, na.rm = T))
FCCell_Cluster <- meanCell_inCluster %>% dplyr::group_by(Var2) %>% 
  dplyr::summarise(FC = meanProt[proteincluster == "ProtCluster1"] - meanProt[proteincluster == "ProtCluster2"])
FCCell_Cluster <- FCCell_Cluster[order(FCCell_Cluster$FC),]
col.order.late <- match(as.character(FCCell_Cluster$Var2),colnames(Clusters.sig.d.1))
Clusters.sig.d.2 <- Clusters.sig.d.1[row.order,col.order.late]
Clusters.sig.d.2.m <- reshape2::melt(as.matrix(Clusters.sig.d.2))
ggplot(data = Clusters.sig.d.2.m,
       aes(y = Var1, x = Var2, fill = (value))) + 
  geom_tile() +
  scale_fill_gradientn(colors = my_col2) + 
  theme_classic() + 
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(), 
        #axis.text.x = element_text(angle = 45),
        title = element_text(size = 16), legend.position = "bottom") + 
  #ggtitle("2-cell stage Differential Proteins") + 
  labs(y = "Protein", x = "Individual Blastomeres", caption = "Late Stages, imputed data", fill = "Log2(Protein Abundance)") + 
  geom_vline(xintercept = (length(col.order.late)/2)+0.5) +
  geom_rect(aes(xmin = 0.5, xmax = (length(col.order.late)/2)+0.5, 
                ymin = length(row.order)+3, ymax = length(row.order)+8), 
            fill = "orchid", alpha = 0.5) + 
  geom_rect(aes(xmin = (length(col.order.late)/2)+0.5, xmax=length(col.order.late)+0.5, 
                ymin = length(row.order)+3, ymax = length(row.order)+8), 
            fill = "darkgoldenrod", alpha = 0.5)


