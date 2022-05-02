###---------------------------------------------------------------------------------------------------------.
###--------------------------------------- Mouse Blastomere Analysis  -------------------------------------------------------
###---------------------------------------------------------------------------------------------------------.

###---------------------------------------------------------------------------------------------------------.
###--------------------------------------- parameter set-up -------------------------------------------------------
###---------------------------------------------------------------------------------------------------------.
WindowsPath <- "G:/My Drive/MS"


# GO term analysis files for mouse
GOterms <- fread("gene_association.mgi.gz")
colnames(GOterms) <- c("DatabaseDes", "MGIMarker", "MouseMarkerSymbol", "NOTdesign", "GOtermID", "MGIRefAccesssionID",
                       "GOEviCode", "InferredFrom", 
                       "Ontology", "MouseMarkerName", "MouseMarkerSynonyms", "MouseMarkerType", "Taxon", "ModDate", 
                       "AssignedBy", "AnnotationExt", "GeneProduct")
# make all letters upercase
GOterms$GeneSymbol <- toupper(GOterms$MouseMarkerSymbol)
GOannotations <- fread("go_terms.mgi.txt", header = F)
GOterms$annotation <- GOannotations$V3[match(GOterms$GOtermID, GOannotations$V2)]

# prot to gene names 
prot2gene <- read.delim("uniprot_Mus_musculus_ProteinsTOGenes.tab")
convertProts <- read.delim("uniprot_Mus_musculus_ProteinsTOGenes.tab")
convertProts$LeadingGene <- sub(" .*", "", convertProts$Gene.names)
convertProts$uppercase <- toupper(convertProts$LeadingGene)


# load in the four current datasets that use within embryo normalization 
x1 <- read.csv("Exp4_early2cell_ProteinsXSamples_WithinEmbryo.csv")
x2 <- read.csv("Exp2_late2cell_ProteinsXSamples_WithinEmbryo.csv")
x3 <- read.csv("Exp3_4c_ProteinsXsamples_WithinEmbryo.csv")
x4 <- read.csv("Exp6_4cMCherry_ProteinsXsamples_WithinEmbryo.csv")
x5 <- read.csv("Exp10_4c_ProteinsXsamples_WithinEmbryo.csv")


# keep names of 2 cell and 4 cell stage embryos in a dataframe
twocell <- as.data.frame(as.matrix(c(colnames(x1)[2:ncol(x1)], colnames(x2)[2:ncol(x2)])))
twocell$type <- ifelse(grepl("406|407|408|409|410|411",twocell$V1) == TRUE, "early_2c", "late_2c")
fourcell <- as.data.frame(as.matrix(c(colnames(x3)[2:ncol(x3)], colnames(x4)[2:ncol(x4)], 
                                      colnames(x5)[2:ncol(x5)])))
fourcell$type <- "_4c"
allcells <- rbind(twocell, fourcell)

# keep all proteins quantified across all datasets
df <- plyr::join_all(list(x1, x2, x3, x4, x5), by = "X", type = "full")
#write.csv(df, "MergedDF_2c_4c_proteinXblastomeres.csv")

# melt the dataframe
df.m <- reshape2::melt(df)

# create a storage df that will contain id information 
dfStore <- data.frame(matrix(ncol = 2, nrow = length(unique(df.m$variable))))
colnames(dfStore) <- c("variable", "label")
dfStore$variable <- unique(df.m$variable)
dfStore$label <- sub(".*\\.", "", dfStore$variable)
dfStore$label2 <- gsub(".*sc_", "", dfStore$label)
dfStore$label3 <- paste(sub("\\..*", "", dfStore$variable),dfStore$label2)
dfStore$type <- allcells$type[match(dfStore$variable, allcells$V1)]
dfStore$cleavage <- sub("_.*", "", dfStore$label2)
dfStore$Raw.file <- sub("\\..*", "", dfStore$variable)
dfStore$embryoNum <- sub(".*_", "", dfStore$variable)
dfStore$embryoID <- paste0(dfStore$Raw.file, "_", dfStore$embryoNum)

# how many proteins per cell are quantified? 
numProts <- df.m %>% dplyr::group_by(variable) %>% dplyr::summarise(numProts = sum(!is.na(value)))
numProts$bin <- floor(round(numProts$numProts, -2))
numProts$type <- dfStore$type[match(numProts$variable, dfStore$variable)]
ggplot(data = numProts, aes(x = type, y = numProts, fill = type)) + geom_boxplot() + theme_bw() +
  geom_beeswarm(size = 2, alpha =0.5)+
  scale_fill_brewer(palette = "OrRd") + 
  theme(axis.text = element_text(size = 22), axis.title =element_text(size = 22), 
        legend.position = "none") + 
  labs(x = "Stage", y = "Number of Proteins per cell") + 
  scale_x_discrete(labels = c("Early 2-cell", "Late 2-cell", "4-cell"))



###---------------------------------------------------------------------------------------------------------.
###---------------------------------- Cell-Cell Correlations & K-means-------------------------------------------------------
###---------------------------------------------------------------------------------------------------------.

colnames(df)[1] <- "Protein"
# missing data filter - let's observe each protein in x amount of cells
df.NA <- df[rowSums(is.na(df)) <= 0.7*ncol(df),]  
# or let's use the proteins found significantly differential in 2-cell stage 
sig2c_new <- read.csv("2-cell-stage-combined-sig-Prots_5percentFDR_2022_02_07.csv")
df.NA <-df %>% dplyr::filter(Protein %in% sig2c_new$Protein)

dist.mat <- cor(df.NA[,2:ncol(df.NA)], use = "pairwise.complete.obs", method = "spearman")

# change the color gradient 
library(circlize)
col_fun = colorRamp2(c(-0.5, 0, 0.5), c("yellow", "white", "purple"))
col_fun(seq(-3, 3))
Heatmap(dist.mat, col = col_fun) 



#assess optimal number of clusters 
library(factoextra)
factoextra::fviz_nbclust(dist.mat, kmeans, method = "silhouette")



###---------------------------------------------------------------------------------------------------------.
###----------------------------- Probability of k-means clustering assignment------------------------------------------------------
###---------------------------------------------------------------------------------------------------------.

colnames(df)[1] <- "Protein"
# missing data filter - let's observe each protein in x amount of cells
df.NA <- df[rowSums(is.na(df)) <= 0.7*ncol(df),]  
dist.mat <- cor(df.NA[,2:ncol(df.NA)], use = "pairwise.complete.obs", method = "spearman")


# kmeans clustering based on similarity metric
seedCHOSE <- 4
set.seed(seedCHOSE)
firTEST <- kmeans(dist.mat, centers = 2, iter.max = 25, nstart = 1) 
firTEST.cluster <- as.data.frame(firTEST$cluster)
colnames(firTEST.cluster) <- "Cluster"
clusplot(dist.mat, firTEST$cluster, main = paste0("Seed = ", seedCHOSE))
# count the number of cells per embryo in this clustering
CluKTEST <- as.data.frame(firTEST$cluster)


# let's implement the loop for 100 iterations
i<-2
CluKTEST <- data.frame()
for(i in 1:200) { 
  seedCHOSE <- i
  set.seed(seedCHOSE)
  firTEST <- kmeans(dist.mat, centers = 2, iter.max = 25, nstart = 1) 
  firTEST.cluster <- as.data.frame(firTEST$cluster)
  colnames(firTEST.cluster) <- "Cluster"
  clusplot(dist.mat, firTEST$cluster, main = paste0("Seed = ", seedCHOSE))
  # count the number of cells per embryo in this clustering
  tempCluKTEST <- as.data.frame(firTEST$cluster)
  if(i == 1) { 
    CluKTEST <- tempCluKTEST
    colnames(CluKTEST) <- paste0("Seed_", seedCHOSE)
  }
  if( i != 1) { 
    CluKTEST <- cbind(CluKTEST, tempCluKTEST)
    colnames(CluKTEST)[i] <- paste0("Seed_", seedCHOSE)
  }
}

corCluK <- cor(CluKTEST)
Heatmap(corCluK)
corCluK.m <- reshape2::melt(corCluK)
ggplot(data = corCluK.m, aes(x = value)) + 
  geom_histogram(aes(y = ..density..), color = "black", fill = "grey50", alpha = 0.3) + 
  geom_density(color = "darkblue", fill = "lightblue", alpha = 0.5) + 
  theme_bw() + 
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) + 
  labs(x = "Pearson Correlation of \nCluster Assignment Vectors")
ggplot(data = corCluK.m, aes(x = value)) + 
  #geom_histogram(color = "black", fill = "grey50", alpha = 0.3, binwidth = 0.05) + 
  geom_histogram(color = "black", fill = "red", alpha = 0.3, binwidth = 0.05) + 
  theme_bw() + 
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) + 
  labs(x = "Pearson Correlations of \nCluster Assignment Vectors")
corCluK.m$absCor <- abs(corCluK.m$value)
ggplot(data = corCluK.m, aes(x = absCor)) + 
  #geom_histogram(color = "black", fill = "grey50", alpha = 0.3, binwidth = 0.05) + 
  geom_histogram(color = "black", fill = "red", alpha = 0.3, binwidth = 0.05) + 
  theme_bw() + 
  theme(axis.text = element_text(size = 30), axis.title = element_text(size = 30)) + 
  labs(x = "Absolute Pearson Correlations of \nCluster Assignment Vectors")
unique(corCluK.m$value)

# calculate the probability of cluster assignment
# cluster assignment is random, so we need to correct for this 
CluKTEST.2 <- c()
i<-2
for(i in 1:ncol(CluKTEST)) { 
  if(i == 1) { CluKTEST.2 <- as.data.frame(CluKTEST[,i]); colnames(CluKTEST.2) <- colnames(CluKTEST)[i]}
  if(i != 2) { 
    temp.cor <- cor(CluKTEST.2[,1], CluKTEST[,i])
    if(temp.cor > 0) { 
      CluKTEST.2 <- cbind(CluKTEST.2,CluKTEST[,i])
      colnames(CluKTEST.2)[i] <- colnames(CluKTEST)[i]
    }
    if(temp.cor < 0) { 
      temp.df <- as.data.frame(CluKTEST[,i])
      temp.df$V2 <- ifelse(temp.df$`CluKTEST[, i]` ==1, 2, 1)
      CluKTEST.2 <- cbind(CluKTEST.2, temp.df$V2)
      colnames(CluKTEST.2)[i] <- colnames(CluKTEST)[i]
    }
  }
}
# now compute the mean across the rows 
CluKTEST.3 <- as.data.frame(rowMeans(CluKTEST.2))
CluKTEST.3$variable <- row.names(as.data.frame(firTEST$cluster))
CluKTEST.4 <- CluKTEST.3 %>% dplyr::filter(`rowMeans(CluKTEST.2)` !=1 )%>% 
  dplyr::filter(`rowMeans(CluKTEST.2)` !=2 )
ggplot(data = CluKTEST.3, aes(x = `rowMeans(CluKTEST.2)`)) + 
  geom_histogram(color = "black", fill = "red", alpha = 0.3, binwidth = 0.05) + 
  theme_bw() + 
  theme(axis.text = element_text(size = 30), axis.title = element_text(size = 30)) + 
  labs(x = "Cluster Assignment", y = "No. of Blastomeres")


i<-12
for(i in 1:nrow(CluKTEST.3)) { 
  if(CluKTEST.3[i,1] == 1) { CluKTEST.3$prob[i] <- 1 }
  if(CluKTEST.3[i,1] == 2) { CluKTEST.3$prob[i] <- 1 }
  if(CluKTEST.3[i,1] != 2 |CluKTEST.3[i,1] != 2) {
    temp.df <- as.data.frame(t(CluKTEST.2[i,]))
    temp.count <- as.data.frame(table(temp.df))
    temp.count$fraction <- temp.count$Freq/200
    CluKTEST.3$prob[i] <- max(temp.count$fraction)
  }
}
ggplot(data = CluKTEST.3, aes(x = prob)) + 
  geom_histogram(color = "black", fill = "red", alpha = 0.3, binwidth = 0.05) + 
  theme_bw() + 
  theme(axis.text = element_text(size = 30, color = "black"), axis.title = element_text(size = 30)) + 
  labs(x = "Probability", y = "No. of Blastomeres")



# count how many times an embryo is in each cluster 
CluK$embryo <- paste(sub("\\..*", "", rownames(CluK)), sub(".*_", "", rownames(CluK)) )
CluK$cleavage <- sub("([^_]*)_([^_]*)_([^_]*)", "\\2", rownames(CluK))
for(i in 1:nrow(CluK)) { 
  if(CluK$cleavage[i] == "v") { CluK$cleavage[i] <- "nid"}
  if(CluK$cleavage[i] == "emv") { CluK$cleavage[i] <- "em"}
  CluK$Raw.file[i] <- str_split(rownames(CluK)[i], fixed("."))[[1]][1]
}
colnames(CluK)[1] <- "cluster"
CluK$clusterName <- paste0("Cluster",CluK$cluster)
countCells <- CluK %>% 
  dplyr::group_by(Raw.file, embryo, cleavage) %>% 
  dplyr::summarise(n_clu1 = sum(cluster == "1"), n_clu2 = sum(cluster == "2"))



# kmeans clustering based on similarity metric 
set.seed(1)
fir <- kmeans(dist.mat, centers = 2, iter.max = 50, nstart = 25) 
fir.cluster <- as.data.frame(fir$cluster)
colnames(fir.cluster) <- "Cluster"
#write.csv(fir.cluster,"2022_01_11_fircluster.csv")
clusplot(dist.mat, fir$cluster)
# count the number of cells per embryo in this clustering
CluK <- as.data.frame(fir$cluster)
# count how many times an embryo is in each cluster 
CluK$embryo <- paste(sub("\\..*", "", rownames(CluK)), sub(".*_", "", rownames(CluK)) )
CluK$cleavage <- sub("([^_]*)_([^_]*)_([^_]*)", "\\2", rownames(CluK))
for(i in 1:nrow(CluK)) { 
  if(CluK$cleavage[i] == "v") { CluK$cleavage[i] <- "nid"}
  if(CluK$cleavage[i] == "emv") { CluK$cleavage[i] <- "em"}
  CluK$Raw.file[i] <- str_split(rownames(CluK)[i], fixed("."))[[1]][1]
}
colnames(CluK)[1] <- "cluster"
CluK$clusterName <- paste0("Cluster",CluK$cluster)
countCells <- CluK %>% 
  dplyr::group_by(Raw.file, embryo, cleavage) %>% 
  dplyr::summarise(n_clu1 = sum(cluster == "1"), n_clu2 = sum(cluster == "2"))


###---------------------------------------------------------------------------------------------------------.
###----------------------------- Significantly differential proteins------------------------------------------------------
###---------------------------------------------------------------------------------------------------------.

# find sig proteins at just the 2-cell stage
colnames(df)[1] <- "Protein"
df.m <- reshape2::melt(df)
df.m$cluster <- fir.cluster$Cluster[match(df.m$variable, rownames(fir.cluster))]
df.m$cluster <- paste0("Cluster",df.m$cluster)
df.m$variable_Cluster <- paste0(df.m$cluster, " ",df.m$variable)
df.m$embryoID <- dfStore$embryoID[match(df.m$variable, dfStore$variable)]
df.m$type <- dfStore$type[match(df.m$variable, dfStore$variable)]
# choose which stage you want to test proteins for? 
df.m.2c <- df.m # %>% dplyr::filter(grepl("_4c", type))
ProteinsAcrossClusters <- data.frame(matrix(nrow = length(unique(df.m.2c$Protein)), ncol = 2))
colnames(ProteinsAcrossClusters) <- c("Protein", "p.value")

for(i in 1:nrow(ProteinsAcrossClusters)) { 
  temp.prot <- df.m.2c %>% dplyr::filter(Protein == unique(df.m.2c$Protein)[i])
  temp.prot.1 <- temp.prot %>% dplyr::filter(cluster == "Cluster1")
  temp.prot.2 <- temp.prot %>% dplyr::filter(cluster == "Cluster2")
  if((sum(!is.na(temp.prot.1$value)) > 2) & (sum(!is.na(temp.prot.2$value)) > 2)  ) { 
    temp.test <- kruskal.test(value ~ cluster, temp.prot)
    ProteinsAcrossClusters$p.value[i] <- temp.test$p.value
  }
  ProteinsAcrossClusters$Protein[i] <- unique(df.m.2c$Protein)[i]
  print(i)
}

ProteinsAcrossClusters.NA <- ProteinsAcrossClusters %>% dplyr::filter(!is.na(p.value))
ProteinsAcrossClusters.NA$q.value <- p.adjust(ProteinsAcrossClusters.NA$p.value, method = "fdr")
ProteinsAcrossClusters.sig <- ProteinsAcrossClusters.NA %>% dplyr::filter(q.value <= 0.05)

# save this list 
df.m.2c.sig <- df.m.2c %>% dplyr::filter(Protein %in% ProteinsAcrossClusters.sig$Protein)
df.m.2c.sig.s <- reshape2::dcast(df.m.2c.sig, Protein ~ variable, value.var = "value")
#write.csv(df.m.2c.sig.s, "Late2-cell-stage-combined-sig-Prots_5percentFDR_2022_03_13.csv", row.names= F)
#write.csv(df.m.2c.sig.s, "Early-cell-stage-combined-sig-Prots_5percentFDR_2022_03_13.csv", row.names= F)
#write.csv(df.m.2c.sig.s, "2-cell-stage-combined-sig-Prots_5percentFDR_2022_02_07.csv", row.names= F)
#write.csv(df.m.2c.sig.s, "4-cell-stage-combined-sig-Prots_5percentFDR_2022_03_13.csv", row.names= F)
#write.csv(df.m.2c.sig.s, "ALL-cell-stage-combined-sig-Prots_5percentFDR_2022_03_13.csv", row.names= F)
#df.m.2c.sig.s$Gene <- convertProts$uppercase[match(df.m.2c.sig.s$Protein, convertProts$Entry)]

# compare the 4-cell stage proteins and 2-cell stage differntial proteins
sig2c_new <- read.csv("2-cell-stage-combined-sig-Prots_5percentFDR_2022_02_07.csv")
sig4c <- read.csv("4-cell-stage-combined-sig-Prots_5percentFDR_2022_03_13.csv")
sigALL <- read.csv("ALL-cell-stage-combined-sig-Prots_5percentFDR_2022_03_13.csv")

# intersect the lists: P61087
# proteins shared in 2c and 4c separate analysis
InFourTwo <- sig2c_new %>% dplyr::filter(Protein %in% sig4c$Protein)
# in 4c but not 2c
InFourNotTwo <- sig4c %>% dplyr::filter(!(Protein %in% sig2c_new$Protein))
# in 4c but not ALL
InFourNotALL <- sig4c %>% dplyr::filter(!(Protein %in% sigALL$Protein))

# in 4c but not 2c but in ALL
InFourInAllNotTwo <- sig4c %>% dplyr::filter(Protein %in% sigALL$Protein) %>% 
  dplyr::filter(!(Protein %in% sig2c_new$Protein))
# in 4c not in 2c or ALL
InFourNotTwoNotAll <- sig4c %>% dplyr::filter(!(Protein %in% sig2c$Protein)) %>% 
  dplyr::filter(!(Protein %in% sigALL$Protein))


# in 2c but not ALL
InTwoNotALL <- sig2c_new %>% dplyr::filter(!(Protein %in% sigALL$Protein))
InTwoInFourNotAll <- sig2c_new %>% dplyr::filter(!(Protein %in% sigALL$Protein)) %>% 
  dplyr::filter(Protein %in% sig4c$Protein)
InTwoNotFourInAll <- sig2c_new %>% dplyr::filter((Protein %in% sigALL$Protein)) %>% 
  dplyr::filter(!(Protein %in% sig4c$Protein))
InTwoNotFourNotAll <- sig2c_new %>% dplyr::filter(!(Protein %in% sigALL$Protein)) %>% 
  dplyr::filter(!(Protein %in% sig4c$Protein))

# make a venn diagram 
# 81 proteins shared between 4-cell and 2-cell 
# 359 in 2-cell
# 157 in 4-cell 
# 455 when all stages 
# 8 in 4-cell not all
INeverything <- Reduce(intersect, (list(sig2c_new$Protein, sig4c$Protein, sigALL$Protein)))

# check which proteins are in the same cluster, but increasing in magnitude across the stage  
# for the 455 proteins 
# set-up dataframe
medianValues1 <- df.m %>% dplyr::group_by(Protein,cluster,type) %>% dplyr::summarise(medianValue = median(value,na.rm =T))
medianValues1$cluster_type <- paste0(medianValues1$cluster, medianValues1$type)

hatter <- data.frame(matrix(nrow=nrow(sigALL), ncol = 4))
colnames(hatter) <- c("Proteins", "SameInCluster", "Magnitude", "QuantifiedIn")
hatter$Proteins <- sigALL$Protein

i<-2
for(i in 1:nrow(hatter)) { 
  temp.prot <- medianValues1 %>% dplyr::filter(Protein == hatter$Proteins[i])
  hatter$QuantifiedIn[i] <- sum(!is.na(temp.prot$medianValue))
  
  if(sum(!is.na(temp.prot$medianValue)) == 6) { 
    if((temp.prot$medianValue[temp.prot$cluster_type == "Cluster1early_2c"]) < 0 &
       (temp.prot$medianValue[temp.prot$cluster_type == "Cluster1late_2c"]) < 0  &
       (temp.prot$medianValue[temp.prot$cluster_type == "Cluster1_4c"]) < 0) { 
      
      hatter$SameInCluster[i] <- "yes"
    }
    if((temp.prot$medianValue[temp.prot$cluster_type == "Cluster1early_2c"]) > 0 &
       (temp.prot$medianValue[temp.prot$cluster_type == "Cluster1late_2c"]) > 0  &
       (temp.prot$medianValue[temp.prot$cluster_type == "Cluster1_4c"]) > 0) { 
      
      hatter$SameInCluster[i] <- "yes"
    }
    if((temp.prot$medianValue[temp.prot$cluster_type == "Cluster1early_2c"]) < 
       (temp.prot$medianValue[temp.prot$cluster_type == "Cluster1late_2c"]) &  
       (temp.prot$medianValue[temp.prot$cluster_type == "Cluster1late_2c"]) <
       (temp.prot$medianValue[temp.prot$cluster_type == "Cluster1_4c"]) ) { 
      
      hatter$Magnitude[i] <- "Up"
    }
    if((temp.prot$medianValue[temp.prot$cluster_type == "Cluster1early_2c"]) > 
       (temp.prot$medianValue[temp.prot$cluster_type == "Cluster1late_2c"]) &  
       (temp.prot$medianValue[temp.prot$cluster_type == "Cluster1late_2c"]) >
       (temp.prot$medianValue[temp.prot$cluster_type == "Cluster1_4c"]) ) { 
      
      hatter$Magnitude[i] <- "Down"
    }
  }
}

hatter.sameincluster1 <- hatter %>% dplyr::filter(SameInCluster == "yes")
hatter.sameincluster2 <- hatter.sameincluster %>% dplyr::filter(!is.na(Magnitude))
hatter.sameincluster3 <- hatter %>% dplyr::filter(QuantifiedIn == 6)
hatter.sameincluster4 <- hatter %>% dplyr::filter(SameInCluster == "yes") %>% dplyr::filter(!is.na(Magnitude))
write.csv( hatter.sameincluster2, "2022_04_20_hatter.csv")
# check the same for the 2-cell stage
hatter2c <- data.frame(matrix(nrow=nrow(sigALL), ncol = 4))
colnames(hatter2c) <- c("Proteins", "SameInCluster", "Magnitude", "QuantifiedIn")
hatter2c$Proteins <- sigALL$Protein

i<-2
for(i in 1:nrow(hatter2c)) { 
  temp.prot <- medianValues1 %>% dplyr::filter(Protein == hatter2c$Proteins[i]) %>% 
    dplyr::filter(type != "_4c")
  hatter2c$QuantifiedIn[i] <- sum(!is.na(temp.prot$medianValue))
  
  if(sum(!is.na(temp.prot$medianValue)) == 4) { 
    if((temp.prot$medianValue[temp.prot$cluster_type == "Cluster1early_2c"]) < 0 &
       (temp.prot$medianValue[temp.prot$cluster_type == "Cluster1late_2c"]) < 0 ) { 
      
      hatter2c$SameInCluster[i] <- "yes"
    }
    if((temp.prot$medianValue[temp.prot$cluster_type == "Cluster1early_2c"]) > 0 &
       (temp.prot$medianValue[temp.prot$cluster_type == "Cluster1late_2c"]) > 0) { 
      
      hatter2c$SameInCluster[i] <- "yes"
    }
    if((temp.prot$medianValue[temp.prot$cluster_type == "Cluster1early_2c"]) < 
       (temp.prot$medianValue[temp.prot$cluster_type == "Cluster1late_2c"])  ) { 
      
      hatter2c$Magnitude[i] <- "Up"
    }
    if((temp.prot$medianValue[temp.prot$cluster_type == "Cluster1early_2c"]) > 
       (temp.prot$medianValue[temp.prot$cluster_type == "Cluster1late_2c"])  ) { 
      
      hatter2c$Magnitude[i] <- "Down"
    }
  }
}

hatter2c.1 <- hatter2c %>% dplyr::filter(SameInCluster == "yes")
hatter2c.2 <- hatter2c.1 %>% dplyr::filter(!is.na(Magnitude))
hatter2c.3 <- hatter2c %>% dplyr::filter(QuantifiedIn == 4)

###---------------------------------------------------------------------------------------------------------.
###----------------------------- 4-cell stage: cleavage pattern ------------------------------------------------------
###---------------------------------------------------------------------------------------------------------.

# using the classification system, need to understand what proteins and processes are differential (if any)

# 01 - calculate the FCs 
mat5a <- df
rownames(mat5a) <- mat5a$Protein
mat5a <- mat5a[,-1]
mat5a.m <- reshape2::melt(as.matrix(mat5a))
# label which cluster each cell belongs to
mat5a.m$cluster <- fir.cluster$Cluster[match(mat5a.m$Var2, rownames(fir.cluster))]
mat5a.m$cluster <- paste0("Cluster",mat5a.m$cluster)
mat5a.m$variable_Cluster <- paste0(mat5a.m$cluster, " ",mat5a.m$Var2)
mat5a.m$embryoID <- dfStore$embryoID[match(mat5a.m$Var2, dfStore$variable)]
mat5a.m$type <- dfStore$type[match(mat5a.m$Var2, dfStore$variable)]
# rename the column names 
colnames(mat5a.m)[1:2] <- c("X", "variable")
# obtain the FC for each protein between sister cells via dplyr
meanValues <- mat5a.m %>% dplyr::group_by(X, embryoID, type) %>%
  dplyr::summarise(FC = value[cluster == "Cluster1"] - value[cluster == "Cluster2"])

# 02 - isolate just the 4 cell stage with defined cleavage pattern 
meanValues$cleavage <- dfStore$cleavage[match(meanValues$embryoID, dfStore$embryoID)]
meanValues.4c <- meanValues %>% dplyr::filter(grepl("em|me|ee", cleavage))
length(unique(meanValues.4c$embryoID))
meanValues.4c$Uniprot <- sub("-.*","",meanValues.4c$X)#meanValues.4c$X#
meanValues.4c$GeneName <- convertProts$uppercase[match(meanValues.4c$Uniprot, convertProts$Entry)]


# 03 - perform protein - level kruskal wallis tests 
cleavProt <- data.frame(matrix(nrow = length(unique(meanValues.4c$X)), ncol = 2))
colnames(cleavProt) <- c("Protein", "p.value")
cleavProt$Protein <- as.character(unique(meanValues.4c$X))
i<-1
for(i in 1:nrow(cleavProt)) { 
  #temp.prot <- meanValues.4c %>% dplyr::filter(X == cleavProt$Protein[i])
  temp.prot <- meanValues.4c[meanValues.4c$X == cleavProt$Protein[i],]
  temp.prot.1 <- temp.prot[temp.prot$cleavage == "ee",]
  temp.prot.2 <- temp.prot[temp.prot$cleavage == "em",]
  temp.prot.3 <- temp.prot[temp.prot$cleavage == "me",]
  if((sum(!is.na(temp.prot.1$FC)) > 2) & (sum(!is.na(temp.prot.2$FC)) > 2) & (sum(!is.na(temp.prot.3$FC)) > 2) ) { 
    #temp.test <- kruskal.test(FC ~ cleavage, temp.prot)
    #cleavProt$p.value[i] <- temp.test$p.value
    temp.test <- aov(FC ~ cleavage, data = temp.prot)
    cleavProt$p.value[i] <- summary(temp.test)[[1]][["Pr(>F)"]][1]
  }
  print(i)
}
# adjust the p values 
cleavProt.NA <- cleavProt %>% dplyr::filter(!is.na(p.value))
cleavProt.NA$q.value <- p.adjust(cleavProt.NA$p.value, method = "fdr")
cleavProt.sig <- cleavProt.NA %>% dplyr::filter(q.value < 0.05)
ggplot(cleavProt.NA, aes(x = -log10(q.value))) + geom_histogram(color = "black", fill = "slateblue") + theme_bw() + 
  geom_vline(xintercept = -log10(0.05), linetype = "dashed") + 
  theme(axis.title = element_text(size = 22), axis.text = element_text(size = 22, color = "black"), 
        plot.title = element_text(size =22)) + 
  labs(title = "Q values of Protein-Level Analaysis")
ggplot(cleavProt.NA, aes(y = -log10(q.value))) + geom_histogram(color = "black", fill = "slateblue") + theme_bw() + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") + 
  theme(axis.title = element_blank(), axis.text = element_text(size = 22, color = "black"), 
        plot.title = element_blank()) + 
  labs(title = "Q values of Protein-Level Analaysis")
# plot the protein distributions that are signficant 
meanValues.4c.sig <- meanValues.4c %>% dplyr::filter(X %in% cleavProt.sig$Protein)
ggplot(data = meanValues.4c.sig, aes(x = GeneName, y = FC, fill = cleavage)) + geom_boxplot() + 
  theme_bw() + theme(legend.position = "top", legend.text = element_text(size = 22),
                     legend.title = element_text(size = 22),
                     legend.key.size = unit(3, "line"),
                     axis.title = element_text(size = 22), 
                     axis.text.y = element_text(size = 22), 
                     axis.text.x = element_text(size = 20, color = "black"), plot.title = element_blank()) + 
  labs(x = "GeneName", y = "Alpha / Beta", title = "Cleavage Patterns, Protein-Level") + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  annotate("rect", xmin = 1.5, xmax = 2.5, 
           ymin = min(meanValues.4c.sig$FC, na.rm = T), ymax = max(meanValues.4c.sig$FC, na.rm = T), 
           alpha = 0.2, fill = "orange") + 
  annotate("rect", xmin = 3.5, xmax = 4.5, 
           ymin = min(meanValues.4c.sig$FC, na.rm = T), ymax = max(meanValues.4c.sig$FC, na.rm = T), 
           alpha = 0.2, fill = "orange")


# 04 - perform PSEA kruskal wallis tests 
GOinDat <- GOterms %>% dplyr::filter(GeneSymbol %in% meanValues.4c$GeneName)
GOCleav <- data.frame(matrix(nrow = length(unique(GOinDat$GOtermID)), ncol = 9))
colnames(GOCleav) <- c("GOterm", "Annotation", "p.value", "Proteins", "NumProteins", "NumProteinsInvolved", 
                       "Median_ee", "Median_em", "Median_me")
GOCleav$GOterm <- unique(GOinDat$GOtermID)
i<-5
#name the amount of proteins that are quantified, and number of proteins associated with GO term
#GOCleavFUNC<-function(GOCleav, GOinDat, GOterms, meanValues.4c){
for(i in 1:nrow(GOCleav)) { 
  temp.GO <- GOinDat %>% dplyr::filter(GOtermID == GOCleav$GOterm[i])
  temp.GO2 <- GOterms %>% dplyr::filter(GOtermID == GOCleav$GOterm[i]) ###
  temp.gene <- temp.GO$GeneSymbol
  #temp.prot <- meanValues.4c %>% dplyr::filter(GeneName %in% temp.gene)###
  temp.prot <- meanValues.4c[meanValues.4c$GeneName %in% temp.gene,]
  #temp.prot.1 <- temp.prot %>% dplyr::filter(cleavage == "ee")
  temp.prot.1 <- temp.prot[temp.prot$cleavage == "ee",]
  #temp.prot.2 <- temp.prot %>% dplyr::filter(cleavage == "em")
  temp.prot.2 <- temp.prot[temp.prot$cleavage == "em",]
  #temp.prot.3 <- temp.prot %>% dplyr::filter(cleavage == "me") 
  temp.prot.3 <- temp.prot[temp.prot$cleavage == "me",]
  if(sum(!is.na(temp.prot.1$FC)) > 2 & sum(!is.na(temp.prot.2$FC)) > 2 & sum(!is.na(temp.prot.3$FC)) > 2 ){
    #temp.test <- kruskal.test(FC ~ cleavage, data = temp.prot)
    #GOCleav$p.value[i] <- temp.test$p.value
    temp.test <- aov(FC ~ cleavage, data = temp.prot)
    GOCleav$p.value[i] <- summary(temp.test)[[1]][["Pr(>F)"]][1]
    GOCleav$Median_ee[i] <- median(temp.prot.1$FC, na.rm = T)
    GOCleav$Median_em[i] <- median(temp.prot.2$FC, na.rm = T)
    GOCleav$Median_me[i] <- median(temp.prot.3$FC, na.rm = T)
  } 
  
  GOCleav$Annotation[i] <- unique(temp.GO$annotation)
  GOCleav$Proteins[i] <- paste(unique(temp.prot$GeneName), collapse = " ")
  GOCleav$NumProteins[i] <- length(unique(temp.prot$GeneName))
  GOCleav$NumProteinsInvolved[i] <- length(unique(temp.GO2$GeneSymbol))
  print(i)
}
#return(GOCleav)
#} 
GOCleav$fraction <- GOCleav$NumProteins / GOCleav$NumProteinsInvolved
GOCleav.NA <- GOCleav %>% dplyr::filter(!is.na(p.value)) %>% dplyr::filter(fraction > 0.5)
GOCleav.NA$q.value <- p.adjust(GOCleav.NA$p.value, method = "fdr")
GOCleav.NA.sig <- GOCleav.NA %>% dplyr::filter(q.value <0.05)
ggplot(GOCleav.NA, aes(x = -log10(q.value))) + geom_histogram(color = "black", fill = "slateblue") + theme_bw() + 
  geom_vline(xintercept = -log10(0.05), linetype = "dashed") + 
  theme(axis.title = element_text(size = 22), axis.text = element_text(size = 22), 
        plot.title = element_text(size =22)) + 
  labs(title = "Q values of GO-Level Analaysis")
ggplot(GOCleav.NA, aes(y = -log10(q.value))) + geom_histogram(color = "black", fill = "slateblue") + theme_bw() + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") + 
  theme(axis.title = element_blank(), axis.text = element_text(size = 22, color = "black"), 
        plot.title = element_blank()) + 
  labs(title = "Q values of  Analaysis")
GOCleav.NA.sig.m <- reshape2::melt(GOCleav.NA.sig[,c(2,7:9)])
ggplot(data = GOCleav.NA.sig.m, aes(fill = Annotation, y = value, x = variable)) + geom_bar(stat = "identity") + 
  theme(legend.position = "bottom", axis.title = element_text(size = 22), axis.text = element_text(size = 22), 
        plot.title = element_text(size =22)) 

# collect the protein FCs for each significant GO term and plot distributions
protdistr <- c()
i<-1
for(i in 1:length(unique(GOCleav.NA.sig$GOterm))) { 
  temp.GO <- GOinDat %>% dplyr::filter(GOtermID == GOCleav.NA.sig$GOterm[i])
  temp.gene <- temp.GO$GeneSymbol
  temp.prot <- meanValues.4c[meanValues.4c$GeneName %in% temp.gene,]
  temp.prot$GOAnnot <- unique(temp.GO$annotation)
  if(i == 1) { protdistr <- temp.prot } 
  if(i != 1) { protdistr <- rbind(protdistr, temp.prot)}
}
ggplot(data = protdistr, aes(y = GOAnnot, x = FC, fill = cleavage)) + geom_boxplot() +  
  theme_bw()  + 
  theme(legend.position = "top", legend.text = element_text(size = 22),
        legend.title = element_text(size =22), 
        legend.key.size = unit(3, "line"),
        axis.title = element_blank(), axis.text.y = element_text(size = 22), 
        axis.text.x = element_text(size = 22)) + 
  labs(y = "Annotation", x = "Alpha / Beta")  + 
  coord_cartesian(xlim = c(-1, 1)) + geom_vline(xintercept = 0, linetype = "dashed")+
  annotate("rect", ymin = c(1.5, 3.5,5.5,7.5, 9.5, 11.5, 13.5), ymax = c(2.5,4.5,6.5,8.5,10.5, 12.5, 14.5), 
           xmin = min(meanValues.4c.sig$FC, na.rm = T), xmax = max(meanValues.4c.sig$FC, na.rm = T), 
           alpha = 0.2, fill = "orange")


###---------------------------------------------------------------------------------------------------------.
###----------------------------- PSEA between clusters ------------------------------------------------------
###---------------------------------------------------------------------------------------------------------.


# find GO terms that are differential between the clusters
# then plot boxplots of specific GO terms between each stage (early 2c, mid/late 2c, 4c)
# code below is taken from Secion 04a of this code 
mat4 <- df
rownames(mat4) <- mat4$Protein
mat4 <- mat4[,-1]

mat4.m <- reshape2::melt(as.matrix(mat4))
# label which cluster each cell belongs to
mat4.m$cluster <- fir.cluster$Cluster[match(mat4.m$Var2, rownames(fir.cluster))]
mat4.m$cluster <- paste0("Cluster",mat4.m$cluster)
# GSEA
mat4.m$Uniprot <- sub("-.*","",mat4.m$Var1)
mat4.m$GeneName <- convertProts$uppercase[match(mat4.m$Uniprot, convertProts$Entry)]
mat4.m$stage <- dfStore$type[match(mat4.m$Var2, dfStore$variable)]
mat4.m <- mat4.m #%>% dplyr::filter(stage == "_4c")
GOinDat <- GOterms %>% dplyr::filter(GeneSymbol %in% mat4.m$GeneName)
GOAcrossClusters <- data.frame(matrix(nrow = length(unique(GOinDat$GOtermID)), ncol = 8))
colnames(GOAcrossClusters) <- c("GOterm", "Annotation", "p.value", "Proteins", "NumProteins", "NumProteinsInvolved", 
                                "median_Cluster1", "median_Cluster2")
GOAcrossClusters$GOterm <- unique(GOinDat$GOtermID)
i<-449
#name the amount of proteins that are quantified, and number of proteins associated with GO term
for(i in 1:nrow(GOAcrossClusters)) { 
  temp.GO <- GOinDat %>% dplyr::filter(GOtermID == GOAcrossClusters$GOterm[i])
  temp.GO2 <- GOterms %>% dplyr::filter(GOtermID == GOAcrossClusters$GOterm[i])
  temp.prot <- mat4.m %>% dplyr::filter(GeneName %in% temp.GO$GeneSymbol)
  temp1 <- temp.prot %>% dplyr::filter( cluster == "Cluster1") 
  temp2 <- temp.prot %>% dplyr::filter( cluster == "Cluster2") 
  if(sum(!is.na(temp1$value)) >3 & sum(!is.na(temp2$value)) >3 ){
    #& length(unique(temp.prot$GeneName)) >= 4 & length(unique(temp.GO2$GeneSymbol)) <= 200) { #
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
    GOAcrossClusters[i,9:(8+ncol(tempmed))] <- tempmed
  } 
  
  GOAcrossClusters$GOterm[i] <- unique(temp.GO$GOtermID)
  GOAcrossClusters$Annotation[i] <- unique(temp.GO$annotation)
  GOAcrossClusters$Proteins[i] <- paste(unique(temp.prot$GeneName), collapse = " ")
  GOAcrossClusters$NumProteins[i] <- length(unique(temp.prot$GeneName))
  GOAcrossClusters$NumProteinsInvolved[i] <- length(unique(temp.GO2$GeneSymbol))
  print(i)
}

#write.csv(GOAcrossClusters, "2021_11_11_GOAcrossCluster_allstages.csv", row.names = FALSE)
#write.csv(GOAcrossClusters, "2021_11_11_GOAcrossCluster_early2c.csv", row.names = FALSE)
#write.csv(GOAcrossClusters, "2021_11_11_GOAcrossCluster_late2c.csv", row.names = FALSE)
#write.csv(GOAcrossClusters, "2021_11_11_GOAcrossCluster_4c.csv", row.names = FALSE)



#................ PLOTTING OPTIONS ARE BELOW ....................................................#
# read in file
#GOAcrossClusters2<- read.csv("GOAcrossClusters_2021_09_29_B_3orMoreObsperGroup_ExtraProtFils.csv")
#GOAcrossClusters <- GOAcrossClusters[,-1]
GOAcrossClusters <- read.csv("2021_11_11_GOAcrossCluster_allstages.csv")
q.value.filter <- 0.05
GOAcrossClusters.NA <- GOAcrossClusters
GOAcrossClusters.NA <- GOAcrossClusters.NA %>% dplyr::filter(!is.na(p.value))
GOAcrossClusters.NA$q.value <- p.adjust(GOAcrossClusters.NA$p.value, method = "fdr")
GOAcrossClusters.NA$fraction <- GOAcrossClusters.NA$NumProteins / GOAcrossClusters.NA$NumProteinsInvolved
GOAcrossClusters.NA$Clus1_vs_Clus2 <- GOAcrossClusters.NA$median_Cluster1 - GOAcrossClusters.NA$median_Cluster2
GOAcrossClusters.NA$group <- ifelse(GOAcrossClusters.NA$median_Cluster1 > GOAcrossClusters.NA$median_Cluster2, 
                                    "Greater in Clus1", "Greater in Clus2")

#GOAcrossClusters.NA <- GOAcrossClusters
#GOAcrossClusters.NA <- GOAcrossClusters.NA %>% dplyr::filter(!is.na(p.value))
# save the GOAcrossClusters
#write.csv(GOAcrossClusters, "GOAcrossClusters_2021_09_28_A_3orMoreObsperGroup.csv") # if(sum(!is.na(temp1$value)) >3 & sum(!is.na(temp2$value)) >3 )
#write.csv(GOAcrossClusters, "GOAcrossClusters_2021_09_29_B_3orMoreObsperGroup_ExtraProtFils.csv")


# take specific GO term themes and try to see patterns that emerge
GOAcrossClusters.sig <- GOAcrossClusters.NA %>% dplyr::filter(q.value <= 0.05)   %>% dplyr::filter(fraction >= 0.3) 
# filter simply on q-value
GOAcrossClusters.sig2 <- GOAcrossClusters.NA %>% dplyr::filter(q.value <= 0.05) 
GOAcrossClusters.sig2.prot <- GOAcrossClusters.sig2 %>% dplyr::filter(NumProteins > 2) 
# look at effect sizes
GOAcrossClusters.sig3 <- GOAcrossClusters.sig[order(GOAcrossClusters.sig$Clus1_vs_Clus2),]
GOAcrossClusters.sig3 <- GOAcrossClusters.NA %>% dplyr::filter(NumProteins > 5)
GOAcrossClusters.sig3$q.value <- p.adjust(GOAcrossClusters.sig3$p.value, method = "fdr")
GOAcrossClusters.sig4 <- GOAcrossClusters.sig3 %>% dplyr::filter(abs(Clus1_vs_Clus2) > 0.2)


bindDF <- GOAcrossClusters.sig2
for(i in 1:nrow(bindDF)) { 
  if(grepl("proteas|ubi|apop|peptidase|autopha|neddyl",bindDF$Annotation[i]) ==TRUE ) { 
    bindDF$type[i] <- "degradation" }
  else if(grepl("transport|vesicle|channel|SNARE|myosin|dyne",bindDF$Annotation[i]) ==TRUE) { 
    bindDF$type[i] <- "transport" }
  else if(grepl("translation|ribosom|endoplasmic reticul",bindDF$Annotation[i]) ==TRUE ) { 
    bindDF$type[i] <- "translation" }
  else if(grepl("splic|snRN",bindDF$Annotation[i]) ==TRUE ) { 
    bindDF$type[i] <- "mRNA splicing" }
  else if(grepl("chromatin|chromosom",bindDF$Annotation[i]) ==TRUE ) { 
    bindDF$type[i] <- "chromatin, chromosome" }
  else if(grepl("kinase|phosph|calcium|calmodulin|migrat|secret|receptor",bindDF$Annotation[i]) ==TRUE  ) { 
    bindDF$type[i] <- "signaling" }
  else if(grepl("mitocho|metabo|electron transport|ATP|catab|glutam|oxid|hydrogenase",bindDF$Annotation[i]) ==TRUE) { 
    bindDF$type[i] <- "mitochondrial|metabolism" }
  else if(grepl("cell cycl|kinetochore|spindle|mitos|meios",bindDF$Annotation[i]) ==TRUE) { 
    bindDF$type[i] <- "cell cycle" }
  else{ bindDF$type[i] <- "No theme"}
  print(i)
}
# plot density plots of themes 
bindDF.plot <- bindDF %>% dplyr::filter(type != "No theme") %>% 
  dplyr::filter(type == "mitochondrial|metabolism" | type == "transport" | type == "translation" | type == "degradation" | type == "mRNA splicing")
ggplot(data = bindDF.plot, aes(x = Clus1_vs_Clus2, fill = type)) + geom_density(alpha = 0.5) + 
  theme_bw() + theme(axis.title = element_text(size = 14), plot.title = element_text(size = 18), 
                     legend.position = "none",#legend.position = c(0.85, 0.1), 
                     legend.text = element_text(size = 14), 
                     legend.title = element_text(size = 14), 
                     strip.text = element_text(size = 16, face = "bold", color = "white"), 
                     strip.background = element_rect(fill = "black")) + 
  scale_fill_discrete(name = "Biological Process Theme") + 
  labs(title = "Median Relative Abundance of GO terms per Cluster", x = "Cluster 1 - Cluster 2") + 
  facet_wrap(~type) + geom_vline(xintercept = 0)


# let's select for a theme and break it down
# for example, degradation- break it down with apoptosis, proteasome, ubiquitin ligases 
bindDF.plot <- bindDF %>% dplyr::filter(type == "degradation") %>% 
  dplyr::mutate(granularType = case_when(#grepl("apop", Annotation) == TRUE ~ "Apoptosis", 
    grepl("proteasom", Annotation) == TRUE & 
      grepl("ubi", Annotation) == FALSE~ "Proteasome", 
    grepl("ubiqui|neddyl", Annotation) == TRUE~ "Ubiquitin", 
    grepl("autop", Annotation) == TRUE ~ "Autophagy")) %>% 
  dplyr::filter(granularType != "NA")
# Version 1 plot 
ggplot(data = bindDF.plot, aes(x = Clus1_vs_Clus2, fill = granularType)) + geom_density(alpha = 0.5) + 
  theme_bw() + theme(axis.title = element_text(size = 24), 
                     plot.title = element_text(size = 24, face = "bold"), 
                     axis.text = element_text(size = 18), 
                     #legend.position = "none",
                     legend.position = c(0.27, 0.77), 
                     legend.background = element_rect(fill = "transparent", color = "transparent"), 
                     legend.text = element_text(size = 18), 
                     legend.title = element_text(size = 24), 
                     strip.text = element_text(size = 16, face = "bold", color = "white"), 
                     strip.background = element_rect(fill = "black")) + 
  coord_cartesian(xlim = c(-0.3, 0.3)) + 
  #facet_wrap(~granularType) + geom_vline(xintercept = 0) + 
  scale_fill_discrete(name = "Protein Degradation \nGO term umbrellas") + 
  #labs(title = "Protein Degradation GO term umbrellas") +
  xlab(bquote("Median"~log[2]~"FC (Alpha - Beta)"))
# Version 2 plot
ggplot(data = bindDF.plot, aes(x = Clus1_vs_Clus2, fill = granularType)) + geom_density(alpha = 0.5) + 
  theme_bw() + theme(axis.title = element_text(size = 24), 
                     plot.title = element_text(size = 24, face = "bold"), 
                     axis.text = element_text(size = 18), 
                     legend.position = c(0.22, 0.85), 
                     #legend.background = element_rect(fill = "transparent", color = "transparent"), 
                     legend.background = element_rect(fill = "transparent"),
                     legend.text = element_text(size = 26, family = "Courier"), 
                     legend.title = element_blank()) + 
  #coord_cartesian(xlim = c(-0.3, 0.3)) + 
  xlab(bquote("Median"~log[2]~"FC (Alpha - Beta)"))
test <- bindDF.plot %>% dplyr::filter(granularType == "Autophagy")


# now let's look at protein transport 
# transport|vesicle|channel|SNARE|myosin|dyne
bindDF.plot <- bindDF %>% dplyr::filter(type == "transport") %>% 
  dplyr::mutate(granularType = case_when(#grepl("voltage", Annotation) == TRUE  &
    grepl("voltage|proton|potassium|calcium|sodium", Annotation) == TRUE  & 
      #!(grepl("negative", Annotation) == TRUE)~ "Ion", 
      !(grepl("negative", Annotation) == TRUE)~ "Channel & Signaling", 
    grepl("vesicle|Golgi|SNARE", Annotation) == TRUE ~ "Vesicle", 
    grepl("dyne|actin", Annotation) == TRUE~ "Molecular Motors" 
    #grepl("ATP", Annotation) == TRUE ~ "ATP"
  )) %>% 
  dplyr::filter(granularType != "NA")
#test <- bindDF.plot %>% dplyr::filter(granularType == "Voltage Channels")
ggplot(data = bindDF.plot, aes(x = Clus1_vs_Clus2, fill = granularType)) + geom_density(alpha = 0.5) + 
  theme_bw() + theme(axis.title = element_text(size = 24), plot.title = element_text(size = 24, face = "bold"), 
                     #legend.position = "none",
                     legend.position = c(0.7, 0.75), 
                     axis.text = element_text(size = 18), 
                     legend.text = element_text(size = 18), 
                     legend.title = element_text(size = 24), 
                     strip.text = element_text(size = 16, face = "bold", color = "white"), 
                     strip.background = element_rect(fill = "indianred")) + 
  coord_cartesian(xlim = c(-0.4, 0.4)) + 
  #labs(title = "Protein Transport GO term umbrellas") + 
  scale_fill_brewer(name = "Protein Transport \nGO term umbrellas", palette = "Set1") + 
  #facet_wrap(~granularType) + 
  #geom_vline(xintercept = 0) + 
  xlab(bquote("Median"~log[2]~"FC (Alpha - Beta)"))
# Version 2 plot 
ggplot(data = bindDF.plot, aes(x = Clus1_vs_Clus2, fill = granularType)) + geom_density(alpha = 0.5) + 
  theme_bw() + theme(axis.title = element_text(size = 24), plot.title = element_text(size = 24, face = "bold"), 
                     legend.position = c(0.73, 0.85), 
                     #legend.position = c(0.7, 0.75), 
                     axis.text = element_text(size = 18), 
                     legend.text = element_text(size = 26), 
                     legend.background = element_rect(fill = "transparent"),
                     legend.title = element_blank()) + 
  coord_cartesian(xlim = c(-0.4, 0.4)) + 
  scale_fill_brewer( palette = "Set1") + 
  xlab(bquote("Median"~log[2]~"FC (Alpha - Beta)"))


test <- bindDF.plot %>% dplyr::filter(granularType == "Ion") %>% dplyr::filter(group == "Greater in Clus2")
test2 <- bindDF.plot %>% dplyr::filter(granularType == "Ion") %>% dplyr::filter(group == "Greater in Clus1")




# now let's look at translation
# translation|ribosom|endoplasmic reticul
bindDF.plot <- bindDF %>% dplyr::filter(type == "translation") %>% 
  dplyr::mutate(granularType = case_when(grepl("initiation", Annotation) == TRUE ~ "Initiation Factors", 
                                         grepl("ribosom", Annotation) == TRUE ~ "Ribosome",
                                         grepl("endoplasmic reticul", Annotation) == TRUE~ "Endoplasmic Reticulum"
  )) %>% 
  dplyr::filter(granularType != "NA")
ggplot(data = bindDF.plot, aes(x = Clus1_vs_Clus2, fill = granularType)) + geom_density(alpha = 0.5) + 
  theme_bw() + theme(axis.title = element_text(size = 20), 
                     #legend.position = "none",
                     legend.position = c(0.23, 0.8), 
                     legend.text = element_text(size = 14), 
                     legend.title = element_text(size = 14), 
                     strip.text = element_text(size = 16, face = "bold", color = "white"), 
                     strip.background = element_rect(fill = "black"), 
                     axis.text = element_text(size = 18), 
                     plot.title = element_text(size = 25)) + 
  #coord_cartesian(xlim = c(-0.2, 0.2)) + 
  labs(title = "Protein Translation Biological Processes", x = "Alpha / Beta") + 
  #facet_wrap(~granularType) + 
  geom_vline(xintercept = 0) + 
  scale_fill_discrete(name = "Themes")
# Version 2 plot 
ggplot(data = bindDF.plot, aes(x = Clus1_vs_Clus2, fill = granularType)) + geom_density(alpha = 0.6) + 
  theme_bw() + theme(axis.title = element_blank(), 
                     legend.position = c(0.25, 0.8), 
                     axis.text = element_text(size = 18, color = "black"), 
                     legend.text = element_text(size = 24), 
                     legend.background = element_rect(fill = "transparent"),
                     legend.title = element_blank()) + 
  #coord_cartesian(xlim = c(-0.4, 0.4)) + 
  scale_fill_brewer(palette = "Spectral",
                    labels = c("Endo. Reticulum", "Initation Factors", "Ribosome")) + 
  xlab(bquote("Median"~log[2]~"FC (Alpha - Beta)"))





###---------------------------------------------------------------------------------------------------------.
###--------------------------------------- Sister Euclidean Distance -------------------------------------------------------
###---------------------------------------------------------------------------------------------------------.

# take the euclidean distance between each sister cell within each embryo 
sigALL <- read.csv("ALL-cell-stage-combined-sig-Prots_5percentFDR_2022_03_13.csv")

df.sister <- df
df.sister[df.sister == 0] <- NA
# use proteins quantified in every cell: #na.omit(df.sister)
r <- na.omit(df.sister)#df.sister %>% dplyr::filter(Protein %in% sigALL$Protein)
df.m <- reshape2::melt(r)
df.m$cluster <- fir.cluster$Cluster[match(df.m$variable, rownames(fir.cluster))]
df.m$cluster <- paste0("Cluster",df.m$cluster)
df.m$variable_Cluster <- paste0(df.m$cluster, " ",df.m$variable)
df.m$embryoID <- dfStore$embryoID[match(df.m$variable, dfStore$variable)]
df.m$type <- dfStore$type[match(df.m$variable, dfStore$variable)]
countembryos <- df.m %>% dplyr::group_by(type) %>% dplyr::summarise(n_embryo = length(unique(embryoID)))

DistanceEmbryo <- as.data.frame(matrix(ncol = 4, nrow = ((21*6)+((15+21))))) #based on embryo comparisons
colnames(DistanceEmbryo) <- c("Embryo", "type", "Comparison", "Distance")

ProteinEmbryo <- as.data.frame(matrix(ncol = 4, nrow = ((21*6)+((15+21)))))
colnames(ProteinEmbryo) <- c("Embryo", "type", "Comparison", "NumProteins")

i<-1
rowCounter <- 1
for(i in 1:length(unique(df.m$embryoID))) { 
  temp.embryo <- df.m %>% dplyr::filter(embryoID == unique(df.m$embryoID)[i])
  temp.embryo.d <- reshape2::dcast(temp.embryo, variable ~ Protein, value.var = "value")
  temp.dis <- dist(temp.embryo.d[,2:ncol(temp.embryo.d)], method = "euclidean")
  
  # count the number of shared observations between protein (i.e how many single cells have both proteins quantified)
  obser <- pairwiseCount(t(temp.embryo.d[,2:ncol(temp.embryo.d)]))
  colnames(obser) <- temp.embryo.d$variable
  rownames(obser) <- temp.embryo.d$variable
  # obtain the upper triangle of shared observation matrix 
  obser[lower.tri(obser, diag = T)] <- 188 
  obser.m <- reshape2::melt(obser)
  obser.m <- obser.m %>% dplyr::filter(value!=188)
  obser.m$Var1Var2 <- paste(obser.m$Var1, ":", obser.m$Var2)
  if(grepl("_2c", unique(temp.embryo$type))) {
    #temp.embryo.d <- reshape2::dcast(temp.embryo, variable ~ Protein, value.var = "value")
    #temp.dis <- dist(temp.embryo.d[,2:ncol(temp.embryo.d)], method = "euclidean")
    DistanceEmbryo$Distance[rowCounter] <- temp.dis
    DistanceEmbryo$Embryo[rowCounter] <- unique(temp.embryo$embryoID)
    DistanceEmbryo$type[rowCounter] <- unique(temp.embryo$type)
    DistanceEmbryo$Comparison[rowCounter] <- paste(unique(temp.embryo$variable)[1], ":", unique(temp.embryo$variable)[2])
    ProteinEmbryo$Comparison[rowCounter] <- paste(unique(temp.embryo.d$variable)[1], ":", unique(temp.embryo$variable)[2])
    ProteinEmbryo$NumProteins[rowCounter] <- obser.m$value[obser.m$Var1Var2 == ProteinEmbryo$Comparison[rowCounter]]
    #ProteinEmbryo$NumProteins[rowCounter] <- rowSums(!is.na(temp.embryo.d[,2:ncol(temp.embryo.d)]))[1]
    rowCounter <- rowCounter + 1 
  }
  else if(unique(temp.embryo$type) == "_4c") {
    
    
    for(j in 1:6) { 
      DistanceEmbryo$Distance[rowCounter] <- temp.dis[j]
      DistanceEmbryo$Embryo[rowCounter] <- unique(temp.embryo$embryoID)
      DistanceEmbryo$type[rowCounter] <- unique(temp.embryo$type)
      
      if(j == 1) { 
        DistanceEmbryo$Comparison[rowCounter] <- paste(unique(temp.embryo.d$variable)[1], ":", unique(temp.embryo$variable)[2])
        ProteinEmbryo$Comparison[rowCounter] <- paste(unique(temp.embryo.d$variable)[1], ":", unique(temp.embryo$variable)[2])
        ProteinEmbryo$NumProteins[rowCounter] <- obser.m$value[obser.m$Var1Var2 == ProteinEmbryo$Comparison[rowCounter]]
        rowCounter <- rowCounter + 1 
      }
      if(j == 2) { 
        DistanceEmbryo$Comparison[rowCounter] <- paste(unique(temp.embryo.d$variable)[1], ":", unique(temp.embryo$variable)[3])
        ProteinEmbryo$Comparison[rowCounter] <- paste(unique(temp.embryo.d$variable)[1], ":", unique(temp.embryo$variable)[3])
        ProteinEmbryo$NumProteins[rowCounter] <- obser.m$value[obser.m$Var1Var2 == ProteinEmbryo$Comparison[rowCounter]]
        rowCounter <- rowCounter + 1 
      }
      if(j == 3) { 
        DistanceEmbryo$Comparison[rowCounter] <- paste(unique(temp.embryo.d$variable)[1], ":", unique(temp.embryo$variable)[4])
        ProteinEmbryo$Comparison[rowCounter] <- paste(unique(temp.embryo.d$variable)[1], ":", unique(temp.embryo$variable)[4])
        ProteinEmbryo$NumProteins[rowCounter] <- obser.m$value[obser.m$Var1Var2 == ProteinEmbryo$Comparison[rowCounter]]
        rowCounter <- rowCounter + 1 
      }
      if(j == 4) { 
        DistanceEmbryo$Comparison[rowCounter] <- paste(unique(temp.embryo.d$variable)[2], ":", unique(temp.embryo$variable)[3])
        ProteinEmbryo$Comparison[rowCounter] <- paste(unique(temp.embryo.d$variable)[2], ":", unique(temp.embryo$variable)[3])
        ProteinEmbryo$NumProteins[rowCounter] <- obser.m$value[obser.m$Var1Var2 == ProteinEmbryo$Comparison[rowCounter]]
        rowCounter <- rowCounter + 1 
      }
      if(j == 5) { 
        DistanceEmbryo$Comparison[rowCounter] <- paste(unique(temp.embryo.d$variable)[2], ":", unique(temp.embryo$variable)[4])
        ProteinEmbryo$Comparison[rowCounter] <- paste(unique(temp.embryo.d$variable)[2], ":", unique(temp.embryo$variable)[4])
        ProteinEmbryo$NumProteins[rowCounter] <- obser.m$value[obser.m$Var1Var2 == ProteinEmbryo$Comparison[rowCounter]]
        rowCounter <- rowCounter + 1 
      }
      if(j == 6) { 
        DistanceEmbryo$Comparison[rowCounter] <- paste(unique(temp.embryo.d$variable)[3], ":", unique(temp.embryo$variable)[4])
        ProteinEmbryo$Comparison[rowCounter] <- paste(unique(temp.embryo.d$variable)[3], ":", unique(temp.embryo$variable)[4])
        ProteinEmbryo$NumProteins[rowCounter] <- obser.m$value[obser.m$Var1Var2 == ProteinEmbryo$Comparison[rowCounter]]
        rowCounter <- rowCounter + 1 
      }
      
      
      
    }
    
  }
}
DistanceEmbryo$type <- factor(DistanceEmbryo$type, levels = c("early_2c", "late_2c", "_4c"))
ggplot(data = DistanceEmbryo, aes(y = Distance, x =type, fill = type)) + geom_boxplot()+ theme_bw() +
  scale_fill_brewer(palette = "OrRd") + 
  #coord_cartesian(ylim = c(3,8)) + 
  theme(axis.text.y  = element_text(size = 25, color = "black"), axis.title =element_blank(), 
        axis.text.x = element_text(size =25, angle = 45, hjust = 1, color = "black"), 
        legend.position = "none") + 
  labs(x = "Stage", y = "Sister Blastomere Distances") + 
  scale_x_discrete(labels = c("Early 2-cell", "Late 2-cell", "4-cell"))
DistanceEmbryo$NumProteins <- ProteinEmbryo$NumProteins
DistanceEmbryo$NormDist <- DistanceEmbryo$Distance / DistanceEmbryo$NumProteins
ggplot(data = DistanceEmbryo, aes(y = log10(NormDist), x =type, fill = type)) + geom_boxplot()+ theme_bw() +
  geom_beeswarm(size = 2, alpha =0.5)+
  scale_fill_brewer(palette = "OrRd") + 
  #coord_cartesian(ylim = c(3,8)) + 
  theme(axis.text = element_text(size = 22), axis.title =element_text(size = 22), 
        legend.position = "none") + 
  labs(x = "Stage", y = "Sister Blastomere Distances") + 
  scale_x_discrete(labels = c("Early 2-cell", "Late 2-cell", "4-cell"))

DistanceEmbryo.1 <- DistanceEmbryo %>% dplyr::filter(type == "early_2c")
DistanceEmbryo.2 <- DistanceEmbryo %>% dplyr::filter(type == "late_2c")
t.test(log10(DistanceEmbryo.1$NormDist), log10(DistanceEmbryo.2$NormDist))
t <- aov(log10(DistanceEmbryo$NormDist) ~ type, data = DistanceEmbryo)
summary(t)[[1]][["Pr(>F)"]][1]








###---------------------------------------------------------------------------------------------------------.
###------------------------------------Protein Fold-change between sisters -------------------------------------------------------
###---------------------------------------------------------------------------------------------------------.

mat5a <- df
rownames(mat5a) <- mat5a$Protein
mat5a <- mat5a[,-1]
mat5a.m <- reshape2::melt(as.matrix(mat5a))
# label which cluster each cell belongs to
mat5a.m$cluster <- fir.cluster$Cluster[match(mat5a.m$Var2, rownames(fir.cluster))]
mat5a.m$cluster <- paste0("Cluster",mat5a.m$cluster)
mat5a.m$variable_Cluster <- paste0(mat5a.m$cluster, " ",mat5a.m$Var2)
mat5a.m$embryoID <- dfStore$embryoID[match(mat5a.m$Var2, dfStore$variable)]
mat5a.m$type <- dfStore$type[match(mat5a.m$Var2, dfStore$variable)]

colnames(mat5a.m)[1:2] <- c("X", "variable")
mat5a.m$X <- as.character(mat5a.m$X)
mat5a.m$variable <- as.character(mat5a.m$variable)
mat5a.m.d <- reshape2::dcast(X~ variable_Cluster, data = mat5a.m, value.var = "value")
# obtain the FC via dplyr
meanValues <- mat5a.m %>% dplyr::group_by(X, embryoID, type) %>%
  dplyr::summarise(FC = value[cluster == "Cluster1"] - value[cluster == "Cluster2"])
ggplot(data = meanValues, aes(x = embryoID, y = FC)) + geom_boxplot()
ggplot(data = meanValues, aes(x = FC)) + geom_histogram(binwidth = 0.2) + coord_cartesian(xlim = c(-2, 2))
ggplot(data = mat5a.m, aes(x = value)) + geom_histogram(binwidth = 0.2) + coord_cartesian(xlim = c(-2, 2))

meanValues.sel <- meanValues %>% dplyr::filter(X %in% sig2c_new$Protein) 
# get the mean fold-change of each protein in each stage 
meanValues.sel.1 <- meanValues.sel %>% dplyr::group_by(type, X) %>% dplyr::summarise(meanFC = mean(FC, na.rm =T))
meanValues.sel.1$type <-factor(meanValues.sel.1$type, levels  = c("early_2c", "late_2c", "_4c"))
ggplot(meanValues.sel.1, aes(x = type, y = meanFC)) + geom_boxplot() + 
  coord_cartesian(ylim = c(-1,1))
ggplot(data = meanValues.sel.1, aes(x = type, y = X, fill = meanFC)) + geom_tile()

#write.csv(meanValues, "2021_12_07_SisterFoldChanges.csv", row.names = FALSE)
meanValues.sel.d <- reshape2::dcast(meanValues.sel, X ~ embryoID, value.var = "FC")
#write.csv(meanValues.sel, "2022_04_20_SisterFC_mouse_sigProteins.csv")


###---------------------------------------------------------------------------------------------------------.
###-------------------------------GO terms that are changing across the stages -------------------------------------------------------
###---------------------------------------------------------------------------------------------------------.
# go through each GO term
# obtain the median level of FCs per stage
# if the level is either increasing or decreasing, test the correlation and obtain pvalue 
FC <- meanValues
FC.1 <- as.data.frame(FC)
FC.2 <- FC.1 %>% dplyr::mutate(timepoint = case_when(type == "early_2c" ~ 0, 
                                                     type == "late_2c" ~1, 
                                                     type == "_4c" ~2))
FC.2$Uniprot <- FC.2$X#sub("-.*","",FC.2$X)
FC.2$GeneName <- convertProts$uppercase[match(FC.2$Uniprot, convertProts$Entry)]

GOinDat <- GOterms %>% dplyr::filter(GeneSymbol %in% FC.2$GeneName)
GOAcrossStages <- data.frame(matrix(nrow = length(unique(GOinDat$GOtermID)), ncol = 11))
colnames(GOAcrossStages) <- c("GOterm", "Annotation", "cor", "p.value", "Proteins", "NumProteins", "NumProteinsInvolved", 
                              "fraction","median_early2c", "median_late2c", "median_4c")
GOAcrossStages$GOterm <- unique(GOinDat$GOtermID)
i<-2465
#name the amount of proteins that are quantified, and number of proteins associated with GO term
for(i in 1:nrow(GOAcrossStages)) { 
  temp.GO <- GOinDat %>% dplyr::filter(GOtermID == GOAcrossStages$GOterm[i])
  temp.GO2 <- GOterms %>% dplyr::filter(GOtermID == GOAcrossStages$GOterm[i])
  temp.gene <- temp.GO$GeneSymbol
  temp.prot <- FC.2 %>% dplyr::filter(GeneName %in% temp.gene)
  temp1 <- temp.prot %>% dplyr::filter( type == "early_2c") 
  temp2 <- temp.prot %>% dplyr::filter( type == "late_2c") 
  temp3 <- temp.prot %>% dplyr::filter( type == "_4c") 
  if(sum(!is.na(temp1$FC)) >5 & sum(!is.na(temp2$FC)) >5 & sum(!is.na(temp3$FC)) >5  #&
     #median(temp1$FC, na.rm = T) < median(temp2$FC, na.rm = T) & 
     #median(temp2$FC, na.rm = T) <  median(temp3$FC, na.rm = T)
  ){
    #& length(unique(temp.prot$GeneName)) >= 4 & length(unique(temp.GO2$GeneSymbol)) <= 200) { #
    temp.cor2 <- cor.test(temp.prot$FC, temp.prot$timepoint, use = "pairwise.complete.obs", method = "spearman")
    GOAcrossStages$p.value[i] <- temp.cor2$p.value
    GOAcrossStages$cor[i] <- temp.cor2$estimate
    GOAcrossStages$median_early2c[i] <- median(temp1$FC, na.rm = T)
    GOAcrossStages$median_late2c[i] <- median(temp2$FC, na.rm = T)
    GOAcrossStages$median_4c[i] <- median(temp3$FC, na.rm = T)
    temp1med <- temp1 %>% dplyr::group_by(embryoID) %>% dplyr::summarise(medLevel = median(FC, na.rm = T))
    temp1med.t <- as.data.frame(as.matrix(t(temp1med[,-1])))
    colnames(temp1med.t) <- temp1med$embryoID
    temp2med <- temp2 %>% dplyr::group_by(embryoID) %>% dplyr::summarise(medLevel = median(FC, na.rm = T))
    temp2med.t <- as.data.frame(as.matrix(t(temp2med[,-1])))
    colnames(temp2med.t) <- temp2med$embryoID
    temp3med <- temp3 %>% dplyr::group_by(embryoID) %>% dplyr::summarise(medLevel = median(FC, na.rm = T))
    temp3med.t <- as.data.frame(as.matrix(t(temp3med[,-1])))
    colnames(temp3med.t) <- temp3med$embryoID
    tempmed <- cbind(temp1med.t, temp2med.t, temp3med.t)
    GOAcrossStages[i,12:(11+ncol(tempmed))] <- tempmed
  } 
  
  GOAcrossStages$Annotation[i] <- unique(temp.GO$annotation)
  GOAcrossStages$Proteins[i] <- paste(unique(temp.prot$GeneName), collapse = " ")
  GOAcrossStages$NumProteins[i] <- length(unique(temp.prot$GeneName))
  GOAcrossStages$NumProteinsInvolved[i] <- length(unique(temp.GO2$GeneSymbol))
  GOAcrossStages$fraction[i] <- GOAcrossStages$NumProteins[i] / GOAcrossStages$NumProteinsInvolved[i]
  print(i)
}


q.value.filter <- 0.05
## OPTION 1 original filter: 
GOAcrossStages.NA.pre <- GOAcrossStages %>% dplyr::filter(!is.na(p.value)) %>% 
  dplyr::filter(NumProteins > 2 & fraction >= 0.5)
# keep only unique protein lists
GOAcrossStages.NA <- GOAcrossStages.NA.pre #%>% dplyr::distinct(Proteins, .keep_all = TRUE)

## OPTION 2 filter on just GO term coverage, and then filter out GO terms with same exact proteins 
GOAcrossStages.NA.pre <- GOAcrossStages %>% dplyr::filter(!is.na(p.value)) %>% 
  dplyr::filter(fraction >= 0.6) 
GOAcrossStages.NA <- GOAcrossStages.NA.pre %>% dplyr::distinct(Proteins, .keep_all = TRUE)
## OPTION 3 filter on just NumProteins before q value correction 
GOAcrossStages.NA.pre <- GOAcrossStages %>% dplyr::filter(!is.na(p.value)) %>% 
  dplyr::filter(NumProteins > 1) 
GOAcrossStages.NA <- GOAcrossStages.NA.pre %>% dplyr::distinct(Proteins, .keep_all = TRUE)
# obtain which rows are duplicated: result <- GOAcrossStages.NA %>% group_by(Proteins) %>% filter(n()>1)
# filter out duplicated/repeated rows
#GOAcrossStages.NA <- GOAcrossStages.NA.pre %>% dplyr::group_by(Proteins) %>% dplyr::mutate(duplicate = n()) %>% 
#  dplyr::filter(duplicate == 1) %>% dplyr::select(-duplicate)

GOAcrossStages.NA$q.value <- p.adjust(GOAcrossStages.NA$p.value, method = "fdr")
GOAcrossStages.NA.sig <- GOAcrossStages.NA %>% dplyr::filter(q.value < q.value.filter) #%>% dplyr::filter(fraction >= 0.5)
colnames(GOAcrossStages.NA.sig)[3] <- "Correlation"



# let's plot some of the GO terms in terms of just the correlation 
Go.plot <- GOAcrossStages.NA.sig %>% #dplyr::filter(fraction >=0.5) %>% 
  #dplyr::filter(Correlation <= -0.2) %>% dplyr::filter(!(grepl("cardiac|viral|sperm|male|glial|dendritic|neuron", Annotation)))
  dplyr::filter(Correlation >= 0.1) %>% dplyr::filter(!(grepl("synaptic|cholesterol", Annotation)))
Go.plot <- Go.plot[order(Go.plot$Correlation),]
Go.plot$Annotation <- factor(Go.plot$Annotation, levels = Go.plot$Annotation)
ggplot(data = Go.plot, aes(y = Annotation, x= Correlation, fill = Annotation)) + 
  scale_fill_brewer(palette = "Blues") + 
  geom_bar(stat = "identity", color = "black") +
  theme_bw() + 
  theme(axis.text = element_text(size = 18, color = "black"), axis.title = element_text(size = 18),
        legend.position = "none", 
        #panel.background = element_rect(fill = 'gray50', color = 'purple'), 
        #panel.grid.major = element_line(color = 'gray87', linetype = 'dotted'), 
        #panel.grid.minor = element_line(color = 'gray87', linetype = 'dotted')
  ) + 
  labs(x = "Spearman Correlation", y = "Annotation") + 
  scale_x_continuous(breaks = c(0, 0.25, 0.5))
# OPTION 1 Plotting 
Go.plot1 <- Go.plot %>% dplyr::filter(grepl("DNA|proteasome|aspartate", Annotation))#[4:6,] #%>% dplyr::filter(Correlation > 0.3)
Go.plot2 <- reshape2::melt(Go.plot1[,c(2,12:68)])
Go.plot2$type <- dfStore$type[match(Go.plot2$variable, dfStore$embryoID)]
Go.plot2$type <- factor(Go.plot2$type, levels = c("early_2c", "late_2c", "_4c"))
ggplot(Go.plot2, aes(x  = type, y = value, fill = Annotation)) + geom_boxplot() + 
  facet_wrap(~Annotation,labeller = label_wrap_gen()) + 
  scale_fill_brewer(palette = "Blues") + 
  #coord_cartesian(ylim = c(-0.65,0.65)) +
  coord_cartesian(ylim = c(-0.3, 0.4)) + 
  theme_bw() + geom_hline(yintercept = 0, linetype = "dashed") + 
  theme(axis.text.x = element_text(size = 22, angle = 45, hjust = 1, color = "black"), #axis.title = element_text(size = 22), 
        axis.title = element_blank(),
        axis.text.y = element_text(size = 22, color = "black"),
        strip.text = element_text(size = 24), legend.position = "none") + 
  labs(y = "Log2(Alpha - Beta)", x = "Developmental Stage") + 
  scale_x_discrete(labels = c("Early 2-c", "Late 2-c", "4-cell"))
# OPTION 2 PLotting 
# for each GO term - need the protein FCs per embryo 
Go.plot1 <- c()
i<-2
for(i in 1:length(unique(Go.plot$Annotation))) {
  temp.GO <- GOinDat %>% dplyr::filter(GOtermID == Go.plot$GOterm[i])
  temp.gene <- temp.GO$GeneSymbol
  temp.prot <- FC.2 %>% dplyr::filter(GeneName %in% temp.gene)
  temp.prot$GO_Annotation <- Go.plot$Annotation[i]
  temp.prot$GO_ID <- Go.plot$GOterm[i]
  if(i == 1) { Go.plot1 <- temp.prot}
  if(i != 1) { Go.plot1 <- rbind(Go.plot1, temp.prot)}
}
Go.plot1$type <- dfStore$type[match(Go.plot1$embryoID, dfStore$embryoID)]
Go.plot1$type <- factor(Go.plot1$type, levels = c("early_2c", "late_2c", "_4c"))
ggplot(Go.plot1, aes(x  = type, y = FC, fill = GO_Annotation)) + geom_boxplot() + 
  facet_wrap(~GO_Annotation,labeller = label_wrap_gen()) + 
  scale_fill_brewer(palette = "Blues") + 
  coord_cartesian(ylim = c(-0.65,0.65)) +
  theme_bw() + geom_hline(yintercept = 0, linetype = "dashed") + 
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1), axis.title = element_text(size = 20), 
        axis.text.y = element_text(size = 20),
        strip.text = element_text(size = 16), legend.position = "none") + 
  labs(y = "Log2(Alpha - Beta)", x = "Developmental Stage")





# FOR NEGATIVE CORRELATION
Go.plot <- GOAcrossStages.NA.sig %>% 
  dplyr::filter(Correlation <= -0.1)
#dplyr::filter(Correlation <= -0.29) %>% 
#dplyr::filter(!(grepl("cardiac|viral|sperm|male|glial|dendritic|neuron|steriod|amyloid", Annotation)))

Go.plot <- Go.plot[order(Go.plot$Correlation),]
Go.plot$Annotation <- factor(Go.plot$Annotation, levels = Go.plot$Annotation)
Go.plot$Annotation_short <- substr(Go.plot$Annotation, 1, 34)
Go.plot$Annotation_short <- factor(Go.plot$Annotation_short, levels = Go.plot$Annotation_short)
ggplot(data = Go.plot, aes(y = Annotation_short, x= Correlation, fill = Annotation)) + 
  geom_bar(stat = "identity", color = "black") + scale_fill_brewer(palette = "Greens") + 
  theme_bw() + 
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18),
        legend.position = "none") + 
  xlab("Spearman Correlation")
Go.plot1 <- Go.plot %>% dplyr::filter(grepl("thioredoxin|dATP", Annotation))#[1:3,] #%>% dplyr::filter(Correlation < -0.32)
Go.plot2 <- reshape2::melt(Go.plot1[,c(2,12:68, 70)])
Go.plot2$type <- dfStore$type[match(Go.plot2$variable, dfStore$embryoID)]
Go.plot2$type <- factor(Go.plot2$type, levels = c("early_2c", "late_2c", "_4c"))
ggplot(Go.plot2, aes(x  = type, y = value, fill = Annotation_short)) + geom_boxplot() + 
  facet_wrap(~Annotation_short,labeller = label_wrap_gen()) + 
  scale_fill_brewer(palette = "Greens") + 
  coord_cartesian(ylim = c(-0.65,0.2)) +
  theme_bw() + geom_hline(yintercept = 0, linetype = "dashed") + 
  theme(axis.text.x = element_text(size = 22, angle = 45, hjust = 1, color = "black"), #axis.title = element_text(size = 22),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 22, color = "black"),legend.position = "none",
        strip.text = element_text(size = 24)) + 
  labs(y = "Log2(Alpha - Beta)", x = "Developmental Stage") + 
  scale_x_discrete(labels = c("Early 2-c", "Late 2-c", "4-cell"))





Go.plot <- Go.plot %>% dplyr::mutate(Class = case_when(grepl("oxi|oxy|aldehyde|electron", Annotation) == TRUE ~ "Oxidation", 
                                                       #grepl("catabolic|ubiquit", Annotation) == TRUE ~ "Degradation", 
                                                       grepl("binding|transferase", Annotation) == TRUE ~ "Modification", 
                                                       grepl("localization|export", Annotation) == TRUE ~ "Protein Movement"))
Go.plot1 <- Go.plot %>% dplyr::filter(!(is.na(Class))) #%>% dplyr::filter(Class != "Modification" )
Go.plot2 <- reshape2::melt(Go.plot1[,c(2,12:68, 70)])
Go.plot2$type <- dfStore$type[match(Go.plot2$variable, dfStore$embryoID)]
Go.plot2$type <- factor(Go.plot2$type, levels = c("early_2c", "late_2c", "_4c"))
Go.plot2$Annotation <- as.character(Go.plot2$Annotation)
#Go.plot2$Annotation <- factor(Go.plot2$Annotation, levels = c(unique(Go.plot2$Annotation)[c(2,3,5)], 
#                                                              unique(Go.plot2$Annotation)[c(1,4,6,7,8)]))
ggplot(Go.plot2, aes(x  = type, y = value, fill = Annotation)) + geom_boxplot() + 
  facet_wrap(~Class,labeller = label_wrap_gen()) + 
  coord_cartesian(ylim = c(-0.65,0.25)) +
  theme_bw() + geom_hline(yintercept = 0, linetype = "dashed") + 
  theme(axis.text.x = element_text(size = 15), axis.title = element_text(size = 20), 
        axis.text.y = element_text(size = 20),
        strip.text = element_text(size = 20), legend.position = "bottom", 
        legend.direction = "vertical", legend.text = element_text(size = 18))
# focus on one class at a time
Go.plot3 <- Go.plot2 %>% dplyr::filter(Class == "Modification")
ggplot(Go.plot3, aes(x  = type, y = value, fill = Annotation)) + geom_boxplot() + 
  coord_cartesian(ylim = c(-0.65,0.25)) +
  theme_bw() + geom_hline(yintercept = 0, linetype = "dashed") + 
  theme(axis.text.x = element_text(size = 18), axis.title = element_text(size = 20), 
        axis.text.y = element_text(size = 20),
        strip.text = element_text(size = 20), 
        #legend.position = "bottom", 
        legend.direction = "vertical", legend.text = element_text(size = 18), 
        legend.position = c(0.5, 0.15)#, 
        #legend.background = element_rect(fill = "transparent", color = "transparent")
  ) +
  ylab(bquote("Median"~log[2]~"FC (Alpha - Beta)")) + 
  labs(x = "Developmental Stage")
Go.plot3 <- Go.plot2 %>% dplyr::filter(Class == "Oxidation")
ggplot(Go.plot3, aes(x  = type, y = value, fill = Annotation)) + geom_boxplot() + 
  coord_cartesian(ylim = c(-0.65,0.25)) +
  theme_bw() + geom_hline(yintercept = 0, linetype = "dashed") + 
  theme(axis.text.x = element_text(size = 18), axis.title = element_text(size = 20), 
        axis.text.y = element_text(size = 20),
        strip.text = element_text(size = 20), 
        #legend.position = "bottom", 
        legend.direction = "vertical", legend.text = element_text(size = 18), 
        legend.position = c(0.41, 0.15), 
        legend.background = element_rect(fill = "transparent", color = "transparent")
  ) +
  ylab(bquote("Median"~log[2]~"FC (Alpha - Beta)")) + 
  labs(x = "Developmental Stage")
