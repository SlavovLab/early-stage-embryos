library(readxl)
library(dplyr)
library(ComplexHeatmap)
library(stringr)
library(circlize)
pathmm <- "G:/.shortcut-targets-by-id/1uQ4exoKlaZAGnOG1iCJPzYN3ooYYZB7g/MS/Collaborators/Zernicka_Goetz"

# load the within_embryo normalization for the early 2c 
c <- read.csv(paste0(pathmm,"/code/CodeForGithub/Exp4_early2cell_ProteinsXSamples_WithinEmbryo.csv"))

# use the same list of proteins used in the split blastomere experiment
useList <- read.csv(paste0(pathmm, "/code/CodeForGithub/2022_10_13_2cellstage_SignificantProteins_kmeansC.csv"))



# experiment map, filtered for Exp04
design <- read_excel(paste0(pathmm,"/dat/EXP_Map.xlsx"), sheet = 1)
design <- design %>% dplyr::filter(grepl("_04", Exp)) 
design <- design[,-1]
design.melt <- reshape2::melt(design, id.vars = "Set")
design.melt$set_RI <- paste(design.melt$Set, design.melt$variable)
design.melt$set_RI_celltype <- paste(design.melt$Set, design.melt$variable, design.melt$value)
design.melt$set_RI_celltype_period <- str_replace_all(design.melt$set_RI_celltype, " ", ".")

# load in polar body information 
pb <- read_excel(paste0(pathmm,"/polarbody.xlsx"), sheet = 3)
pb <- pb[,-1]
pb.m <- reshape2::melt(pb, id.vars = "Set")
pb.m$set_RI <- paste(pb.m$Set, pb.m$variable)

# match to have the correct naming
pb.m$set_RI_celltype_period <- design.melt$set_RI_celltype_period[match(pb.m$set_RI, 
                                                                        design.melt$set_RI)]

# filter protein data and calculate correlation
ab <- c %>% dplyr::filter(X %in% useList$Protein)
ab.1 <- ab
rownames(ab.1) <- ab.1$X
ab.1 <- ab.1[,-1]
abcor <- cor(ab.1,  use = "pairwise.complete.obs")
Heatmap(abcor)

# now calculate alpha-beta character for each blastomere
ab.2 <- ab
ab.2$proteinType <- useList$GreaterIN[match(ab.2$X, useList$Protein)]
ab.2.m <- reshape2::melt(ab.2)

sigprotm.1 <- ab.2.m %>% 
  dplyr::group_by(variable, proteinType) %>% 
  dplyr::summarise(medianProtein = median(value,na.rm = T))
sigprotm.2 <- sigprotm.1 %>% 
  dplyr::group_by(variable) %>% 
  dplyr::summarise(Alpha_Beta = 
                     medianProtein[proteinType == "Cluster1_Alpha"] 
                   - medianProtein[proteinType == "Cluster2_Beta"], 
                   Alpha = medianProtein[proteinType == "Cluster1_Alpha"], 
                   Beta = medianProtein[proteinType == "Cluster2_Beta"], 
                   Alpha_Beta_corrected = 
                     medianProtein[proteinType == "Cluster2_Beta"] 
                   - medianProtein[proteinType == "Cluster1_Alpha"])

threshold_ab <- 0
sigprotm.2 <- sigprotm.2 %>% 
  dplyr::mutate(alpha_beta_character = case_when(Alpha_Beta_corrected > threshold_ab ~ "Alpha", 
                                                 Alpha_Beta_corrected < -threshold_ab ~ "Beta"))



column_anno <- data.frame(colnames(abcor))

column_anno$alpha_beata <- sigprotm.2$alpha_beta_character[match(column_anno$colnames.abcor., 
                                                                 sigprotm.2$variable)]
column_anno$PolarBody <- pb.m$value[match(column_anno$colnames.abcor., 
                                          pb.m$set_RI_celltype_period)]


column_ha <- HeatmapAnnotation(`Alpha or Beta` = column_anno$alpha_beata,
                               `Polar Body` = column_anno$PolarBody,
                               col = list(`Alpha or Beta` = c(Alpha = "#D60093", Beta ="#E8AC04"),
                                          `Polar Body` = c("x" = "bisque", "p" = "cadetblue1")))
colfun = colorRamp2(c(-1,0,1), c("purple", "white", "yellow"), space = 'LAB')

htlist <- Heatmap(abcor, show_row_names = FALSE, show_column_names = FALSE,#col = colfun,
                  bottom_annotation = column_ha, 
                  heatmap_legend_param = list(direction = 'vertical', title = "Correlation"))
draw(htlist, annotation_legend_side = 'bottom', heatmap_legend_side = 'right')

# count the number of polar bodies per cell type
count.1 <- column_anno %>% dplyr::filter(PolarBody == "p" & alpha_beata == "Beta")
count.2 <- column_anno %>% dplyr::filter(PolarBody == "p" & alpha_beata == "Alpha")
count.3 <- column_anno %>% dplyr::filter(alpha_beata == "Beta")
count.4 <- column_anno %>% dplyr::filter(alpha_beata == "Alpha")

dhyper(nrow(count.1), nrow(count.3), 
       nrow(count.4), nrow(count.1)+nrow(count.2))

dhyper(nrow(count.1), nrow(count.3), 
       6, 6)
