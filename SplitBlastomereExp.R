WindowsPath <- "G:/.shortcut-targets-by-id/1uQ4exoKlaZAGnOG1iCJPzYN3ooYYZB7g/MS/Archive"
source(paste0(WindowsPath, "/Ribo/Aleks_LungCancer/Code/LibsFunsAP.R"))

reference_channel = 2

# batch2 data files
design2 <- read.csv(paste0(WindowsPath, "/Collaborators/Zernicka_Goetz/dat/Experiment18_newbatch_split/embryo_pSCoPE_scSets/annotation_splitBlastomeres_part2.csv"))
# batch 1 related files 
design1 <- read.csv(paste0(WindowsPath, "/Collaborators/Zernicka_Goetz/dat/Experiment16_single_blastomeres/annotation_splitBlastomeres.csv"))
# combine the design files 
design <- rbind(design1,design2)
design.melt <- reshape2::melt(design, id.vars = "Set")
design.melt$set_RI <- paste(design.melt$Set, design.melt$variable)
design.melt$set_RI_celltype <- paste(design.melt$Set, design.melt$variable, design.melt$value)

# dart-id ran on all sets and both batchee 
ev.load <- read.delim(paste0(WindowsPath,"/Collaborators/Zernicka_Goetz/dat/DART_ID_output/06_exp18_exp16_combined_split/dartid_output/ev_updated.txt"))
ev <- ev.load %>% dplyr::filter(grepl("105|106|107|108|162|163|164|165|166", Raw.file))


######################################################################################################.
################################### Peptide Data Prep (Filter) #################################################
######################################################################################################.

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
length(unique(ev$Leading.razor.protein))
# filter for PEP and PIF 
ev.1 <- ev %>% dplyr::filter(dart_PEP <= 0.03) %>% dplyr::filter(PIF >=0.8)
length(unique(ev.1$Leading.razor.protein))
# select the most necessary columns
ev.2 <- ev.1 %>% dplyr::select(Raw.file, Leading.razor.protein, SeqCharge, contains("intensity.corrected"))
# change O's to NA's 
ev.2[ev.2 == 0] <- NA 



######################################################################################################.
################################### Peptide Data Prep (CV Filter) #################################################
######################################################################################################.
ev.2.ref <- ev.2
RI.ref <- paste0("Reporter.intensity.corrected.", reference_channel)
ev.2.ref[,4:19] <- ev.2.ref[,4:19] / ev.2.ref[,which(colnames(ev.2.ref) == RI.ref)]

# melt the dataframe
ev.2.melt <- reshape2::melt(ev.2.ref)
ev.2.melt$RI <- paste0("RI",gsub("(?:[^.]+\\.){3}([^.]+).*", "\\1", ev.2.melt$variable))
ev.2.melt$Raw.file_RI <- paste(ev.2.melt$Raw.file, ev.2.melt$RI)

# now we can match the annotation file with this dataframe
ev.2.melt$celltype <- design.melt$value[match(ev.2.melt$Raw.file_RI, design.melt$set_RI)]
ev.2.melt$Raw.file_RI_celltype <- paste(ev.2.melt$Raw.file_RI, ev.2.melt$celltype)
# filter out the carrier and the unused channels 
ev.2.melt.quant <- ev.2.melt %>% dplyr::filter(!grepl("carrier|unused|reference|sc_x|sc_pos|sc_c", Raw.file_RI_celltype))

# take the median of repeat peptide observations 
ev.2.melt.quant <- ev.2.melt.quant %>% dplyr::group_by(Raw.file_RI_celltype, Leading.razor.protein, SeqCharge) %>% 
  dplyr::summarise(medPep = median(value, na.rm = T))
# unmelt so we can do relative quant 
ev.2.melt.quant.d <- reshape2::dcast(ev.2.melt.quant, Leading.razor.protein + SeqCharge  ~ Raw.file_RI_celltype, 
                                     value.var = "medPep")

# normalize for sample loading
ev.2.melt.quant.d <- columned(ev.2.melt.quant.d, start.num = 3, end.num = ncol(ev.2.melt.quant.d))
# normalize each row by the median 
ev.2.melt.quant.d$rowMedians <- apply(ev.2.melt.quant.d[,3:ncol(ev.2.melt.quant.d)],1, median, na.rm = TRUE)
ev.2.melt.quant.d <- rowed(ev.2.melt.quant.d, start.num = 3, end.num = (ncol(ev.2.melt.quant.d)-1), static.num = ncol(ev.2.melt.quant.d))
#remove the column with the rowMedians
ev.2.melt.quant.d <- ev.2.melt.quant.d[,-ncol(ev.2.melt.quant.d)]
# melt the row normalized peptide data 
sc.CV <- reshape2::melt(ev.2.melt.quant.d)

sc.CV.1 <- sc.CV %>% dplyr::group_by(variable, Leading.razor.protein) %>% 
  dplyr::summarise(CV = sd(value, na.rm = T) / mean(value, na.rm = T))
SC.CV.median <- sc.CV.1 %>% dplyr::group_by(variable) %>% dplyr::summarise(medCV = median(CV, na.rm = T))
SC.CV.median$type <- ifelse(grepl("control",SC.CV.median$variable), "control", "sc")



cells2keep <- SC.CV.median %>% dplyr::filter(medCV < 0.4)
cells2keep$embryo <- sub(".*_", "", cells2keep$variable)



ggplot(data = SC.CV.median, aes(x = medCV, fill = type)) + geom_histogram(color = "black") + 
  theme_bw() + theme(axis.title = element_text(size = 16), axis.text = element_text(size = 16), title = element_text(size = 16)) + 
  labs(title = "Median CV of each single cell in late 2-cell stage\nSecond Batch", 
       subtitle = paste0("# of cells kept = ", nrow(cells2keep))) + 
  geom_vline(xintercept = 0.4, linetype = "dashed")


######################################################################################################.
################################### Peptides Normalization #################################################
######################################################################################################.
# normalize to the reference 
ev.3 <- ev.2
RI.ref <- paste0("Reporter.intensity.corrected.", reference_channel)
ev.3[,4:19] <- ev.3[,4:19] / ev.3[,which(colnames(ev.3) == RI.ref)]

# melt 
ev.3.melt <- reshape2::melt(ev.3, id.vars = c("Raw.file", "Leading.razor.protein", "SeqCharge"))
ev.3.melt$RI <- paste0("RI",gsub("(?:[^.]+\\.){3}([^.]+).*", "\\1", ev.3.melt$variable))
ev.3.melt$Raw.file_RI <- paste(ev.3.melt$Raw.file, ev.3.melt$RI)

# now we can match the annotation file with this dataframe
ev.3.melt$celltype <- design.melt$value[match(ev.3.melt$Raw.file_RI, design.melt$set_RI)]
ev.3.melt$Raw.file_RI_celltype <- paste(ev.3.melt$Raw.file_RI, ev.3.melt$celltype)
# filter out the carrier and the ununsed channels 
ev.3.melt.quant <- ev.3.melt %>% dplyr::filter(!grepl("carrier|unused|reference|control|sc_x|sc_c", Raw.file_RI_celltype)) %>% 
  dplyr::filter(Raw.file_RI_celltype %in% cells2keep$variable)


# take the median of repeat peptide observations 
ev.3.melt.quant.2 <- ev.3.melt.quant %>% dplyr::group_by(Raw.file_RI_celltype, Leading.razor.protein, SeqCharge) %>% 
  dplyr::summarise(medPep = median(value, na.rm = T))
# unmelt so we can do relative quant
ev.4 <- reshape2::dcast(ev.3.melt.quant.2, Leading.razor.protein + SeqCharge  ~ Raw.file_RI_celltype, 
                        value.var = "medPep")

# count the number of peptides per blastomere 
numpeptides <- as.data.frame(colSums(!is.na(ev.4[,3:ncol(ev.4)])))

# normalize for sample loading
ev.6 <- columned(ev.4, start.num = 3, end.num = ncol(ev.4))
# normalize each row by the median 
ev.6$rowMedians <- apply(ev.6[,3:ncol(ev.6)],1, median, na.rm = TRUE)
ev.7a <- rowed(ev.6, start.num = 3, end.num = (ncol(ev.6)-1), static.num = ncol(ev.6))
#remove the column with the rowMedians
ev.7b <- ev.7a[,-ncol(ev.7a)]
# melt the row normalized peptide data 
ev.7b.melt <- reshape2::melt(ev.7b)

# obtain number of peptides 
numberPeps <- ev.7b.melt %>% 
  dplyr::group_by(variable) %>% 
  dplyr::summarise(numberOFpeptides = sum(!is.na(value)))

length(unique(ev.7b.melt$Leading.razor.protein))

######################################################################################################.
################# Peptides --> Proteins  ################################################################
######################################################################################################.

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
#write.csv(sc.collapse.d, "Exp2_late2cell_CollapsedPeptides.csv")
pct.r = 0.9
sc.collapse.1.missing <- c()
for(k in 1:nrow(sc.collapse.d)){
  
  pct.na<-length(which(is.na(sc.collapse.d[k,]))) / length(sc.collapse.d[k,])
  if(pct.na <= pct.r){ sc.collapse.1.missing<-rbind(sc.collapse.1.missing, sc.collapse.d[k,])}
  
}

# obtain number of proteins 
numberProteins <- data.frame(colSums(!is.na(sc.collapse.1.missing[,2:ncol(sc.collapse.1.missing )])) )

# column normalize the proteins x samples matrix 
sc.collapse.1 <- columned_log2(sc.collapse.1.missing, start.num = 2, end.num = ncol(sc.collapse.1.missing))
sc.collapse.1.m <- reshape2::melt(sc.collapse.1)
sc.collapse.1.m$embryo <- sub(".* ", "", sc.collapse.1.m$variable)
ggplot(data = sc.collapse.1.m,  aes(x = variable, y = value)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# row normalize the proteins x samples matrix 
sc.collapse.1$rowMeans <- rowMeans(sc.collapse.1[,2:ncol(sc.collapse.1)], na.rm = T)
sc.collapse.2 <- rowed_log2(sc.collapse.1, start.num = 2, end.num = ncol(sc.collapse.1), static.num = ncol(sc.collapse.1))
sc.splitExperiment <- sc.collapse.2[,-ncol(sc.collapse.2)]

#write.csv(sc.splitExperiment, "C:/Users/apete/Desktop/Embryo_Project_2022/Split_Experiment_ProteinsXsamples.csv")




######################################################################################################.
################# Prepare cell count data ################################################################
######################################################################################################.

# prepare the biological cell count data 
LineageCounts <- 
  read_xlsx(paste0(WindowsPath,
                   "/Collaborators/Zernicka_Goetz/dat/Experiment18_16_combined_dataFiles/LineageCounts_Alpha_Beta_toDrive_batch1AND2.xlsx"),
                           sheet = 6)

LineageCounts <- LineageCounts %>% dplyr::mutate(PlateNumber = case_when(trial == 1 & plate == 1 ~ 1, 
                                                                         trial == 1 & plate == 2 ~ 2, 
                                                                         trial == 2 & plate == 1 ~ 3, 
                                                                         trial == 2 & plate == 2 ~ 4,
                                                                         trial == 3 & plate == 1 ~ 5, 
                                                                         trial == 4 & plate == 1 ~ 6))
LineageCounts$PlateNumber_RowLetter <- paste0(LineageCounts$PlateNumber, "_", LineageCounts$`Row letter`)
LineageCounts$PlateNumber_RowLetter_Column <- paste0(LineageCounts$PlateNumber_RowLetter, LineageCounts$Column)

PlateNumber <- read.csv(paste0(WindowsPath,
                               "/Collaborators/Zernicka_Goetz/dat/Experiment18_16_combined_dataFiles/PlateNumber_bothBatches.csv"))
RowLetter <- read.csv(paste0(WindowsPath,
                             "/Collaborators/Zernicka_Goetz/dat/Experiment18_16_combined_dataFiles/annotation_splitBlastomeres_BothBatches.csv"))

# manipulate plate number and row letter information 
PlateNumber.m <- reshape2::melt(PlateNumber, id.vars = "Set")
PlateNumber.m$Set_RI <- paste(PlateNumber.m$Set, PlateNumber.m$variable)
RowLetter.m <- reshape2::melt(RowLetter, id.vars = "Set")
RowLetter.m$Set_RI <- paste(RowLetter.m$Set, PlateNumber.m$variable)

PlateNumber.m$RowLetter <- RowLetter.m$value[match(PlateNumber.m$Set_RI, RowLetter.m$Set_RI)]
PlateNumber.m$PlateNumber_RowLetter <- paste0(PlateNumber.m$value,"_", PlateNumber.m$RowLetter)

# match back the single blastomere id 
LineageCounts$SingleCellId <- PlateNumber.m$Set_RI[match(LineageCounts$PlateNumber_RowLetter_Column,PlateNumber.m$PlateNumber_RowLetter)]

# Normalized cell counts
LineageCounts$normEpi <- LineageCounts$epi / LineageCounts$total
LineageCounts$normTE <- LineageCounts$te / LineageCounts$total
LineageCounts$normPe <- LineageCounts$pe / LineageCounts$total

# find the number of lineages that are present in each embryo 
LineageCounts <- LineageCounts %>% 
  dplyr::mutate(Num_Lineages_Present = case_when(te>0 & epi>0 & pe>0  ~ 3, 
                                        te==0 & epi==0 & pe==0 ~ 0, 
                                        
                                        te==0 & epi>0 & pe>0  ~ 2, 
                                        te>0 & epi==0 & pe>0  ~ 2, 
                                        te>0 & epi>0 & pe==0  ~ 2,
                                        
                                        te==0 & epi==0 & pe>0  ~ 1, 
                                        te==0 & epi>0 & pe==0  ~ 1, 
                                        te>0 & epi==0 & pe==0  ~ 1))

LineageCounts$Num_Lineages_Present_char <- as.character(LineageCounts$Num_Lineages_Present)

# combined cell counts and fractions
LineageCounts$te_epi <- LineageCounts$te + LineageCounts$epi
LineageCounts$te_pe <- LineageCounts$te + LineageCounts$pe
LineageCounts$pe_epi <- LineageCounts$pe + LineageCounts$epi

LineageCounts$te_epi_fraction <- (LineageCounts$te + LineageCounts$epi)/ LineageCounts$total
LineageCounts$te_pe_fraction <- (LineageCounts$te + LineageCounts$pe)/ LineageCounts$total
LineageCounts$pe_epi_fraction <- (LineageCounts$pe + LineageCounts$epi)/ LineageCounts$total


# download alpha-beta Proteins
useList <- read.csv(paste0(WindowsPath, "/Collaborators/Zernicka_Goetz/code/CodeForGithub/2022_10_13_2cellstage_SignificantProteins_kmeansC__.csv"))
test <- useList %>% dplyr::filter(GreaterIN == "Cluster2_Beta")
# intersect current data with this list 
sigProt <-sc.splitExperiment %>% dplyr::filter(Leading.razor.protein %in% useList$Protein)
sigprotm <- reshape2::melt(sigProt)
sigprotm$proteinType <- useList$GreaterIN[match(sigprotm$Leading.razor.protein, useList$Protein)]
# plot distributions of each blastomere 
ggplot(data = sigprotm, aes(x = variable, color = proteinType, y = value)) + geom_boxplot()

#write.csv(sigProt, paste0(Sys.Date(), "SigProteins_splitblastExp.csv"))





######################################################################################################.
################# Protein Abundance ################################################################
######################################################################################################.


# MEDIAN ABUNDANCE OF ALPHA AND BETA PROTEINS 
# obtain median value of alpha proteins and beta proteins for each blastomere
sigprotm.1 <- sigprotm %>% 
  dplyr::group_by(variable, proteinType) %>% 
  dplyr::summarise(medianProtein = median(value,na.rm = T))
# obtain the number of proteins quantified 
sigprotm.1a <- sigprotm %>% 
  dplyr::group_by(variable, proteinType) %>% 
  dplyr::summarise( numProteins = sum(!is.na(value)))
ggplot(data = sigprotm.1a, aes(x = numProteins, fill = proteinType)) + geom_histogram() + theme_bw() + 
  theme(text = element_text(size = 18), legend.position = "bottom")
# obtain the fold-change between the alpha and beta protein medians 
# what was originally used: 
sigprotm.2 <- sigprotm.1 %>% 
  dplyr::group_by(variable) %>% 
  dplyr::summarise(Alpha_Beta_ = medianProtein[proteinType == "Cluster2"] - medianProtein[proteinType == "Cluster1"])


# remove the "b" in each single cell id 
sigprotm.2$singleCellID <- sub(" b.*", "", sigprotm.2$variable)
ggplot(data = sigprotm.2, aes(x = Alpha_Beta)) + geom_histogram() + geom_density() + 
  theme_classic() + 
  theme(text = element_text(size =25))




# now see if there is any relation to biological output that we have 
LineageCounts$Alpha_Beta <- sigprotm.2$Alpha_Beta[match(LineageCounts$SingleCellId, sigprotm.2$singleCellID)]
LineageCounts$Alpha <- sigprotm.2$Alpha[match(LineageCounts$SingleCellId, sigprotm.2$singleCellID)]
LineageCounts$Beta <- sigprotm.2$Beta[match(LineageCounts$SingleCellId, sigprotm.2$singleCellID)]
LineageCounts$Alpha_Beta_ <- sigprotm.2$Alpha_Beta_[match(LineageCounts$SingleCellId, sigprotm.2$singleCellID)]
# save csv and remove some columns after
#write.csv(LineageCounts, "LineageCounts.csv")


LineageCounts.plot <- LineageCounts %>% dplyr::filter(!is.na(Alpha_Beta))

LineageCounts.plot$pe_char <- as.character(LineageCounts.plot$pe)

##### FOCUS ON PROTEIN ABUNDANCES PLOTTING ######



# remove blastomeres with sister that have 10 or less cells overall or 2 lineages absent 
LineageCounts.plot.remove <- LineageCounts.plot %>% dplyr::filter(total > 10) %>% dplyr::filter(Num_Lineages_Present != 1) 
round(cor(LineageCounts.plot.remove$normEpi, LineageCounts.plot.remove$Alpha_Beta, use = 'pairwise.complete.obs'),digits =3)

p<-cor.test(LineageCounts.plot.remove$normEpi, LineageCounts.plot.remove$Alpha_Beta_, use = 'pairwise.complete.obs')$p.value
filename<- paste(gsub(":", "-", Sys.Date()),"_file.csv",sep="")
#write.csv(LineageCounts.plot.remove, paste0("C:/Users/apete/Desktop/Embryo_Project_2022/", filename))


# PLOT BELOW
ggplot(data = LineageCounts.plot.remove, 
       aes(y = normEpi, x = Alpha_Beta_, size = total, color = Num_Lineages_Present_char)) + 
  geom_point(alpha = 0.65) +  
  theme_bw() + 
  theme(axis.title = element_text(size = 25), 
        legend.position = "bottom",legend.box = "vertical", 
        legend.text = element_text(size = 19), 
        plot.title = element_text(size = 25), 
        axis.text = element_text(size = 19), 
        text = element_text(size = 25)) + 
  scale_color_manual(values = c("#CD6600", "#9A32CD"), name = "# of Lineages Present") + 
  labs(#x = "Epiblast Cells / Total Cell Count", y = "Alpha / Beta", 
    y = "Proportion of epiblast cells in\nblastocyst from cultured cell",
    x = "Alpha-Beta Fold Change of analyzed cell",
    size = "Total Cell Count",
    title = paste0("Correlation = ", 
                   round(cor(LineageCounts.plot.remove$normEpi, LineageCounts.plot.remove$Alpha_Beta_, 
                             use = 'pairwise.complete.obs'),digits =3),", ",
                   paste0("p value = ", round(p, digits = 4)))) + 
  
  guides(color = guide_legend(override.aes = list(size = 10)))


# take a look at raw epiblast cell counts 
LineageCounts.plot$epi_char <- as.character(LineageCounts.plot$epi)
LineageCounts.plot.remove <- LineageCounts.plot %>% dplyr::filter(total > 10) %>% dplyr::filter(Num_Lineages_Present != 1) 


# correct plot below 
LineageCounts.plot.remove$EPI_Health <- ifelse(LineageCounts.plot.remove$epi > 3, "Healthy", "Less Healthy")

LineageCounts.plot.remove.fil <- LineageCounts.plot.remove %>% dplyr::filter(EPI_Health == "Healthy")
LineageCounts.plot.remove.fil2 <- LineageCounts.plot.remove %>% dplyr::filter(EPI_Health != "Healthy")
epip = t.test(LineageCounts.plot.remove.fil$Alpha_Beta_, LineageCounts.plot.remove.fil2$Alpha_Beta_)

unique(design2$Set)
LineageCounts.plot.remove.2exp <- LineageCounts.plot.remove %>% dplyr::filter(grepl("AP162|AP163|AP164|AP165|AP166",SingleCellId))
LineageCounts.plot.remove.2exp1 <- LineageCounts.plot.remove.2exp %>% dplyr::filter(EPI_Health == "Healthy")
LineageCounts.plot.remove.2exp2 <- LineageCounts.plot.remove.2exp %>% dplyr::filter(EPI_Health != "Healthy")
epip.2 = t.test(LineageCounts.plot.remove.2exp1$Alpha_Beta_, LineageCounts.plot.remove.2exp2$Alpha_Beta_)
epip.2$p.value

ggplot(data = LineageCounts.plot.remove, aes(x = EPI_Health, y = Alpha_Beta_, fill = EPI_Health)) + 
  #geom_boxplot(fill = 'gray95') +
  geom_violin() +
  #geom_point(color  = 'red', alpha = 0.67, aes(size = total, shape = Num_Lineages_Present_char)) + 
  geom_jitter(aes(size = total), position = position_jitter( 0.15), alpha = 0.67)+
  theme_bw() + 
  theme(text = element_text(size = 25), legend.position = "bottom",legend.box = "vertical") + 
  guides(shape = guide_legend(override.aes = list(size = 7))) + 
  scale_fill_brewer(palette = "Blues") +
  scale_x_discrete(labels = c(expression("epi">="4"), expression("epi"<"4"))) +
  labs(x = "No. of EPI cells", 
       y = "Alpha-Beta Fold Change of\n analyzed cell", 
       shape = "# of Lineages Present", 
       size = "Total Cell Count", 
       fill  = 'Epiblast Health', 
       title = paste0('p-value = ', round(epip$p.value, digits = 6)) )  







LineageCounts.plot.remove <- LineageCounts.plot %>% dplyr::filter(total > 10) %>% dplyr::filter(Num_Lineages_Present != 1) %>% dplyr::filter(normPe != 0)
round(cor(LineageCounts.plot.remove$normPe, LineageCounts.plot.remove$Alpha_Beta, use = 'pairwise.complete.obs'),digits =3)
p<-cor.test(LineageCounts.plot.remove$normPe, LineageCounts.plot.remove$Alpha_Beta, use = 'pairwise.complete.obs')$p.value
p


# PE plot: 
p<-cor.test(LineageCounts.plot.remove$normPe, LineageCounts.plot.remove$Alpha_Beta_, use = 'pairwise.complete.obs')$p.value

ggplot(data = LineageCounts.plot.remove, 
       aes(y = normPe, x = Alpha_Beta_, size = total, color = Num_Lineages_Present_char)) + 
  geom_point(alpha = 0.65) +  
  theme_bw() + 
  theme(text = element_text(size = 25), 
        legend.position = "bottom",legend.box = "vertical") + 
  scale_color_manual(values = c("#CD6600", "#9A32CD"), name = "# of Lineages Present") + 
  guides(color = guide_legend(override.aes = list(size = 10))) + 
  labs(#x = "Epiblast Cells / Total Cell Count", y = "Alpha / Beta", 
    y = "Proportion of PE cells in\nblastocyst from cultured cell",
    x = "Alpha-Beta Fold Change of analyzed cell",
    size = "Total Cell Count",
    title = paste0("Correlation = ", 
                   round(cor(LineageCounts.plot.remove$normPe, LineageCounts.plot.remove$Alpha_Beta_, use = 'pairwise.complete.obs'),digits =3)), 
    subtitle = paste0("P value = ", round(p, digits = 4)))






# corrected plot below
p<-cor.test(LineageCounts.plot.remove$total, LineageCounts.plot.remove$Alpha_Beta_, use = 'pairwise.complete.obs')$p.value
ggplot(data = LineageCounts.plot.remove, 
       aes(y = total, x = Alpha_Beta_, size = total, color = Num_Lineages_Present_char)) + 
  geom_point(alpha = 0.65) +  
  theme_bw() + 
  theme(text = element_text(size = 25), 
        legend.position = "bottom",legend.box = "vertical") + 
  scale_color_manual(values = c("#CD6600", "#9A32CD"), name = "# of Lineages Present") + 
  guides(color = guide_legend(override.aes = list(size = 10))) + 
  labs(#x = "Epiblast Cells / Total Cell Count", y = "Alpha / Beta", 
    y = "Total cells in blastocyst from\ncultured cell",
    x = "Alpha-Beta Fold Change of analyzed cell",
    size = "Total Cell Count",
    title = paste0("Correlation = ", 
                   round(cor(LineageCounts.plot.remove$total, LineageCounts.plot.remove$Alpha_Beta_, 
                             use = 'pairwise.complete.obs'),digits =3)), 
    subtitle = paste0("P value = ", round(p, digits = 4)))














