# understand the vegetal cell analysis more by including more blastomeres 
#########----------------------------------------------------##############.
##### paths and report file ####
#########----------------------------------------------------##############.
pathMS <- "G:/.shortcut-targets-by-id/1uQ4exoKlaZAGnOG1iCJPzYN3ooYYZB7g/MS/Archive"
pathMS_nA <- "G:/.shortcut-targets-by-id/1uQ4exoKlaZAGnOG1iCJPzYN3ooYYZB7g/MS"
pathRibO <- "G:/.shortcut-targets-by-id/18ukgK8yNyl2m8EsUhorWmPGgqYV81-hl/Ribo"
source(paste0(pathMS_nA, "/Ribo/Aleks_LungCancer/Code/LibsFunsAP.R"))

ev <- read.delim(paste0(pathMS,"/Collaborators/Zernicka_Goetz/dat/Experiment20_4cell.02.2024/01_DIANN_mtraq_NEUlib/report.tsv"))

#########----------------------------------------------------##############.
##### meta file load ####
#########----------------------------------------------------##############.
# load meta files
meta_embryo <- read_excel(paste0(pathMS,"/Collaborators/Zernicka_Goetz/dat/4CS_Sample_from_caltech_withembryo.xlsx"), sheet = 7) 
meta_embryo$PlateID_row <- paste0(meta_embryo$Plate, "_", meta_embryo$`Plate Row`)
meta <- read_excel(paste0(pathMS,"/Collaborators/Zernicka_Goetz/dat/4CS_Sample_from_caltech.xlsx"), sheet = 6) 
meta$injectionWellCaps <- toupper(meta$`Injection Well`)
meta$injectionWellCaps_label <- paste0(meta$injectionWellCaps, "_", meta$Label)



# connect blastomere information on whether a vegetal cell was identified 
meta$row <- toupper(substr(meta$`Cell Well`,1,1))
meta$PlateID_row <- paste0(meta$Plate, "_", meta$row)
meta$vegetal_ided <- meta_embryo$`Vegetal Celll`[match(meta$PlateID_row, meta_embryo$PlateID_row)]

# download alpha-beta Proteins
useList <- read.csv(paste0(pathMS, "/Collaborators/Zernicka_Goetz/code/CodeForGithub/2022_10_13_2cellstage_SignificantProteins_kmeansC__.csv"))

#########----------------------------------------------------##############.
##### isotopic correction #####
#########----------------------------------------------------##############.
pD_channel <- function(df){
  df <- df %>% dplyr::mutate("channel_name" = ifelse(grepl("-0|mTRAQ0", Modified.Sequence), "mTRAQ0",
                                                     ifelse(grepl("-4|mTRAQ4", Modified.Sequence), "mTRAQ4", "mTRAQ8")))
  return(df)
}

pD_seqcharge <- function(df){
  df$seqcharge <- paste0(df$Modified.Sequence, df$Precursor.Charge)
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ0\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ-K-0\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ-n-0\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ4\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ-K-4\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ-n-4\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ8\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ-K-8\\)")
  df$seqcharge <- str_remove_all(df$seqcharge, "\\(mTRAQ-n-8\\)")
  return(df)
}

df <- ev.fil
#test <- df %>% dplyr::filter(grepl("ETFIQNDXDAFENLLISK", Precursor.Id))
pD_isotopicCO_modified <- function(df,ppm=5,monoisotopic_only=T){
  cnames <- colnames(df)
  df <- pD_seqcharge(df)
  df$seqcharge_run <- paste0(df$seqcharge, df$Run)
  mTRAQ0 <- list(C=7, H=12, N=2, O=1)
  mTRAQ4 <- list(C=4, H=12, N=1, O=1)
  mTRAQ8 <- list(C=1, H=12, N=0, O=1)
  # count number of each mTRAQ label type per precursor:
  df$d0 <- str_count(df$Modified.Sequence, "mTRAQ0|mTRAQ-K-0|mTRAQ-n-0")
  df$d4 <- str_count(df$Modified.Sequence, "mTRAQ4|mTRAQ-K-4|mTRAQ-n-4")
  df$d8 <- str_count(df$Modified.Sequence, "mTRAQ8|mTRAQ-K-8|mTRAQ-n-8")
  uni_prec <- df %>% dplyr::distinct(Modified.Sequence, .keep_all=T)
  deamid <- list(N=-1, H=-1, O=1)
  ox <- list(O=1)
  carba <- list(C=2, H=3, N=1, O=1)
  precs <- vector(mode = "list", length = nrow(uni_prec))
  #i<-13305
  # get chemical formula for each peptide
  for(i in 1:nrow(uni_prec)){
    tempseq <- paste0(uni_prec$Stripped.Sequence[i])
    modseq <- paste0(uni_prec$Modified.Sequence[i])
    numLab_d0 <- uni_prec$d0[i]
    numLab_d4 <- uni_prec$d4[i]
    numLab_d8 <- uni_prec$d8[i]
    ## other modifications
    d <- str_count(uni_prec$Precursor.Id[i], "UniMod:7")
    o <- str_count(uni_prec$Precursor.Id[i], "UniMod:35")
    c <- str_count(uni_prec$Precursor.Id[i], "UniMod:4")
    ch <- as.numeric(paste0(uni_prec$Precursor.Charge[i]))
    el_temp <- ConvertPeptide(tempseq)
    el_temp$C <- el_temp$C+(mTRAQ0$C*numLab_d0)+(mTRAQ4$C*numLab_d4)+(mTRAQ8$C*numLab_d8)+(carba$C*c)
    el_temp$H <- el_temp$H+(mTRAQ0$H*numLab_d0)+(mTRAQ4$H*numLab_d4)+(mTRAQ8$H*numLab_d8)+(carba$H*c)
    el_temp$N <- el_temp$N+(mTRAQ0$N*numLab_d0)+(mTRAQ4$N*numLab_d4)+(mTRAQ8$N*numLab_d8)+(carba$N*c)
    el_temp$O <- el_temp$O+(mTRAQ0$O*numLab_d0)+(mTRAQ4$O*numLab_d4)+(mTRAQ8$O*numLab_d8)+(carba$O*c)
    el_temp_df <- data.frame(el_temp)
    el_temp_df$charge <- paste0(ch)
    rownames(el_temp_df) <- modseq
    precs[[i]] <- el_temp_df
  }
  precs_df <- do.call("rbind", precs)
  precs_df <- precs_df[!grepl("K",row.names(precs_df)),] #remove K peptides because they have an extra label which will shift them far enough away to not be problematic (isotopically)
  precs_df$charge <- as.numeric(precs_df$charge)
  #precs_df$form <- paste0("C",precs_df$C,"H",precs_df$H,"N",precs_df$N,"O",precs_df$O,"S",precs_df$S)
  precs_df <- precs_df %>% dplyr::mutate(form = ifelse(S==0, paste0("C",precs_df$C,"H",precs_df$H,"N",precs_df$N,"O",precs_df$O),
                                                       paste0("C",precs_df$C,"H",precs_df$H,"N",precs_df$N,"O",precs_df$O,"S",precs_df$S)))
  #precs_df <- precs_df %>% dplyr::mutate(form = paste0("C",precs_df$C))
  data(isotopes)
  iso <- isopattern(isotopes, precs_df$form, threshold=0.001, plotit=FALSE,
                    algo=1,rel_to=0,verbose=TRUE, return_iso_calc_amount=FALSE, charge=precs_df$charge)
  # iso_1 <- iso[1:5]
  # profiles<-envelope( iso_1, ppm=F, dmz="get", frac=1, env="Gaussian", resolution=1E6, plotit=F)
  # v <- vdetect(profiles,detect="centroid",plotit=F,verbose=F)
  iso_env_list <- list()
  for(i in 1:length(iso)){
    expected_mz <- as.numeric(iso[[i]][1,1])
    iso_temp <- data.frame(iso[[i]])
    charge <- as.numeric(precs_df$charge[i])
    expected_mzs <- c(expected_mz, expected_mz+(1.0034)/charge, expected_mz+(2*1.0034)/charge, expected_mz+(3*1.0034)/charge, expected_mz+(4*1.0034)/charge, expected_mz+(5*1.0034)/charge)
    keep <-data.frame("X1" = sum(iso_temp$abundance[which(iso_temp$'m.z'>(expected_mzs[1]-expected_mzs[1]*ppm*0.000001)&iso_temp$'m.z'<(expected_mzs[1]+expected_mzs[1]*ppm*0.000001))]),
                      "X2" = sum(iso_temp$abundance[which(iso_temp$'m.z'>(expected_mzs[2]-expected_mzs[2]*ppm*0.000001)&iso_temp$'m.z'<(expected_mzs[2]+expected_mzs[2]*ppm*0.000001))]),
                      "X3" = sum(iso_temp$abundance[which(iso_temp$'m.z'>(expected_mzs[3]-expected_mzs[3]*ppm*0.000001)&iso_temp$'m.z'<(expected_mzs[3]+expected_mzs[3]*ppm*0.000001))]),
                      "X4" = sum(iso_temp$abundance[which(iso_temp$'m.z'>(expected_mzs[4]-expected_mzs[4]*ppm*0.000001)&iso_temp$'m.z'<(expected_mzs[4]+expected_mzs[4]*ppm*0.000001))]),
                      "X5" = sum(iso_temp$abundance[which(iso_temp$'m.z'>(expected_mzs[5]-expected_mzs[5]*ppm*0.000001)&iso_temp$'m.z'<(expected_mzs[5]+expected_mzs[5]*ppm*0.000001))]),
                      "X6" = sum(iso_temp$abundance[which(iso_temp$'m.z'>(expected_mzs[6]-expected_mzs[6]*ppm*0.000001)&iso_temp$'m.z'<(expected_mzs[6]+expected_mzs[6]*ppm*0.000001))]))
    row.names(keep) <- row.names(precs_df)[i]
    iso_env_list[[i]] <- keep
  }
  iso_env_df <- do.call("rbind", iso_env_list)
  iso_env_df$modseq <- row.names(iso_env_df)
  # Now join with the quantitative data, and correct for isotopic carryover
  dat <- df %>% left_join(iso_env_df, by =c("Modified.Sequence" = "modseq"))
  dat$num_labs <- dat$d0+dat$d4+dat$d8
  dat <- pD_channel(dat)
  if(monoisotopic_only==T){
    dat0 <- dat %>% dplyr::filter(channel_name=="mTRAQ0") %>% dplyr::mutate("y0" = (X5)/(X1+X5))
    dat4 <- dat %>% dplyr::filter(channel_name=="mTRAQ4") %>% dplyr::mutate("y4" = (X5)/(X1+X5))
    dat8 <- dat %>% dplyr::filter(channel_name=="mTRAQ8") %>% dplyr::mutate("y8" = (X5)/(X1+X5))
    dat_fin <- rbind.fill(dat0,dat4,dat8)
  } else{
    dat0 <- dat %>% dplyr::filter(channel_name=="mTRAQ0") %>% dplyr::mutate("y0" = (X5+X6)/(X1+X2+X5+X6))
    dat4 <- dat %>% dplyr::filter(channel_name=="mTRAQ4") %>% dplyr::mutate("y4" = (X5+X6)/(X1+X2+X5+X6))
    dat8 <- dat %>% dplyr::filter(channel_name=="mTRAQ8") %>% dplyr::mutate("y8" = (X5+X6)/(X1+X2+X5+X6))
    dat_fin <- rbind.fill(dat0,dat4,dat8)
  }
  if(monoisotopic_only==T){
    dat <- dat %>% dplyr::mutate("y" = (X5)/(X1+X5))
  } else{
    dat <- dat %>% dplyr::mutate("y" = (X5+X6)/(X1+X2+X5+X6))
  }
  dat_d <- dcast(dat_fin, seqcharge_run~channel_name,value.var="Ms1.Area")
  dat_d[is.na(dat_d)] <- 0
  dat_d_y <- dcast(dat, seqcharge_run~channel_name,value.var="y")
  dat_d_y <- dat_d_y %>% dplyr::rename("y0"="mTRAQ0","y4"="mTRAQ4","y8"="mTRAQ8")
  all_dat <- dat_d %>% full_join(dat_d_y, by =c("seqcharge_run"="seqcharge_run"))
  all_dat[is.na(all_dat)]<-0
  all_dat$d0_iso <- all_dat$mTRAQ0+all_dat$mTRAQ0*all_dat$y0
  all_dat$d4_iso <- all_dat$mTRAQ4 - all_dat$mTRAQ0*all_dat$y0 + all_dat$mTRAQ4*all_dat$y4 #remove what was from d0, and add what bled into d8
  all_dat$d8_iso <- all_dat$mTRAQ8 - all_dat$mTRAQ4*all_dat$y4 + all_dat$mTRAQ8*all_dat$y8 #remove what was from d4, and add what bled into d12
  #only do the isotopic correction for precursors with < 2 labels.. otherwise, the difference in mass is enough to avoid impactful isotopic carryover (such is the case with K-containing peptides).
  all_dat$d0_iso <- ifelse(is.na(all_dat$d0_iso), all_dat$mTRAQ0, all_dat$d0_iso)
  all_dat$d4_iso <- ifelse(is.na(all_dat$d4_iso), all_dat$mTRAQ4, all_dat$d4_iso)
  all_dat$d8_iso <- ifelse(is.na(all_dat$d8_iso), all_dat$mTRAQ8, all_dat$d8_iso)
  all_dat[all_dat<0]<-0 #replace negatives with 0
  all_dat_d0 <- all_dat %>% dplyr::select("seqcharge_run", "d0_iso") %>% dplyr::rename("Ms1.Area_iso" = "d0_iso") %>% na.omit()
  all_dat_d4 <- all_dat %>% dplyr::select("seqcharge_run", "d4_iso")%>% dplyr::rename("Ms1.Area_iso" = "d4_iso") %>% na.omit()
  all_dat_d8 <- all_dat %>% dplyr::select("seqcharge_run", "d8_iso")%>% dplyr::rename("Ms1.Area_iso" = "d8_iso") %>% na.omit()
  d0_dat <- dat[grepl("mTRAQ0|mTRAQ-K-0|mTRAQ-n-0", dat$Modified.Sequence),]
  d4_dat <- dat[grepl("mTRAQ4|mTRAQ-K-4|mTRAQ-n-4", dat$Modified.Sequence),]
  d8_dat <- dat[grepl("mTRAQ8|mTRAQ-K-8|mTRAQ-n-8", dat$Modified.Sequence),]
  d0_dat <- d0_dat %>% left_join(all_dat_d0, by =c("seqcharge_run" = "seqcharge_run"))
  d4_dat <- d4_dat %>% left_join(all_dat_d4, by =c("seqcharge_run" = "seqcharge_run"))
  d8_dat <- d8_dat %>% left_join(all_dat_d8, by =c("seqcharge_run" = "seqcharge_run"))
  dat <- rbind(d0_dat, d4_dat, d8_dat)
  dat$Ms1.Area_iso[is.na(dat$Ms1.Area_iso)] <-0
  dat <- dat %>% dplyr::select(all_of(cnames), "Ms1.Area_iso")
  rm(df)
  return(dat)
}



#install.packages("OrgMassSpecR", repos="http://R-Forge.R-project.org")
library('OrgMassSpecR')

#install.packages('enviPat')
library('enviPat')

ev.fil <- ev %>% 
  #dplyr::filter(Lib.PG.Q.Value <= 0.01) %>% 
  #dplyr::filter(Ms1.Area >0 )  %>% 
  dplyr::filter(!grepl("X", Precursor.Id))
test <- ev.fil %>% dplyr::filter(grepl("ETFIQNDXDAFENLLISK", Precursor.Id))

isoCorrected <- pD_isotopicCO_modified(ev.fil)

#########----------------------------------------------------##############.
##### filters on Q value #####
#########----------------------------------------------------##############.
#CHOOSE - non-isotopically corrected 
x.fil <- ev %>% 
  dplyr::filter(Lib.Q.Value <= 0.01 & Lib.PG.Q.Value <= 0.01 & PG.Q.Value <= 0.01) %>% 
  dplyr::filter(Ms1.Area >0 ) 

# or isotopically corrected data
x.fil <- isoCorrected %>% 
  dplyr::filter(Ms1.Area_iso >0 ) %>%
  dplyr::filter(Lib.PG.Q.Value <= 0.05)
  

# label the unique precursors
x.fil$label <- sub( ".*\\(","", (sub(").*", "", x.fil$Modified.Sequence)))
x.fil$Run_Label <- paste0(x.fil$Run, "_", x.fil$label )
unique(x.fil$label)[1]
# shorten the label to d0, d4, d8
x.fil$label_short <- paste0("d", sub(".*-", "", x.fil$label))
unique(x.fil$label_short)
# paste together the Run and the mtraq label used 
x.fil$Run_Label <- paste0(x.fil$Run, "_", x.fil$label_short )
unique(x.fil$Run_Label)

x.fil$SeqCharge <- paste0(x.fil$Stripped.Sequence, x.fil$Precursor.Charge)

x.fil$InjectionWell <- gsub("_.*", "", gsub("_.*","",gsub(".*-","",x.fil$Run)) )
test <- x.fil %>% dplyr::filter(!(InjectionWell %in% meta$injectionWellCaps))
unique(test$InjectionWell)

#########----------------------------------------------------##############.
##### more filters: MS1 Area ####
#########----------------------------------------------------##############.

# can filter here for a particular plate or something like this
x.fil1 <- x.fil %>% 
  dplyr::filter(InjectionWell %in% meta$injectionWellCaps)
length(unique(x.fil$InjectionWell))
length(unique(x.fil1$InjectionWell)) 
x.fil1$injectionwell_Label <- paste0(x.fil1$InjectionWell, "_", x.fil1$label_short)
x.fil1$CellType <- meta$type[match(x.fil1$injectionwell_Label, meta$injectionWellCaps_label)]
unique(x.fil1$CellType)
test <- x.fil1 %>% dplyr::filter(is.na(CellType))
unique(test$InjectionWell) # J8 set did not have d8 labeled sample, which is why we have NAs for CellType

# see the amount of log10 signal of MS1 level for each label 
ms1sig <- x.fil1 %>% dplyr::group_by(injectionwell_Label, CellType) %>% 
  dplyr::summarise(medms1 = median(Ms1.Area_iso))
  #dplyr::summarise(medms1 = median(Ms1.Area))
ms1sig$log10_medianMS1 <- log10(ms1sig$medms1)
ms1sig$embryo <- meta$Embryo[match(ms1sig$injectionwell_Label, meta$injectionWellCaps_label)]

threshold_ms1 <- 4
ggplot(ms1sig, aes(x = log10_medianMS1, fill = CellType )) + 
  geom_histogram(binwidth = 0.01, color = 'black',  alpha = 0.5) + 
  geom_vline(xintercept = threshold_ms1, linetype = 'dotted', size = 1.5, color = 'blue')

ms1sig$CellType_Vegetal <- meta$Vegetal[match(ms1sig$injectionwell_Label, meta$injectionWellCaps_label)]

# based on this metric, we are ready to filter out cells that did not work 
ms1sig.fil <- ms1sig %>% dplyr::filter(log10_medianMS1 > threshold_ms1) %>% dplyr::filter(CellType == 'cell')
ms1sig.cell <- ms1sig %>% dplyr::filter(CellType == 'cell') 
length(unique(ms1sig.fil$injectionwell_Label))
length(unique(ms1sig.cell$injectionwell_Label))
# 76 total blastomeres collected. With the median MS1 filter, we lose 10 blastomeres, leaving us with 66
# we are able to keep 11/12 of the vegetal cells with this filter 

#########----------------------------------------------------##############.
##### more filters: Precursors Amount ####
#########----------------------------------------------------##############.

protein.count <- x.fil1 %>% dplyr::group_by(injectionwell_Label, CellType)%>% 
  dplyr::summarise(n_Genes = length(unique(Genes)),  n_PG = length(unique(Protein.Group)),
                   n_Precursors = length(unique(SeqCharge)))
protein.count$embryo <- meta$Embryo[match(protein.count$injectionwell_Label, meta$injectionWellCaps_label)]
protein.count$CellType_Vegetal <- meta$Vegetal[match(protein.count$injectionwell_Label, meta$injectionWellCaps_label)]

ggplot(protein.count, aes(x = n_PG, fill = CellType )) + 
  geom_histogram(color = 'black',  alpha = 0.5)
protein.count.cell <- protein.count %>% dplyr::filter(CellType == 'cell') 

ggplot(protein.count.cell, aes(x = CellType_Vegetal, y = n_PG)) + 
  geom_boxplot() + geom_point()
ggplot(protein.count.cell, aes(x = CellType_Vegetal, y = n_Precursors)) + 
  geom_boxplot() + geom_point()


protein.count.fil <- protein.count %>% dplyr::filter(CellType == 'cell' & n_PG > 400 & n_Precursors > 2500)

length(unique(protein.count.fil$injectionwell_Label))
length(unique(protein.count.cell$injectionwell_Label))
# 76 total blastomeres collected. With the protein/precursor cutoff, we lose 8 blastomeres, leaving us with 68
# we are able to keep 10/12 of the vegetal cells with this filter 


#########----------------------------------------------------##############.
##### more filters: Combine MS1 + Precursors Amount ####
#########----------------------------------------------------##############.
meta_cells <- meta %>% dplyr::filter(type == 'cell')
meta_cells$log10_medianMS1 <- ms1sig$log10_medianMS1[match(meta_cells$injectionWellCaps_label, ms1sig$injectionwell_Label)]
meta_cells$n_Precursors <- protein.count$n_Precursors[match(meta_cells$injectionWellCaps_label, ms1sig$injectionwell_Label)]
meta_cells$n_PG <- protein.count$n_PG[match(meta_cells$injectionWellCaps_label, ms1sig$injectionwell_Label)]

threshold_PG <- 400
threshold_Precursors <- 2000
meta_cells.fil <- meta_cells %>% 
  dplyr::filter(n_PG > threshold_PG & n_Precursors > threshold_Precursors & log10_medianMS1 > threshold_ms1)
# out of 76 cells, we are able to keep 59 blastomeres
# this includes 9/12 vegetal cells 

# let's go back and look at vegeatl cells
meta_cells.vegetal <- meta_cells %>% dplyr::filter(Vegetal == "V")
ggplot(data = meta_cells.vegetal, aes(x = n_PG, y = log10_medianMS1, color = Vegetal)) + geom_point() + 
  geom_vline(xintercept = threshold_PG, linetype = 'dotted', size = 1) + 
  geom_hline(yintercept = threshold_ms1,, linetype = 'dotted', size = 1)
ggplot(data = meta_cells.vegetal, aes(x = n_Precursors, y = log10_medianMS1, color = Vegetal)) + geom_point() + 
  geom_vline(xintercept = threshold_Precursors, linetype = 'dotted', size = 1) + 
  geom_hline(yintercept = threshold_ms1, linetype = 'dotted', size = 1)


#########----------------------------------------------------##############.
##### more filters: all 4 cells in embryo ####
#########----------------------------------------------------##############.

#choose what you want to filter on
uchoose <- meta_cells %>% 
  dplyr::filter(n_PG > threshold_PG & n_Precursors > threshold_Precursors & log10_medianMS1 > threshold_ms1)
length(unique(uchoose$Embryo))
length(unique(meta_cells$Embryo))
embryonot <- meta_cells %>% dplyr::filter(!(Embryo %in% uchoose$Embryo))

# count the number of cells kept per embryo 
uchoose.count <- uchoose %>% 
  dplyr::group_by(Embryo) %>% 
  dplyr::summarise(n_cells = length(unique(injectionWellCaps_label)))


wholeEmbryos <- uchoose.count %>% dplyr::filter(n_cells == 4)
uchooseWhole <- uchoose %>% dplyr::filter(Embryo %in% wholeEmbryos$Embryo)
length(unique(uchooseWhole$Embryo))

#########----------------------------------------------------##############.
##### filter the data ####
#########----------------------------------------------------##############.


x.fil1$embryo <- meta$Embryo[match(x.fil1$injectionwell_Label, meta$injectionWellCaps_label)]
# filter to have all successfully quantified blastomeres, regardless of whole embryo set
x.fil2.all <- x.fil1 %>% dplyr::filter(injectionwell_Label %in% uchoose$injectionWellCaps_label)

# filter to have all successfully quantified blastomeres part of whole embryo set
x.fil2.whole <- x.fil1 %>% dplyr::filter(injectionwell_Label %in% uchooseWhole$injectionWellCaps_label)


#########----------------------------------------------------##############.
##### normalize cells ####
#########----------------------------------------------------##############.

# CHOOSE
useDF <- x.fil2.whole 
#filter on amount of peptides seen? 
chooseProtein <- 'Y'

test <- useDF %>% dplyr::filter(grepl(";", Protein.Group))
unique(test$Protein.Group)


# quantify the number of times the precursor has been quantified per plexDIA run
precursorcount <- useDF %>% 
  dplyr::group_by(SeqCharge) %>% 
  dplyr::summarise(numRuns = length(unique(Run)))
ggplot(data = precursorcount, aes(x = numRuns)) + geom_histogram()
# filter for precursors in at least 50% of the runs
precursorcount.fil <- precursorcount %>% 
  dplyr::filter(numRuns >= length(unique(useDF$Run))/2)

# filter the data for only these precursors
x.fil50 <- useDF %>% 
  dplyr::filter(SeqCharge %in% precursorcount.fil$SeqCharge)


# Normalize for sample loading 
x.fil3 <- useDF %>% 
  dplyr::group_by(injectionwell_Label) %>% 
  dplyr::mutate(q1 = Ms1.Area_iso / median(Ms1.Area_iso, na.rm = T))
  #dplyr::mutate(q1 = Ms1.Area / median(Ms1.Area, na.rm = T))
ggplot(data = x.fil3, aes(x = injectionwell_Label, y = log2(q1))) + geom_boxplot()
test <- x.fil3 %>% dplyr::filter(grepl("E9Q414",Protein.Group))

# Normalize each peptide by the mean across all single cells 
x.fil4 <- x.fil3 %>% 
  dplyr::group_by(SeqCharge) %>% 
  dplyr::mutate(q2 = q1 / mean(q1, na.rm = T))

# collapse peptide by median and then normalize across all blastomeres
x.fil5 <- x.fil4 %>% 
  dplyr::group_by(injectionwell_Label, Protein.Group, Genes, embryo) %>% 
  dplyr::summarise(medPep = median(q2, na.rm = T))



# quantify the number of times the precursor has been quantified per plexDIA run
proteincountRUN <- x.fil5 %>% 
  dplyr::group_by(Protein.Group) %>% 
  dplyr::summarise(numCells = length(unique(injectionwell_Label)))
ggplot(proteincountRUN, aes(numCells)) + geom_histogram()
# filter for precursors in at least x% of the cells
proteincountRUN.fil <- proteincountRUN %>% 
  dplyr::filter(numCells >= 0.1*(length(unique(x.fil5$injectionwell_Label))))

if(chooseProtein == "Y") { 
  x.fil5.protein <- x.fil5 %>% dplyr::filter(Protein.Group %in% proteincountRUN.fil$Protein.Group)
}

if(chooseProtein == "N") { 
  x.fil5.protein <- x.fil5 
}



x.fil6 <- x.fil5.protein %>% 
  dplyr::group_by(injectionwell_Label) %>% 
  dplyr::mutate(q3 = medPep / median(medPep, na.rm = T))
ggplot(data = x.fil6, aes(x = injectionwell_Label, y = log2(q3))) + geom_boxplot()

x.fil.across <- x.fil6 %>% 
  dplyr::group_by(Protein.Group) %>% 
  dplyr::mutate(q4 = q3 / mean(q3, na.rm = T))
length(unique(x.fil5.protein$Protein.Group))


x.fil.within <- x.fil6 %>% 
  dplyr::group_by(embryo, Protein.Group) %>% 
  dplyr::mutate(q4 = q3 / mean(q3, na.rm = T))
length(unique(x.fil5.protein$Protein.Group))
length(unique(x.fil.across$injectionwell_Label))


test <- x.fil.within %>% dplyr::filter(embryo ==12)
unique(test$injectionwell_Label)
test.d <- reshape2::dcast(test, Protein.Group ~ injectionwell_Label, value.var = 'q4')

x.fil.within.noOnes <- x.fil.within %>% dplyr::filter(q4 !=1)
test <- x.fil.within.noOnes %>% dplyr::filter(embryo == unique(x.fil.within.noOnes$embryo)[6])
unique(test$injectionwell_Label)
test.d <- reshape2::dcast(test, Protein.Group ~ injectionwell_Label, value.var = 'q4')
test.d3 <- reshape2::dcast(test, Protein.Group ~ injectionwell_Label, value.var = 'q3')
test.d3$rowMean <- rowMeans(test.d3[,2:5], na.rm = T)




#########----------------------------------------------------##############.
##### imputation per embryo + Combat ####
#########----------------------------------------------------##############.
x.use <- x.fil.within.noOnes
i<-2
k.t <-2
for(i in 1:length(unique(x.use$embryo))) { 
  temp.embryo <- x.use %>% dplyr::filter(embryo == unique(x.use$embryo)[i] )
  temp.embryo.d <- reshape2::dcast(temp.embryo, Protein.Group ~ injectionwell_Label, value.var = 'q4')
  temp.embryo.d.log2 <- temp.embryo.d
  temp.embryo.d.log2[,2:5] <- log2(temp.embryo.d.log2[,2:5])
  rownames(temp.embryo.d.log2) <- temp.embryo.d.log2$Protein.Group
  temp.embryo.d.log2 <- temp.embryo.d.log2 %>% dplyr::select(-Protein.Group)
  temp.sc.imp <- as.data.frame(hknn(as.matrix(temp.embryo.d.log2), k.t))
  temp.sc.imp$Protein.Group <- rownames(temp.sc.imp)
  temp.sc.imp.m <- reshape2::melt(temp.sc.imp)
  temp.sc.imp.m$Genes <- temp.embryo$Genes[match(temp.sc.imp.m$Protein.Group, temp.embryo$Protein.Group)]
  temp.sc.imp.m$embryo <- unique(x.use$embryo)[i]
  if(i == 1) { 
    embyro.imp <- temp.sc.imp.m
  }
  if(i !=1) { 
    embyro.imp <- rbind(embyro.imp, temp.sc.imp.m)
    }
  
}

colnames(embyro.imp) <- c("Protein.Group", "injectionwell_Label", "q4", "Genes", "embryo")


x.combat <- embyro.imp
x.combat.d <- reshape2::dcast(x.combat, Protein.Group ~ injectionwell_Label, value.var = 'q4')
rownames(x.combat.d) <- x.combat.d$Protein.Group
x.combat.d <- x.combat.d %>% dplyr::select(-Protein.Group)

batch.cov <- data.frame(unique(colnames(x.combat.d)))
batch.cov$label <- gsub(".*_", "", batch.cov$unique.colnames.x.combat.d..)
batch.cov <- batch.cov %>% dplyr::mutate(label_batch = case_when(label == "d0" ~ 1, 
                                                                 label == "d4" ~ 2, 
                                                                 label == "d8" ~ 3))
batch.cov.c <- batch.cov$label_batch

## Impute single celldata
imp.input<-as.matrix((x.combat.d))
k.t <- 3
sc.imp <- hknn(imp.input, k.t)

Heatmap(sc.imp)

# Complete data
sc.imp.complete <- na.omit((x.combat.d))


# CHOOSE what do combat on 
tobeCOMBAT <- sc.imp.complete
library(sva)
x.combat.d.c <- ComBat(dat = tobeCOMBAT, batch = batch.cov.c)



#########----------------------------------------------------##############.
##### alpha beta with COMBAT ####
#########----------------------------------------------------##############.
library(circlize)
colfun = colorRamp2(c(-1,0,1), c("blue", "white", "red"))
# CHOOSE
useDF.ab <- reshape2::melt((x.combat.d.c))
colnames(useDF.ab) <- c("Protein.Group", "injectionwell_Label", "q4")


ab <- useDF.ab %>% dplyr::filter(Protein.Group %in% useList$Protein)
length(unique(ab$Protein.Group))
length(unique(ab$injectionwell_Label))
ab$proteinType <- useList$GreaterIN[match(ab$Protein.Group, useList$Protein)]
# plot distributions of each blastomere 
ggplot(data = ab, aes(x = injectionwell_Label, color = proteinType, y = q4)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust =1)) +
  coord_cartesian(ylim = c(-1.2, 1.2))




sigprotm.1 <- ab %>% 
  dplyr::group_by(injectionwell_Label, proteinType) %>% 
  dplyr::summarise(medianProtein = median(q4,na.rm = T))
# obtain the number of proteins quantified 
sigprotm.1a <- ab %>% 
  dplyr::group_by(injectionwell_Label, proteinType) %>% 
  dplyr::summarise( numProteins = sum(!is.na(q4)))
ggplot(data = sigprotm.1a, aes(x = numProteins, fill = proteinType)) + geom_histogram() + theme_bw() + 
  theme(text = element_text(size = 18), legend.position = "bottom")
# obtain the fold-change between the alpha and beta protein medians 
# what was originally used: 
sigprotm.2 <- sigprotm.1 %>% 
  dplyr::group_by(injectionwell_Label) %>% 
  dplyr::summarise( Alpha_Beta_ = medianProtein[proteinType == "Cluster2"] - medianProtein[proteinType == "Cluster1"])
ggplot(data = sigprotm.2, aes(x = Alpha_Beta_)) + geom_histogram() + geom_density() + 
  theme_classic() + 
  theme(text = element_text(size =25))


#connect which cell is a vegetal cell 
sigprotm.2$CellType <- meta$Vegetal[match(sigprotm.2$injectionwell_Label, meta$injectionWellCaps_label)]

sigprotm.2$vegetalcell_ided <- meta$vegetal_ided[match(sigprotm.2$injectionwell_Label, meta$injectionWellCaps_label)]
sigprotm.2.v <- sigprotm.2 %>% dplyr::filter(vegetalcell_ided == 'yes')
sigprotm.2.v$embryo <- meta$Embryo[match(sigprotm.2.v$injectionwell_Label, meta$injectionWellCaps_label)]
length(unique(sigprotm.2.v$injectionwell_Label))


ggplot(data = sigprotm.2.v, aes(x = CellType, y = Alpha_Beta_)) + geom_boxplot() + geom_point() + theme_bw() +
  theme(text = element_text(size = 20)) + 
  labs(x = "Vegetal?", y = "Alpha - Beta Fold Change")

sigprotm.2.v.fil.v <- sigprotm.2.v %>% dplyr::filter(CellType=="V")
sigprotm.2.v.fil.nv <- sigprotm.2.v %>% dplyr::filter(CellType=="x")
t.test(sigprotm.2.v.fil.v$Alpha_Beta_, sigprotm.2.v.fil.nv$Alpha_Beta_)





threshold_ab <- 0.02

sigprotm.2 <- sigprotm.2 %>% dplyr::mutate(alpha_beta_character = case_when(Alpha_Beta_ > threshold_ab ~ "Alpha", 
                                                                            Alpha_Beta_ < -threshold_ab ~ "Beta", 
                                                                     TRUE ~ "Ill-Defined"))
# obtain the median value of each protein for each cell type, disregarding Ill-Definded + vegetal
ab$alpha_beta_character <- sigprotm.2$alpha_beta_character[match(ab$injectionwell_Label, 
                                                                 sigprotm.2$injectionwell_Label)]
ab$vegetal <- meta_cells$Vegetal[match(ab$injectionwell_Label, meta_cells$injectionWellCaps_label)]
ab$vegetal_ided <- meta_cells$vegetal_ided[match(ab$injectionwell_Label, meta_cells$injectionWellCaps_label)]

ab.fil <- ab %>% dplyr::filter(vegetal_ided == "yes") %>% 
  dplyr::filter(alpha_beta_character != "Ill-Defined") %>% 
  dplyr::filter(vegetal != "V")

# calculate median protein level
ab.fil.median <- ab.fil %>% dplyr::group_by(Protein.Group) %>% 
  dplyr::summarise(ALPHA = median(q4[alpha_beta_character == "Alpha"]), 
                   BETA = median(q4[alpha_beta_character == "Beta"]))
ab.fil.median$FC_alpha_beta <- ab.fil.median$ALPHA - ab.fil.median$BETA

threshold_FC <- 0.1
ab.fil.median.fil <- ab.fil.median %>% dplyr::filter(abs(FC_alpha_beta) >threshold_FC)


# filter on the highest FC proteins and then recalculate alpha-beta polarization
ab.fillerup <- ab %>% 
  dplyr::filter(Protein.Group %in% ab.fil.median.fil$Protein.Group) #%>% 
  #dplyr::filter(alpha_beta_character == "Alpha" | alpha_beta_character == "Beta" | vegetal == "V")
length(unique(ab$injectionwell_Label))
length(unique(ab.fillerup$injectionwell_Label))
sigprotm.1.filler <- ab.fillerup %>% 
  dplyr::group_by(injectionwell_Label, proteinType) %>% 
  dplyr::summarise(medianProtein = median(q4,na.rm = T))
sigprotm.2.filler <- sigprotm.1.filler %>% 
  dplyr::group_by(injectionwell_Label) %>% 
  dplyr::summarise( Alpha_Beta_ = medianProtein[proteinType == "Cluster2"] - medianProtein[proteinType == "Cluster1_"])








threshold_ab2 <- 0.05
sigprotm.2.filler <- sigprotm.2.filler %>% dplyr::mutate(alpha_beta_character = 
                                                           case_when(Alpha_Beta_ > threshold_ab2 ~ "Alpha", 
                                                                            Alpha_Beta_ < -threshold_ab2 ~ "Beta", 
                                                                            TRUE ~ "Ill-Defined"))
length(unique(sigprotm.2.filler$injectionwell_Label))
ab$alpha_beta_character2 <- sigprotm.2.filler$alpha_beta_character[match(ab$injectionwell_Label, 
                                                                         sigprotm.2.filler$injectionwell_Label)]


test <-  ab %>% 
  dplyr::filter(Protein.Group %in% ab.fil.median.fil$Protein.Group)
unique(test$injectionwell_Label)
ab.fillerup2 <- ab %>% 
  dplyr::filter(Protein.Group %in% ab.fil.median.fil$Protein.Group) %>% 
  dplyr::filter(alpha_beta_character == "Alpha" | alpha_beta_character == "Beta" | vegetal == "V") 
length(unique(ab.fillerup2$injectionwell_Label))
cc <- reshape2::dcast(ab.fillerup2, Protein.Group ~ injectionwell_Label, value.var = 'q4')
rownames(cc) <- cc$Protein.Group
cc <- cc[,-1]

cc.cor <- cor(cc)

column_anno <- data.frame(colnames(cc.cor))
column_anno$VegetalCell <- meta$Vegetal[match(column_anno$colnames.cc.cor., 
                                              meta$injectionWellCaps_label)]

column_anno$VegetalCell_ID <- meta$vegetal_ided[match(column_anno$colnames.cc.cor., 
                                                      meta$injectionWellCaps_label)]
column_anno$alpha_beata2 <- sigprotm.2.filler$alpha_beta_character[match(column_anno$colnames.cc.cor., 
                                                                         sigprotm.2.filler$injectionwell_Label)]
column_anno$alpha_beata <- sigprotm.2$alpha_beta_character[match(column_anno$colnames.cc.cor., 
                                                                        sigprotm.2$injectionwell_Label)]






# reorder the heatmap  
colum_anno_reorder <- column_anno
colum_anno_reorder$ab_val <- sigprotm.2.filler$Alpha_Beta_[match(colum_anno_reorder$colnames.cc.cor., 
                                                                          sigprotm.2.filler$injectionwell_Label)]
colum_anno_reorder$ab_val <- sigprotm.2$Alpha_Beta_[match(colum_anno_reorder$colnames.cc.cor., 
                                                                          sigprotm.2$injectionwell_Label)]
colum_anno_reorder <- colum_anno_reorder[order(colum_anno_reorder$ab_val),]

colum_anno_reorder <- colum_anno_reorder %>% 
  dplyr::mutate(alpha_betaMINUSv = case_when(alpha_beata == "Alpha" & VegetalCell == "x" ~ "Alpha", 
                                             alpha_beata == "Beta" & VegetalCell == "x" ~ "Beta",
                                             VegetalCell == "V" ~ "Vegetal"))


colum_anno_reorder <- colum_anno_reorder %>% 
  dplyr::mutate(alpha_betaMINUSv = case_when(ab_val > 0  ~ "Alpha", 
                                             ab_val < 0  ~ "Beta"))


cc.cor2 <- cc.cor
cc.cor2 <- cc.cor2[colum_anno_reorder$colnames.cc.cor.,]
cc.cor3 <- cc.cor2[,colum_anno_reorder$colnames.cc.cor.]
column_ha <- HeatmapAnnotation(Vegetal = colum_anno_reorder$VegetalCell, 
                               #VegetalEmbryo = colum_anno_reorder$VegetalCell_ID, 
                               `Alpha or Beta` = colum_anno_reorder$alpha_betaMINUSv,
                               col = list(Vegetal = c("V" = "cyan", "x" = "bisque"), 
                                          #VegetalEmbryo = c("no" = "bisque", "yes" = "cadetblue1"), 
                                          `Alpha or Beta` = c(Alpha = "#D60093", Beta ="#E8AC04", Vegetal = 'grey')))

colfun = colorRamp2(c(-1,0,1), c("purple", "white", "yellow"), space = 'LAB')
colfun(seq(-3,3))
htlist <- Heatmap(cc.cor3, show_row_names = FALSE, show_column_names = FALSE,col = colfun,
                  column_order = colum_anno_reorder$colnames.cc.cor.,
                  row_order = colum_anno_reorder$colnames.cc.cor.,
                  bottom_annotation = column_ha, 
                  heatmap_legend_param = list(direction = 'vertical', title = "Correlation"))
draw(htlist, annotation_legend_side = 'bottom', heatmap_legend_side = 'right')









