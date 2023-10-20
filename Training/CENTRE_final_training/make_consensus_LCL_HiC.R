##########################################################################################
#
# make_consensus_LcL_HiC: add the HiC features to GM12878 datasets
#
##########################################################################################

#Add the HiC features to go into training

#GM BENGI datasets
CHiC <- read.table("/project/CRUP_scores/toSara/BENGI_processed_datasets/GM12878.CHiC-Benchmark.v38.txt", header = T)

CTCF <- read.table("/project/CRUP_scores/toSara/BENGI_processed_datasets/GM12878.CTCF-ChIAPET-Benchmark.v38.txt", header = T)

GEUVADIS<- read.table("/project/CRUP_scores/toSara/BENGI_processed_datasets/GM12878.GEUVADIS-Benchmark.v38.txt", header = T)

RNAPII <- read.table("/project/CRUP_scores/toSara/BENGI_processed_datasets/GM12878.RNAPII-ChIAPET-Benchmark.v38.txt", header = T)

HiC <- read.table("/project/CRUP_scores/toSara/BENGI_processed_datasets/GM12878_HiC-Benchmark.v38.txt", header = T)

#MinIns datasets
CHiC_MI <- read.table("/project/CRUP_scores/CENTRE_HiC/WilcoxonTestResults/MinInsulationScore/minInsulationWilcoxtestGM12878.CHiC.csv", 
                      header = T,
                      sep = ",")

CTCF_MI <- read.table("/project/CRUP_scores/CENTRE_HiC/WilcoxonTestResults/MinInsulationScore/minInsulationWilcoxtestGM12878.CTCF-ChIAPET.csv",
                      header = T,
                      sep = ",")

GEUVADIS_MI <- read.table("/project/CRUP_scores/CENTRE_HiC/WilcoxonTestResults/MinInsulationScore/minInsulationWilcoxtestGM12878.GEUVADIS.csv", 
                          header = T,
                          sep = ",")

RNAPII_MI <- read.table("/project/CRUP_scores/CENTRE_HiC/WilcoxonTestResults/MinInsulationScore/minInsulationWilcoxtestGM12878.RNAPII-ChIAPET.csv", 
                        header = T,
                        sep = ",")

HiC_MI <- read.table("/project/CRUP_scores/CENTRE_HiC/WilcoxonTestResults/MinInsulationScore/minInsulationWilcoxtestGM12878.HiC.csv", 
                     header = T,
                     sep = ",")

#MeanSwitch Intensity

CHiC_MSI <- read.table("/project/CRUP_scores/CENTRE_HiC/WilcoxonTestResults/MeanSwitchIntensity/meanSwitchIntensityWilcoxtestGM12878.CHiC.csv", header = T)

CTCF_MSI <- read.table("/project/CRUP_scores/CENTRE_HiC/WilcoxonTestResults/MeanSwitchIntensity/meanSwitchIntensityWilcoxtestGM12878.CTCF-ChIAPET.csv", header = T)

GEUVADIS_MSI <- read.table("/project/CRUP_scores/CENTRE_HiC/WilcoxonTestResults/MeanSwitchIntensity/meanSwitchIntensityWilcoxtestGM12878.GEUVADIS.csv", header = T)
  
RNAPII_MSI <- read.table("/project/CRUP_scores/CENTRE_HiC/WilcoxonTestResults/MeanSwitchIntensity/meanSwitchIntensityWilcoxtestGM12878.RNAPII-ChIAPET.csv", header = T)
 
HiC_MSI <- read.table("/project/CRUP_scores/CENTRE_HiC/WilcoxonTestResults/MeanSwitchIntensity/meanSwitchIntensityWilcoxtestGM12878.HiC.csv", header = T)
 
consensus <- read.table("/project/CRUP_scores/CENTRE_HiC/Training/CENTRE_final_training/consensus_LCLs.txt", header = T)
#### The merge is first with RNAPII then CTCF then CHiC then CTCF then HiC

##RNAPII
RNAPII.consensus.subset <- subset(consensus, consensus$pair %in% RNAPII$pair)
RNAPII_MI$pair <- paste(RNAPII_MI$symbol38, RNAPII_MI$gene_id1, sep = "_")
RNAPII_MI <- RNAPII_MI[c("pair", "score_min")]
colnames(RNAPII_MI) <- c("pair", "min_insulation")
RNAPII.consensus.MI <- merge(RNAPII.consensus.subset, RNAPII_MI,
                             by.x = "pair", by.y = "pair", all.x = T)

RNAPII.consensus.MI.MSI <- merge(RNAPII.consensus.MI, RNAPII_MSI[c("pair", "mean_switch_intensity")],
                             by.x = "pair", by.y = "pair", all.x = T)

##CTCF
consensus.subset <- subset(consensus, !(consensus$pair %in% RNAPII$pair))
CTCF.consensus.subset <- subset(consensus.subset, consensus.subset$pair %in% CTCF$pair)
CTCF_MI$pair <- paste(CTCF_MI$symbol38, CTCF_MI$gene_id1, sep = "_")
CTCF_MI <- CTCF_MI[c("pair", "score_min")]
colnames(CTCF_MI) <- c("pair", "min_insulation")
CTCF.consensus.MI <- merge(CTCF.consensus.subset, CTCF_MI,
                           by.x = "pair", by.y = "pair", all.x = T)

CTCF.consensus.MI.MSI <- merge(CTCF.consensus.MI, CTCF_MSI[c("pair", "mean_switch_intensity")],
                               by.x = "pair", by.y = "pair", all.x = T)
##cHIC
consensus.remainder <- subset(consensus.subset, !(consensus.subset$pair %in% CTCF$pair)) 
CHiC.consensus.subset <- subset(consensus.remainder, consensus.remainder$pair %in% CHiC$pair)
CHiC_MI$pair <- paste(CHiC_MI$symbol38, CHiC_MI$gene_id1, sep = "_")
CHiC_MI <- CHiC_MI[c("pair", "score_min")]
colnames(CHiC_MI) <- c("pair", "min_insulation")
CHiC.consensus.MI <- merge(CHiC.consensus.subset, CHiC_MI,
                             by.x = "pair", by.y = "pair", all.x = T)

CHiC.consensus.MI.MSI <- merge(CHiC.consensus.MI, CHiC_MSI[c("pair", "mean_switch_intensity")],
                                 by.x = "pair", by.y = "pair", all.x = T)


##now Hic 
consensus.rest <- subset(consensus.remainder, !(consensus.remainder$pair %in% CHiC$pair))
HiC.consensus.subset <- subset(consensus.rest, consensus.rest$pair %in% HiC$pair)
HiC_MI$pair <- paste(HiC_MI$symbol38, HiC_MI$gene_id1, sep = "_")
head(HiC_MI)
HiC_MI <- HiC_MI[c("pair", "score_min")]
colnames(HiC_MI) <- c("pair", "min_insulation")
HiC.consensus.MI <- merge(HiC.consensus.subset, HiC_MI,
                           by.x = "pair", by.y = "pair", all.x = T)

HiC.consensus.MI.MSI <- merge(HiC.consensus.MI, HiC_MSI[c("pair", "mean_switch_intensity")],
                               by.x = "pair", by.y = "pair", all.x = T)

##GEUVADIS
consensus.rest2 <- subset(consensus.rest, !(consensus.rest$pair %in% HiC$pair))
GEUVADIS.consensus.subset <- subset(consensus.rest2, consensus.rest2$pair %in% GEUVADIS$pair)
GEUVADIS_MI$pair <- paste(GEUVADIS_MI$symbol38, GEUVADIS_MI$gene_id1, sep = "_")
head(GEUVADIS_MI)
GEUVADIS_MI <- GEUVADIS_MI[c("pair", "score_min")]
colnames(GEUVADIS_MI) <- c("pair", "min_insulation")
GEUVADIS.consensus.MI <- merge(GEUVADIS.consensus.subset, GEUVADIS_MI,
                          by.x = "pair", by.y = "pair", all.x = T)

GEUVADIS.consensus.MI.MSI <- merge(GEUVADIS.consensus.MI, GEUVADIS_MSI[c("pair", "mean_switch_intensity")],
                              by.x = "pair", by.y = "pair", all.x = T)
nrow(CTCF.consensus.MI.MSI) + nrow(RNAPII.consensus.MI.MSI) + nrow(CHiC.consensus.MI.MSI)+nrow(HiC.consensus.MI.MSI)+ nrow(GEUVADIS.consensus.MI.MSI)

consensus.HIC <- rbind(CTCF.consensus.MI.MSI, RNAPII.consensus.MI.MSI, CHiC.consensus.MI.MSI,
                       HiC.consensus.MI.MSI, GEUVADIS.consensus.MI.MSI)

consensusLcL <- merge(consensus, consensus.HIC[c("pair","min_insulation","mean_switch_intensity")],
                      by.x = "pair", by.y = "pair", all.x = T)
head(consensus)

RNAPII.gm.MI <- merge(RNAPII, RNAPII_MI,
                               by.x = "pair", by.y = "pair", all.x = T)

RNAPII.GM.MI.MSI <- merge(RNAPII.gm.MI, RNAPII_MSI[c("pair", "mean_switch_intensity")],
                                   by.x = "pair", by.y = "pair", all.x = T)


min_ins <- RNAPII.GM.MI.MSI$min_insulation[!is.infinite(RNAPII.GM.MI.MSI$min_insulation)]
min_ins <- min_ins[!is.na(min_ins)]
max(min_ins)
RNAPII.GM.MI.MSI[mapply(is.infinite, RNAPII.GM.MI.MSI)] <- 6.180041
write.csv(RNAPII.GM.MI.MSI, "/project/CRUP_scores/CENTRE_HiC/Training/CENTRE_final_training/GM12878.RNAPII-ChIAPET-Benchmark.MI.MSI.v38.txt",
          row.names = F)
?write.csv
colnames(RNAPII.GM.MI.MSI)
indx <- apply(RNAPII, 2, function(x) any(is.na(x) | is.infinite(x)))
colnames[indx]

DF <- lapply(RNAPII.GM.MI.MSI, is.infinite)
length(RNAPII.GM.MI.MSI$min_insulation[is.infinite(RNAPII.GM.MI.MSI$min_insulation)])
