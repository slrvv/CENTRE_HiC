#-------------------------------------------------------------------------------
# Author:  Sara Lopez Ruiz de Vargas
# Email:  lopez_s@molgen.mpg.de
#
# Date:    2024-04-04
#
# Script Name:  /project/CRUP_scores/CENTRE_HiC/Training/CENTRE_final_training/accrossCTMIMSIaccuracyPlots.R
#
# Script Description: Plot the accuracy of samples with MI MSI values coming from
# other samples
#
#-------------------------------------------------------------------------------

library(ggplot2)

# The F1 scores taking the consensus LcL CENTRE MIMSI model for the bengi celltypes

f1scoresSameSample <- read.table("/project/CRUP_scores/CENTRE_HiC/Training/CENTRE_final_training/F1_scores/F1_scores_CENTREMIMSI_AllBENGI.csv",
                                 sep = ",", header = T) 
f1scoresdiffSample <- read.table("/project/CRUP_scores/CENTRE_HiC/Training/CENTRE_final_training/F1scoresAcrossCT/F1scoresAllCT.csv",
                                 sep = ",", header = T)
f1scoresdiffSample
selectlist <- function(x){
  if(x[2] == "ChIAPET"){
    return(x[3:length(x)])
  }
  return(x[2:length(x)])
}

##Making both datasets match for rbind
head(f1scoresdiffSample) #error in the naming in python. There are 
# empty columns just remove them

f1scoresdiffSample$X <- NULL
f1scoresdiffSample$SampleBENGI <- NULL
f1scoresdiffSample$SampleMSIMI<- NULL
unique(f1scoresdiffSample$Sample)
samenames <- c("Ovary.GTEx-ovary", "Pancreas.GTEx-pancreas", "Stomach.GTEx-stomach")
f1scoresdiffSample <- f1scoresdiffSample[!f1scoresdiffSample$Sample %in% samenames, ]
## This part is a bit dirty bc I chose a bad string joiner character
splitnames <- strsplit(f1scoresdiffSample$Sample, "-")
f1scoresdiffSample$SampleBENGI <- sapply(splitnames, "[[", 1)
f1scoresdiffSample$SampleBENGI[f1scoresdiffSample$SampleBENGI == "GM12878.CTCF"]  <- "GM12878.CTCF-ChIAPET"
f1scoresdiffSample$SampleBENGI[f1scoresdiffSample$SampleBENGI == "GM12878.RNAPII"]  <- "GM12878.RNAPII-ChIAPET"
listhic <- sapply(splitnames, selectlist)
length(listhic)
flatlisthic <- c()
for (element in listhic){
  if (length(element) == 2){
    flatlisthic <- c(flatlisthic,
                     paste(element[1], element[2], sep = "-"))
  } else {
  flatlisthic <- c(flatlisthic, element[1])
  }
}
f1scoresdiffSample$SampleMSIMI <- flatlisthic
f1scoresdiffSample$Sample <- NULL
f1scoresdiffSample

f1scoresdiffSample <- f1scoresdiffSample[,c("SampleBENGI", "SampleMSIMI", 
                                            "F1_Score_MIMSI", 
                                            "NumRows")]


f1scoresSameSample$SampleBENGI <- f1scoresSameSample$Sample
head(f1scoresSameSample)
f1scoresSameSample <- f1scoresSameSample[,c("SampleBENGI", "Sample", 
                                            "F1_Score_MIMSI", 
                                            "NumRows")]
colnames(f1scoresSameSample) <- c("SampleBENGI", "SampleMSIMI", 
                                  "F1_Score_MIMSI", 
                                  "NumRows")
unique(f1scoresSameSample$Sample)
f1scoresSameSample$SampleMSIMI <- c(rep("sigmoid-colon", 12),
                                    rep("GM12878", 12*5),
                                    rep("IMR-90", 12),
                                    rep("K562",2*12),
                                    rep("ovary", 12),
                                    rep("pancreas", 12),
                                    rep("stomach", 12))


samples <- unique(f1scoresdiffSample$SampleBENGI)
f1scoresSameSample <- f1scoresSameSample[f1scoresSameSample$SampleBENGI %in%
                                         samples,]
                   
f1scoresAll <- rbind(f1scoresSameSample, f1scoresdiffSample)

f1scoresAll$Type <- c(rep("Same", nrow(f1scoresSameSample)), rep("Diff", nrow(f1scoresdiffSample)))
unique(f1scoresAll$SampleBENGI)
tissues <- c("Colon.GTEx", "Ovary.GTEx", "Pancreas.GTEx", "Stomach.GTEx")

f1scoresAllTissues <- f1scoresAll[f1scoresAll$SampleBENGI %in% tissues,]
f1scoresAllTissues
f1scoresAllmean <- aggregate(x=f1scoresAll$F1_Score_MIMSI,
                             # Specify group indicator
                             by = list(paste(f1scoresAll$SampleBENGI, f1scoresAll$SampleMSIMI)),      
                             # Specify function (i.e. mean)
                             FUN = mean)
names <- strsplit(f1scoresAllmean$Group.1, " ")

f1scoresAllmean$BengiBenchmark <- sapply(names, "[[", 1)

f1scoresAllmean$HicSample <- sapply(names, "[[", 2)

f1scoresAllmean <- f1scoresAllmean[,c("BengiBenchmark", "HicSample", "x")]
colnames(f1scoresAllmean) <- c("BengiBenchmark", "HicSample", "F1scoreMean")

write.csv(f1scoresAllmean, 
          "/project/CRUP_scores/CENTRE_HiC/AcrossCTComparisonMSIMI/F1scoresAllMean.csv",
          row.names = F,
          quote = F)

p <- ggplot(data = f1scoresAllTissues, mapping = aes(x=SampleMSIMI, y = F1_Score_MIMSI, fill = Type)) + 
  geom_boxplot() + labs(x = "Tissue of MI MSI features", y = "F1 score") +
  facet_grid(SampleBENGI~.) 


p

ct <- c("GM12878.CHiC", "GM12878.CTCF-ChIAPET", "GM12878.GEUVADIS",      
        "GM12878.HiC", "GM12878.RNAPII-ChIAPET","IMR90.HiC", 
        "K562.CRISPR", "K562.HiC")

f1scoresAllct <- f1scoresAll[f1scoresAll$SampleBENGI %in% ct,]

p <- ggplot(data = f1scoresAllct, mapping = aes(x=SampleMSIMI, y = F1_Score_MIMSI, fill = Type)) + 
  geom_boxplot() + labs(x = "CT of MI MSI features", y = "F1 score") +
  facet_grid(SampleBENGI~.) 


p

