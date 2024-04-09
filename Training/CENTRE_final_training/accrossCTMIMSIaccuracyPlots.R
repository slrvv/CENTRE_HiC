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
f1scoresdiffSample <- read.table("/project/CRUP_scores/CENTRE_HiC/Training/CENTRE_final_training/F1scoresAcrossCT/F1scoresTissues.csv",
                                 sep = ",", header = T)


##Making both datasets match for rbind
f1scoresdiffSample$X <- NULL

f1scoresSameSample$SampleBENGI <- f1scoresSameSample$Sample
f1scoresSameSample <- f1scoresSameSample[,c("SampleBENGI", "Sample", 
                                            "F1_Score_MIMSI", 
                                            "NumRows")]
colnames(f1scoresSameSample) <- c("SampleBENGI", "SampleMSIMI", 
                                  "F1_Score_MIMSI", 
                                  "NumRows")



samples <- unique(f1scoresdiffSample$SampleBENGI)
f1scoresSameSample <- f1scoresSameSample[f1scoresSameSample$SampleBENGI %in%
                                           samples,]

                      
f1scoresAll <- rbind(f1scoresSameSample, f1scoresdiffSample)

f1scoresAll$Type <- c(rep("Same", nrow(f1scoresSameSample)), rep("Diff", nrow(f1scoresdiffSample)))

f1scoresAllmean <- aggregate(x=f1scoresAll$F1_Score_MIMSI,
                             # Specify group indicator
                             by = list(paste(f1scoresAll$SampleBENGI, f1scoresAll$SampleMSIMI)),      
                             # Specify function (i.e. mean)
                             FUN = mean)

library(stringr)

f1scorestext <- cbind(str_split_fixed(f1scoresAllmean$Group.1, " ", 2), 
                      as.numeric(f1scoresAllmean$x))
f1scorestext <- as.data.frame(f1scorestext)
f1scorestext$V3 <- as.numeric(f1scorestext$V3)
colnames(f1scorestext) <- c("SampleBENGI", "SampleMSIMI", 
                            "F1_Score_MIMSI")
f1scoresAll$helper <- paste(f1scoresAll$SampleBENGI, f1scoresAll$SampleMSIMI)
f1scorestext$helper <- paste(f1scorestext$SampleBENGI, f1scorestext$SampleMSIMI)

f1scoreAllhelp <- aggregate(f1scoresAll$NumRows,list(f1scoresAll$helper),unique)

f1scorestext$NumRows <- f1scoreAllhelp$x
f1scorestext
f1scorestext$Type <- c("Same", rep("Diff", 8), "Same", rep("Diff", 4), 
                       "Same", rep("Diff", 4), "Same")
f1scorestext$F1_Score_MIMSI <- rep(1.0, nrow(f1scorestext))

p <- ggplot(data = f1scoresAll, mapping = aes(x=SampleMSIMI, y = F1_Score_MIMSI, fill = Type)) + 
  geom_boxplot() + labs(x = "Tissue of MI MSI features", y = "F1 score") +
  facet_nested("Tissue of remaining features"*SampleBENGI~.) + 
  geom_text(data = f1scorestext, aes(x = SampleMSIMI,
                                     y = F1_Score_MIMSI,
                                     label = NumRows))


p
