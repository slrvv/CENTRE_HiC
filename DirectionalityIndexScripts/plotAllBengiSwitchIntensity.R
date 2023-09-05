################################################################################
#
# Mean Switch Intensity
#
################################################################################

# So far we did Wilcoxon tests on all of the BENGI data sets we had put them 
# all in one facetted plot

library(ggplot2)
library(ggpubr)
root <- "/project/CRUP_scores/CENTRE_HiC/WilcoxonTestResults/MeanSwitchIntensity/"
setwd(root)


paths <- c("meanSwitchIntensityWilcoxtestGM12878.RNAPII-ChIAPET.csv",
           "meanSwitchIntensityWilcoxtestGM12878.GEUVADIS.csv",
           "meanSwitchIntensityWilcoxtestGM12878_HiC.csv",
           "meanSwitchIntensityWilcoxtestGM12878.CTCF-ChIAPET.csv",
           "meanSwitchIntensityWilcoxtestColon.GTEx.csv",
           "meanSwitchIntensityWilcoxtestGM12878.CHiC.csv",
           "meanSwitchIntensityWilcoxtestIMR90.HiC.csv",
           "meanSwitchIntensityWilcoxtestK562.CRISPR.csv",
           "meanSwitchIntensityWilcoxtestK562.HiC.csv",
           "meanSwitchIntensityWilcoxtestOvary.GTEx.csv",
           "meanSwitchIntensityWilcoxtestPancreas.GTEx.csv",
           "meanSwitchIntensityWilcoxtestStomach.GTEx.csv")

samplenames <- c("GM12878.RNAPII.ChIAPET",
                 "GM12878.GEUVADIS",
                 "GM12878.HiC",
                 "GM12878.CTCF.ChIAPET",
                 "Colon.GTEx",
                 "GM12878.cHiC",
                 "IMR90.HiC",
                 "K562.CRISPR", 
                 "K562.HiC",
                 "Ovary.GTEx",
                 "Pancreas.GTEx",
                 "Stomach.GTEx")
pathdata <- data.frame(samplenames, paths)
for (i in seq_len(nrow(pathdata))){
  
  data <- read.table(pathdata[i, 2], header=T, sep =" ")
  data <- data[,c("label", "mean_switch_intensity")]
  namesample <- rep(pathdata[i,1], nrow(data))
  data$namesample <- namesample
  assign(pathdata[i,1], data)
}


allBengi <- rbind(GM12878.RNAPII.ChIAPET,
                  GM12878.GEUVADIS,
                  GM12878.HiC,
                  GM12878.CTCF.ChIAPET,
                  Colon.GTEx,
                  GM12878.cHiC,
                  IMR90.HiC,
                  K562.CRISPR, 
                  K562.HiC,
                  Ovary.GTEx,
                  Pancreas.GTEx,
                  Stomach.GTEx)


allBengi$label_s <- factor(allBengi$label==1, labels = c("non-interacting pairs", 
                                                "interacting pairs"))

stat.test <- compare_means(mean_switch_intensity~label_s, allBengi, 
                           group.by = "namesample", method = "wilcox.test",
                           alternative = "l")

stat.test
data_text <- data.frame(namesample = stat.test$namesample, 
                        mean_switch_intensity = rep(225, 12),
                        lab = stat.test$p.format)

p <- ggplot(allBengi, aes(x = namesample, y = mean_switch_intensity)) +
  geom_boxplot(aes(color = label_s),outlier.shape = NA, na.rm = T) +
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(limits = c(NA, 225)) +
  theme(
    legend.title = element_blank(),
    axis.title.x = element_blank(), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)
  ) +
  labs(y= "mean of switch intensities over enhancer target window")
p + geom_text(data = data_text,
              mapping = aes(x = namesample, y = mean_switch_intensity, label = lab),
              size = 4)
