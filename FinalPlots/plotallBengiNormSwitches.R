################################################################################
#
# Plot of all Wilcoxtest Normalized count switches
#
################################################################################

#Plot in one plot the normalized Dir Count Switches

library(ggplot2)
library(ggpubr)
root <- "/project/CRUP_scores/CENTRE_HiC/WilcoxonTestResults/NormCountsDirIndex/"
setwd(root)


paths <- c("normCountSwitchesWilcoxtestGM12878.RNAPII-ChIAPET.csv",
           "normCountSwitchesWilcoxtestGM12878GEUVADIS.csv",
           "normCountSwitchesWilcoxtestGM12878_HiC-Benchmark.v38.txt.csv",
           "normCountSwitchesWilcoxtestCTCF-ChIAPET.csv",
           "normCountSwitchesWilcoxtestColon.GTEx.csv",
           "normCountSwitchesWilcoxtestGMCHIC.png.csv",
           "normCountSwitchesWilcoxtestIMR90.HiC.csv",
           "normCountSwitchesWilcoxtestK562.CRISPR.csv",
           "normCountSwitchesWilcoxtestK562.HiC.csv",
           "normCountSwitchesWilcoxtestOvary.GTEx.csv",
           "normCountSwitchesWilcoxtestPancreas.GTEx.csv",
           "normCountSwitchesWilcoxtestStomach.GTEx.csv")

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
  
  data <- read.table(pathdata[i, 2], header=T, sep =",")
  head(data)
  data <- data[c("label", "score_min")]
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

stat.test <- compare_means(score_min~label_s, allBengi, 
                           group.by = "namesample", method = "wilcox.test")

data_text <- data.frame(namesample = stat.test$namesample, 
                        score_min = rep(3, 12),
                        lab = stat.test$p.format)

p <- ggplot(allBengi, aes(x = namesample, y = score_min)) +
  geom_boxplot(aes(color = label_s),outlier.shape = NA, na.rm = T) +
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(limits = c(-8, 3)) +
  theme(
    legend.title = element_blank(),
    axis.title.x = element_blank(), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)
  ) +
  labs(y= "minimum insulation score over enhancer target window")
p + geom_text(data = data_text,
              mapping = aes(x = namesample, y = score_min, label = lab),
              size = 4)
