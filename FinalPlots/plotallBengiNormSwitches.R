################################################################################
#
# Plot of all Wilcoxtest Normalized count switches
#
################################################################################

#Plot in one plot the normalized Dir Count Switches

library(ggplot2)
library(ggpubr)
root <- "/project/CRUP_scores/CENTRE_HiC/Features/NormCountsDirIndex/"
setwd(root)


paths <- c("normCountSwitchesWilcoxtestGM12878.RNAPII-ChIAPET.csv",
           "normCountSwitchesWilcoxtestGM12878.GEUVADIS.csv",
           "normCountSwitchesWilcoxtestGM12878.HiC.csv",
           "normCountSwitchesWilcoxtestGM12878.CTCF-ChIAPET.csv",
           "normCountSwitchesWilcoxtestColon.GTEx.csv",
           "normCountSwitchesWilcoxtestGM12878.CHiC.csv",
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
  data <- data[c("label", "norm_switch", "Freq.x")]
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

allBengi[allBengi$namesample == "GM12878.CTCF.ChIAPET", 4] <- "GM12878.CTCF-\nChIAPET"
allBengi[allBengi$namesample == "GM12878.RNAPII.ChIAPET", 4] <- "GM12878.\nRNAPII-ChIAPET"
stat.test <- compare_means(norm_switch~label_s, allBengi, 
                           group.by = "namesample", method = "wilcox.test")

data_text <- data.frame(namesample = stat.test$namesample, 
                        score_min = rep(15, 12),
                        lab = stat.test$p.signif)

p <- ggplot(allBengi, aes(x = namesample, y = norm_switch)) +
  geom_boxplot(aes(color = label_s),outlier.shape = NA, na.rm = T) +
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(limits = c(NA, 16)) +
  labs(y= "normalized switch counts over enhancer target window")

p + geom_text(data = data_text,
              mapping = aes(x = namesample, y = score_min, label = lab),
              size = 6)+theme(
              legend.title = element_blank(),
              legend.text = element_text(size = 12),
              legend.key.height= unit(2, 'cm'),
              legend.key.width= unit(2, 'cm'),
              axis.title.x = element_blank(), 
              axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5, size = 12),
              axis.text.y = element_text(size = 12),
              axis.title.y = element_text(size = 14)
            ) 


stat.test <- compare_means(Freq.x~label_s, allBengi, 
                           group.by = "namesample", method = "wilcox.test")

data_text <- data.frame(namesample = stat.test$namesample, 
                        Freq.x = rep(47, 12),
                        lab = stat.test$p.signif)
stat.test
p <- ggplot(allBengi, aes(x = namesample, y = Freq.x)) +
  geom_boxplot(aes(color = label_s),outlier.shape = NA, na.rm = T) +
  scale_color_brewer(palette = "Dark2") + 
  labs(y= "absolute switch counts over enhancer target window")
p
p + geom_text(data = data_text,
              mapping = aes(x = namesample, y = Freq.x, label = lab),
              size = 6)+theme(
                legend.title = element_blank(),
                legend.text = element_text(size = 12),
                legend.key.height= unit(2, 'cm'),
                legend.key.width= unit(2, 'cm'),
                axis.title.x = element_blank(), 
                axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5, size = 12),
                axis.text.y = element_text(size = 12),
                axis.title.y = element_text(size = 14)
              ) 
