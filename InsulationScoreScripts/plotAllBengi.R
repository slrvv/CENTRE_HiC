################################################################################
#
# Min Insulation : make a plot with facets with all the wilcoxon tests
#
################################################################################

# So far we did Wilcoxon tests on all of the BENGI data sets we had put them 
# all in one facetted plot

library(ggplot2)
library(ggpubr)
root <- "/project/CRUP_scores/CENTRE_HiC/"
setwd(root)


paths <- c("minInsulationGM12878.RNAPII-ChIAPET.csv",
           "minInsulationGM12878GEUVADIS.csv",
           "minInsulationGM12878_HiC-Benchmark.v38.txt.csv",
           "minInsulationWilcoxtestCTCF-ChIAPET.csv",
           "minInsulationWilcoxtestColon.GTEx.csv",
           "minInsulationWilcoxtestGMCHIC.png.csv",
           "minInsulationWilcoxtestIMR90.HiC.csv",
           "minInsulationWilcoxtestK562.CRISPR.csv",
           "minInsulationWilcoxtestK562.HiC.csv",
           "minInsulationWilcoxtestOvary.GTEx.csv",
           "minInsulationWilcoxtestPancreas.GTEx.csv",
           "minInsulationWilcoxtestStomach.GTEx.csv")

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
                           group.by = "namesample", method = "wilcox.test", 
                           alternative = "g")

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
