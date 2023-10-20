#Plot F1 scores after 12 CV for all BENGI and plot correlation heatmap

library(ggplot2)
library(ggpubr)
library(ggcorrplot)

#::install_github("kassambara/ggcorrplot")

##Load the data with all of the F1 scores 
file <- "/project/CRUP_scores/CENTRE_HiC/Training/CENTRE_final_training/f1_allBENGI_10KB.csv"
f1dat <- read.table(file, header = T, sep = ",")

## T.test to compare CENTRE to CENTRE.MSI.MI
stat.test <- compare_means(f1~Model, f1dat, 
                           group.by = "SampleName", method = "t.test")

data_text <- data.frame(SampleName = stat.test$SampleName, 
                        f1 = rep(0.75, 12),
                        lab = stat.test$p.format)

p <- ggplot(f1dat, aes(x = SampleName, y = f1)) +
  geom_boxplot(aes(color = Model),) +
  scale_color_brewer(palette = "Dark2") +
  theme(
    legend.title = element_blank(),
    axis.title.x = element_blank(), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)
  ) +
  labs(y= "F1 score")
p + geom_text(data = data_text,
              mapping = aes(x = SampleName, y = f1, label = lab),
              size = 4)

##Feature correlation heatmap for GM12878.RNAPII-ChIAPET
file <- "/project/CRUP_scores/CENTRE_HiC/Training/BENGI_MSI_MI_datasets/GM12878.RNAPII-ChIAPET-Benchmark.MI.MSI.v38.csv"
featdat <- read.table(file, header = T, sep=",")
featdat[is.na(featdat)] <- 0
head(featdat)

featdat$mean_EP_prob_enh <- rowMeans(featdat[, c("EP_prob_enh.1", 
                                           "EP_prob_enh.2",
                                           "EP_prob_enh.3", 
                                           "EP_prob_enh.4",
                                           "EP_prob_enh.5")])



featdat$mean_EP_prob_gene <- rowMeans(featdat[, c("EP_prob_gene.1", 
                                            "EP_prob_gene.2",
                                            "EP_prob_gene.3", 
                                            "EP_prob_gene.4",
                                            "EP_prob_gene.5")])

featdat$mean_PP_prob_enh <-rowMeans(featdat[, c("PP_prob_enh.1", 
                                          "PP_prob_enh.2",
                                          "PP_prob_enh.3", 
                                          "PP_prob_enh.4",
                                          "PP_prob_enh.5")])

featdat$mean_PP_prob_gene <- rowMeans(featdat[, c("PP_prob_gene.1", 
                                            "PP_prob_gene.2",
                                            "PP_prob_gene.3", 
                                            "PP_prob_gene.4",
                                            "PP_prob_gene.5")])
correl <- round(cor(featdat[,29:42]), 1)
ggcorrplot(correl,
           type = "lower",
           lab = TRUE)


