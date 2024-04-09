#Plot F1 scores after 12 CV for all BENGI and plot correlation heatmap

library(ggplot2)
library(ggpubr)
library(ggcorrplot)

#::install_github("kassambara/ggcorrplot")

##Load the data with all of the F1 scores 
file <- "/project/CRUP_scores/CENTRE_HiC/Training/CENTRE_final_training/F1_scores/f1_allBENGI_10KB.csv"
f1dat <- read.table(file, header = T, sep = ",")

f1dat[f1dat$SampleName == "GM12878.CTCF-ChIAPET", 2] <- "GM12878.CTCF-\nChIAPET"

f1dat[f1dat$SampleName == "GM12878.RNAPII-ChIAPET",2] <- "GM12878.\nRNAPII-ChIAPET"
## T.test to compare CENTRE to CENTRE.MSI.MI
?compare_means
stat.test <- compare_means(f1~Model, f1dat, 
                           group.by = "SampleName", method = "wilcox.test")

stat.test
data_text <- data.frame(SampleName = stat.test$SampleName, 
                        f1 = rep(0.75, 12),
                        lab = stat.test$p.signif)

p <- ggplot(f1dat, aes(x = SampleName, y = f1)) +
  geom_boxplot(aes(color = Model),) +
  scale_color_brewer(palette = "Dark2") +
  labs(y= "F1 score")
p + geom_text(data = data_text,
              mapping = aes(x = SampleName, y = f1, label = lab),
              size = 4) +theme(
                legend.title = element_blank(),
                legend.text = element_text(size = 12),
                legend.key.height= unit(2, 'cm'),
                legend.key.width= unit(2, 'cm'),
                axis.title.x = element_blank(), 
                axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5, size = 12),
                axis.text.y = element_text(size = 12),
                axis.title.y = element_text(size = 14)
              ) 

##Feature correlation heatmap for GM12878.RNAPII-ChIAPET
file <- "/project/CRUP_scores/CENTRE_HiC/Training/BENGI_MSI_MI_datasets/10Kb/GM12878.RNAPII-ChIAPET-Benchmark.MI.MSI.v38.csv"
featdat <- read.table(file, header = T, sep=",")
featdat[is.na(featdat)] <- 0
head(featdat)
featdat$distance <- abs(featdat$distance)

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



file <- "/project/CRUP_scores/CENTRE_HiC/Training/BENGI_MSI_MI_datasets/consensusLcL.MI.MSI.v38.10kb30kb.csv"
featdat <- read.table(file, header = T, sep=",")
featdat[is.na(featdat)] <- 0
head(featdat)
featdat$distance <- abs(featdat$distance)

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


