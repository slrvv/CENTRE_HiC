
#### Script to plot the F1 scores of the resolution optimized datasets

library(ggplot2)
library(ggpubr)
###-------------------FUNCTIONS-----------------------------------------------##

multicomparison_plot <- function(sample, centref1, centremsimif1, 
                                 maxf1, wilcox, compnum){
  
  ## Formatting both F1 score dasets accordingly
  centref1 <- centref1[,c("SampleName", "Model", "f1")] 
  centref1 <- centref1[centref1$Model == "CENTRE", ]
  centref1 <- centref1[centref1$SampleName == sample, ]

  
  centremsimif1$Model <- paste("CENTRE.MI.MSI",centremsimif1$Resolution, centremsimif1$WindowSize, sep = "_")
  centremsimif1 <- centremsimif1[centremsimif1$SampleName == sample, ]
  centremsimif1 <- centremsimif1[,c("SampleName","Model", "f1")]
  
  ## putting both datasets into one
  f1dat <- rbind(centref1, centremsimif1)
  
  if (wilcox){
    stat.test <- compare_means(f1~Model, f1dat, ref.group = "CENTRE")
    data_text <- data.frame(Model = stat.test$group2, 
                            f1 = rep(maxf1, compnum),
                            lab = stat.test$p.adj)
    p <- ggplot(f1dat, aes(x = Model, y = f1)) +
      geom_boxplot(aes(color = Model),) +
      scale_color_brewer(palette = "Dark2") +
      theme(
        legend.title = element_blank(),
        axis.title.x = element_blank(), 
        legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)
      ) +
      labs(y= "F1 score")
    p + geom_text(data = data_text,
                  mapping = aes(x = Model, y = f1, label = lab),
                  size = 4) + ggtitle(paste("Window optimization on", sample, "dataset", sep = " " ))
    
  } else {
    p <- ggplot(f1dat, aes(x = Model, y = f1)) +
      geom_boxplot(aes(color = Model),) +
      scale_color_brewer(palette = "Dark2") +
      theme(
        legend.title = element_blank(),
        axis.title.x = element_blank(), 
        legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5, size = 20),
        axis.text.y = element_text(size = 16)
      ) +
      labs(y= "F1 score")
    p + ggtitle(paste("Window optimization on", sample, "dataset", sep = " " ))
    
  }
  
  
  
}

###---------------------K562.HiC----------------------------------------------##

## K562 F1 boxplots diff window sizes compared to 10kb30KB

file10Kb <- "/project/CRUP_scores/CENTRE_HiC/Training/CENTRE_final_training/F1_scores/f1_allBENGI_10KB.csv"
f1dat <- read.table(file10Kb, header = T, sep = ",")
file5kb <- "/project/CRUP_scores/CENTRE_HiC/Training/CENTRE_final_training/F1_scores/f1_K562.HiC_5Kb.csv"
f1dat5kb <- read.table(file5kb, header = T, sep = "," )

multicomparison_plot("K562.HiC", f1dat, f1dat5kb, 0.4, F, 5)
multicomparison_plot("K562.HiC", f1dat, f1dat5kb, 0.4, T, 5)

######-----------------IMR-90.HiC---------------------------------------------##

file10Kb <- "/project/CRUP_scores/CENTRE_HiC/Training/CENTRE_final_training/F1_scores/f1_allBENGI_10KB.csv"
f1dat <- read.table(file10Kb, header = T, sep = ",")
file5kb <- "/project/CRUP_scores/CENTRE_HiC/Training/CENTRE_final_training/F1_scores/f1_K562.HiCandIMR90_5KbandK562CRISPR.csv"
f1dat5kb <- read.table(file5kb, header = T, sep = ",")
head(f1dat5kb)
multicomparison_plot("IMR90.HiC", f1dat, f1dat5kb, 0.5, F, 5)
multicomparison_plot("IMR90.HiC", f1dat, f1dat5kb, 0.5, T, 5)

head(f1dat10kb)
f1dat5kb <- f1dat5kb[f1dat5kb$SampleName == "IMR90.HiC", ]
head(f1dat5kb)
f1dat5kb$Model <- paste("CENTRE.MI.MSI",f1dat5kb$Resolution, f1dat5kb$WindowSize, sep = "_")
head(f1dat5kb)
f1dat5kb <- f1dat5kb[,c("SampleName","Model", "f1")]
f1dat10kb <- f1dat10kb[,c("SampleName","Model", "f1")]

f1dat <- rbind(f1dat10kb, f1dat5kb)
f1dat
?compare_means
stat.test <- compare_means(f1~Model, f1dat, ref.group = "CENTRE",
                           method = "t.test")

stat.test
data_text <- data.frame(Model = stat.test$group2, 
                        f1 = rep(0.5, 5),
                        lab = stat.test$p.signif)
p <- ggplot(f1dat, aes(x = Model, y = f1)) +
  geom_boxplot(aes(color = Model),) +
  scale_color_brewer(palette = "Dark2") +
  theme(
    legend.title = element_blank(),
    axis.title.x = element_blank(), 
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)
  ) +
  labs(y= "F1 score")
p
p + geom_text(data = data_text,
              mapping = aes(x = Model, y = f1, label = lab),
              size = 4)

###-------------------------K562.CRISPR--------------------------------------###

file10Kb <- "/project/CRUP_scores/CENTRE_HiC/Training/CENTRE_final_training/F1_scores/f1_allBENGI_10KB.csv"
f1dat10kb <- read.table(file10Kb, header = T, sep = ",")
file5kb <- "/project/CRUP_scores/CENTRE_HiC/Training/CENTRE_final_training/F1_scores/f1_K562.HiCandIMR90_5KbandK562CRISPR.csv"
f1dat5kb <- read.table(file5kb, header = T, sep = ",")


multicomparison_plot("K562.CRISPR", f1dat10kb, f1dat5kb, 0.5, F, 5)
multicomparison_plot("K562.CRISPR", f1dat10kb, f1dat5kb, 0.5, T, 5)

head(f1dat10kb)
f1dat5kb <- f1dat5kb[f1dat5kb$SampleName == "K562.CRISPR", ]
head(f1dat5kb)
f1dat5kb$Model <- paste("CENTRE.MI.MSI",f1dat5kb$Resolution, f1dat5kb$WindowSize, sep = "_")
head(f1dat5kb)
f1dat5kb <- f1dat5kb[,c("SampleName","Model", "f1")]
f1dat10kb <- f1dat10kb[,c("SampleName","Model", "f1")]

f1dat <- rbind(f1dat10kb, f1dat5kb)
f1dat
?compare_means
stat.test <- compare_means(f1~Model, f1dat, ref.group = "CENTRE",
                           method = "t.test")

stat.test
data_text <- data.frame(Model = stat.test$group2, 
                        f1 = rep(0.8, 5),
                        lab = stat.test$p.adj)
p <- ggplot(f1dat, aes(x = Model, y = f1)) +
  geom_boxplot(aes(color = Model),) +
  scale_color_brewer(palette = "Dark2") +
  theme(
    legend.title = element_blank(),
    axis.title.x = element_blank(), 
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)
  ) +
  labs(y= "F1 score")
p
p + geom_text(data = data_text,
              mapping = aes(x = Model, y = f1, label = lab),
              size = 4)

###--------------------GM12878 RNAPII ---------------------------------------###

### Load F1 scores for CENTRE 
file10Kb <- "/project/CRUP_scores/CENTRE_HiC/Training/CENTRE_final_training/F1_scores/f1_allBENGI_10KB.csv"
f1dat <- read.table(file10Kb, header = T, sep = ",")

file5kb <- "/project/CRUP_scores/CENTRE_HiC/Training/CENTRE_final_training/f1_GMRNAPII1kb.csv"

f1dat5kb <- read.table(file5kb, header = T, sep = ",")

multicomparison_plot("GM12878.RNAPII-ChIAPET", f1dat, f1dat5kb, 0.8, F)
multicomparison_plot("GM12878.RNAPII-ChIAPET", f1dat, f1dat5kb, 0.8, T, 3)


###--------------------GM12878 HIC---------------------------------------###

### Load F1 scores for CENTRE 
file10Kb <- "/project/CRUP_scores/CENTRE_HiC/Training/CENTRE_final_training/F1_scores/f1_allBENGI_10KB.csv"
f1dat <- read.table(file10Kb, header = T, sep = ",")

file5kb <- "/project/CRUP_scores/CENTRE_HiC/Training/CENTRE_final_training/f1_GMHiC1kb.csv"

f1dat5kb <- read.table(file5kb, header = T, sep = ",")

multicomparison_plot("GM12878.HiC", f1dat, f1dat5kb, 0.8, F)
multicomparison_plot("GM12878.RNAPII-ChIAPET", f1dat, f1dat5kb, 0.8, T, 3)





###-------------------------Dunnet procedure---------------------------------###
f1datMSI <- "/project/CRUP_scores/CENTRE_HiC/Training/CENTRE_final_training/f1_K562.HiCMSI.csv"
f1datMSI <- read.table(f1datMSI, header = T, sep = ",")
f1datMSI$Model <- paste("CENTRE.MSI",f1datMSI$Resolution, f1datMSI$WindowSize, sep = "_")
f1datMSI <- f1datMSI[,c("SampleName","Model", "f1")]

f1datMI <- "/project/CRUP_scores/CENTRE_HiC/Training/CENTRE_final_training/f1_K562.HiCMI.csv"
f1datMI <- read.table(f1datMI, header = T, sep = ",")
f1datMI$Model <- paste("CENTRE.MI",f1datMI$Resolution, f1datMI$WindowSize, sep = "_")
f1datMI <- f1datMI[,c("SampleName","Model", "f1")]
f1dat <- rbind(f1dat10kb, f1dat5kb, f1datMSI, f1datMI)
f1dat$ModelName <- c(rep("CENTRE", times = 12),rep("CENTRE.MI.MSI", times = 60),
                     rep("CENTRE.MSI", times = 60), rep("CENTRE.MI", times = 60))
f1dat
?emmeans
fit.f1dat <- aov(f1 ~ Model, data = f1dat)
f1.emm <- emmeans(fit.f1dat, specs = ~ Model)
cont <- contrast(f1.emm, method = "dunnett")
cont <- as.data.frame(cont)

stat.test <- compare_means(f1~Model, f1dat, ref.group = "CENTRE",
                           method = "t.test", p.adjust.method = "holm")

stat.test
models <- stringr::str_remove(cont$contrast, " -.*")
?strsplit
models
data_text <- data.frame(Model = as.factor(models), 
                        f1 = rep(0.4, 5),
                        lab = round(cont$p.value,3))
p <- ggplot(f1dat, aes(x = Model, y = f1)) +
  geom_boxplot(aes(color = ModelName),) +
  scale_color_brewer(palette = "Dark2") +
  theme(
    legend.title = element_blank(),
    axis.title.x = element_blank(), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)
  ) +
  labs(y= "F1 score")
p
p + geom_text(data = data_text,
              mapping = aes(x = Model, y = f1, label = lab),
              size = 4)