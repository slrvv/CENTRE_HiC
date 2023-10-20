
## K562 F1 boxplots diff window sizes compared to 10kb30KB

file10Kb <- "/project/CRUP_scores/CENTRE_HiC/Training/CENTRE_final_training/f1_allBENGI_10KB.csv"
f1dat <- read.table(file, header = T, sep = ",")
file5kb <- "/project/CRUP_scores/CENTRE_HiC/Training/CENTRE_final_training/f1_K562.HiC_5Kb.csv"
#f1dat5kb <- read.table(file5kb, header = T, se 
head(f1dat5kb)
f1dat10kb <- f1dat[,c("SampleName", "Model", "f1")] 
f1dat10kb <- f1dat10kb[f1dat10kb$Model == "CENTRE", ]
f1dat10kb <- f1dat10kb[f1dat10kb$SampleName == "K562.HiC", ]

head(f1dat10kb)
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
                        f1 = rep(0.4, 5),
                        lab = stat.test$p.format)
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


### imr90

file10Kb <- "/project/CRUP_scores/CENTRE_HiC/Training/CENTRE_final_training/f1_allBENGI_10KB.csv"
f1dat <- read.table(file, header = T, sep = ",")
file5kb <- "/project/CRUP_scores/CENTRE_HiC/Training/CENTRE_final_training/f1_K562.HiCandIMR90_5Kb.csv"
f1dat5kb <- read.table(file5kb, header = T, sep = ",")
head(f1dat5kb)
f1dat10kb <- f1dat[,c("SampleName", "Model", "f1")] 
f1dat10kb <- f1dat10kb[f1dat10kb$Model == "CENTRE", ]
f1dat10kb <- f1dat10kb[f1dat10kb$SampleName == "IMR90.HiC", ]

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

file10Kb <- "/project/CRUP_scores/CENTRE_HiC/Training/CENTRE_final_training/f1_allBENGI_10KB.csv"
f1dat <- read.table(file10Kb, header = T, sep = ",")
file5kb <- "/project/CRUP_scores/CENTRE_HiC/Training/CENTRE_final_training/f1_K562.HiCandIMR90_5KbandK562CRISPR.csv"
f1dat5kb <- read.table(file5kb, header = T, sep = ",")
head(f1dat5kb)
f1dat10kb <- f1dat[,c("SampleName", "Model", "f1")] 
f1dat10kb <- f1dat10kb[f1dat10kb$Model == "CENTRE", ]
f1dat10kb <- f1dat10kb[f1dat10kb$SampleName == "K562.CRISPR", ]

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


file10Kb <- "/project/CRUP_scores/CENTRE_HiC/Training/CENTRE_final_training/f1_allBENGI_10KB.csv"
f1dat <- read.table(file10Kb, header = T, sep = ",")
file5kb <- "/project/CRUP_scores/CENTRE_HiC/Training/CENTRE_final_training/f1_K562.HiCandIMR90_5KbandK562CRISPR.csv"
f1dat5kb <- read.table(file5kb, header = T, sep = ",")
head(f1dat5kb)
f1dat10kb <- f1dat[,c("SampleName", "Model", "f1")] 
f1dat10kb <- f1dat10kb[f1dat10kb$Model == "CENTRE", ]
f1dat10kb <- f1dat10kb[f1dat10kb$SampleName == "K562.HiC", ]

head(f1dat10kb)
f1dat5kb <- f1dat5kb[f1dat5kb$SampleName == "K562.HiC", ]
head(f1dat5kb)
f1dat5kb$Model <- paste("CENTRE.MI.MSI",f1dat5kb$Resolution, f1dat5kb$WindowSize, sep = "_")
head(f1dat5kb)
f1dat5kb <- f1dat5kb[,c("SampleName","Model", "f1")]
f1dat10kb <- f1dat10kb[,c("SampleName","Model", "f1")]

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