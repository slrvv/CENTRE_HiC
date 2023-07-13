################################################################################
#                                                                              #
#   Wilcox_test_insscore: A script to test whether the minimum insulation score#                                                                         #
#                         is a good metric                                     #
################################################################################

# We designed the following feature: minimum of the insulation scores overlapping
# with the ET window
# To test whether it is different between active et pairs and non-active ET pairs
# we do a Wilcoxon test comparing the group of active and inactive pairs for
# each BENGI dataset.
# The insulation scores are calculated on a 10kb window
# We expect the minimum score to be greater when the pair is active

#------------------------Libraries and functions-------------------------------#

library(CENTRE)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(ggprism)

readInsulation <- function(path){
  ##### A function to read in the insulation scores and make a dataframe #####
  insulation <- read.table(path, sep = "\t")
  insulation[!is.finite(insulation$V5),] <- NA ##remove infinite values 
  insulation <- na.omit(insulation) #remove NA values
  
  #NA values produced by insulation scores are bc.
  #In the insulation score calculation, if the insulation window is covered by
  #more than 50% of unmappable regions, the score will be NaN 
  return(insulation)
}

makeInsRanges <- function(insulation){
  ##### Turn insulation dataframe into GRanges object ######
  
  insulationRange <- with(insulation,
                          GenomicRanges::GRanges(V1,
                                                 IRanges::IRanges(start = V2,
                                                                  end = V3),
                                                 mcols = insulation))
  return(insulationRange)
}

readBENGI <- function(path){
  BENGIdf <- read.table(path,
                             header = T,
                             sep = "\t")
  print(head(BENGIdf))
  
  BENGIdf <- BENGIdf[,c("gene_id1",
                        "symbol38",
                        "label",
                        "distance")]
  ##Get the TTS and middle point
  conn <- RSQLite::dbConnect(RSQLite::SQLite(),
                             system.file("extdata",
                                         "Annotation.db",
                                         package = "CENTRE"))
  #get chromosome   and tts of our genes
  query <- paste("SELECT  gene_id1, chr, transcription_start FROM gencode WHERE gene_id1 in (",
                 paste0(sprintf("'%s'", BENGIdf$gene_id1), collapse = ", "),")",sep="" )
  gene <- RSQLite::dbGetQuery(conn, query)
  

  
  
  query_enh <-  paste("SELECT V1, V5, middle_point FROM ccres_enhancer WHERE V5 in (",
                      paste0(sprintf("'%s'", BENGIdf$symbol38), collapse = ", "),")",sep="" )
  geneMiddlePoint <- RSQLite::dbGetQuery(conn, query_enh)
  
  RSQLite::dbDisconnect(conn)
  
  BENGIdf <- merge(BENGIdf,
                             geneMiddlePoint[, c("V5", "middle_point")],
                             by.x = "symbol38",
                             by.y = "V5") #change V1 and V5 to more meaningful names
  #Getting the chrosomes and transcription start sites for the provided genes
  BENGIdf <- merge(BENGIdf,
                             gene[, c("gene_id1", "transcription_start", "chr")],
                             by.x = "gene_id1",
                             by.y = "gene_id1")
  return(BENGIdf)
}

makeETRanges <- function(BENGIdf){
  # ordering if start is at enhancer or promoter
  start <- c()
  end <- c()
  for (i in 1:nrow(BENGIdf)){
    if (BENGIdf$transcription_start[i] <= BENGIdf$middle_point[i]){
      start <- c(start,BENGIdf$transcription_start[i])
      end <- c(end, BENGIdf$middle_point[i])
    } else {
      start <- c(start, BENGIdf$middle_point[i])
      end <- c(end, BENGIdf$transcription_start[i] )
    }
  }
  BENGIdf$start <- start
  BENGIdf$end <- end
  BENGIdf$pair <- paste(BENGIdf$symbol38, BENGIdf$gene_id1, sep ="_")
  
  ETRanges<-  with(BENGIdf,
                          GenomicRanges::GRanges( chr,
                                                  IRanges::IRanges(start = start,
                                                                   end = end),
                                                  mcols = pair))
  return(ETRanges)
}

findETScoreOverlaps <- function(ETRanges, insulationRange, BENGIdf, insulation){
  overlaps <- GenomicRanges::findOverlaps(ETRanges,
                                          insulationRange, 
                                          ignore.strand = TRUE)
  score_overlapping <-data.frame(ET = overlaps@from, score = overlaps@to)
  score_overlapping$pair <- BENGIdf$pair[score_overlapping$ET]
  score_overlapping$insulation <- insulation$V5[score_overlapping$score]
  
  
  
  return(score_overlapping)
}

computeMin <- function(BENGIdf, score_overlapping){
  score_min <- c()
  for (i in 1:nrow(BENGIdf)){
    score <- min(score_overlapping[score_overlapping$ET == i, "insulation"])
    score_min <- c(score_min, score)
  }
  
  BENGIdf$score_min <- score_min
  return(BENGIdf)
}


#---------------------------Script start --------------------------------------#

args = commandArgs(trailingOnly=TRUE)
cat("Reading Data in")
GM12878BENGI <- readBENGI(args[1])
insulationscores <- readInsulation(args[2])

GM12878Ranges <- makeETRanges(GM12878BENGI)
insulationRanges <- makeInsRanges(insulationscores)

cat("Find overlaps")
overlaps <- findETScoreOverlaps(GM12878Ranges,
                                insulationRanges,
                                GM12878BENGI,
                                insulationscores)
cat("Compute Insulation minimum")
overlapsmin <- computeMin(GM12878BENGI, overlaps)


cat("Wilcoxon test")
## Wilcoxon test
x <- overlapsmin[overlapsmin$label == 1, ]
y <- overlapsmin[overlapsmin$label == 0, ]
stat.test <- wilcox.test(x$score_min, y$score_min, alternative = "g")
overlapsmin$label <- as.factor(overlapsmin$label)

cat("Box Plot")
df_p_val <- data.frame(
  group1 = "0",
  group2 = "1",
  label = stat.test$p.value,
  y.position = 3
)
png(paste0("/project/CRUP_scores/CENTRE_HiC/", args[3]))
p <- ggplot(overlapsmin, aes(x = label, y = score_min)) +
  geom_boxplot()

p+add_pvalue(df_p_val,
             xmin = "group1",
             xmax = "group2",
             label = "label",
             y.position = "y.position") + 
  labs(title=paste0("Wilcoxon test on BENGI", args[4]),
       x ="Active ET pair", y = "Minimum insulation score")
dev.off()
write.csv(overlapsmin, paste0("/project/CRUP_scores/CENTRE_HiC/minInsulation",
                              args[3], ".csv"),
          row.names = F)
write.csv(overlaps, paste0("/project/CRUP_scores/CENTRE_HiC/overlaps",
                           args[3], ".csv"),
          row.names = F)
