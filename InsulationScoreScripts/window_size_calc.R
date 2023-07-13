################################################################################
#                                                                              #
#   window_size_calc: Check how many times the Insulation score falls inside   #                                                                            #
#                     ET pair and how many times outside.                      #
################################################################################

# Issue: How fine grained does our insulation score need to be? 
#
# What we would ideally want: 
# Insulation score                  score
# window:           -------------#++++++++++++#-------------
# ET pair window : ---------###+++++++++++++++++######--------
#                           enhancer            promoter
#
# What we maybe don't want: 
# Insulation score                  score
# window:           ---#+++++++++++++++++++++++++++++++++++++#------
# ET pair window : ---------###+++++++++++++++++######--------
#                           enhancer            promoter
#
#
# EXTRA: Also question whether it matters at all.


#--------------------------Libraries and functions-----------------------------#

library(CENTRE)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(ggprism)
#install.packages("ggprism")
#-------------------------Script start ----------------------------------------#

# HiC window of 10kb, insulation score window of 30kb.
insulation10kb <- read.table("/project/CRUP_scores/CENTRE_HiC/InsulationScores/GM12878/ENCFF237QCN.hic@10kb.insulation_30kb.bed", 
                             sep = "\t")
insulation10kb[!is.finite(insulation10kb$V5),] <- NA ##remove infinite values 
insulation10kb <- na.omit(insulation10kb) #remove NA values

#NA values produced by insulation scores are bc.
#In the insulation score calculation, if the insulation window is covered by
#more than 50% of unmappable regions, the score will be NaN 

insulation10kb <- insulation10kb[insulation10kb$V1 == "chr20",]

## BENGI data 
GM12878BENGI <- read.table("/project/CRUP_scores/toSara/BENGI_processed_datasets/GM12878.RNAPII-ChIAPET-Benchmark.v38.txt",
                           header = T,
                           sep = "\t")

##Take only BENGI data from chr20

conn <- RSQLite::dbConnect(RSQLite::SQLite(),
                           system.file("extdata",
                                       "Annotation.db",
                                       package = "CENTRE"))
#get chromosome   and tts of our genes
query <- paste("SELECT  gene_id1, chr, transcription_start FROM gencode WHERE chr = 'chr20'")
gene <- RSQLite::dbGetQuery(conn, query)


GM12878BENGIchr20 <- GM12878BENGI[GM12878BENGI$gene_id1 %in% gene$gene_id1,] 

GM12878BENGIchr20 <- GM12878BENGIchr20[,c("gene_id1",
                                          "symbol38",
                                          "label",
                                          "distance")]

## Add TTS and middle point


query_enh <-  paste("SELECT  V5, middle_point FROM ccres_enhancer WHERE V1 = 'chr20'")
geneMiddlePoint <- RSQLite::dbGetQuery(conn, query_enh)

RSQLite::dbDisconnect(conn)

GM12878BENGIchr20 <- merge(GM12878BENGIchr20,
                geneMiddlePoint[, c("V5", "middle_point")],
                by.x = "symbol38",
                by.y = "V5") #change V1 and V5 to more meaningful names
#Getting the chrosomes and transcription start sites for the provided genes
GM12878BENGIchr20 <- merge(GM12878BENGIchr20,
                gene[, c("gene_id1", "transcription_start")],
                by.x = "gene_id1",
                by.y = "gene_id1")


##MAKE IT INTO Granges object
head(insulation10kb)

insulationRange <- with(insulation10kb,
                   GenomicRanges::GRanges(V1,
                                          IRanges::IRanges(start = V2,
                                                           end = V3),
                                          mcols = insulation10kb$V5))

# ordering if start is at enhancer or promoter
start <- c()
end <- c()
for (i in 1:nrow(GM12878BENGIchr20)){
  if (GM12878BENGIchr20$transcription_start[i] <= GM12878BENGIchr20$middle_point[i]){
    start <- c(start,GM12878BENGIchr20$transcription_start[i])
    end <- c(end, GM12878BENGIchr20$middle_point[i])
  } else {
    start <- c(start, GM12878BENGIchr20$middle_point[i])
    end <- c(end, GM12878BENGIchr20$transcription_start[i] )
  }
}
GM12878BENGIchr20$start <- start
GM12878BENGIchr20$end <- end
GM12878BENGIchr20$chr <- rep("chr20", times = nrow(GM12878BENGIchr20))
GM12878BENGIchr20$pair <- paste(GM12878BENGIchr20$symbol38, GM12878BENGIchr20$gene_id1, sep ="_")

GM12878ETRanges<-  with(GM12878BENGIchr20,
                      GenomicRanges::GRanges( chr,
                                             IRanges::IRanges(start = start,
                                                              end = end),
                                             mcols = pair))


# Find the overlaps in the windows
overlaps <- GenomicRanges::findOverlaps(GM12878ETRanges,
                                        insulationRange, 
                                        ignore.strand = TRUE)
score_overlapping <-data.frame(ET = overlaps@from, score = overlaps@to)
score_overlapping$pair <- GM12878BENGIchr20$pair[score_overlapping$ET]
score_overlapping$insulation <- insulation10kb$V5[score_overlapping$score]

length(unique(score_overlapping$pair)) == nrow(GM12878BENGIchr20) ## This is true 
## meaning each pair overlaps with at least one window score.
repeats <- table(score_overlapping$pair)
table(score_overlapping$score) ###the same score also usually maps to multiple pairs

outside <- c() #is the score window bigger than the pair window
inside <- c() #is the score window smaller than the pair window
for (i in 1:nrow(score_overlapping)){
  startET <- GM12878BENGIchr20$start[score_overlapping$ET[i]]
  endET <- GM12878BENGIchr20$end[score_overlapping$ET[i]]
  
  startIns <- insulation10kb$V2[score_overlapping$score[i]]
  endIns <- insulation10kb$V3[score_overlapping$score[i]]
  if (startET > startIns && endET < endIns){
    outside <- c(outside, T)
    inside <- c(inside, F)
  } else {
    outside <- c(outside, F)
    inside <- c(inside, T)
  }
}

score_overlapping$inside <- inside
score_overlapping$outside <- outside
outside <- score_overlapping[score_overlapping$outside == T, ]
inside <- score_overlapping[score_overlapping$outside == F, ]
cat(paste("There are ",
         length(unique(outside$pair)),
               " ET pairs where the score window falls outside"))
cat(paste("There are ",
          length(unique(inside$pair)),
          " ET pairs where the score window falls inside"))

score_min <- c()
for (i in 1:nrow(GM12878BENGIchr20)){
  score <- min(score_overlapping[score_overlapping$ET == i, "insulation"])
  score_min <- c(score_min, score)
}

GM12878BENGIchr20$score_min <- score_min

x <- GM12878BENGIchr20[GM12878BENGIchr20$label == 1, ]
y <- GM12878BENGIchr20[GM12878BENGIchr20$label == 0, ]
stat.test <- wilcox.test(x$score_min, y$score_min, alternative = "greater") 
GM12878BENGIchr20$label <- as.factor(GM12878BENGIchr20$label)

df_p_val <- data.frame(
  group1 = "0",
  group2 = "1",
  label = stat.test$p.value,
  y.position = 2
)
p <- ggplot(GM12878BENGIchr20, aes(x = label, y = score_min)) +
  geom_boxplot()

p+add_pvalue(df_p_val,
             xmin = "group1",
             xmax = "group2",
             label = "label",
             y.position = "y.position") + 
      labs(title="Wilcoxon test on BENGI GM12878.RNAPII-ChIAPET chr 20",
       x ="Active ET pair", y = "Minimum insulation score")

###Now expand this to all BENGI datasets