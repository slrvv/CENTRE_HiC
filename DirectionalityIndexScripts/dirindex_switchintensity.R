################################################################################
#
# mean of Switch intensity
#
################################################################################

# Finding a way to aggregate the dir index values over ET window
# Idea: Switch intensity as defined by the absolute value of the substraction of
# the switch 
# Then map ET window over switch intensity and take the mean if multiple values
# map  
# The idea is the higher the switch the more likely a TAD boundary and the less
# likely a Enhancer target link 

library(CENTRE)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(ggprism)
library(dplyr)

readDirectionality <- function(path){
  ##### A function to read in the directionality scores and make a dataframe #####
  directionality <- read.table(path, sep = "\t")
  directionality[!is.finite(directionality$V5),] <- NA ##remove infinite values 
  directionality <- na.omit(directionality) #remove NA values
  
  return(directionality)
}

computeDirSwitches <- function(directionality){
  ###Compute where the sign switch happens###
  #dirswitch 
  dirswitch <- data.frame(matrix(ncol = 6, nrow = 0))
  colnames(dirswitch) <- c("chr", "start_switch_neg", "end_switch_neg",
                           "satrt_switch_pos", "end_switch_pos", "switch_intensity")
  
  for (i in 2:nrow(directionality)){
    if (directionality[i-1,5] < 0 & directionality[i,5] > 0 & directionality[i,2] == directionality[i-1,3]+1){
      row <- c(directionality[i-1, 1],
               directionality[i-1, 2],
               directionality[i-1, 3],
               directionality[i, 2],
               directionality[i, 3],
               abs(directionality[i-1,5]-directionality[i,5]))
      dirswitch <- rbind(dirswitch, row)  
    }
  }
  return(dirswitch)
}

makeDirSwitchRanges <- function(directionality){
  ##### Turn directionality switches dataframe into GRanges object ######
  
  directionalityRange <- with(directionality,
                              GenomicRanges::GRanges(chr,
                                                     IRanges::IRanges(start = end_switch_neg,
                                                                      end = end_switch_pos),
                                                     mcols = switch_intensity))
  return(directionalityRange)
}

readBENGI <- function(path){
  ##### Read in BENGI benchmark and get the ET start and end positions ######
  BENGIdf <- read.table(path,
                        header = T,
                        sep = "\t")
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
  ### Make the ET GRanges object
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

findETScoreOverlaps <- function(ETRanges, directionalityRange, BENGIdf, directionality){
  ### overlap the BENGI ranges and the directonality ranges ####
  overlaps <- GenomicRanges::findOverlaps(ETRanges,
                                          directionalityRange, 
                                          ignore.strand = TRUE)
  score_overlapping <-data.frame(ET = overlaps@from, score = overlaps@to)
  score_overlapping$pair <- BENGIdf$pair[score_overlapping$ET]
  score_overlapping$directionality <- directionality$switch_intensity[score_overlapping$score]
  
  
  
  return(score_overlapping)
}


computeMeanSwitchIntensity <- function(overlaps_switches,
                                       directionality,
                                       bengidf){
  mean_intensity <- overlaps_switches %>% 
    group_by(ET) %>% 
    summarise(mean_switch_intensity = mean(directionality))
  mean_intensity$pair <- bengidf$pair[mean_intensity$ET]
  mean_intensity <- as.data.frame(mean_intensity)
  
  meandf <- merge(bengidf[,c("pair", "label")],
                mean_intensity[,c("pair", "mean_switch_intensity")],
                by.x = "pair", by.y = "pair")
  
  meandf[is.na(meandf)] <- 0
  return(meandf)
}

##----------------------------Start script------------------------------------##

args = commandArgs(trailingOnly=TRUE)
cat("Reading Data in\n")
cat(paste0("Path to BENGI is ",args[1], "\n"))
cat(paste0("Path to Directionality index is ", args[2], "\n"))
cat(paste0("Name of file is", args[3], "\n"))

resol <- args[4]
windowsize <- args[5]

cat("\n")

cat("Reading in all data ...")

dirIndex <- readDirectionality(args[2])
switches <- computeDirSwitches(dirIndex)

colnames(switches) <- c("chr", "start_switch_neg", "end_switch_neg",
                         "start_switch_pos", "end_switch_pos", "switch_intensity")

#making sure these columns are numeric
numcols <- c("start_switch_neg", "end_switch_neg",
             "start_switch_pos", "end_switch_pos", "switch_intensity")
switches[numcols] <- sapply(switches[numcols],as.numeric)



switchesRanges <- makeDirSwitchRanges(switches)

bengidf <- readBENGI(args[1])

bengidfRanges <- makeETRanges(bengidf)

cat("Computing overlaps\n")
overlaps_switches <- findETScoreOverlaps(bengidfRanges,
                                         switchesRanges,
                                         bengidf,
                                         switches)

bengidf$pair <- paste(bengidf$symbol38,
                      bengidf$gene_id1,
                      sep = "_")
cat("Computing mean Switch")
meandf <- computeMeanSwitchIntensity(overlaps_switches, switches, bengidf)

cat("Made it to write table")
rt <- "/project/CRUP_scores/CENTRE_HiC/Features/MeanSwitchIntensity/"
namefile <- paste0(rt, "meanSwitchIntensity", args[3], resol, windowsize, ".csv") 
write.table(meandf, 
            namefile, 
            row.names = F)
