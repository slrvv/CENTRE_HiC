################################################################################
#                                                                              #
# wilcoxtest_dirindex; Using the directionality index to figure a way          #
# of aggregating it over the ET window.                                        #
#                                                                              #
################################################################################

# We designed the minimum insulation score feature that takes the minimum of the
# directionality scores over the Enhancer Target window at 10kb resolution.
# We find out that there is a significant difference the min. directionality score
# distribution of active ET pairs vs. inactive ET pairs for most of the BENGI
# datasets.

#Now, we want to do the same for Directionality Index. The first idea that comes
#to mind is that we do not need all of the directionality indices. we just need
#to know where they swtich from negative to positive (right boundary of the TAD)
#so we should save that.
#Next, we should think on what to do if there are multiple switches in an ET window.

#------------------------Libraries and functions-------------------------------#

library(CENTRE)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(ggprism)

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
  dirswitch <- data.frame()
  for (i in 2:nrow(directionality)){
    if (directionality[i-1,5] < 0 & directionality[i,5] > 0 & directionality[i,2] == directionality[i-1,3]+1){
      dirswitch <- rbind(dirswitch, directionality[i-1,]) 
      dirswitch <- rbind(dirswitch, directionality[i,]) ## does this make sense??
    }
  }
  return(dirswitch)
}

makeDirSwitchRanges <- function(directionality){
  ##### Turn directionality switches dataframe into GRanges object ######
  
  directionalityRange <- with(directionality,
                          GenomicRanges::GRanges(V1,
                                                 IRanges::IRanges(start = V2,
                                                                  end = V3),
                                                 mcols = directionality))
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
  score_overlapping$directionality <- directionality$V5[score_overlapping$score]
  
  
  
  return(score_overlapping)
}

### mEH THINK ABOUT IT ONCE YOU ARENT SUPER HUNGRY
computeDirSwitchesCounts<- function(overlaps, BENGIdf){
  bins <- as.data.frame(table(cres_EP$between))
  bins$ID
  all_bins <- merge(bins, BENGIdf, by.x = "Var1", by.y = "Var1", all.x = TRUE)
  all_bins[is.na(all_bins)] <- 0
  
}

#-----------------------------Script Start-------------------------------------#
args = commandArgs(trailingOnly=TRUE)
cat("Reading Data in\n")
cat(paste0("Path to BENGI is ",args[1], "\n"))
cat(paste0("Path to Insulation score is ", args[2], "\n"))


#directionality index
directionalitydf <- readDirectionality(args[2])
dirdf_chr20 <- directionalitydf
dirswitches <- computeDirSwitches(dirdf_chr20)

dirswitchesRanges <- makeDirSwitchRanges(dirswitches)
dirRanges <- makeDirSwitchRanges(dirdf_chr20)

## Bengi datac
GM12878df <- readBENGI(args[1])
GM12878df_chr20 <- GM12878df
GM12878Ranges <- makeETRanges(GM12878df_chr20)

##Finding overlaps
overlaps_switches <- findETScoreOverlaps(GM12878Ranges,
                                dirswitchesRanges,
                                GM12878df_chr20,
                                dirswitches)

overlaps_dir <- findETScoreOverlaps(GM12878Ranges,
                                    dirRanges,
                                    GM12878df_chr20,
                                    dirdf_chr20)

GM12878df_chr20$pair <- paste(GM12878df_chr20$symbol38,
                              GM12878df_chr20$gene_id1,
                              sep = "_")
###Count number of switches
bins_switch <- as.data.frame(table(overlaps_switches$ET))
bins_switch$ID <- GM12878df_chr20[bins_switch$Var1,]$pair
all_bins_switch <- merge(bins_switch, GM12878df_chr20, by.x = "ID", by.y = "pair", all.x = TRUE)
all_bins_switch[is.na(all_bins_switch)] <- 0
head(all_bins_switch)


###Count total number of directionality maps
bins_total <- as.data.frame(table(overlaps_dir$ET))
bins_total$ID <- GM12878df_chr20[bins_total$Var1,]$pair
all_bins_total <- merge(bins_total, GM12878df_chr20, by.x = "ID", by.y = "pair", all.x = TRUE)
all_bins_total[is.na(all_bins_total)] <- 0
head(all_bins_total)

all_bins <- merge(all_bins_switch, all_bins_total[c("ID", "Freq")], by.x = "ID", by.y = "ID", all.x = TRUE)
head(all_bins)


all_bins$norm_switch <- all_bins$Freq.x/all_bins$Freq.y
write.csv(all_bins, paste0("/project/CRUP_scores/CENTRE_HiC/Features/NormCountsDirIndex/normCountSwitches",
                           args[3], ".csv"),
          row.names = F)

