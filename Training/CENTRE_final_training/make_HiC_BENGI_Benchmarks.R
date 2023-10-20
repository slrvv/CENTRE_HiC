#########################################################################
#
# make_HiC_BENGI_Benchmarks: Add the MI and MSI to the BENGI datasets
#
########################################################################


# We add all of the computed features to the BENGI benchmarked datasets

merge_MI_MSI <- function(bengipath, MIpath, MSIpath, savepath){
  bengi <- read.table(bengipath, header = T)
  correct_rows <- nrow(bengi)
  MI <- read.table(MIpath,
                   header = T,
                   sep = ",")

  MSI <- read.table(MSIpath, header = T)

  MI$pair <- paste(MI$symbol38, MI$gene_id1, sep = "_")
  MI <- MI[c("pair", "score_min")]
  colnames(MI) <- c("pair", "min_insulation")
  BENGI.MI <- merge(bengi, MI,
                    by.x = "pair", by.y = "pair", all.x = T)

  BENGI.MI.MSI <- merge(BENGI.MI, MSI[c("pair", "mean_switch_intensity")],
                              by.x = "pair", by.y = "pair", all.x = T)
  BENGI.MI.MSI$min_insulation[is.na(BENGI.MI.MSI$min_insulation)] <- 0
  BENGI.MI.MSI$min_insulation[is.infinite(BENGI.MI.MSI$min_insulation)] <- max(BENGI.MI.MSI$min_insulation[!is.infinite(BENGI.MI.MSI$min_insulation)])
  write.csv(BENGI.MI.MSI, savepath, row.names = F)
  print(paste0("Number of rows is correct: ", correct_rows == nrow(BENGI.MI.MSI)))
  return(0)
}


###### Argument passing ##

args = commandArgs(trailingOnly=TRUE)
cat("Reading Data in\n")
cat(paste0("Path to BENGI is ",args[1], "\n"))
cat(paste0("Path to MI is", args[2], "\n"))
cat(paste0("Path to MSI is", args[3], "\n"))
cat(paste0("Path to save file is", args[4], "\n"))
merge_MI_MSI(args[1], args[2], args[3], args[4])
