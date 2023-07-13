################################################################################
#                                                                              
# Pair length distribution BENGI: Calculate the pair size of active BENGI ETs
#
################################################################################

# In order to find an appropiate window size for the indices we want to plot
# the distribution of the length of positive pairs of BENGI annotated datasets

#---------------Libraries and functions----------------------------------------#

library(ggplot2)

producePlot <- function(data, dataSetName){
  png(paste0("/project/CRUP_scores/CENTRE_HiC/DataAggregation/PairDistanceDistribution/",
             dataSetName, ".png"))
  print({
  ggplot(data, aes(x=abs(distance))) + geom_histogram()
  })
  dev.off()
  return()
}

bengiDistribution <- function(directory, dataSetName){
  bengi_pairs <- read.table(directory, 
                            header = T)
  #take only BENGI positive pairs
  bengi_pos <- bengi_pairs[bengi_pairs$label == 1, ]
  producePlot(bengi_pos, dataSetName)
  cat(paste0("Mean distance is ", mean(abs(bengi_pos$distance))))
  cat("\n")
  return()
}

#------------------------------------------------------------------------------#

#------------------------Script start------------------------------------------#


path <- "/project/CRUP_scores/toSara/BENGI_processed_datasets/"
bash_comm <- paste("ls", path, sep =" ")
samples <- system(bash_comm, intern=T)
for (i in 1:length(samples)){
  cat(samples[i])
  cat("\n")
  directory <- paste0(path, samples[i])
  bengiDistribution(directory, samples[i])
}



#------------------------------------------------------------------------------#

