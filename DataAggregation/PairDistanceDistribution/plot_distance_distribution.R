#########################################################################
#
# plot_distance_distribution: Use ggplot to plot the distance dist.
# between the pairs over all chromosomes and separately for each chrom
#
########################################################################

#---------------------------Libraries and functions--------------------#
library(ggplot2)

#----------------------------------------------------------------------#

#-----------------------------Script-----------------------------------#

#load the distance data 

data <- read.csv("/project/CRUP_scores/CENTRE_HiC/distance_distribution.tsv",
		head = T)
png("/project/CRUP_scores/CENTRE_HiC/distances_all_chrom_dist.png")
ggplot(data, aes(x=abs(distance)) + geom_histogram(binwidth=10)
dev.off()
#----------------------------------------------------------------------#
