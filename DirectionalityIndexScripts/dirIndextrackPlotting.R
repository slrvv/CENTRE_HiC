################################################################################
#
# Example to visualize Directionality switches
#
################################################################################

# Try to do a visualization of the MALAT1 gene that includes the following tracks
# 1st a track with the chromosome
# 2nd a track with the genes and enhancers
# 3rd a track with the HiC map
# 4th a track of the directionality index
# 5th a track with the annotated switches in directionality

install.packages("ggcoverage")
library("rtracklayer")
library("ggcoverage")
library("ggpattern")
library("idr2d")

hicpath <- "/project/CRUP_scores/CENTRE_HiC/DirectionalityIndexScripts/hicmap.txt"
hicmap <- read.table(file = hicpath, sep = " ")
hicmap <- as.matrix(hicmap)

ggcoverage() + geom_tad(matrix = hicmap, value.cut = 0.99,
         color.palette = "viridis",
         top = FALSE, show.rect = TRUE)
