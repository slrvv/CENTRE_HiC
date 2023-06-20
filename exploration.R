###############################################################
# 
# Exploration: Script to check out diff packages to compute
# insulation score and directionality index
#
###############################################################

# First package I'm looking at is GENOVA. Looks promising 
# but is not on Bioconductor or CRAN :C

# install.packages("remotes")
# Install from github
#remotes::install_github("robinweide/GENOVA")

# Loading cooler format (.hic) contact matrix

GM12878_contacts <- GENOVA::load_contacts(
	signal_path = "/project/CRUP_scores/CENTRE_HiC/HiCmaps/GM12878/ENCFF256UOW.hic",
	sample_name = "GM12878",
	resolution = 10e3, 
	colour = black) 
)

## Making a matrix plot of chr11
png(file = "/project/CRUP_scores/CENTRE_HiC/GM12878map.png", 
	width = 600,
	height = 350 )
GENOVA::hic_matrixplot(exp1 = GM12878,
		chrom = 'chr11',
		start = 25e6,
		end=30e6,
		cut.off = 50)
dev.off()
