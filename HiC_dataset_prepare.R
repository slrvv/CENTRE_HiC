################################################################################
##  Hi-C data set preparation                                                                          ##
################################################################################
library(stringr)

# 1. Download metadata for all available HiC datasets in encode.
# 2. Compare with the 66 cell lines for which we have RNA-seq and Histone ChIP seq
# 3. Download the HiC data for the intersecting cell lines --> This is done
# with a python script 


HicMeta <- read.csv("/project/CRUP_scores/CENTRE_HiC/HiCmaps/all_HiC_encode_datasets.tsv",
                       sep = "\t", header = T) ## All of the HiC available in ENCODE

RNAseqSet <- read.csv("/project/CRUP_scores/CENTRE_HiC/Rna_seq_metadata.csv",
                       sep = "\t", header = T)
##RNaseqSet is the metadata of the samples for which we have RNAseq and Histone
##ChIP seq

## Formating Biosample name so that we can compare with RNAseqSet
head(HicMeta)

HicMeta$Biosample.term.name <- str_replace_all(HicMeta$Biosample.term.name, " ", "-") 
HicMeta$Biosample.term.name <- str_replace_all(HicMeta$Biosample.term.name, ",", "") 

HiC_samples <- HicMeta$Biosample.term.name
HiC_samples <- str_replace_all(HiC_samples , "-", "")
HiC_samples <- tolower(HiC_samples)

HicMeta$Formatted.name <- HiC_samples

RNAseqSet_samples <- str_replace_all(RNAseqSet$Name.of.sample, "-", "") 
RNAseqSet_samples <- tolower(RNAseqSet_samples)

overlap <- unique(HiC_samples[HiC_samples %in% RNAseqSet_samples])

length(overlap)
## Of all the HiC data in ENCODE 37 correspond to samples we have RNA-seq and 
## Histone ChIP seq on 

subsetted_HiC <- HicMeta[HicMeta$Formatted.name %in% overlap,]
subsetted_HiC <- subsetted_HiC[!duplicated(subsetted_HiC$Formatted.name),]
### Save the data
head(subsetted_HiC)
subsetted_HiC <- subsetted_HiC[, c("Accession", "Assay.title",
                                   "Biosample.classification",
                                   "Biosample.summary",
                                   "Biosample.term.name",
                                   "Lab")]
write.csv(subsetted_HiC, "/project/CRUP_scores/CENTRE_HiC/subsetHiC_summary.tsv",
          sep = "\t", row.names = F)


