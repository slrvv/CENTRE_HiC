#-------------------------------------------------------------------------------
# Author:  Sara Lopez Ruiz de Vargas
# Email:  lopez_s@molgen.mpg.de
#
# Date:    2024-04-05
#
# Script Name:  
#
# Script Description: We want to calculate MI and MSI for the BENGI ET pairs 
# using Insulation scores and Directionality indices from other Samples. For this
# we want to create a file for which samples we want to compare
#
#-------------------------------------------------------------------------------


# It doesn't make sense to compare cell types and tissues.
pathIS <- "/project/CRUP_scores/CENTRE_HiC/InsulationScores/"
pathDI <- "/project/CRUP_scores/CENTRE_HiC/DirectionalityIndex/"
pathBENGI <- "/project/CRUP_scores/CENTRE_HiC/Training/BENGI_MSI_MI_datasets/10Kb/"

tissuesHiC <- c("Colon", "Ovary", "Pancreas", "Stomach")
tissuesIS <- c("Sigmoid-Colon", "Ovary", "Pancreas", "Stomach")
tissueslower <- tolower(tissuesIS)
#tissues IS and tissues lower are just used to search through DI and IS
#directories
tissuesBengi <- c("Colon.GTEx", 
                  "Ovary.GTEx",
                  "Pancreas.GTEx",
                  "Stomach.GTEx")

##for each Bengi dataset we want to use DI and IS coming from different tissues/celltypes


BengiSample <- c()
HiCSample <- c()
ConjointName <- c()

for (bengi in tissuesBengi){
  for (hic in tissueslower) {
    namebengi <- strsplit(bengi, ".", fixed = T)
    namebengi <- sapply(namebengi,"[[",1)
    if (namebengi != hic){
      BengiSample <- c(BengiSample, bengi)
      HiCSample <- c(HiCSample, hic)
    }
  }
}

datacompare <- as.data.frame(cbind(BengiSample, HiCSample))
datacompare$JoinedName <- paste(datacompare$BengiSample, 
                                datacompare$HiCSample, 
                                sep = "-")

datacompare
##get the path for the IS file and DI for each tissue
fileIS <- c()
fileDI <- c()
for (tissue in datacompare$HiCSample){
  pathtofile <- paste0(pathIS, tissue, "/")
  file <- system(paste0('find ', pathtofile, ' -name "*10kb.insulation_30kb.bed"'), intern = T)
  fileIS <- c(fileIS, file)
  pathtofile <- paste0(pathDI, tissue, "/")
  file <- system(paste0('find ', pathtofile, ' -name "*directionality_30kb.bed"'), intern = T)
  fileDI <- c(fileDI, file)
}

datacompare$FileIS <- fileIS

datacompare$FileDI <- fileDI

##Now the same thing for cell types and then we cbind it

celltypeHiC<- c("GM12878", "IMR-90", "K562")

celltypeBengi <- c("GM12878.CHiC", 
                   "GM12878.CTCF-ChIAPET",
                   "GM12878.GEUVADIS",       
                   "GM12878.HiC",            
                   "GM12878.RNAPII-ChIAPET",
                   "IMR90.HiC",
                   "K562.CRISPR",
                   "K562.HiC")

BengiSample <- c()
HiCSample <- c()
ConjointName <- c()

for (bengi in celltypeBengi){
  for (hic in celltypeHiC) {
    namebengi <- strsplit(bengi, ".", fixed = T)
    namebengi <- sapply(namebengi,"[[",1)
    if (namebengi != hic){
      BengiSample <- c(BengiSample, bengi)
      HiCSample <- c(HiCSample, hic)
    }
  }
}

datacompare1 <- as.data.frame(cbind(BengiSample, HiCSample))
datacompare1$JoinedName <- paste(datacompare1$BengiSample, 
                                datacompare1$HiCSample, 
                                sep = "-")


##get the path for the IS file and DI for each tissue
fileIS <- c()
fileDI <- c()
for (cell in datacompare1$HiCSample){
  pathtofile <- paste0(pathIS, cell, "/")
  file <- system(paste0('find ', pathtofile, ' -name "*10kb.insulation_30kb.bed"'), intern = T)
  fileIS <- c(fileIS, file[1])
  pathtofile <- paste0(pathDI, cell, "/")
  file <- system(paste0('find ', pathtofile, ' -name "*directionality_30kb.bed"'), intern = T)
  fileDI <- c(fileDI, file[1])
}


datacompare1$FileIS <- fileIS

datacompare1$FileDI <- fileDI

datacompare <- rbind(datacompare, datacompare1)

write.csv(datacompare, 
          "/project/CRUP_scores/CENTRE_HiC/AcrossCTComparisonMSIMI/AcrossCTMetadata.csv",
          quote = F ,
          row.names = F)
