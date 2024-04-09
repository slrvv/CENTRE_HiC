#-------------------------------------------------------------------------------
# Author:  Sara Lopez Ruiz de Vargas
# Email:  lopez_s@molgen.mpg.de
#
# Date:    2024-04-03
#
# Script Name:/project/CRUP_scores/CENTRE_HiC/Training/CENTRE_final_training/acrossCTcomparisonGetChrom.R
#
# Script Description: Comparison with MSIMI accross CT. Extract the chr for
# each of the enhancers
#-------------------------------------------------------------------------------

library(devtools)
load_all("/project/CRUP_scores/CENTRE/")
samplelisthic <- c("Colon.GTEx", "GM12878.CHiC", "GM12878.CTCF-ChIAPET",
                  "GM12878.GEUVADIS", "GM12878.HiC", "GM12878.RNAPII-ChIAPET",
                  "IMR90.HiC", "K562.CRISPR", "K562.HiC", "Ovary.GTEx",
                  "Pancreas.GTEx", "Stomach.GTEx")

samplelistbengi <- c("Colon.GTEx", "GM12878.CHiC", "GM12878.CTCF-ChIAPET", 
                    "GM12878.GEUVADIS", "GM12878.HiC", "GM12878.RNAPII-ChIAPET",
                    "IMR90.HiC", "K562.CRISPR", "K562.HiC", "Ovary.GTEx",
                    "Pancreas.GTEx", "Stomach.GTEx", "Liver.GTEx", "NHEK.HiC", 
                    "HeLa.CTCF-ChIAPET", "HeLa.HiC", "HeLa.RNAPII-ChIAPET")

for (samplebengi in samplelistbengi){
  for (samplehic in samplelisthic){
    if (samplebengi != samplehic){
      path <- paste0('/project/CRUP_scores/CENTRE_HiC/Training/BENGI_MSI_MI_acrossCT/',
                     samplebengi,
                     'BENGI_',
                     samplehic,
                     'MIMSI.csv')
      data <- read.table(path, header = T, sep = " ")
      print(path)
      print(head(data))
      
      conn <- RSQLite::dbConnect(RSQLite::SQLite(),
                                 system.file("extdata",
                                             "Annotation.db",
                                             package = "CENTRE"))
      #get chromosome for the genes we have
      query <- paste("SELECT  gene_id1, chr FROM gencode WHERE gene_id1 in (",
                     paste0(sprintf("'%s'", data$gene_id1), collapse = ", "),")",sep="" )
      gencode <- RSQLite::dbGetQuery(conn, query)
      RSQLite::dbDisconnect(conn)
      #Get the chr gene_id and transcription_start from gencode annotation
      #Getting the chrosomes and the middle points for the provided enhancers
      
      #Getting the chrosomes and transcription start sites for the provided genes
      result <- merge(data,
                      gencode[, c("chr", "gene_id1")],
                      by.x = "gene_id1",
                      by.y = "gene_id1")
      write.table(result, path, row.names=F, quote = F)
    }
  }
}
