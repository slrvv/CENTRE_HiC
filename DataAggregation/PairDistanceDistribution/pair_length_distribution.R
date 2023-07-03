################################################################################
#                                                                              
# Pair length distribution: Calculate the pair size of all possible length using 
# CENTRE::pairs()
#
################################################################################

# In order to find an appropiate window size for the indices we want to plot
# the distribution of the length of the pairs created by CENTRE::pairs()

#---------------Libraries and functions----------------------------------------#

library(CENTRE)


#function to compute the distance from pairs dataframe
computeDistances <- function(x) {
  ## Computing the distance features

  
  colnames(x) <- c("gene_id", "enhancer_id")
  # connect to annotation dataBase
  conn <- RSQLite::dbConnect(RSQLite::SQLite(),
                             system.file("extdata",
                                         "Annotation.db",
                                         package = "CENTRE"))
  #get chromosome and tts of our genes
  query <- paste("SELECT  gene_id1, chr, transcription_start FROM gencode WHERE gene_id1 in (",
                 paste0(sprintf("'%s'", x$gene_id), collapse = ", "),")",sep="" )
  gencode <- RSQLite::dbGetQuery(conn, query)
  
  #get chr and middle point of enhancers
  query_enh <-  paste("SELECT  V5, V1, middle_point FROM ccres_enhancer WHERE V5 in (",
                      paste0(sprintf("'%s'", x$enhancer_id), collapse = ", "),")",sep="" )
  ccres_enhancer <- RSQLite::dbGetQuery(conn, query_enh)
  RSQLite::dbDisconnect(conn)
  
  #Get the chr gene_id and transcription_start from gencode annotation
  #Getting the chrosomes and the middle points for the provided enhancers
  result <- merge(x,
                  ccres_enhancer[, c("V1", "V5", "middle_point")],
                  by.x = "enhancer_id",
                  by.y = "V5") #change V1 and V5 to more meaningful names
  #Getting the chrosomes and transcription start sites for the provided genes
  
  result <- merge(result,
                  gencode[, c("chr", "gene_id1", "transcription_start")],
                  by.x = "gene_id",
                  by.y = "gene_id1")
  print(head(result))
  
  cat("Removing all gene enhancer pairs that are not in the same chromosome.\n")
  
  result <- result[(result$V1 == result$chr), ]
  
  result$distance <- result$middle_point - result$transcription_start
  
  return(result)
}
#------------------------------------------------------------------------------#

#------------------------Script start------------------------------------------#

#connect to the gencode annotation database in the CENTRE package.
conn <- RSQLite::dbConnect(RSQLite::SQLite(),
                           system.file("extdata",
                                       "Annotation.db",
                                       package = "CENTRE"))
#get ALL gene_ids and chr

query <- paste("SELECT  gene_id1, chr FROM gencode")

gene_data <- RSQLite::dbGetQuery(conn, query)

gene_data <- gene_data[!gene_data$chr == "chrM",] ##REMOVE mitochondrial dna


genes <- as.data.frame(gene_data$gene_id1)

##compute the possible pairs for all gencode genes
pairs <- CENTRE::createPairs(genes)
print(head(pairs))
##compute the distances using the distance function (using the one CENTRE 
##internally uses but is not user facing)
distances <- computeDistances(pairs)

print(head(distances))

distance_distribution <- merge(distances, 
                               gene_data,
                               by.x = "gene_id",
                               by.y = "gene_id1")

print(head(distance_distribution))

write.csv(distance_distribution,
          "/project/CRUP_scores/CENTRE_HiC/distance_distribution.tsv",
          sep = "\t",
          row.names = F)
#------------------------------------------------------------------------------#

