################################################################################
#
# Looking at MI for different pairs
#
################################################################################

## Genome track figure to compare MI in positive vs. negative case.

imr90 <- "/project/CRUP_scores/CENTRE_HiC/Training/BENGI_MSI_MI_datasets/10Kb/IMR90.HiC-Benchmark.MI.MSI.v38.csv"
imr90dat <- read.table(imr90, header = T, sep=",")
imr90dat[is.na(imr90dat)] <- 0 
gm <- "/project/CRUP_scores/CENTRE_HiC/Training/BENGI_MSI_MI_datasets/10Kb/GM12878.RNAPII-ChIAPET-Benchmark.MI.MSI.v38.csv"
gmdat <- read.table(gm, header = T, sep=",")
gmdat[is.na(gmdat)] <- 0
gmdatpos <- gmdat[gmdat$label == 1,] 

imr90datpos <- imr90dat[imr90dat$label == 0,] 

overlapids <- merge(gmdatpos, 
                    imr90datpos, 
                    by = "pair")
overlapids_MI <- overlapids[overlapids$min_insulation.x > overlapids$min_insulation.y,]


overlapids_MI
overlapids_MSI <- overlapids[overlapids$mean_switch_intensity.x < overlapids$mean_switch_intensity.y,]
overlapids_MSI

overlapids_MSI <- overlapids_MSI[overlapids$min_insulation.x > overlapids$min_insulation.y,]
overlapids_MSI
overlapids <- overlapids[,c("pair", 
                            "gene_id.x", 
                            "gene_id.y", 
                            "symbol38.x",
                            "symbol38.y",
                            "min_insulation.x",
                            "min_insulation.y",
                            "mean_switch_intensity.x",
                            "mean_switch_intensity.y",
                            "RNA_seq.x",
                            "RNA_seq.y",)]

conn <- RSQLite::dbConnect(RSQLite::SQLite(), "/project/CRUP_scores/CENTRE/inst/extdata/Annotation.db")

#get chromosome and tts of our genes
query <- paste("SELECT  gene_id1, chr, start, end FROM gencode WHERE gene_id1 in (",
               paste0(sprintf("'%s'", overlapids_MI$gene_id1.x), collapse = ", "),")",sep="" )

gene <- RSQLite::dbGetQuery(conn, query)
gene <- gene[,c("chr", "start", "end")]
write.table(gene,
            "/project/CRUP_scores/CENTRE_HiC/FinalPlots/GenomeTrack/gene_candidates.db",
            quote = F,
            row.names = F)

#Select all of the annotation for ccres

query <- paste("SELECT  V5, V1, V2, V3 FROM ccres_enhancer WHERE V5 in (",
               paste0(sprintf("'%s'", overlapids_MI$symbol38.x), collapse = ", "),")",sep="" )
ccres_enhancer <- RSQLite::dbGetQuery
ccres_enhancer <- ccres_enhancer[, c("V1", "V2", "V3", "V5")] 
ccres_enhancer
write.table(ccres_enhancer,
            "/project/CRUP_scores/CENTRE_HiC/FinalPlots/GenomeTrack/enhancer_candidates.db",
            quote = F,
            row.names = F)


load("/project/CRUP_scores/CENTRE//sysdata.rda")

###we have the candidate pairs 
### For the figure we need :
# 1. Genomic positions of the pairs
# 2. RNA-seq
# 3. min_insulation scores in bigWig file format

