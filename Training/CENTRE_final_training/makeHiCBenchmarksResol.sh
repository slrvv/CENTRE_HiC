###############################################################################
#
# make Hic benchmarks for 1kb 5kb resol
#
###############################################################################

BENGIROOT=/project/CRUP_scores/toSara/BENGI_processed_datasets
MIROOT=/project/CRUP_scores/CENTRE_HiC/WilcoxonTestResults/MinInsulationScore
MSIROOT=/project/CRUP_scores/CENTRE_HiC/WilcoxonTestResults/MeanSwitchIntensity
SAVEROOT=/project/CRUP_scores/CENTRE_HiC/Training/BENGI_MSI_MI_datasets


# for WINDOW in 3kb 5kb 7kb 10kb 15kb
# do
# 
# #GMCTCF
#   Rscript make_HiC_BENGI_Benchmarks.R $BENGIROOT/GM12878.CTCF-ChIAPET-Benchmark.v38.txt \
#   $MIROOT/minInsulationWilcoxtestGM12878.CTCF-ChIAPET1kb$WINDOW.csv $MSIROOT/meanSwitchIntensityWilcoxtestGM12878.CTCF-ChIAPET1kb$WINDOW.csv \
#   $SAVEROOT/GM12878.CTCF-ChIAPET-Benchmark.MI.MSI.v38.1kb$WINDOW.csv
# 
# #GMGEUVADIS
#   Rscript make_HiC_BENGI_Benchmarks.R $BENGIROOT/GM12878.GEUVADIS-Benchmark.v38.txt \
#   $MIROOT/minInsulationWilcoxtestGM12878.GEUVADIS1kb$WINDOW.csv \
#   $MSIROOT/meanSwitchIntensityWilcoxtestGM12878.GEUVADIS1kb$WINDOW.csv \
#   $SAVEROOT/GM12878.GEUVADIS-Benchmark.MI.MSI.v38.1kb$WINDOW.csv
# #GMHiC
#   Rscript make_HiC_BENGI_Benchmarks.R $BENGIROOT/GM12878_HiC-Benchmark.v38.txt \
#   $MIROOT/minInsulationWilcoxtestGM12878.HiC1kb$WINDOW.csv $MSIROOT/meanSwitchIntensityWilcoxtestGM12878.HiC1kb$WINDOW.csv \
#   $SAVEROOT/GM12878.HiC-Benchmark.MI.MSI.v38.1kb$WINDOW.csv
# #GMRNAPII
#   
#   Rscript make_HiC_BENGI_Benchmarks.R $BENGIROOT/GM12878.RNAPII-ChIAPET-Benchmark.v38.txt \
#   $MIROOT/minInsulationWilcoxtestGM12878.RNAPII-ChIAPET1kb$WINDOW.csv \
#   $MSIROOT/meanSwitchIntensityWilcoxtestGM12878.RNAPII-ChIAPET1kb$WINDOW.csv \
#   $SAVEROOT/GM12878.RNAPII-ChIAPET-Benchmark.MI.MSI.v38.1kb$WINDOW.csv
# #GMcHiC
#   Rscript make_HiC_BENGI_Benchmarks.R $BENGIROOT/GM12878.CHiC-Benchmark.v38.txt \
#   $MIROOT/minInsulationWilcoxtestGM12878.CHiC1kb$WINDOW.csv $MSIROOT/meanSwitchIntensityWilcoxtestGM12878.CHiC1kb$WINDOW.csv \
#   $SAVEROOT/GM12878.CHiC-Benchmark.MI.MSI.v38.1kb$WINDOW.csv
#   
# done 

for WINDOW in 15kb 25kb 35kb 50kb 75kb
do
  #IMR90
  # Rscript make_HiC_BENGI_Benchmarks.R $BENGIROOT/IMR90.HiC-Benchmark.v38.txt \
  # $MIROOT/minInsulationWilcoxtestIMR90.HiC5kb$WINDOW.csv $MSIROOT/meanSwitchIntensityWilcoxtestIMR90.HiC5kb$WINDOW.csv \
  # $SAVEROOT/IMR90.HiC-Benchmark.MI.MSI.v38.5kb$WINDOW.csv
  # 
  # #K562
  # Rscript make_HiC_BENGI_Benchmarks.R $BENGIROOT/K562.HiC-Benchmark.v38.txt \
  # $MIROOT/minInsulationWilcoxtestK562.HiC5kb$WINDOW.csv $MSIROOT/meanSwitchIntensityWilcoxtestK562.HiC5kb$WINDOW.csv \
  # $SAVEROOT/K562.HiC-Benchmark.MI.MSI.v38.5kb$WINDOW.csv
  
  Rscript make_HiC_BENGI_Benchmarks.R $BENGIROOT/K562.CRISPR-Benchmark.v38.txt \
  $MIROOT/minInsulationWilcoxtestK562.CRISPR5kb$WINDOW.csv $MSIROOT/meanSwitchIntensityWilcoxtestK562.CRISPR5kb$WINDOW.csv \
  $SAVEROOT/K562.CRISPR-Benchmark.MI.MSI.v38.5kb$WINDOW.csv
done

