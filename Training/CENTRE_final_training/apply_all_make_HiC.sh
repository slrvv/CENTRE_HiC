#### Apply the merging script to all BENGI

BENGIROOT=/project/CRUP_scores/toSara/BENGI_processed_datasets
MIROOT=/project/CRUP_scores/CENTRE_HiC/WilcoxonTestResults/MinInsulationScore
MSIROOT=/project/CRUP_scores/CENTRE_HiC/WilcoxonTestResults/MeanSwitchIntensity
SAVEROOT=/project/CRUP_scores/CENTRE_HiC/Training/BENGI_MSI_MI_datasets

#Colon
#Rscript make_HiC_BENGI_Benchmarks.R $BENGIROOT/Colon.GTEx-Benchmark.v38.txt \
#$MIROOT/minInsulationWilcoxtestColon.GTEx.csv $MSIROOT/meanSwitchIntensityWilcoxtestColon.GTEx.csv \
#$SAVEROOT/Colon.GTEx-Benchmark.MI.MSI.v38.csv

#GMCTCF
#Rscript make_HiC_BENGI_Benchmarks.R $BENGIROOT/GM12878.CTCF-ChIAPET-Benchmark.v38.txt \
#$MIROOT/minInsulationWilcoxtestGM12878.CTCF-ChIAPET10kb.csv $MSIROOT/meanSwitchIntensityWilcoxtestGM12878.CTCF-ChIAPET.csv \
#$SAVEROOT/GM12878.CTCF-ChIAPET-Benchmark.MI.MSI.v38.csv

#GMGEUVADIS
Rscript make_HiC_BENGI_Benchmarks.R $BENGIROOT/GM12878.GEUVADIS-Benchmark.v38.txt \
$MIROOT/minInsulationWilcoxtestGM12878.GEUVADIS10kb.csv \
$MSIROOT/meanSwitchIntensityWilcoxtestGM12878.GEUVADIS.csv \
$SAVEROOT/GM12878.GEUVADIS-Benchmark.MI.MSI.v38.csv
#GMHiC
#Rscript make_HiC_BENGI_Benchmarks.R $BENGIROOT/GM12878_HiC-Benchmark.v38.txt \
#$MIROOT/minInsulationWilcoxtestGM12878.HiC10kb.csv $MSIROOT/meanSwitchIntensityWilcoxtestGM12878.HiC.csv \
#$SAVEROOT/GM12878.HiC-Benchmark.MI.MSI.v38.csv
#GMRNAPII
Rscript make_HiC_BENGI_Benchmarks.R $BENGIROOT/GM12878.RNAPII-ChIAPET-Benchmark.v38.txt \
$MIROOT/minInsulationWilcoxtestGM12878.RNAPII-ChIAPET10kb.csv \
$MSIROOT/meanSwitchIntensityWilcoxtestGM12878.RNAPII-ChIAPET.csv \
$SAVEROOT/GM12878.RNAPII-ChIAPET-Benchmark.MI.MSI.v38.csv
#GMcHiC
#Rscript make_HiC_BENGI_Benchmarks.R $BENGIROOT/GM12878.CHiC-Benchmark.v38.txt \
#$MIROOT/minInsulationWilcoxtestGM12878.CHiC10kb.csv $MSIROOT/meanSwitchIntensityWilcoxtestGM12878.CHiC.csv \
#$SAVEROOT/GM12878.CHiC-Benchmark.MI.MSI.v38.csv
#IMR90
#Rscript make_HiC_BENGI_Benchmarks.R $BENGIROOT/IMR90.HiC-Benchmark.v38.txt \
#$MIROOT/minInsulationWilcoxtestIMR90.HiC.csv $MSIROOT/meanSwitchIntensityWilcoxtestIMR90.HiC.csv \
#$SAVEROOT/IMR90.HiC-Benchmark.MI.MSI.v38.csv
#K562
#Rscript make_HiC_BENGI_Benchmarks.R $BENGIROOT/K562.HiC-Benchmark.v38.txt \
#$MIROOT/minInsulationWilcoxtestK562.HiC.csv $MSIROOT/meanSwitchIntensityWilcoxtestK562.HiC.csv \
#$SAVEROOT/K562.HiC-Benchmark.MI.MSI.v38.csv

#Rscript make_HiC_BENGI_Benchmarks.R $BENGIROOT/K562.CRISPR-Benchmark.v38.txt \
#$MIROOT/minInsulationWilcoxtestK562.CRISPR.csv $MSIROOT/meanSwitchIntensityWilcoxtestK562.CRISPR.csv \
#$SAVEROOT/K562.CRISPR-Benchmark.MI.MSI.v38.csv

#Ovary
#Rscript make_HiC_BENGI_Benchmarks.R $BENGIROOT/Ovary.GTEx-Benchmark.v38.txt \
#$MIROOT/minInsulationWilcoxtestOvary.GTEx.csv $MSIROOT/meanSwitchIntensityWilcoxtestOvary.GTEx.csv \
#$SAVEROOT/Ovary.GTEx-Benchmark.MI.MSI.v38.csv

#Pancreas
#Rscript make_HiC_BENGI_Benchmarks.R $BENGIROOT/Pancreas.GTEx-Benchmark.v38.txt \
#$MIROOT/minInsulationWilcoxtestPancreas.GTEx.csv $MSIROOT/meanSwitchIntensityWilcoxtestPancreas.GTEx.csv \
#$SAVEROOT/Pancreas.GTEx-Benchmark.MI.MSI.v38.csv

#Stomach
#Rscript make_HiC_BENGI_Benchmarks.R $BENGIROOT/Stomach.GTEx-Benchmark.v38.txt \
#$MIROOT/minInsulationWilcoxtestStomach.GTEx.csv $MSIROOT/meanSwitchIntensityWilcoxtestStomach.GTEx.csv \
#$SAVEROOT/Stomach.GTEx-Benchmark.MI.MSI.v38.csv
