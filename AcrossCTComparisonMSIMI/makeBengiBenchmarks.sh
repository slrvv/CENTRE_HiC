#-------------------------------------------------------------------------------
# Author:  Sara Lopez Ruiz de Vargas
# Email:  lopez_s@molgen.mpg.de
#
# Date:    2024-04-08
#
# Script Name:  makeBengiBenchmarks.sh
#
# Script Description: make the bengi Benchmarks for the datasets with MI MSI 
# computed across dataset
#
#-------------------------------------------------------------------------------

#---------------------------------Paths----------------------------------------#
BENGIPATH="/project/CRUP_scores/toSara/BENGI_processed_datasets"
BENGIEND="-Benchmark.v38.txt"
SUMMARY="/project/CRUP_scores/CENTRE_HiC/AcrossCTComparisonMSIMI/AcrossCTMetadata.csv"
MSIPATH="/project/CRUP_scores/CENTRE_HiC/Features/MeanSwitchIntensityDiffCT"
MIPATH="/project/CRUP_scores/CENTRE_HiC/Features/MinInsulationScoreDiffCT"
IS=/project/CRUP_scores/CENTRE_HiC/InsulationScoreScripts/wilcoxtest_insscore.R
DI=/project/CRUP_scores/CENTRE_HiC/DirectionalityIndexScripts/dirindex_switchintensity.R
MERGE=/project/CRUP_scores/CENTRE_HiC/Training/CENTRE_final_training/make_HiC_BENGI_Benchmarks.R
MIMSIEND=.MIMSI.csv
SAVEPATH=/project/CRUP_scores/CENTRE_HiC/AcrossCTComparisonMSIMI/BENGI_MI_MSI_acrossCT

#summary file 
#BengiSample	HiCSample	JoinedName	FileIS	FileDI

while read line ; do
     set $line
     IFS=$','; split=($line); unset IFS;
     
     Rscript $MERGE $BENGIPATH/${split[0]}$BENGIEND \
     $MIPATH/minInsulation${split[2]}.csv $MSIPATH/meanSwitchIntensity${split[2]}.csv\
     $SAVEPATH/${split[2]}$MIMSIEND

done < <(tail -n +2 $SUMMARY)