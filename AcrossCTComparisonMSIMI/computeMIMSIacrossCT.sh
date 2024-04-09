#-------------------------------------------------------------------------------
# Author:  Sara Lopez Ruiz de Vargas
# Email:  lopez_s@molgen.mpg.de
#
# Date:    2024-04-05
#
# Script Name:  computeMIMSIacrossCT.sh
#
# Script Description: Compute the MI and MSI for BENGI sets of diff CT
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

#-----------------------------Script-------------------------------------------#

#summary file 
#BengiSample	HiCSample	JoinedName	FileIS	FileDI

while read line ; do
     set $line
     IFS=$','; split=($line); unset IFS;
     echo "Computing MI"
     Rscript $IS $BENGIPATH/${split[0]}$BENGIEND \
     ${split[3]} $MIPATH/minInsulation${split[2]}.csv 10kb 30kb
     
     echo "computing MSI"
     Rscript $DI $BENGIPATH/${split[0]}$BENGIEND \
     ${split[4]} $MSIPATH/meanSwitchIntensity${split[2]}.csv 10kb 30kb
     
     echo " "
done < <(tail -n +2 $SUMMARY)