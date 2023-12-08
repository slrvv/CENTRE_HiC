# Script to compute the MI and MSI with the 1kb 5kb resol

BENGIPATH="/project/CRUP_scores/toSara/BENGI_processed_datasets"
BENGIEND="-Benchmark.v38.txt"

ISPATH="/project/CRUP_scores/CENTRE_HiC/InsulationScores"

DIPATH="/project/CRUP_scores/CENTRE_HiC/DirectionalityIndex"

SAVEPATH="/project/CRUP_scores/CENTRE_HiC/WilcoxonTestResults"

#Rscript /project/CRUP_scores/CENTRE_HiC/DirectionalityIndexScripts/dirindex_switchintensity.R $BENGIPATH/K562.CRISPR$BENGIEND $DIPATH/K562/ENCFF080DPJ@5kb.directionality_15kb.bed \
#.CRISPR 5kb 15kb

# for WINDOW in  "25kb" "35kb" "50kb" "75kb"
# do
#   Rscript wilcoxtest_insscore.R $BENGIPATH/IMR90.HiC$BENGIEND $ISPATH/IMR-90/ENCFF636CEM@5kb.insulation_$WINDOW.bed \
#   WilcoxtestIMR90.HiC 5kb $WINDOW
# 
#   Rscript wilcoxtest_insscore.R $BENGIPATH/K562.CRISPR$BENGIEND $ISPATH/K562/ENCFF080DPJ@5kb.insulation_$WINDOW.bed \
#   WilcoxtestK562.CRISPR 5kb $WINDOW
# 
#   Rscript wilcoxtest_insscore.R $BENGIPATH/K562.HiC$BENGIEND $ISPATH/K562/ENCFF080DPJ@5kb.insulation_$WINDOW.bed \
#   WilcoxtestK562.HiC 5kb $WINDOW
# 
#   Rscript /project/CRUP_scores/CENTRE_HiC/DirectionalityIndexScripts/dirindex_switchintensity.R $BENGIPATH/IMR90.HiC$BENGIEND $DIPATH/IMR-90/ENCFF636CEM@5kb.directionality_$WINDOW.bed \
#   WilcoxtestIMR90.HiC 5kb $WINDOW
# 
#   Rscript /project/CRUP_scores/CENTRE_HiC/DirectionalityIndexScripts/dirindex_switchintensity.R $BENGIPATH/K562.CRISPR$BENGIEND $DIPATH/K562/ENCFF080DPJ@5kb.directionality_$WINDOW.bed \
#   WilcoxtestK562.CRISPR 5kb $WINDOW
# 
#   Rscript /project/CRUP_scores/CENTRE_HiC/DirectionalityIndexScripts/dirindex_switchintensity.R $BENGIPATH/K562.HiC$BENGIEND $DIPATH/K562/ENCFF080DPJ@5kb.directionality_$WINDOW.bed \
#   WilcoxtestK562.HiC 5kb $WINDOW
# 
# done

for SAMPLE in "GM12878.HiC"
do
  for WINDOW in "3kb" "5kb" "7kb"
  do
    FILE=$SAVEPATH/MinInsulationScore/minInsulationWilcoxtest$SAMPLE"1kb"$WINDOW.csv
    if [ -f "$FILE" ]; then
        echo "$FILE exists."
    else 
        Rscript wilcoxtest_insscore.R $BENGIPATH/$SAMPLE$BENGIEND \
        $ISPATH/GM12878/ENCFF237QCN@1kb.insulation_$WINDOW.bed \
        Wilcoxtest$SAMPLE 1kb $WINDOW
    fi
    
    FILE=$SAVEPATH/MeanSwitchIntensity/meanSwitchIntensityWilcoxtest$SAMPLE"1kb"$WINDOW.csv
    if [ -f "$FILE" ]; then
        echo "$FILE exists."
    else 
        Rscript /project/CRUP_scores/CENTRE_HiC/DirectionalityIndexScripts/dirindex_switchintensity.R \
        $BENGIPATH/$SAMPLE$BENGIEND $DIPATH/GM12878/ENCFF237QCN@1kb.directionality_$WINDOW.bed \
        Wilcoxtest$SAMPLE 1kb $WINDOW
    fi
    
  done

done
