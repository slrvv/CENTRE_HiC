BENGIPATH="/project/CRUP_scores/toSara/BENGI_processed_datasets"
BENGIEND="-Benchmark.v38.txt"

ISPATH="/project/CRUP_scores/CENTRE_HiC/InsulationScores"

DIPATH="/project/CRUP_scores/CENTRE_HiC/DirectionalityIndex"

SAVEPATH="/project/CRUP_scores/CENTRE_HiC/WilcoxonTestResults"
for SAMPLE in "GM12878.HiC"
do
for WINDOW in "3kb" "5kb" "7kb"
do
FILE=$SAVEPATH/MinInsulationScore/minInsulationWilcoxtest$SAMPLE"1kb"$WINDOW.csv
if [ -f "$FILE" ]; then
echo "$FILE exists."
else 
  Rscript /project/CRUP_scores/CENTRE_HiC/InsulationScoreScripts/wilcoxtest_insscore.R $BENGIPATH/$SAMPLE$BENGIEND \
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