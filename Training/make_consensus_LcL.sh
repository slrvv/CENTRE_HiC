##################################################################################################
#                                                                                                #
# make_consensus_LcL: Calculate the minimum insulation score and DI measures on the consensus    #
# LcL dataset.                                                                                   #
#                                                                                                #
##################################################################################################

# Calculate the newly designed features for the consesus Lcl dataset in order to build the training
# dataset.

DIRLCL=/project/CRUP_scores/CENTRE_HiC/Training/CENTRE_final_training/consensus_LCLs.txt
INS=/project/CRUP_scores/CENTRE_HiC/InsulationScoreScripts/wilcoxtest_insscore.R
DI=/project/CRUP_scores/CENTRE_HiC/DirectionalityIndexScripts/wilcoxtest_dirindex.R
DII=/project/CRUP_scores/CENTRE_HiC/DirectionalityIndexScripts/dirindex_switchintensity.R

FILEINS=/project/CRUP_scores/CENTRE_HiC/InsulationScores/GM12878/ENCFF237QCN.hic@10kb.insulation_30kb.bed
FILEDI=/project/CRUP_scores/CENTRE_HiC/DirectionalityIndex/GM12878/ENCFF237QCN.hic@10kb.directionality_30kb.bed

MIROOT=/project/CRUP_scores/CENTRE_HiC/WilcoxonTestResults/MinInsulationScore
MSIROOT=/project/CRUP_scores/CENTRE_HiC/WilcoxonTestResults/MeanSwitchIntensity
SAVEROOT=/project/CRUP_scores/CENTRE_HiC/Training/BENGI_MSI_MI_datasets

#echo "Minimum Insulation Score"
#Rscript $INS $DIRLCL $FILEINS consensusLcL 10kb 30kb

#echo "Mean Switch Intensity  Directionality Index"
#Rscript $DII $DIRLCL $FILEDI consensusLcL 10kb 30kb


echo "Merging all files"
Rscript /project/CRUP_scores/CENTRE_HiC/Training/CENTRE_final_training/make_HiC_BENGI_Benchmarks.R $DIRLCL \
$MIROOT/minInsulationconsensusLcL10kb30kb.csv \
$MSIROOT/meanSwitchIntensityconsensusLcL10kb30kb.csv \
$SAVEROOT/consensusLcL.MI.MSI.v38.10kb30kb.csv

