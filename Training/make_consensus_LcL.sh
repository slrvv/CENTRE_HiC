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

FILEINS=/project/CRUP_scores/CENTRE_HiC/InsulationScores/ENCFF237QCN.hic@1kb.insulation_3kb.bed
FILEDI=/project/CRUP_scores/CENTRE_HiC/DirectionalityIndex/ENCFF237QCN.hic@1kb.directionality_3kb.bed

echo "Minimum Insulation Score"
Rscript $INS $DIRLCL $FILEINS consensusLcL " consensus LcL"

echo "Norm Counts Directionality Index"
Rscript $DI $DIRLCL $FILEDI consensusLcL " consensus LcL"

echo "Mean Switch Intensity  Directionality Index"
Rscript $DII $DIRLCL $FILEDI consensusLcL " consensus LcL"

