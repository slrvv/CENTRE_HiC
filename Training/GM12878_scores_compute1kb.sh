##################################################################################################
#                                                                                                #
#  GM12878_scores_compute1kb.sh: Compute the insulation and directionality for GM12878 on a      #
#  window of 1kb (original window of the experiments)                                            #
#                                                                                                #
##################################################################################################

# For now we used the same window size for all experiments (10kb) we want to use the window size
# corresponding to original resolution eeach of the HiC experiments. In GM12878  it is 1kb.

DIRMAPS=/project/CRUP_scores/CENTRE_HiC/HiCmaps/GM12878/thresholded_matrix/ENCFF237QCN.hic
DIRINS=/project/CRUP_scores/CENTRE_HiC/InsulationScores/GM12878
DIRDI=/project/CRUP_scores/CENTRE_HiC/DirectionalityIndex/GM12878

echo "Fanc insulation"
fanc insulation $DIRMAPS@1kb \
$DIRINS@1kb.insulation -o bed

echo "Fanc directionality"
fanc directionality $DIRMAPS@1kb \
$DIRDI@1kb.directionality -o bed

