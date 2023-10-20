########################################################################
#
# HiC insulation scores & dir index compute
#
########################################################################


#Script to compute the insulation scores of GM12878 IMR-90 and K562.
#Adjust to original resolution of experiments 
DIRFILES=/project/CRUP_scores/CENTRE_HiC/HiCmaps
DIRINS=/project/CRUP_scores/CENTRE_HiC/InsulationScores
DIRDI=/project/CRUP_scores/CENTRE_HiC/DirectionalityIndex

echo "GM12878"

fanc insulation $DIRFILES/GM12878/thresholded_matrix/ENCFF237QCN.hic@1kb \
	$DIRINS/GM12878/ENCFF237QCN@1kb.insulation -o bed 
	
fanc directionality $DIRFILES/GM12878/thresholded_matrix/ENCFF237QCN.hic@1kb \
	$DIRDI/GM12878/ENCFF237QCN@1kb.directionality -o bed 

echo "IMR-90"

fanc insulation $DIRFILES/IMR-90/thresholded_matrix/ENCFF636CEM.hic@5kb \
	$DIRINS/IMR-90/ENCFF636CEM@5kb.insulation -o bed

fanc directionality $DIRFILES/IMR-90/thresholded_matrix/ENCFF636CEM.hic@5kb \
	$DIRDI/IMR-90/ENCFF636CEM@5kb.directionality -o bed

echo "K562"

fanc insulation $DIRFILES/K562/thresholded_matrix/ENCFF080DPJ.hic@5kb \
	$DIRINS/K562/ENCFF080DPJ@5kb.insulation -o bed

fanc directionality $DIRFILES/K562/thresholded_matrix/ENCFF080DPJ.hic@5kb \
	$DIRINS/K562/ENCFF080DPJ@5kb.directionality -o bed

	
