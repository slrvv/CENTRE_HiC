########################################################################
#
# HiC insulation scores & dir index compute bigWig
#
########################################################################


#Script to compute the insulation scores of GM12878 IMR-90 and K562.
#Adjust to original resolution of experiments 
DIRFILES=/project/CRUP_scores/CENTRE_HiC/HiCmaps
DIRINS=/project/CRUP_scores/CENTRE_HiC/InsulationScores
DIRDI=/project/CRUP_scores/CENTRE_HiC/DirectionalityIndex

#echo "GM12878"

#fanc insulation $DIRFILES/GM12878/thresholded_matrix/ENCFF237QCN.hic@10kb \
#	/project/CRUP_scores/CENTRE_HiC/FinalPlots/GenomeTrack/GM12878.ENCFF237QCN@10kb.insulation \
#	-w 30kb -o bigwig
	
#fanc directionality $DIRFILES/GM12878/thresholded_matrix/ENCFF237QCN.hic@10kb \
#	/project/CRUP_scores/CENTRE_HiC/FinalPlots/GenomeTrack/GM12878ENCFF237QCN@10kb.directionality \
#	-w 30kb -o bigwig 

echo "K562"

fanc insulation $DIRFILES/K562/thresholded_matrix/ENCFF080DPJ.hic@10kb \
	/project/CRUP_scores/CENTRE_HiC/FinalPlots/GenomeTrack/K562.ENCFF080DPJ@10kb.insulation -w 30kb -o bigwig

fanc directionality $DIRFILES/K562/thresholded_matrix/ENCFF080DPJ.hic@10kb \
	/project/CRUP_scores/CENTRE_HiC/FinalPlots/GenomeTrack/K562.ENCFF080DPJ@10kb.directionality -w 30kb -o bigwig
