########################################################################
#
# HiC insulation scores compute
#
########################################################################


# Script to compute the insulation scores of the 33 ENCODE samples
# we downloaded

DIRFILES=/project/CRUP_scores/CENTRE_HiC/HiCmaps
DIRDI=/project/CRUP_scores/CENTRE_HiC/InsulationScores

for SAMPLE in $(ls -F $DIRFILES | grep / )
do
    ###Make the directories to store the results
    echo $SAMPLE
    T=thresholded_matrix
    echo $DIRFILES/$SAMPLE$T
    ls $DIRFILES/$SAMPLE$T
    for FILE in $(ls $DIRFILES/$SAMPLE$T | grep .hic)
    do
	echo $FILE
	echo "Fanc insulation"
	fanc insulation $DIRFILES/$SAMPLE$T/$FILE \
	$DIRDI/$SAMPLE/$FILE.insulation -o bed
    done
done
