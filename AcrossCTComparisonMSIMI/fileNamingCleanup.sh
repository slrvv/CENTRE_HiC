#-------------------------------------------------------------------------------
# Author:  Sara Lopez Ruiz de Vargas
# Email:  lopez_s@molgen.mpg.de
#
# Date:    2024-04-08
#
# Script Name:  fileNamingCleanup.sh
#
# Script Description: small script to clean up the names of the Feature files.
# Remove 10kb30kb.csv from the end of files
#
#-------------------------------------------------------------------------------



MIPATH=/project/CRUP_scores/CENTRE_HiC/Features/MinInsulationScoreDiffCT/
MSIPATH=/project/CRUP_scores/CENTRE_HiC/Features/MeanSwitchIntensityDiffCT/

cd $MIPATH  
rename 10kb30kb.csv "" *

cd $MSIPATH
rename 10kb30kb.csv "" *
