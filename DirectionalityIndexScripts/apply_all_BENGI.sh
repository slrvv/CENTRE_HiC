#######################################################################
# apply_all_BENGI: Summary on how we applied wilcoxtest to all BENGI
#######################################################################

# In order to find out whether using for each ET pair the minimum over 
# the directionality 10kb overlapping window, we apply the procedure
# of wilcoxtest_dirindex script to all of the BENGI processed datasets

BENGIPATH="/project/CRUP_scores/toSara/BENGI_processed_datasets"
BENGIEND="-Benchmark.v38.txt"

SCOREPATH="/project/CRUP_scores/CENTRE_HiC/DirectionalityIndex"

SCOREEND=".hic@10kb.directionality_30kb.bed"

#PROGRAM=wilcoxtest_dirindex.R

PROGRAM=dirindex_switchintensity.R

Rscript $PROGRAM $BENGIPATH/Colon.GTEx$BENGIEND $SCOREPATH/transverse-colon/ENCFF858ZGP$SCOREEND \
WilcoxtestColon.GTEx " Colon.GTEx"

Rscript $PROGRAM $BENGIPATH/IMR90.HiC$BENGIEND $SCOREPATH/IMR-90/ENCFF636CEM$SCOREEND \
WilcoxtestIMR90.HiC " IMR90.HiC"

Rscript $PROGRAM $BENGIPATH/K562.CRISPR$BENGIEND $SCOREPATH/K562/ENCFF080DPJ$SCOREEND \
WilcoxtestK562.CRISPR " K562.CRISPR"

Rscript $PROGRAM $BENGIPATH/K562.HiC$BENGIEND $SCOREPATH/K562/ENCFF080DPJ$SCOREEND \
WilcoxtestK562.HiC " K562.HiC"

Rscript $PROGRAM $BENGIPATH/Ovary.GTEx$BENGIEND $SCOREPATH/ovary/ENCFF700CYI$SCOREEND \
WilcoxtestOvary.GTEx " Ovary.GTEx"

Rscript $PROGRAM $BENGIPATH/Pancreas.GTEx$BENGIEND $SCOREPATH/pancreas/ENCFF586MQY$SCOREEND \
WilcoxtestPancreas.GTEx " Pancreas.GTEx"

Rscript $PROGRAM $BENGIPATH/Stomach.GTEx$BENGIEND $SCOREPATH/stomach/ENCFF883XXW$SCOREEND \
WilcoxtestStomach.GTEx " Stomach.GTEx"

Rscript $PROGRAM $BENGIPATH/GM12878.RNAPII-ChIAPET$BENGIEND $SCOREPATH/GM12878/ENCFF237QCN$SCOREEND \
WilcoxtestGM12878.RNAPII-ChIAPET.png " GM12878.RNAPII-ChIAPET"

Rscript $PROGRAM $BENGIPATH/GM12878.RNAPII-ChIAPET$BENGIEND $SCOREPATH/GM12878/ENCFF237QCN$SCOREEND \
WilcoxtestGM12878.RNAPII-ChIAPET.png " GM12878.RNAPII-ChIAPET"

Rscript $PROGRAM $BENGIPATH/GM12878.CHiC$BENGIEND $SCOREPATH/GM12878/ENCFF237QCN$SCOREEND \
WilcoxtestGM12878.CHiC.png " GM12878.CHiC"

Rscript $PROGRAM $BENGIPATH/GM12878.CTCF-ChIAPET$BENGIEND $SCOREPATH/GM12878/ENCFF237QCN$SCOREEND \
WilcoxtestGM12878.CTCF-ChIAPET.png " GM12878.CTCF-ChIAPET"

Rscript $PROGRAM $BENGIPATH/GM12878.GEUVADIS$BENGIEND $SCOREPATH/GM12878/ENCFF237QCN$SCOREEND \
WilcoxtestGM12878.GEUVADIS.png " GM12878.GEUVADIS"

Rscript $PROGRAM $BENGIPATH/GM12878_HiC$BENGIEND $SCOREPATH/GM12878/ENCFF237QCN$SCOREEND \
WilcoxtestGM12878_HiC.png " GM12878.HiC"
