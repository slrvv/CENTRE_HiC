#######################################################################
# apply_all_BENGI: Summary on how we applied wilcoxtest to all BENGI
#######################################################################

# In order to find out whether using for each ET pair the minimum over 
# the Insulation Score 10kb overlapping window, we apply the procedure
# of wilcoxtest_insscore script to all of the BENGI processed datasets

BENGIPATH="/project/CRUP_scores/toSara/BENGI_processed_datasets"
BENGIEND="-Benchmark.v38.txt"

SCOREPATH="/project/CRUP_scores/CENTRE_HiC/InsulationScores"

SCOREEND=".hic@10kb.insulation_30kb.bed"

#Colon.GTEx IMR90.HiC K562.CRISPR K562.HiC Liver.GTEx
#NHEK.HiC Ovary.GTEx HeLa.CTCF-ChIAPET	  Pancreas.GTEx
#HeLa.HiC		  Stomach.GTEx HeLa.RNAPII-ChIAPET

#For now we dont have data for HeLa, NHEK or liver
#Liver --> right lobe liver available
#HeLa --> only raw sequences available


Rscript wilcoxtest_insscore.R $BENGIPATH/Colon.GTEx$BENGIEND $SCOREPATH/transverse-colon/ENCFF858ZGP$SCOREEND \
WilcoxtestColon.GTEx " Colon.GTEx"

Rscript wilcoxtest_insscore.R $BENGIPATH/IMR90.HiC$BENGIEND $SCOREPATH/IMR-90/ENCFF636CEM$SCOREEND \
WilcoxtestIMR90.HiC " IMR90.HiC"

Rscript wilcoxtest_insscore.R $BENGIPATH/K562.CRISPR$BENGIEND $SCOREPATH/K562/ENCFF080DPJ$SCOREEND \
WilcoxtestK562.CRISPR " K562.CRISPR"

Rscript wilcoxtest_insscore.R $BENGIPATH/K562.HiC$BENGIEND $SCOREPATH/K562/ENCFF080DPJ$SCOREEND \
WilcoxtestK562.HiC " K562.HiC"

Rscript wilcoxtest_insscore.R $BENGIPATH/Ovary.GTEx$BENGIEND $SCOREPATH/ovary/ENCFF700CYI$SCOREEND \
WilcoxtestOvary.GTEx " Ovary.GTEx"

Rscript wilcoxtest_insscore.R $BENGIPATH/Pancreas.GTEx$BENGIEND $SCOREPATH/pancreas/ENCFF586MQY$SCOREEND \
WilcoxtestPancreas.GTEx " Pancreas.GTEx"

Rscript wilcoxtest_insscore.R $BENGIPATH/Stomach.GTEx$BENGIEND $SCOREPATH/stomach/ENCFF883XXW$SCOREEND \
WilcoxtestStomach.GTEx " Stomach.GTEx"

Rscript wilcoxtest_insscore.R $BENGIPATH/GM12878.RNAPII-ChIAPET$BENGIEND $SCOREPATH/GM12878/ENCFF237QCN$SCOREEND \
WilcoxtestGM12878.RNAPII-ChIAPET.png " GM12878.RNAPII-ChIAPET"
