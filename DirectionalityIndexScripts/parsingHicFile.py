import hicstraw
import numpy as np

path = "/project/CRUP_scores/CENTRE_HiC/HiCmaps/GM12878/thresholded_matrix/ENCFF237QCN.hic"
hic = hicstraw.HiCFile(path)

mzd = hic.getMatrixZoomData('chr11', 'chr11', "observed", "NONE", "BP", 10000)
print(mzd)
numpy_matrix = mzd.getRecordsAsMatrix(60000000, 66000000, 60000000, 66000000)
print(numpy_matrix)
np.savetxt('/project/CRUP_scores/CENTRE_HiC/DirectionalityIndexScripts/hicmap.txt', numpy_matrix)
