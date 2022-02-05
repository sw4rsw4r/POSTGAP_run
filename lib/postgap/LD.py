import pandas as pd
import postgap.Globals
import numpy

def get_pairwise_ld(ld_snps, population='EUR'):
	"""

		For large numbers of SNPs, best to specify SNP region with chrom:to-from, e.g. 1:7654947-8155562

		Args:
		* [ SNP ], SNPs of interest
		* string, population name
		Returntype: [String (rsID)], Numpy.Array (2D)

	"""
        r2_array = pd.read_csv(postgap.Globals.ld_file,header=None,sep=' ').to_numpy()
        return [x.rsID for x in ld_snps], r2_array


