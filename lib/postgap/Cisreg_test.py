

class S_GTEx(Cisreg_source):
	display_name = "S_GTEx"
	hash = dict()
	with open(postgap.Globals.summary_stats_eqtl) as f:
		for line in f:
			Chromosome,Position,rsID,Effect_allele,Non_Effect_allele,Beta,SE,Pvalue,z_score,MAF = line.strip().split('\t')
			if Chromosome =='Chromosome':
				continue
			hash[rsID] = dict()
			hash[rsID]['Pvalue']= Pvalue
			hash[rsID]['Beta'] = Beta

	def run(self, snps):
		res = concatenate(map(self.snp, snps))
		return res

	def snp(self, snp):
		Pvalue = float(self.hash[snp.rsID]['Pvalue'])
		Beta   = float(self.hash[snp.rsID]['Beta'])
		z_score = postgap.FinemapIntegration.z_score_from_pvalue(Pvalue, Beta)
		res=[ Cisregulatory_Evidence(
				snp = snp,
				gene= postgap.Ensembl_lookup.get_ensembl_gene('ENSG00000078900'),
				tissue = 'TEST_tissue',
				score = 1 - Pvalue,
				source = self.display_name,
				study = None,
				info = None,
				z_score = z_score,
				pvalue = Pvalue,
				beta = Beta
				)
			]
		return res

sources = Cisreg_source.__subclasses__()

