
import collections

Gene = collections.namedtuple("Gene", ['name', 'id', 'chrom', 'tss', 'biotype'])
SNP = collections.namedtuple(
	"SNP", 
	[
		'rsID',
		'chrom',
		'pos',
		'approximated_zscore'
	]
)
GWAS_Association = collections.namedtuple(
	'GWAS_Association', 
	[
		'snp',
		'disease',
		'reported_trait',
		'pvalue',
		'pvalue_description',
		'sample_size',
		'source',
		'publication',
		'study',
		'odds_ratio',
		'odds_ratio_ci_start',
		'odds_ratio_ci_end',
		'beta_coefficient',
		'beta_coefficient_unit',
		'beta_coefficient_direction',
		# rest_hash is no longer needed, risk_alleles_present_in_reference 
		# gets set when the object is created and that is all that is ever 
		# needed.
		'rest_hash',
		'risk_alleles_present_in_reference'
	]
)

GWAS_SNP = collections.namedtuple(
	'GWAS_SNP', 
	[
		'snp',
		'pvalue',
		'z_score',
		'evidence',
		'beta'
	]
)


GWAS_Cluster = collections.namedtuple(
	'GWAS_Cluster', 
	[
		'gwas_snps',
		'ld_snps',
		'ld_matrix',
		'z_scores',
		'betas',
		'mafs',
		'annotations',
		'gwas_configuration_posteriors'
	]
)



GWAS_Cluster_with_lambdas = collections.namedtuple(
        'GWAS_Cluster_with_lambdas',
        [
                'gwas_snps',
                'ld_snps',
                'ld_matrix',
                'z_scores',
                'betas',
                'mafs',
                'annotations',
                'gwas_configuration_posteriors',
                'lambdas'
        ]
)

Cisregulatory_Evidence = collections.namedtuple(
	'Cisregulatory_Evidence', 
	[
		'snp',
		'gene',
		'score',
		'source',
		'study',
		'tissue',
		'info',
		'z_score',
		'pvalue',
		'beta'
	]
)
Regulatory_Evidence = collections.namedtuple('Regulatory_Evidence', ['snp','score','source','study','tissue','info'])
GeneSNP_Association = collections.namedtuple('GeneSNP_Association', ['gene', 'snp', 'score', 'rank', 'intermediary_scores', 'cisregulatory_evidence', 'regulatory_evidence'])

GeneCluster_Association = collections.namedtuple('GeneCluster_Association', ['gene', 'cluster', 'score', 'collocation_posterior', 'evidence','r2'])

