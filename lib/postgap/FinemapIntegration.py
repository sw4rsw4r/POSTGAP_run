
from postgap.DataModel import *
import cPickle as pickle
import postgap.Globals
import postgap.LD
import numpy
import postgap.Finemap



def compute_gwas_posteriors(cluster_associations):
    prepped_clusters = []
    for cluster, associations in cluster_associations:
        if (len(cluster.ld_snps) < 3):  # hsw
            continue
        prepped_clusters.append((prepare_cluster_for_finemap(
            cluster, associations), associations))
    # MLE calculation
    prepped_clusters = [(postgap.Finemap.mk_modified_clusters(
        cluster), associations) for cluster, associations in prepped_clusters]
    return [(finemap_gwas_cluster(cluster), associations) for cluster, associations in prepped_clusters]



def prepare_cluster_for_finemap(cluster, associations):
    ld_snps, ld_matrix, z_scores, betas = extract_z_scores_from_file(
            cluster)

    cluster_f = GWAS_Cluster(cluster.gwas_snps, ld_snps,
                             ld_matrix, z_scores, betas, None, None, None)

    mafs = extract_snp_mafs(cluster_f, associations)
    mafs = mafs.astype(float)
    annotations = (extract_snp_annotations(
        cluster_f, associations) > 0.).astype('float')

    assert len(ld_snps) == ld_matrix.shape[0]
    assert len(ld_snps) == ld_matrix.shape[1]
    return GWAS_Cluster(cluster.gwas_snps, ld_snps, ld_matrix, z_scores, betas, mafs, annotations, None)


def extract_snp_mafs(cluster, associations):
    maf_hash = dict()
    for line in open(postgap.Globals.GWAS_SUMMARY_STATS_FILE):
        chromosome, position, rsID, effect_allele, non_effect_allele, beta, se, pvalue, z_score, MAF = line.rstrip().split('\t')
        if chromosome == 'Chromosome':
            continue
        rsID = rsID.strip()  # hsw
        maf_hash[rsID] = float(MAF)
    return numpy.array([maf_hash[snp.rsID] for snp in cluster.ld_snps])



def extract_snp_annotations(cluster, associations):
    annotation_hash = collections.defaultdict(
        lambda: collections.defaultdict(float))
    for association in associations:
        for evidence in association.regulatory_evidence:
            annotation_hash[evidence.source][evidence.snp.rsID] = evidence.score
    return numpy.array([[annotation_hash[annotation][snp.rsID] for snp in cluster.ld_snps] for annotation in postgap.Globals.source_lst])


def extract_z_scores_from_file(cluster):
    ld_snp_hash = dict((ld_snp.rsID, ld_snp)
                       for index, ld_snp in enumerate(cluster.ld_snps))
    ld_snp_results = dict()
    missing = len(ld_snp_hash)
    found_ld_snps = []
    file = open(postgap.Globals.GWAS_SUMMARY_STATS_FILE)
    for line in file:
        # Chromosome	Position	MarkerName	Effect_allele	Non_Effect_allele	Beta	SE	Pvalue
        # 1	751343	rs28544273	A	T	-0.0146	0.0338	0.6651
        chromosome, position, rsID, effect_allele, non_effect_allele, beta, se, pvalue, z_score, MAF = line.rstrip().split('\t')
        rsID = rsID.strip()  # hsw
        if rsID in ld_snp_hash:
            ld_snp_results[rsID] = (float(pvalue), float(beta))
            found_ld_snps.append(ld_snp_hash[rsID])
            missing -= 1
            if missing == 0:
                break

    # Compute ld_matrix
    ld_snp_ids, ld_matrix = postgap.LD.get_pairwise_ld(found_ld_snps)

    # Update list of LD SNPs
    ld_snps = [ld_snp_hash[rsID] for rsID in ld_snp_ids]
    z_scores = [z_score_from_pvalue(
        ld_snp_results[rsID][0], ld_snp_results[rsID][1]) for rsID in ld_snp_ids]
    betas = [ld_snp_results[rsID][1] for rsID in ld_snp_ids]
    assert len(ld_snps) == ld_matrix.shape[0]
    assert len(ld_snps) == ld_matrix.shape[1]
    return ld_snps, ld_matrix, z_scores, betas


def finemap_gwas_cluster(cluster):
    '''

            Enriches GWAS clusters with z-scores and GWAS posteriors
            Arg1: GWAS_Cluster
            Arg2: mafs (Numpy vector)
            Arg3: annotations (2D Numpy array)
            Returntype: GWAS_Cluster

    '''
    # Define experiment label (serves for debugging logs)
    chrom = cluster.ld_snps[0].chrom
    start = min(ld_snp.pos for ld_snp in cluster.ld_snps)
    end = max(ld_snp.pos for ld_snp in cluster.ld_snps)
    sample_label = 'GWAS_Cluster_%s:%i-%i' % (chrom, start, end)

    # Define LD SNP labels (serves for debugging logs)
    ld_snp_ids = [ld_snp.rsID for ld_snp in cluster.ld_snps]

    # Define sample size: mean of max for each SNP
    sample_sizes = map(lambda gwas_snp: max(
        gwas_association.sample_size for gwas_association in gwas_snp.evidence), cluster.gwas_snps)
    sample_size = sum(sample_sizes) / len(sample_sizes)

    # hwangse ====
    with open(postgap.Globals.OUTPUT+'_GWAS_lambdas.txt', 'a') as fw1:
        for idx, L in enumerate(cluster.lambdas):
            fw1.write(
                '\t'.join(map(str, [sample_label, postgap.Globals.source_lst[idx], L]))+'\n')

    # ======

    # Compute posterior
    if postgap.Globals.TYPE == 'binom' or postgap.Globals.TYPE == 'ML':
        configuration_posteriors = postgap.Finemap.finemap_v1(
        z_scores=numpy.array(cluster.z_scores),
        beta_scores=numpy.array(cluster.betas),
        cov_matrix=cluster.ld_matrix,
        n=sample_size,
        labels=ld_snp_ids,
        sample_label=sample_label,
        lambdas=cluster.lambdas,
        mafs=cluster.mafs,
        annotations=cluster.annotations,
        kmax=postgap.Globals.kmax_gwas,
        isGWAS=True
        )
    elif postgap.Globals.TYPE == 'EM' or postgap.Globals.TYPE == 'ML_EM':
        configuration_posteriors = postgap.Finemap.finemap_v2(
            z_scores=numpy.array(cluster.z_scores),
            beta_scores=numpy.array(cluster.betas),
            cov_matrix=cluster.ld_matrix,
            n=sample_size,
            labels=ld_snp_ids,
            sample_label=sample_label,
            lambdas=cluster.lambdas,
            mafs=cluster.mafs,
            annotations=cluster.annotations,
            kmax=postgap.Globals.kmax_gwas,
            isGWAS=True
        )
    print postgap.Globals.TYPE

    return GWAS_Cluster_with_lambdas(cluster.gwas_snps, cluster.ld_snps, cluster.ld_matrix, cluster.z_scores, cluster.betas, cluster.mafs, cluster.annotations, configuration_posteriors, cluster.lambdas)


def compute_joint_posterior(cluster, associations):
    """
            Compute collocation posterior of gene expression and GWAS phenotype at the specified cluster and tissue
            Arg1: GWAS_Cluster
            Arg2: [GeneSNP_Association]
            Arg3: String
            Returntype: Hash of hashes: Gene => Tissue => Float
    """
    print "Computing finemap on new cluster"
    # print cluster

    assert len(cluster.ld_snps) == cluster.ld_matrix.shape[0], (len(
        cluster.ld_snps), cluster.ld_matrix.shape[0], cluster.ld_matrix.shape[1])
    assert len(cluster.ld_snps) == cluster.ld_matrix.shape[1], (len(
        cluster.ld_snps), cluster.ld_matrix.shape[0], cluster.ld_matrix.shape[1])
    gene_tissue_snp_eQTL_hash = organise_eQTL_data(associations)

    return dict((gene, compute_gene_joint_posterior(cluster, gene, gene_tissue_snp_eQTL_hash[gene], cluster.gwas_configuration_posteriors, cluster.mafs, cluster.annotations)) for gene in gene_tissue_snp_eQTL_hash)



def compute_gene_joint_posterior(cluster, gene, tissue_snp_eQTL_hash, gwas_configuration_posteriors, mafs, annotations):
    """
            Compute collocation posterior of gene expression and GWAS phenotype at the specified cluster and tissue
            Arg1: GWAS_Cluster
            Arg2: Gene
            Arg3: Hash of hashes: Tissue => SNP => (Float, Float)
            Arg4: Hash of hashes: configuration => posterior
            Arg5: Numpy Vector
            Arg6: Numpy 2D Array
            Returntype: Hash of hashes: Tissue => Float
    """
    assert len(cluster.ld_snps) == cluster.ld_matrix.shape[0], (len(
        cluster.ld_snps), cluster.ld_matrix.shape[0], cluster.ld_matrix.shape[1])
    assert len(cluster.ld_snps) == cluster.ld_matrix.shape[1], (len(
        cluster.ld_snps), cluster.ld_matrix.shape[0], cluster.ld_matrix.shape[1])
    return dict((tissue, compute_gene_tissue_joint_posterior(cluster, gene, tissue, tissue_snp_eQTL_hash[tissue], gwas_configuration_posteriors, mafs, annotations)) for tissue in tissue_snp_eQTL_hash)


def compute_gene_tissue_joint_posterior(cluster, gene, tissue, eQTL_snp_hash, gwas_configuration_posteriors, mafs, annotations):
    """
            Compute posterior of gene expression regulation at the specified cluster and tissue
            Arg1: GWAS_Cluster
            Arg2: Tissue (string)
            Arg3: Gene
            Arg4: Hash string (rsID) => (Float (z-score), Float (beta))
            Arg5: Hash of hashes: configuration => posterior
            Arg6: Numpy Vector
            Arg7: Numpy 2D Array
            Returntype: Float
    """
    # eQTL posteriors
    eQTL_configuration_posteriors = compute_eqtl_posteriors(
        cluster, tissue, gene, eQTL_snp_hash, mafs, annotations)
    # Joint posterior
    joint_out = eQTL_configuration_posteriors.joint_posterior(
        gwas_configuration_posteriors)
    pickle.dump(joint_out[0], open(postgap.Globals.OUTPUT+'_'+str(tissue) +
                                   '_'+str(gene.name)+'_snp_clpp.pkl', "w"))  # DEBUG remove hard coded path
    pickle.dump(joint_out[1], open(postgap.Globals.OUTPUT +
                                   '_'+str(tissue)+'_'+str(gene.name)+'_eqtl_PIP.pkl', "w"))
    return joint_out[2], joint_out[3]



def compute_eqtl_posteriors(cluster, tissue, gene, eQTL_snp_hash, mafs, annotations):
    """
            Compute posterior of gene expression regulation at the specified cluster and tissue
            Arg1: GWAS_Cluster
            Arg2: Tissue (string)
            Arg3: Gene
            Arg4: Hash string (rsID) => (Float (z-score), Float (beta))
            Arg5: Numpy Vector
            Arg6: Numpy 2D Array
            Returntype: Float
    """
    assert len(cluster.ld_snps) == cluster.ld_matrix.shape[0], (len(
        cluster.ld_snps), cluster.ld_matrix.shape[0], cluster.ld_matrix.shape[1])
    assert len(cluster.ld_snps) == cluster.ld_matrix.shape[1], (len(
        cluster.ld_snps), cluster.ld_matrix.shape[0], cluster.ld_matrix.shape[1])
    # Determine which SNPs are missing values
    missing_indices = numpy.array([index for index, ld_snp in enumerate(
        cluster.ld_snps) if ld_snp.rsID not in eQTL_snp_hash]).astype(int)
    known_z_scores = numpy.array([eQTL_snp_hash[ld_snp.rsID][0]
                                  for ld_snp in cluster.ld_snps if ld_snp.rsID in eQTL_snp_hash])
    known_betas = numpy.array([eQTL_snp_hash[ld_snp.rsID][1]
                               for ld_snp in cluster.ld_snps if ld_snp.rsID in eQTL_snp_hash])
    assert len(known_z_scores) > 0
    assert len(missing_indices) != len(
        cluster.ld_snps), (missing_indices, known_z_scores)
    assert len(cluster.ld_snps) == cluster.ld_matrix.shape[0], (len(
        cluster.ld_snps), cluster.ld_matrix.shape[0], cluster.ld_matrix.shape[1])
    assert len(cluster.ld_snps) == cluster.ld_matrix.shape[1], (len(
        cluster.ld_snps), cluster.ld_matrix.shape[0], cluster.ld_matrix.shape[1])

    if len(missing_indices) > 0:
        # Generate LD matrix of known values
        ld_matrix_known = numpy.delete(
            cluster.ld_matrix, missing_indices, axis=1)
        ld_matrix_known = numpy.delete(
            ld_matrix_known, missing_indices, axis=0)
        assert ld_matrix_known.size > 0, (missing_indices,
                                          cluster.ld_matrix, ld_matrix_known)

        # Generate LD matrix of known SNPs to missing SNPs
        ld_matrix_k2m = cluster.ld_matrix[missing_indices, :]
        ld_matrix_k2m = numpy.delete(ld_matrix_k2m, missing_indices, axis=1)

        # Imputation
        shrink_lambda = 0.1  # shrinkage factor, magic number
        ld_matrix_known_shrink = shrink_lambda * numpy.diag(numpy.ones(
            ld_matrix_known.shape[0])) + (1-shrink_lambda) * ld_matrix_known
        assert ld_matrix_known_shrink.size > 0, (
            missing_indices, ld_matrix, ld_matrix_known)
        ld_matrix_k2m_shrink = (1-shrink_lambda) * ld_matrix_k2m
        z_shrink_imputed = numpy.dot(numpy.dot(ld_matrix_k2m_shrink, numpy.linalg.pinv(
            ld_matrix_known_shrink, 0.0001)), known_z_scores)
        beta_shrink_imputed = numpy.dot(numpy.dot(ld_matrix_k2m_shrink, numpy.linalg.pinv(
            ld_matrix_known_shrink, 0.0001)), known_betas)

        # Aggregate z_scores into a single vector
        z_scores = []
        betas = []
        for i in cluster.ld_snps:
            z_scores.append(0)
            betas.append(0)
        for index, z_score, beta in zip(missing_indices, z_shrink_imputed, beta_shrink_imputed):
            z_scores[index] = z_score
            betas[index] = beta
        for index, ld_snp in enumerate(cluster.ld_snps):
            if ld_snp.rsID in eQTL_snp_hash:
                z_scores[index] = eQTL_snp_hash[ld_snp.rsID][0]
                betas[index] = eQTL_snp_hash[ld_snp.rsID][1]
    else:
        z_scores = known_z_scores
        betas = known_betas

    # Define experiment label (serves for debugging logs)
    chrom = cluster.ld_snps[0].chrom
    start = min(ld_snp.pos for ld_snp in cluster.ld_snps)
    end = max(ld_snp.pos for ld_snp in cluster.ld_snps)
    sample_label = 'eQTL_Cluster_%s:%i-%i_%s' % (chrom, start, end, gene)

    # hwangse ====
    cluster_label = 'Cluster_%s:%i-%i' % (chrom, start, end)
    lambdas = postgap.Finemap.mk_eqtl_lambdas(cluster, numpy.array(z_scores))
    with open(postgap.Globals.OUTPUT+'_eQTL_lambdas.txt', 'a') as fw2:
        for idx, L in enumerate(lambdas):
            if idx < len(postgap.Globals.source_lst):
                fw2.write('\t'.join(map(str, [
                          cluster_label, tissue, gene.name, postgap.Globals.source_lst[idx], L]))+'\n')
            else:
                fw2.write(
                    '\t'.join(map(str, [cluster_label, tissue, gene.name, idx, L]))+'\n')
    # ======

    # Compute posterior
    return postgap.Finemap.finemap_v1(
        z_scores=numpy.array(z_scores),
        beta_scores=numpy.array(betas),
        cov_matrix=cluster.ld_matrix,
        n=1000,  # TODO extract eQTL sample sizes
        labels=[ld_snp.rsID for ld_snp in cluster.ld_snps],
        sample_label=sample_label,
        lambdas=lambdas,  # hwangsw cluster.lambdas,
        mafs=mafs,
        annotations=annotations,
        kmax=postgap.Globals.kmax_eqtl,  # eQTL_kmax
        kstart=1,  # eQTL_kstart
        isGWAS=False
    )



def organise_eQTL_data(associations):
    """
            Organise unsorted eQTL data into easily read hash of hashes:
            Arg1: [GeneSNP_Association] 
            Returntype: Hash of hashes: Gene => Tissue => SNP => Float
    """
    res = collections.defaultdict(lambda: collections.defaultdict(
        lambda: collections.defaultdict(float)))
    for association in associations:
        for evidence in association.cisregulatory_evidence:
            # if evidence.source == 'GTEx':
            if evidence.source == 'S_GTEx':
                res[association.gene][evidence.tissue][association.snp.rsID] = (
                    evidence.z_score, evidence.beta)
    return res


def sign(number):
	"""
		Returns the sign of the number (-1, 0 or 1)
		Arg1: float
		Returntype: int
	"""
	if number > 0:
		return 1
	elif number < 0:
		return -1
	else:
		return 0


def z_score_from_pvalue(p_value, direction):
    """
            Estimates z-score from p-value and effect direction
            Arg1: float
            Arg2: float
            Returntype: float
    """
    import scipy
    import scipy.stats
    if p_value == 0:
        p_value = 4.2e-317
    return -scipy.stats.norm.ppf(p_value/2) * sign(direction)

 