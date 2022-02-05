import postgap.Globals
import postgap.GWAS
import postgap.FinemapIntegration

import postgap.Ensembl_lookup
from postgap.DataModel import *
from postgap.Utils import *
 

def scan_disease_databases():
    if postgap.Globals.GWAS_SUMMARY_STATS_FILE is not None:
        gwas_associations = postgap.GWAS.GWAS_File().run()

    associations_by_snp = dict()
    for gwas_association in gwas_associations:
        if gwas_association.snp not in associations_by_snp or associations_by_snp[gwas_association.snp].pvalue < gwas_association.pvalue:
            associations_by_snp[gwas_association.snp] = GWAS_SNP(
                snp=SNP(rsID=gwas_association.snp,
                        chrom=None,
                        pos=None,
                        approximated_zscore=None),
                pvalue=gwas_association.pvalue,
                evidence=[gwas_association],
                z_score=None,
                beta=None
            )

    gwas_snps = associations_by_snp.values()
    return gwas_snps


def get_gwas_snp_locations(gwas_snps):
    original_gwas_snp = dict((gwas_snp.snp.rsID, gwas_snp)
                             for gwas_snp in gwas_snps)
    mapped_snps = postgap.Ensembl_lookup.get_snp_locations(
        original_gwas_snp.keys())
    return [
        GWAS_SNP(
            snp=mapped_snp,
            pvalue=original_gwas_snp[mapped_snp.rsID].pvalue,
            evidence=original_gwas_snp[mapped_snp.rsID].evidence,
            z_score=original_gwas_snp[mapped_snp.rsID].z_score,
            beta=original_gwas_snp[mapped_snp.rsID].beta
        )
        for mapped_snp in mapped_snps
        if mapped_snp.rsID in original_gwas_snp
    ]


def gwas_snp_to_precluster(gwas_snp):
    mapped_ld_snps = []
    for line in open(postgap.Globals.GWAS_SUMMARY_STATS_FILE):
        chromosome, position, rsID, effect_allele, non_effect_allele, beta, se, pvalue, z_score, MAF = line.rstrip().split('\t')
        if chromosome == 'Chromosome':
            continue
        mapped_ld_snps.append(SNP(
            rsID=rsID,
            chrom=chromosome,
            pos=int(position),
            approximated_zscore=None
        ))

    return GWAS_Cluster(
        gwas_snps=[gwas_snp],
        ld_snps=mapped_ld_snps,
        ld_matrix=None,
        z_scores=None,
        betas=None,
        mafs=None,
        annotations=None,
        gwas_configuration_posteriors=None
    )


def cisregulatory_evidence(ld_snps):
    evidence = concatenate(source().run(ld_snps)
                           for source in postgap.Cisreg.sources)
    flaten_evidence = []
    for evi in evidence:
        if type(evi) == list:
            for e in evi:
                flaten_evidence.append(e)
        else:
            flaten_evidence.append(evi)

    res = collections.defaultdict(lambda: collections.defaultdict(list))
    for association in flaten_evidence:
        res[association.snp][association.gene].append(association)
    return res


def regulatory_evidence(snps):
    res = concatenate(source().run(snps)
                      for source in postgap.Reg.get_filtered_subclasses(postgap.Globals.Reg_adaptors))
    hash = collections.defaultdict(list)
    for hit in res:
        hash[hit.snp].append(hit)
    return hash


def compute_v2g_scores(reg, cisreg):
    intermediary_scores = dict()
    gene_scores = dict()
    for gene in cisreg:
        intermediary_scores[gene] = collections.defaultdict(int)
        seen = set()
        for evidence in cisreg[gene] + reg:
            if evidence.source not in seen or float(evidence.score) > intermediary_scores[gene][evidence.source]:
                intermediary_scores[gene][evidence.source] = float(
                    evidence.score)
                seen.add(evidence.source)
        gene_scores[gene] = sum(intermediary_scores[gene][source] * postgap.Globals.EVIDENCE_WEIGHTS[source]
                                for source in intermediary_scores[gene] if source in postgap.Globals.EVIDENCE_WEIGHTS)
    return intermediary_scores, gene_scores


def create_SNP_GeneSNP_Associations(snp, reg, cisreg):
    intermediary_scores, gene_scores = compute_v2g_scores(reg, cisreg)
    rank = dict((score, index) for index, score in enumerate(
        sorted(gene_scores.values(), reverse=True)))

    return [GeneSNP_Association(
        gene=gene,
        snp=snp,
        cisregulatory_evidence=cisreg[gene],
        regulatory_evidence=reg,
        intermediary_scores=intermediary_scores[gene],
        score=gene_scores[gene],
        rank=rank[gene_scores[gene]] + 1)
        for gene in cisreg]


def ld_snps_to_genes(ld_snps):
    cisreg = cisregulatory_evidence(ld_snps)
    reg = regulatory_evidence(ld_snps)
    return concatenate((create_SNP_GeneSNP_Associations(snp, reg[snp], cisreg[snp]) for snp in ld_snps))


def cluster_to_genes(cluster, associations):
    """

        Associated Genes to a cluster of gwas_snps
        Args:
        * Cluster 
        * [ GeneSNP_Associations ] 
        * { tissue_name: scalar (weights) }
        * { population_name: scalar (weight) }
        Returntype: [ GeneCluster_Association ]

    """
    assert len(cluster.ld_snps) == cluster.ld_matrix.shape[0], (len(
        cluster.ld_snps), cluster.ld_matrix.shape[0], cluster.ld_matrix.shape[1])
    assert len(cluster.ld_snps) == cluster.ld_matrix.shape[1], (len(
        cluster.ld_snps), cluster.ld_matrix.shape[0], cluster.ld_matrix.shape[1])
    gene_tissue_posteriors = postgap.FinemapIntegration.compute_joint_posterior(
        cluster, associations)
    res = [
        GeneCluster_Association(
            gene=gene,
            score=None,
            collocation_posterior=gene_tissue_posteriors[gene],
            cluster=cluster,
            # This is a [ GeneSNP_Association ]
            evidence=filter(lambda X: X.gene == gene, associations),
            r2=None
        )
        for gene in gene_tissue_posteriors
    ]
    # Pick the association with the highest score
    return sorted(res, key=lambda X: X.score)


def clusters_to_genes(clusters):
    cluster_associations = [(cluster, ld_snps_to_genes(
        cluster.ld_snps)) for cluster in clusters]
    # hwangse
    with open(postgap.Globals.OUTPUT+'_GWAS_lambdas.txt', 'w') as fw1:
        fw1.write('source\tlambdas\n')
    with open(postgap.Globals.OUTPUT+'_eQTL_lambdas.txt', 'w') as fw2:
        fw2.write('cluster\ttissue\tgene\tlambdas\n')
    # ===

    if postgap.Globals.PERFORM_BAYESIAN:
        cluster_associations = postgap.FinemapIntegration.compute_gwas_posteriors(
            cluster_associations)

    # Perform cluster by cluster finemapping
    res = concatenate(cluster_to_genes(cluster, associations)
                      for cluster, associations in cluster_associations)

    if len(res) == 0:
        return []
    else:
        return sorted(res, key=lambda X: X.score)


def gwas_snps_to_genes():
    import postgap.Cisreg
    import postgap.Reg

    gwas_snps = scan_disease_databases()
    min_p = 1.
    selected = ['']
    for gwas_snp in gwas_snps:
        if gwas_snp.pvalue < min_p:
            min_p = gwas_snp.pvalue
            selected[0] = gwas_snp

    gwas_snp_locations = get_gwas_snp_locations(selected)
    clusters = filter(lambda X: X is not None, [gwas_snp_to_precluster(
        gwas_snp_location) for gwas_snp_location in gwas_snp_locations])

    return clusters_to_genes(clusters)
