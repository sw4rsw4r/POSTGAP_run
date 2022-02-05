import postgap.Globals
from postgap.DataModel import *
GRCH37_ENSEMBL_REST_SERVER = "http://grch37.rest.ensembl.org"
known_genes = {}
def get_snp_locations(rsIDs):
    hash = dict()
    with open(postgap.Globals.GWAS_SUMMARY_STATS_FILE) as f:
        for line in f:
            Chromosome, Position, rsID, Effect_allele, Non_Effect_allele, Beta, SE, Pvalue, z_score, MAF = line.strip().split('\t')
            hash[rsID] = dict()
            hash[rsID]['seq_region_name'] = Chromosome
            hash[rsID]['start'] = Position
            hash[rsID]['end'] = Position

    results = []
    for rsID in rsIDs:
        mapping = hash[rsID]
        snp = SNP(
            rsID=rsID,
            chrom=mapping['seq_region_name'],
            pos=(int(mapping['start']) + int(mapping['end'])) / 2,
            approximated_zscore=None
        )
        results.append(snp)
    return results


def get_ensembl_gene(ensembl_id, ENSEMBL_REST_SERVER=GRCH37_ENSEMBL_REST_SERVER):
    """

            Get gene details from name
            * string
            Returntype: Gene

    """
    key = (ensembl_id, ENSEMBL_REST_SERVER)
    if key not in known_genes:
        known_genes[key] = Gene(
            name=ensembl_id,
            id=ensembl_id,
            chrom=None,
            tss=None,
            biotype=None
        )
    return known_genes[key]
