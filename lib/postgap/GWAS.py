import postgap.Globals
import postgap.Integration
from postgap.DataModel import *


class GWAS_source(object):
    def run(self, diseases, iris):
        """

                Returns all GWAS SNPs associated to a disease in this source
                Args:
                * [ string ] (trait descriptions)
                * [ string ] (trait Ontology IRIs)
                Returntype: [ GWAS_Association ]

        """
        assert False, "This stub should be defined"


class GWAS_File(GWAS_source):
    display_name = "GWAS File"

    def create_gwas_association_collector(self):
        class gwas_association_collector:
            def __init__(self):
                self.found_list = []

            def add_to_found_list(self, gwas_association):
                self.found_list.append(gwas_association)

            def get_found_list(self):
                return self.found_list
        return gwas_association_collector()

    def create_pvalue_filter(self, pvalue_threshold):
        def filter_for_pvalues_smaller_than(pvalue):
            if float(pvalue) < pvalue_threshold:
                return True
            return False
        return filter_for_pvalues_smaller_than

    def parse_gwas_data_file(
        self,
        gwas_data_file,
        callback,
        want_this_gwas_association_filter,
        column_labels=[
            "Chromosome",
            "Position",
            "MarkerName",
            "Effect_allele",
            "Non_Effect_allele",
            "Beta",
            "SE",
            "Pvalue",
            "z_score",
            "MAF"
        ]
    ):

        file = open(gwas_data_file)

        for line in file:
            # Skip the line with headers
            if line.startswith("Chromosome"):
                continue

            items = line.rstrip().split('\t')
            parsed = dict()
            for column_index in range(len(column_labels)):
                column_label = column_labels[column_index]
                parsed[column_label] = items[column_index]

            if not want_this_gwas_association_filter(parsed["Pvalue"]):
                continue
 
            try:
                gwas_association = GWAS_Association(
                    pvalue=float(parsed["Pvalue"]),
                    pvalue_description='Manual',
                    snp=parsed["MarkerName"],
                    disease='None',
                    reported_trait="Manual",
                    source="Manual",
                    publication="PMID000",
                    study="Manual",
                    sample_size=postgap.Globals.gwas_sample_size,
                    odds_ratio="Manual",
                    odds_ratio_ci_start=None,
                    odds_ratio_ci_end=None,
                    beta_coefficient=float(parsed["Beta"]),
                    beta_coefficient_unit="Manual",
                    beta_coefficient_direction="Manual",
                    rest_hash=None,
                    risk_alleles_present_in_reference=None,
                )
            except ValueError:
                continue
            callback(gwas_association)

    def run(self):
        gwas_data_file = postgap.Globals.GWAS_SUMMARY_STATS_FILE
        pvalue_filtered_gwas_associations = self.create_gwas_association_collector()
        pvalue_filter = self.create_pvalue_filter(
            pvalue_threshold=postgap.Globals.GWAS_PVALUE_CUTOFF)

        self.parse_gwas_data_file(
            gwas_data_file=gwas_data_file,
            want_this_gwas_association_filter=pvalue_filter,
            callback=pvalue_filtered_gwas_associations.add_to_found_list
        )
        return pvalue_filtered_gwas_associations.get_found_list()
