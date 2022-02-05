import postgap.Integration

class GWAS_File(GWAS_source):
    display_name = "GWAS File"

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

                
