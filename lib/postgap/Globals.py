#! /usr/bin/env python

"""

Copyright [1999-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License")
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

		 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

"""

"""

	Please email comments or questions to the public Ensembl
	developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

	Questions may also be sent to the Ensembl help desk at
	<http://www.ensembl.org/Help/Contact>.

"""

DATABASES_DIR = None
GWAS_SUMMARY_STATS_FILE = None
PERFORM_BAYESIAN = False
TYPE = None
kmax_gwas = 1
kmax_eqtl = 1
OUTPUT = None

EVIDENCE_WEIGHTS = {
        'S_DHS': 1,
        'S_TSS': 1,
        'S_coding': 1
}
GWAS_adaptors = None
Cisreg_adaptors = None
Reg_adaptors = None

# summary_stats_eqtl = 'abc'
ld_file = None
GTEx_path = None

# GWAS_PVALUE_CUTOFF = 1e-4
GWAS_PVALUE_CUTOFF = 99.

# GWAS sample size
gwas_sample_size = 100
 


ALL_TISSUES=[]

#hwangse
source_lst=['S_DHS', 'S_TSS', 'S_coding']
