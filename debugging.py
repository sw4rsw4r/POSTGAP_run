from postgap.FinemapIntegration import *
from postgap.Utils import *
from postgap.DataModel import *
import postgap.Ensembl_lookup
import postgap.FinemapIntegration
import postgap.Cisreg
import postgap.Reg
import postgap.GWAS
import os
import sys
import argparse
import postgap.Globals
import postgap.Integration
import postgap.GWAS
import collections
import postgap.DataModel
import cPickle as pickle
from argparse import RawTextHelpFormatter
from postgap.Finemap import *

# EM debugging
postgap.Globals.PERFORM_BAYESIAN = True
postgap.Globals.GWAS_SUMMARY_STATS_FILE = '/home/seongwonhwang/Desktop/projects/POSTGAP/p_postgap/postgap/generate_synthetic_1000g/final_simulated_input/100samples/causal_110_0.03_str_high/summary-data/summary_stats_gwas_tss0'
# postgap.Globals.GWAS_SUMMARY_STATS_FILE = '/Users/msung/projects/p_postgap/postgap/generate_synthetic_1000g/final_simulated_input/100samples/causal_110_0.03_str_high/summary-data/summary_stats_gwas_tss0'
postgap.Globals.summary_stats_eqtl = '/home/seongwonhwang/Desktop/projects/POSTGAP/p_postgap/postgap/generate_synthetic_1000g/final_simulated_input/100samples/causal_110_0.03_str_high/summary-data/summary_stats_eqtl_tss0'
# postgap.Globals.summary_stats_eqtl = '/Users/msung/projects/p_postgap/postgap/generate_synthetic_1000g/final_simulated_input/100samples/causal_110_0.03_str_high/summary-data/summary_stats_eqtl_tss0'
postgap.Globals.ld_file = '/home/seongwonhwang/Desktop/projects/POSTGAP/p_postgap/postgap/generate_synthetic_1000g/final_simulated_input/100samples/causal_110_0.03_str_high/summary-data/ld_gwas_0'
# postgap.Globals.ld_file = '/Users/msung/projects/p_postgap/postgap/generate_synthetic_1000g/final_simulated_input/100samples/causal_110_0.03_str_high/summary-data/ld_gwas_0'
postgap.Globals.Reg_adaptors = ['S_DHS', 'S_TSS', 'S_coding']
# postgap.Globals.TYPE = 'EM'
postgap.Globals.TYPE = 'ML'
postgap.Globals.kmax_gwas = 1
postgap.Globals.kmax_eqtl = 1
postgap.Globals.DATABASES_DIR = 'databases'
# postgap.Globals.OUTPUT = 'EM/causal_110_0.03_str_high/res0'
postgap.Globals.OUTPUT = 'ML/causal_110_0.03_str_high/res0'
from postgap.Integration import * 

# res = postgap.Integration.gwas_snps_to_genes()
#################################################
# Integration
#################################################

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
cluster_associations = [(cluster, ld_snps_to_genes(
    cluster.ld_snps)) for cluster in clusters]

with open(postgap.Globals.OUTPUT+'_GWAS_lambdas.txt', 'w') as fw1:
    fw1.write('source\tlambdas\n')

with open(postgap.Globals.OUTPUT+'_eQTL_lambdas.txt', 'w') as fw2:
    fw2.write('cluster\ttissue\tgene\tlambdas\n')

#cluster_associations = postgap.FinemapIntegration.compute_gwas_posteriors(cluster_associations)
#################################################
# FinemapIntegration
#################################################

# cluster_associations = postgap.FinemapIntegration.compute_gwas_posteriors(
    # cluster_associations)

############ MLE calculation. lambda was added #############
# Perform cluster by cluster finemapping
cluster, associations = cluster_associations[0]
prepped_clusters = []
prepped_clusters.append((prepare_cluster_for_finemap(cluster, associations), associations))
cluster, associations = prepped_clusters[0]

### Finemapping.mk_modified_clusters
import pandas as pd

p_cluster = cluster 
W=[0.01, 0.1, 0.5]
pi = numpy.exp(calc_logbinom(1, postgap.Globals.kmax_gwas, len(p_cluster.z_scores)))
initial_lambdas = [0.0] * p_cluster.annotations.shape[0]
z_scores = p_cluster.z_scores
MAFs = map(float, p_cluster.mafs)
sample_size = p_cluster.gwas_snps[0].evidence[0].sample_size
approx_v = [calc_approx_v(maf, sample_size) for maf in MAFs]
logBFs = [calc_logBF_ML(approx_v[i], W, z_scores[i])for i in range(len(z_scores))]
pd_logBFs = pd.DataFrame(logBFs)
pd_logBFs.describe()



# gene_tissue_snp_eQTL_hash = organise_eQTL_data(associations)












### for entire debugging

cluster_associations = [(cluster, ld_snps_to_genes(cluster.ld_snps)) for cluster in clusters]
cluster_associations = postgap.FinemapIntegration.compute_gwas_posteriors(cluster_associations)
res = concatenate(cluster_to_genes(cluster, associations) for cluster, associations in cluster_associations)