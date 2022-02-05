#!/usr/src/Python-2.7.16 python
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


def get_options():
    parser = argparse.ArgumentParser(description="""
    """, formatter_class=RawTextHelpFormatter
                                     )

    Reg_options = ["S_DHS", "S_TSS", "S_coding"]
    TYPE_options = ['binom', 'ML', 'EM', 'ML_EM']

    parser.add_argument('--bayesian', action='store_true')
    parser.add_argument('--summary_stats')
    parser.add_argument('--summary_stats_eqtl')
    parser.add_argument('--Reg', default=None,
                        nargs='*', choices=(Reg_options))
    parser.add_argument('--TYPE', default='binom',
                        type=str, choices=(TYPE_options))
    parser.add_argument('--kmax_gwas', type=int)
    parser.add_argument('--kmax_eqtl', type=int)
    parser.add_argument('--ld_file')
    parser.add_argument('--database_dir', dest='databases',
                        default='databases')
    parser.add_argument('--output')
    options = parser.parse_args()

    postgap.Globals.PERFORM_BAYESIAN = options.bayesian
    postgap.Globals.GWAS_SUMMARY_STATS_FILE = options.summary_stats
    postgap.Globals.summary_stats_eqtl = options.summary_stats_eqtl
    postgap.Globals.Reg_adaptors = options.Reg
    postgap.Globals.TYPE = options.TYPE
    postgap.Globals.kmax_gwas = options.kmax_gwas
    postgap.Globals.kmax_eqtl = options.kmax_eqtl
    postgap.Globals.ld_file = options.ld_file
    postgap.Globals.DATABASES_DIR = options.databases
    postgap.Globals.OUTPUT = options.output


    # EM debugging
    # postgap.Globals.PERFORM_BAYESIAN = True
    # postgap.Globals.GWAS_SUMMARY_STATS_FILE = '/home/seongwonhwang/Desktop/projects/POSTGAP/p_postgap/postgap/generate_synthetic_1000g/final_simulated_input/100samples/causal_110_0.03_str_high/summary-data/summary_stats_gwas_tss0'
    # postgap.Globals.summary_stats_eqtl = '/home/seongwonhwang/Desktop/projects/POSTGAP/p_postgap/postgap/generate_synthetic_1000g/final_simulated_input/100samples/causal_110_0.03_str_high/summary-data/summary_stats_eqtl_tss0'
    # postgap.Globals.Reg_adaptors = ['S_DHS','S_TSS','S_coding']
    # postgap.Globals.TYPE = 'EM'
    # postgap.Globals.kmax_gwas = 1
    # postgap.Globals.kmax_eqtl = 1
    # postgap.Globals.ld_file = '/home/seongwonhwang/Desktop/projects/POSTGAP/p_postgap/postgap/generate_synthetic_1000g/final_simulated_input/100samples/causal_110_0.03_str_high/summary-data/ld_gwas_0'
    # postgap.Globals.DATABASES_DIR = 'databases'
    # postgap.Globals.OUTPUT = 'EM/causal_110_0.03_str_high/res0'

    assert postgap.Globals.DATABASES_DIR is not None
    if not os.path.exists(os.path.dirname(postgap.Globals.OUTPUT)):
        os.makedirs(os.path.dirname(postgap.Globals.OUTPUT))
    return options


def main():
    options = get_options()
    res = postgap.Integration.gwas_snps_to_genes()
    pickle.dump(res, open(postgap.Globals.OUTPUT + '.pkl', "w"))


if __name__ == "__main__":
    main()