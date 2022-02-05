#!/bin/bash
date

dirname=$1
outdir=$2
kmax_gwas=$3
kmax_eqtl=$4
type=$5

inputdir='/home/seongwonhwang/Desktop/projects/POSTGAP/p_postgap/postgap/generate_synthetic_1000g/final_simulated_input/'$dirname'/summary-data'
#inputdir='/Users/msung/projects/p_postgap/postgap/generate_synthetic_1000g/final_simulated_input/'$dirname'/summary-data'

if [[ -d $outdir ]]; then
	rm -r $outdir
fi
if [[ ! -d $outdir ]]; then
        mkdir -p $outdir
fi

ls $inputdir/summary_stats_gwas_tss* |while read FILE_GWAS
 do
	IDX=$(basename $FILE_GWAS)
	IDX=${IDX#summary_stats_gwas_tss}

 python POSTGAP.py --bayesian --summary_stats $FILE_GWAS \
                   --summary_stats_eqtl $inputdir/summary_stats_eqtl_tss$IDX \
                   --Reg S_DHS S_TSS S_coding \
                   --TYPE $type \
                   --kmax_gwas $kmax_gwas \
                   --kmax_eqtl $kmax_eqtl \
                   --ld_file $inputdir/ld_gwas_$IDX \
                   --database_dir 'databases' \
                   --output $outdir/res$IDX 

done

date
