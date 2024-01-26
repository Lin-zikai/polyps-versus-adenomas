#!/bin/bash

export PATH=~/User/lin/projects/bin:$PATH
N_THREADS=14

OTTERS_DIR=/GPUFS/gyfyy_jxhe_1/User/lin/OTTERS
SDPR_DIR=/GPUFS/gyfyy_jxhe_1/User/lin/SDPR

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${SDPR_DIR}/MKL/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${SDPR_DIR}/gsl/lib
export PATH=~/User/lin/projects/bin:$PATH

export NUMEXPR_MAX_THREADS=64

export MKL_NUM_THREADS=$N_THREADS
export NUMEXPR_NUM_THREADS=$N_THREADS
export OMP_NUM_THREADS=$N_THREADS


exp_anno=Colon_hg38_anno.txt
sst_file=Colon_hg38.txt

eqtl_dir=colon

Exp_geno_chr=hg38/1000G_hg38_chr

clump_r2=0.99

for chr in {1..22}; do
  out_dir=${eqtl_dir}/Results_chr${chr}
  geno_dir=${Exp_geno_chr}${chr}

  python3 ${OTTERS_DIR}/training.py \
  --OTTERS_dir=${OTTERS_DIR} \
  --SDPR_dir=${SDPR_DIR} \
  --anno_dir=${exp_anno} \
  --geno_dir=${geno_dir} \
  --sst_file=${sst_file} \
  --out_dir=${out_dir} \
  --chrom=${chr} \
  --r2=${clump_r2} \
  --models=PT,lassosum,SDPR,PRScs \
  --thread=$N_THREADS;
done

