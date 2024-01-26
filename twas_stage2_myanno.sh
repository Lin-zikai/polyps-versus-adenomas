#!/bin/bash

conda activate lin
export PATH=/GPUFS/gyfyy_jxhe_1/User/lin/projects/bin:$PATH 



## gwas文件名称
file=crcgwas.txt


mkdir /GPUFS/gyfyy_jxhe_1/User/zyt/0821stage2_myanno/${file}


##核数
N_THREADS=10

##修改为实际程序文件夹
OTTERS_DIR=/GPUFS/gyfyy_jxhe_1/User/lin/OTTERS
SDPR_DIR=/GPUFS/gyfyy_jxhe_1/User/lin/SDPR


#不用改====================================================
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${SDPR_DIR}/MKL/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${SDPR_DIR}/gsl/lib
export PATH=/GPUFS/gyfyy_jxhe_1/User/lin/projects/bin:$PATH 

export MKL_NUM_THREADS=$N_THREADS
export NUMEXPR_NUM_THREADS=$N_THREADS
export OMP_NUM_THREADS=$N_THREADS


#======================================

## 修改文件存放文件夹===================

exp_anno=/GPUFS/gyfyy_jxhe_1/User/lin/crc_hg38_anno.txt

gwas_sst_file=/GPUFS/gyfyy_jxhe_1/User/lin/${file}

Exp_geno_chr=/GPUFS/gyfyy_jxhe_1/User/lin/hg38/1000G_hg38_chr

##================================
clump_r2=0.99


## stage1文件夹路径
out_dir=/GPUFS/gyfyy_jxhe_1/User/zyt/colon/stage1

for chr in {1..22}; do
  
  geno_dir=${Exp_geno_chr}${chr}

  twas_dir=/GPUFS/gyfyy_jxhe_1/User/zyt/0821stage2_myanno/${file}/TWAS_chr${chr}

  python3 ${OTTERS_DIR}/testing.py \
  --OTTERS_dir=${OTTERS_DIR} \
  --weight_dir=${out_dir} \
  --models=P0.001,P0.05,lassosum,SDPR,PRScs \
  --anno_dir=${exp_anno} \
  --geno_dir=${geno_dir} \
  --out_dir=${twas_dir} \
  --gwas_file=${gwas_sst_file} \
  --chrom=${chr} \
  --thread=$N_THREADS;
done

