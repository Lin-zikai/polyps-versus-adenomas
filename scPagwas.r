#################在服务器运行之前命令行运行export OPENBLAS_NUM_THREADS=1
library(ggplot2)
library(Seurat)
library(dplyr)
library(scCustomize)
library(ComplexHeatmap)
library(scPagwas)
Pagwas<-scPagwas_main(Pagwas =NULL,
                      gwas_data="colongwas.txt",
                      Single_data ='merge_seurat_verify_Ident.rds',
                      output.prefix="colon",
                      output.dirs="merge_seurat_verify",
                      Pathway_list=Genes_by_pathway_kegg,
                      n.cores=2,
                      assay="RNA",
                      singlecell=T, 
                      iters_singlecell = 100,
                      celltype=T,
                      block_annotation = block_annotation_hg37,
                      chrom_ld = chrom_ld)
save(Pagwas,file="./merge_seurat_verify_scPagwas.RData")



