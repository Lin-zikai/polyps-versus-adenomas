source("S:/twas/5.8.2cyclemr.r")
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(KEGGREST)
library(stringr)
library(pheatmap)
library(scPagwas)
library(ggvenn)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(Seurat)
library(tidyr)
library(GSVA)
library(xlsx)
library(futile.logger)
library(forcats)
library(AnnotationDbi)
library(enrichplot)
library(R.utils)
library(readxl)
library(limma)
library(mice)
library(ggpubr)
library(gridExtra)
library(future)
library(cli)
library(mr.raps)
library(TwoSampleMR)
library(fs) 
library(lubridate) 
library(ieugwasr)
library(plinkbinr) 
library(progress)
library(MRPRESSO)
library(parallel)
library(foreach)
library(doParallel)
library(pbapply)
library(cause)
library(vroom)
library(MungeSumstats)
library(meta)
library(readr)
library(forestploter)
library(ldscr)
library(MRlap)
library(meta)
library(forestplot)
library(coloc)
library(ACAT)
library(biomaRt)
library(Rfast)
suppressMessages(if(!require(CellChat))devtools::install_github("sqjin/CellChat"))
suppressMessages(if(!require(patchwork))install.packages('patchwork') )
suppressMessages(if(!require(ggalluvial))install.packages('ggalluvial'))
suppressMessages(if(!require(igraph))install.packages('igraph'))
suppressMessages(options(stringsAsFactors = FALSE))
suppressMessages(options(futrue.globlas.Maxsize=2*1024**3))



get_otters<-function(dir,chr=c(1:22),eqtlgen=F,
                     method=c("P0.001","P0.05","lassosum","PRScs","SDPR")){
  setwd(dir)
  dir_r<-paste0(dir,'/TWAS_chr',chr)
  dir_e<-lapply(dir_r,dir_info)
  names(dir_e)<-c(1:length(dir_e))
  dir_e<-clean_list(dir_e,0)
  chr<-names(dir_e)
  dir_r<-paste0(dir,'/TWAS_chr',chr)
  all_m<-c("P0.001","P0.05","lassosum","PRScs","SDPR")
  tofill<-all_m[!(all_m%in%method)]
  if(length(tofill)>0)cli_progress_step("填充空白的，没有运行方法的文件")
  fill_file<-function(tofill,dir_r){
    
    x<-data.frame(CHROM=NA,	GeneStart=NA,	GeneEnd=NA,	TargetID=NA,
                  n_snps=NA,	FUSION_Z=NA,	FUSION_PVAL=NA)
    
    wfile<-paste0(dir_r,"/",tofill,".txt")
    lapply(wfile,x,FUN=function(x,y)fwrite(y,x,sep="\t"))
  }
  lapply(tofill,dir_r,FUN=fill_file)
  
  read_result<-function(dir_r,eqtlgen){
    setwd(dir_r)
    a<-fread('lassosum.txt')
    b<-fread('P0.001.txt')
    c<-fread('P0.05.txt')
    d<-fread('PRScs.txt')
    e<-fread('SDPR.txt')
    if(eqtlgen){
      names(a)<-c('TargetID','CHROM','GeneStart','GeneEnd','n_snps',   'FUSION_Z', 'FUSION_PVAL')
      names(b)<-c('TargetID','CHROM','GeneStart','GeneEnd','n_snps',   'FUSION_Z', 'FUSION_PVAL')
      names(c)<-c('TargetID','CHROM','GeneStart','GeneEnd','n_snps',   'FUSION_Z', 'FUSION_PVAL')
      names(d)<-c('TargetID','CHROM','GeneStart','GeneEnd','n_snps',   'FUSION_Z', 'FUSION_PVAL')
      names(e)<-c('TargetID','CHROM','GeneStart','GeneEnd','n_snps',   'FUSION_Z', 'FUSION_PVAL')
    }
    l<-which(a$FUSION_PVAL=="NaN" | a$FUSION_PVAL=="nan" | a$FUSION_PVAL=='0.e+')
    a$FUSION_PVAL[l]<-NA
    a$FUSION_Z[l]<-NA
    
    l<-which(b$FUSION_PVAL=="NaN" | b$FUSION_PVAL=="nan" | b$FUSION_PVAL=='0.e+')
    b$FUSION_PVAL[l]<-NA
    b$FUSION_Z[l]<-NA
    
    l<-which(c$FUSION_PVAL=="NaN" | c$FUSION_PVAL=="nan"| c$FUSION_PVAL=='0.e+')
    c$FUSION_PVAL[l]<-NA
    c$FUSION_Z[l]<-NA
    
    l<-which(d$FUSION_PVAL=="NaN" | d$FUSION_PVAL=="nan" | d$FUSION_PVAL=='0.e+')
    d$FUSION_PVAL[l]<-NA
    d$FUSION_Z[l]<-NA
    
    l<-which(e$FUSION_PVAL=="NaN" | e$FUSION_PVAL=="nan" | e$FUSION_PVAL=='0.e+')
    e$FUSION_PVAL[l]<-NA
    e$FUSION_Z[l]<-NA
    
    
    a<-data.frame(Gene=a$TargetID,lassosum_z=a$FUSION_Z,lassosum_p=a$FUSION_PVAL)
    b<-data.frame(Gene=b$TargetID,PT_0.001_z=b$FUSION_Z,PT_0.001_p=b$FUSION_PVAL)
    c<-data.frame(Gene=c$TargetID,PT_0.05_z=c$FUSION_Z,PT_0.05_p=c$FUSION_PVAL)
    d<-data.frame(Gene=d$TargetID,PRScs_z=d$FUSION_Z,PRScs_p=d$FUSION_PVAL)
    e<-data.frame(Gene=e$TargetID,SDPR_z=e$FUSION_Z,SDPR_p=e$FUSION_PVAL)
    
    
    result<-merge(a,b,by='Gene',all=T)
    result<-merge(result,c,by='Gene',all=T)
    result<-merge(result,d,by='Gene',all=T)
    result<-merge(result,e,by='Gene',all=T)
    result<-result[!(is.na(result$Gene)),]
    sum_na<-t(apply(result,1,is.na))
    sum_T<-data.frame(t(apply(sum_na,1,table)))
    nl<-which(sum_T$`TRUE.`>9)
    if(length(nl)>0)result<-result[-nl,]
    
    #colnames(result)[3:4]<-c('P+T(0.001)','P+T(0.05)')
    for( i in 2:11){result[,i]<-as.numeric(result[,i])}
    
    p<-result%>%select(lassosum_p,PT_0.001_p,PT_0.05_p,PRScs_p,SDPR_p)
    p<-split(p,seq(nrow(p)))
    l<-unlist(lapply(p,function(x){
      which.min(x)
    }))
    for(i in 1:nrow(result)){
      if(!is.na(l[i]) && l[i] == 1){
        result$z_minp[i]<-result$lassosum_z[i]
        result$z_from[i]<-"lassosum"
      }
      if(!is.na(l[i]) && l[i] == 2){
        result$z_minp[i]<-result$PT_0.001_z[i]
        result$z_from[i]<-"P+T(0.001)"
      }
      if(!is.na(l[i]) && l[i] == 3){
        result$z_minp[i]<-result$PT_0.05_z[i]
        result$z_from[i]<-"P+T(0.05)"
      }
      if(!is.na(l[i]) && l[i] == 4){
        result$z_minp[i]<-result$PRScs_z[i]
        result$z_from[i]<-"PRScs"
      }
      if(!is.na(l[i]) && l[i] == 5){
        result$z_minp[i]<-result$SDPR_z[i]
        result$z_from[i]<-"SDPR"
      }
    }
    cli_progress_step(paste0("完成第",dir_r,"的读取"))
    return(result)
  }
  cli_progress_step(paste0("正在读取",dir,"的OTTERS结果文件\n"))
  
  result<-pblapply(dir_r,eqtlgen, FUN=read_result)
  acat<-function(result){
    id<-result$Gene
    for_acat<-split(result,seq(nrow(result)))
    names(for_acat)<-id
    for_acat<-lapply(for_acat,FUN=function(x)unlist(c(x[,c(3,5,7,9,11)])))
    for_acat<-lapply(for_acat,na.omit)
    res<-lapply(for_acat,ACAT)
    res<-data.frame(Gene=id,OTTERS_p=unlist(res))
    res<-merge(res,result,by="Gene")
    res$OTTERS_fdr<-p.adjust(res$OTTERS_p,'fdr')
    res$OTTERS_bonferroni<-p.adjust(res$OTTERS_p,'bonferroni')
    return(res)
  }
  cli_progress_step(paste0("\n正在对",dir,"运行ACAT-O\n"))
  
  res<-pblapply(result,acat)
  
  anno<-function(dir_r,eqtlgen){
    setwd(dir_r)
    a<-fread('lassosum.txt')
    l<-which(a$FUSION_PVAL=="NaN")
    a$FUSION_PVAL[l]<-NA
    b<-fread('P0.001.txt')
    c<-fread('P0.05.txt')
    d<-fread('PRScs.txt')
    e<-fread('SDPR.txt')
    a<-a[,1:4]
    b<-b[,1:4]
    c<-c[,1:4]
    d<-d[,1:4]
    e<-e[,1:4]
    toint<-function(x){
      x$CHROM<-as.character( x$CHROM)
      x$GeneStart<-as.character( x$GeneStart)
      x$GeneEnd<-as.character( x$GeneEnd)
      x$TargetID<-as.character(x$TargetID)
      return(x)
    }
    
    
    a<-toint(a)
    b<-toint(b)
    c<-toint(c)
    d<-toint(d)
    e<-toint(e)
    
    if(eqtlgen==T){
      names(a)<-c('TargetID','CHROM','GeneStart','GeneEnd')
      names(b)<-c('TargetID','CHROM','GeneStart','GeneEnd')
      names(c)<-c('TargetID','CHROM','GeneStart','GeneEnd')
      names(d)<-c('TargetID','CHROM','GeneStart','GeneEnd')
      names(e)<-c('TargetID','CHROM','GeneStart','GeneEnd')
    }
    
    
    dat<-merge(a,b,by=c('CHROM','GeneStart','GeneEnd','TargetID'),all=T)
    dat<-merge(dat,c,by=c('CHROM','GeneStart','GeneEnd','TargetID'),all=T)
    dat<-merge(dat,d,by=c('CHROM','GeneStart','GeneEnd','TargetID'),all=T)
    dat<-merge(dat,e,by=c('CHROM','GeneStart','GeneEnd','TargetID'),all=T)
    colnames(dat)<-c('CHROM','GeneStart','GeneEnd','Gene')
    return(dat)
  }
  cli_progress_step(paste0("\n正在对",dir,"进行处理\n"))
  
  res_anno<-pblapply(dir_r,eqtlgen,FUN=anno)
  res_all<-data.frame()
  for(i in 1:length(chr)){
    res1<-merge(res_anno[[i]],res[[i]],by="Gene")
    res_all<-rbind(res_all,res1)
  }
  if(nrow(res_all)>0)res_all$trait<-basename(dir)
  
  cli_progress_step("完成")
  return(res_all)
}

get_pathway<-function(x){
  org <- keggList('organism') 
  hsa_path <- keggLink("pathway",x) 
  meta= unique(hsa_path)[grepl(paste(x,"0",sep = ""),unique(hsa_path))]
  hsa_info <- lapply(meta, keggGet)  
  nm=unlist(lapply( hsa_info , function(x) x[[1]]$NAME))
  genes <- unlist(lapply(hsa_info, function(x) {
    if (length(x[[1]]$GENE) <= 1) {
      return(NA)  # 如果没有基因信息，返回NA
    } else {
      g <- x[[1]]$GENE
      if (length(g) > 1) {  # 如果有足够的基因信息，继续处理
        # 请注意，这里我们假设'g'的结构是每两项一组，第一项是基因名，第二项是基因描述
        gene_info <- g[seq(2, length(g), by = 2)]
        # 分隔基因信息，并只获取基因名
        gene_names <- sapply(strsplit(gene_info, ';'), `[`, 1)
        # 合并基因名
        return(paste(gene_names, collapse = ';'))
      } else {
        return(NA)  # 如果'g'的结构不是预期的，返回NA
      }
    }
  }))
  df =data.frame(
    pathway= str_extract(meta, "(?<=path:)hsa\\d+"),
    description=nm,
    genes =genes
  )
}


get_genes_from_pathway<-function(df){
  cc<-data.frame()
  for (i in 1:nrow(df)){
    aa<-str_split(df$genes[i],";")
    genename<-aa[[1]]
    pathway<-df$pathway[i]
    bb<-data.frame(genename=genename,pathway=pathway)
    cc<-rbind(cc,bb)
  }
  return(cc)
}

perform_gseKEGG <- function(df, id_col, species = "hsa") {
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
  library(KEGGREST)
  library(stringr)
  # 确定OrgDb和toType
  if(species == "hsa") {
    OrgDb <- "org.Hs.eg.db"
    organism<-"hsa"
  } else if(species == "mmu") {
    OrgDb <- "org.Mm.eg.db"
    organism<-"mmu"
  } else {
    stop("Invalid species. Please choose 'hsa' or 'mmu'.")
  }
  
  # 确定ID类型
  if(all(str_detect(df[[id_col]], "ENSG"))) {
    id_type <- "ENSEMBL"
  } else if(all(str_detect(df[[id_col]], "^\\d+$"))) {
    id_type <- "ENTREZID"
  } else {
    id_type <- "SYMBOL"
  }
  
  # 如果不是ENTREZID，进行转换
  if(id_type != "ENTREZID") {
    df_id <- bitr(df[[id_col]], 
                  fromType = id_type,
                  toType = "ENTREZID",
                  OrgDb = OrgDb)
    
    if(id_type == "SYMBOL") {
      df <- merge(df_id,df, by.x ="SYMBOL" , by.y =id_col , all.x = TRUE)
    } else if (id_type == "ENSEMBL") {
      df <- merge(df_id,df, by.x ="ENSEMBL" , by.y =id_col , all.x = TRUE)
    }
  }
  
  gene_list <- df$ENTREZID
  
  kegg <- enrichKEGG(gene_list, 
                     organism=organism, 
                     pvalueCutoff=0.05, 
                     #pAdjustMethod="BH",
                     keyType="kegg") #pvaluecutoff 是pvalue的阈值，显著富集性要<0.01
  kegg <- setReadable(kegg, OrgDb=OrgDb,keyType = 'ENTREZID')
  
  return(kegg@result)
}

prepare_GSEA_data <- function(df, id_col, species = "hsa") {
  # 确定OrgDb和toType
  if(species == "hsa") {
    OrgDb <- "org.Hs.eg.db"
  } else if(species == "mmu") {
    OrgDb <- "org.Mm.eg.db"
  } else {
    stop("Invalid species. Please choose 'hsa' or 'mmu'.")
  }
  
  # 确定ID类型
  id_type <- ifelse(str_detect(id_col, "ENSG"), "ENSEMBL", 
                    ifelse(str_detect(id_col, "^\\d+$"), "ENTREZID", "SYMBOL"))
  
  # 如果不是ENTREZID，进行转换
  if(id_type != "ENTREZID") {
    df_id <- bitr(df[[id_col]], 
                  fromType = id_type,
                  toType = "ENTREZID",
                  OrgDb = OrgDb)
    df <- merge(df, df_id, by = id_col, all = FALSE)
  }
  
  # 按logFC进行排序
  df_sort <- df[order(df$log2FC, decreasing = TRUE),]
  gene_fc <- df_sort$log2FC
  names(gene_fc) <- df_sort$ENTREZID
  
  return(gene_fc)
}

perform_GSEA <- function(df, logFC_col, gene_col, species = "hsa") {
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(stringr)
  library(org.Mm.eg.db)
  library(clusterProfiler)
  library(enrichplot)
  library(R.utils)
  # 确定OrgDb和toType
  if(species == "hsa") {
    OrgDb <- org.Hs.eg.db
    organism<-"hsa"
  } else if(species == "mmu") {
    OrgDb <- org.Mm.eg.db
    organism<-"mmu"
  } else {
    stop("Invalid species. Please choose 'hsa' or 'mmu'.")
  }
  
  # 确定ID类型
  if(all(str_detect(df[[gene_col]], "ENSG"))) {
    id_type <- "ENSEMBL"
  } else if(all(str_detect(df[[gene_col]], "^\\d+$"))) {
    id_type <- "ENTREZID"
  } else {
    id_type <- "SYMBOL"
  }
  
  # 如果不是ENTREZID，进行转换
  if(id_type != "ENTREZID") {
    df_id <- bitr(df[[gene_col]], 
                  fromType = id_type,
                  toType = "ENTREZID",
                  OrgDb = OrgDb)
    
    if(id_type == "SYMBOL") {
      df <- merge(df_id,df, by.x ="SYMBOL" , by.y =gene_col , all.x = TRUE)
    } else if (id_type == "ENSEMBL") {
      df <- merge(df_id,df, by.x ="ENSEMBL" , by.y =gene_col , all.x = TRUE)
    }
  }
  
  # 提取基因列表和logFC
  gene_list <- df[[logFC_col]]
  names(gene_list) <- df$ENTREZID
  
  # 对基因列表按logFC排序
  gene_list_sorted <- sort(gene_list, decreasing = TRUE)
  
  # 进行gseKEGG分析
  kegg <- gseKEGG(geneList     = gene_list_sorted,
                  organism     = organism,
                  nPerm        = 1000,    # 可以根据需要调整排列次数
                  minGSSize    = 10,      # 富集分析的最小基因集大小
                  maxGSSize    = 500,     # 富集分析的最大基因集大小
                  pvalueCutoff = 1,    # p-value的阈值
                  verbose      = FALSE)
  
  # 设置可读的基因名
  kegg <- setReadable(kegg, OrgDb = OrgDb, keyType = 'ENTREZID')
  
  return(kegg)
}

df<-get_pathway("hsa")
genelist <- list()

for (i in seq(length(unique(df$pathway)))){
  
  genelist[i] <- df$genes[i]
  
}

names(genelist) <- paste0("gs",c(1:length(unique(df$pathway))))
genelist<- lapply(genelist,function(x) unlist(strsplit(x,split = ";")))

mid<-data.frame()
for (i in 1:nrow(df)){
  aa<-str_split(df$genes[i],";")
  genename<-aa[[1]]
  pathway<-df$pathway[i]
  bb<-data.frame(genename=genename,pathway=pathway)
  mid<-rbind(mid,bb)
}

