source("prepare.R")
########################################
setwd("S:/twas/twascrc/")
seurat_object<-readRDS("merge_seurat.rds")
Idents(seurat_object)<-seurat_object$Cell_Type
seurat_object$disease_status <- ifelse(seurat_object$disease %in% c("colon sessile serrated adenoma/polyp"),"SS", 
                                       ifelse(seurat_object$disease %in% c("tubulovillous adenoma", "tubular adenoma"), "AP", "Other"))
ABS<-subset(seurat_object,Cell_Type=="ABS")
ABS$disease_status <- ifelse(ABS$disease %in%c("colon sessile serrated adenoma/polyp") , "SS", 
                             ifelse(ABS$disease %in% c("tubulovillous adenoma", "tubular adenoma"), "AP", "Other"))
Idents(ABS)<-ABS$Cell_Type

#####################################################
load("S:/twas/twascrc/colon_scPagwas.RData")
testa<-as.data.frame(Pagwas@assays$scPagwasPaHeritability@data)
APcell<-subset(ABS,disease_status=="AP")
APcell<-colnames(APcell@assays$RNA@data)
SScell<-subset(ABS,disease_status=="SS")
SScell<-colnames(SScell@assays$RNA@data)
AP_paher<-testa[,APcell]
SS_paher<-testa[,SScell]
AP_paher1<-data.frame(mean_AP=rowSums(AP_paher))
SS_paher1<-data.frame(mean_SS=rowSums(SS_paher))
pheatplot<-cbind(AP_paher1,SS_paher1)


group<-rbind(data.frame(cell=SScell,group="SS"),data.frame(cell=APcell,group="AP"))
data<-cbind(AP_paher,SS_paher)
library(limma)
Group=factor(group$group,levels = c("SS","AP"))
design=model.matrix(~0+Group)
colnames(design)<-c("SS","AP")
design
fit<-lmFit(data,design)
contrast.matrix<-makeContrasts(SS-AP,levels = design)

fit2<-contrasts.fit(fit,contrast.matrix)
fit2<-eBayes(fit2)
diff=topTable(fit2,adjust.method = "fdr",p.value=0.05,number=5000,,sort.by = "logFC")
top5 <- diff[order(diff$logFC, decreasing = TRUE)[1:5], ]
bottom5 <- diff[order(diff$logFC, decreasing = FALSE)[1:5], ]
diff <- rbind(top5, bottom5)
diff$pathway<-rownames(diff)
for (i in 1:nrow(df)){
  df$pathwaydiff[i]<- paste0("gs",i)
}
diff<-diff%>%left_join(df,by="pathway")
diff<-diff%>%dplyr::select(pathway,description,logFC:B)


result<-rbind(top5, bottom5)
pheatplot2<-pheatplot[rownames(pheatplot)%in%rownames(result),]
setwd("S:/twas/twascrc/figure/正文")
pdf("Figure5 heatmap_scPagwas_SSCvsASC.pdf",width = 8,height = 6)
pheatmap(pheatplot2,cluster_rows = F,cluster_cols = F,angle_col = 0,labels_col = c("AP","SS"))
dev.off()

test<-df[df$pathway%in%rownames(pheatplot2),]
ABS_GSEA_result2<-ABS_GSEA_result[rownames(ABS_GSEA_result)%in%rownames(pheatplot2),]
data.frame(test$description)

test_gene<-get_genes_from_pathway(test)
test_gene<-subset(test_gene,test_gene$genename%in%twasgene$SYMBOL)
intersect(test_gene$genename,twasgene$SYMBOL)
DotPlot(ABS,features = intersect(test_gene$genename,twasgene$SYMBOL))


#####################################GSEA

setwd("S:/twas/twascrc/")
seurat_object<-readRDS("merge_seurat.rds")
Idents(seurat_object)<-seurat_object$Cell_Type
seurat_object$disease_status <- ifelse(seurat_object$disease %in% c("colon sessile serrated adenoma/polyp"),"SS", 
                                       ifelse(seurat_object$disease %in% c("tubulovillous adenoma", "tubular adenoma"), "AP", "Other"))
ABS<-subset(seurat_object,Cell_Type=="ABS")
ABS$disease_status <- ifelse(ABS$disease %in%c("colon sessile serrated adenoma/polyp") , "SS", 
                             ifelse(ABS$disease %in% c("tubulovillous adenoma", "tubular adenoma"), "AP", "Other"))
Idents(ABS)<-ABS$disease_status
##############################hsa00190,hsa04714,hsa05208, gseaplot2
de_results <- FindMarkers(ABS, ident.1 = "SS", ident.2 = "AP", min.pct = 0.25)
de_results$gene<-rownames(de_results)
ABS_SS_GSEA<-perform_GSEA(df = de_results,"avg_log2FC","gene",species = "hsa")
setwd("S:/twas/twascrc/figure/正文")
pdf("Figure5 GSEA_hsa04714_SS.pdf",width = 8,height = 6)
gseaplot2(ABS_SS_GSEA, "hsa04714", color = "firebrick", rel_heights=c(1, .2, .6),pvalue_table = F,"Thermogenesis")
dev.off()
#############################################################
de_results <- FindMarkers(ABS, ident.1 = "AP", ident.2 = "SS", min.pct = 0.25)
de_results$gene<-rownames(de_results)
ABS_AP_GSEA<-perform_GSEA(df = de_results,"avg_log2FC","gene",species = "hsa")
setwd("S:/twas/twascrc/figure/正文")
pdf("Figure5 GSEA_hsa04714_AP.pdf",width = 8,height = 6)
gseaplot2(ABS_AP_GSEA, "hsa04714", color = "firebrick", rel_heights=c(1, .2, .6),pvalue_table = F,"Thermogenesis")
dev.off()


#########################################################################
# 
ABS$Cell_Type[WhichCells(ABS, expression = disease_status == "AP")] <- "ABS.AP"
#
ABS$Cell_Type[WhichCells(ABS, expression = disease_status == "SS")] <- "ABS.SS"
ABS$Cell_Type[WhichCells(ABS, expression = disease_status %in%"Other" )] <- "ABS.Other"
unique(ABS$Cell_Type)



############################################
APcell<-subset(ABS,disease_status=="AP")
APcell<-colnames(APcell@assays$RNA@data)
SScell<-subset(ABS,disease_status=="SS")
SScell<-colnames(SScell@assays$RNA@data)
gene_data <- FetchData(ABS,c("GUCA2A","GUCY2C","GRN","SORT1","PRKAA1","COX5B","COX14","ATG10"))
#gene_data <- FetchData(ABS,intersect(test_gene$genename,twasgene$SYMBOL))
testa<-as.data.frame(Pagwas@assays$scPagwasPaHeritability@data)
df1<-df%>%filter(pathway%in%diff$pathway)
aaa<-as.data.frame(t(testa%>%filter(rownames(testa)%in%df1$pathway)))
scPagwas.result.ABS.SS<-subset(aaa,rownames(aaa)%in%SScell)
scPagwas.result.ABS.AP<-subset(aaa,rownames(aaa)%in%APcell)
test1<-subset(gene_data,rownames(gene_data)%in%rownames(scPagwas.result.ABS.AP))
test2<-subset(gene_data,rownames(gene_data)%in%rownames(scPagwas.result.ABS.SS))

# 
results_AP <- matrix(nrow = ncol(test1), ncol = ncol(scPagwas.result.ABS.AP))
p_values_AP <- matrix(nrow = ncol(test1), ncol = ncol(scPagwas.result.ABS.AP))
for (i in 1:ncol(test1)) {  
  for (j in 1:ncol(scPagwas.result.ABS.AP)) {  
    test_result <- cor.test(test1[, i], scPagwas.result.ABS.AP[, j], method = "pearson")
    results_AP[i, j] <- test_result$estimate
    p_values_AP[i, j] <- test_result$p.value  
  }  
}
results_AP_df <- as.data.frame(results_AP)
p_values_AP_df <- as.data.frame(p_values_AP)
rownames(results_AP_df) <- colnames(test1)
colnames(results_AP_df) <- df1$pathway
rownames(p_values_AP_df) <- colnames(test1)
colnames(p_values_AP_df) <- df1$pathway
results_AP_df$gene <- rownames(results_AP_df)
long_AP <- pivot_longer(results_AP_df, cols = -gene, names_to = "pathway", values_to = "correlation")
long_AP$group <- "AP"
long_AP$p_value <- c(t(p_values_AP_df))
long_AP$significance <- ifelse(long_AP$p_value < 0.05, "*", "")

results_SS <- matrix(nrow = ncol(test2), ncol = ncol(scPagwas.result.ABS.SS))
p_values_SS <- matrix(nrow = ncol(test2), ncol = ncol(scPagwas.result.ABS.SS))
for (i in 1:ncol(test2)) {  
  for (j in 1:ncol(scPagwas.result.ABS.SS)) {  
    test_result <- cor.test(test2[, i], scPagwas.result.ABS.SS[, j], method = "pearson")
    results_SS[i, j] <- test_result$estimate
    p_values_SS[i, j] <- test_result$p.value  
  }  
}
results_SS_df <- as.data.frame(results_SS)
p_values_SS_df <- as.data.frame(p_values_SS)
rownames(results_SS_df) <- colnames(test2)
colnames(results_SS_df) <- df1$pathway
rownames(p_values_SS_df) <- colnames(test2)
colnames(p_values_SS_df) <- df1$pathway
results_SS_df$gene <- rownames(results_SS_df)
long_SS <- pivot_longer(results_SS_df, cols = -gene, names_to = "pathway", values_to = "correlation")
long_SS$group <- "SS"
long_SS$p_value <- c(t(p_values_SS_df))
long_SS$significance <- ifelse(long_SS$p_value < 0.05, "*", "")

color_scheme <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
long_data <- rbind(long_AP, long_SS)
long_data<-subset(long_data,long_data$pathway%in%c("hsa00190","hsa04714","hsa05208"))
setwd("S:/twas/twascrc/figure/正文")
pdf("Figure5 heatmap_scPagwasgene.pdf",width = 8,height = 6)
ggplot(long_data, aes(x = pathway, y = gene, fill = correlation)) +
  geom_tile() +
  geom_text(aes(label = significance), color = "black", size = 3) +
  facet_wrap(~ group, ncol = 1) +
  scale_fill_gradientn(colors = rainbow(5), limits = c(-0.55, 0.3)) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


#############################################MR analysis
ao <- readRDS("S:/twas/twascrc/ao.rds")
a<-read.csv("S:/twas/GSVA/mr_all.csv")
a<-subset(a,a$method=="Inverse variance weighted"&a$pval<0.05)
a<-a%>%left_join(ao,by=c("id.outcome"="id"))%>%dplyr::select(id.outcome,b:gene,trait)


dataset<-read.csv("S:/twas/GSVA/mr_all.csv")
dataset<-subset(dataset,dataset$method=="Inverse variance weighted"&dataset$pval<0.05&dataset$gene!="EGFR")
dataset<-dataset%>%dplyr::select(Group=gene,exposure=id.outcome,b,nsnp,se,pval,Gene=gene)
dataset<-subset(dataset,dataset$exposure!="prot-a-2024")

level<-dataset$exposure
dataset$exposure<-factor(dataset$exposure)
dataset$exposure<-fct_inorder(dataset$exposure)
dataset$OR<-exp(dataset$b)
dataset$Upper<-exp(dataset$b+1.96*dataset$se)
dataset$Lower<-exp(dataset$b-1.96*dataset$se)
dataset$Group<-factor(dataset$Group)
dataset$pval<-round(dataset$pval,2)
dataset$pval[dataset$pval>=0.05]<-"< 0.05"
annotation <- data.frame(x=c(rep(-1,20),rep(-0.4,20)),
                         y=c(rep(c(19.44,19:1),2)),
                         label=c(c("nSNP",dataset$nsnp[1:8],"5,11",dataset$nsnp[11:20]),
                                 c('P Value', dataset$pval[1:8],"0.02,0.02",dataset$pval[11:20])))
setwd("S:/twas/twascrc/figure/正文")
pdf("Figure6 MR_Forestplot.pdf",width = 16,height =12)
p <- ggplot(dataset, aes(OR, exposure))+ # 
 geom_point(size=3.6,aes(col=Group)) +
  geom_errorbarh(aes(xmax =Upper, xmin = Lower,color=Group), height = 0.4) +
  scale_x_continuous(limits= c(-1,6), breaks= seq(-1, 6, 1))+
  geom_vline(aes(xintercept = 1)) +
  geom_text(data=annotation,aes(x=x,y=y,label=label))+ 
  theme_bw() + theme(legend.position ="top") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        plot.background = element_rect(fill ='#F7F7F7'),
        axis.text=element_text(size=10, face = "bold"),
        legend.text=element_text(size=11))+xlab('Odds Ratio ') + ylab("ID")
p
dev.off()