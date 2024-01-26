source("prepare.R")

setwd("S:/twas/twascrc/")
seurat_object<-readRDS("merge_seurat.rds")
Idents(seurat_object)<-seurat_object$Cell_Type
seurat_object$disease_status <- ifelse(seurat_object$disease %in% c("colon sessile serrated adenoma/polyp"),"SS", 
                 
                                                             ifelse(seurat_object$disease %in% c("tubulovillous adenoma", "tubular adenoma"), "AP", "Other"))
###############################################################################################
setwd("S:/twas/twascrc/figure/正文/")
pdf("Figure2 UMAP.pdf",width = 8,height = 6)
umap =seurat_object@reductions$pca@cell.embeddings %>% as.data.frame()
umap$cell<-rownames(umap)
cell_type<-data.frame(cell_type=seurat_object$Cell_Type,cell=names(seurat_object$Cell_Type))
umap<-merge(umap,cell_type,by="cell")
colour=c("#9370DB","#20B2AA","#7CFC00","#1E90FF","#DC143C","#F08080","#575554","#40E0D0","#FFA500")
p <- ggplot(umap,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_type)) +  geom_point(size = 1 , alpha =1 )  +  scale_color_manual(values = colour)
p2 <- p  +
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"))
p3 <- p2 +         
  theme(
    legend.title = element_blank(), #去掉legend.title 
    legend.key=element_rect(fill='white'), #
    legend.text = element_text(size=20), #设置legend标签的大小
    legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5))) #设置legend中 点的大小 
offset_x <- min(umap$UMAP_1) - (range(umap$UMAP_1)[2] - range(umap$UMAP_1)[1]) * 0.08
offset_y <- min(umap$UMAP_2) - (range(umap$UMAP_2)[2] - range(umap$UMAP_2)[1]) * 0.08
p4 <- p3 + 
  geom_segment(aes(x = offset_x, y = offset_y,
                   xend = offset_x +2, yend = offset_y),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = offset_x, y = offset_y,
                   xend = offset_x, yend = offset_y + 2),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = offset_x +1, y = offset_y -1, label = "UMAP_1",
           color="black",size = 3, fontface="bold" ) + 
  annotate("text", x = offset_x -1, y = offset_y + 1, label = "UMAP_2",
           color="black",size = 3, fontface="bold" ,angle=90)
cell_type_med <- umap %>%
  group_by(cell_type) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )
library(ggrepel)
p4 +
  geom_label_repel(aes(label=cell_type), fontface="bold",data = cell_type_med,
                   point.padding=unit(0.5, "lines")) +
  theme(legend.position = "none")
rm(p1)
rm(p2)
rm(p3)
dev.off()

###########################################################################Bootstrap
load("S:/twas/twascrc/coloncancer_scPagwas/coloncancer_scPagwas.RData")
setwd("S:/twas/twascrc/figure/正文")
pdf("Figure2 Bootstrap.pdf",width = 8,height = 6)
Bootstrap_estimate_Plot(bootstrap_results=Pagwas@misc$bootstrap_results,
                        figurenames = "estimateplot.pdf",
                        width = 9,
                        height = 7,
                        do_plot=T)
dev.off()

###########################################GSEA
group_of_cluster<-FindMarkers(object = seurat_object,ident.1 = "SSC",ident.2 = "ASC",test.use = "wilcox", only.pos = FALSE,min.pct = 0.1)
group_of_cluster$gene<-rownames(group_of_cluster)
SSC<-perform_GSEA(df = group_of_cluster,"avg_log2FC","gene",species = "hsa")
SSC_GSEA_result<-SSC@result
SSC_GSEA_result<-SSC_GSEA_result%>%filter(pvalue<0.05)
SSC_GSEA_result%>%filter(ID=="hsa04350")
mid<-data.frame()
for (i in 1:nrow(SSC_GSEA_result)){
  aa<-str_split(SSC_GSEA_result$core_enrichment[i],"/")
  genename<-aa[[1]]
  pathway<-SSC_GSEA_result$ID[i]
  bb<-data.frame(genename=genename,pathway=pathway)
  mid<-rbind(mid,bb)
}
mid<-subset(mid,mid$genename%in%twasgene$SYMBOL)

setwd("S:/twas/twascrc/figure/正文")
pdf("Figure2 GSEA_SSC_TGFβ.pdf",width = 8,height = 6)
gseaplot2(SSC, "hsa04080",color = "firebrick", rel_heights=c(1, .2, .6),pvalue_table = T,title="TGF-beta signaling pathway in SSC")
dev.off()

group_of_cluster<-FindMarkers(object = seurat_object,ident.1 = "ASC",ident.2 = "SSC",test.use = "wilcox", only.pos = FALSE,min.pct = 0.1)
group_of_cluster$gene<-rownames(group_of_cluster)
ASC<-perform_GSEA(df = group_of_cluster,"avg_log2FC","gene",species = "hsa")
pdf("Figure2 GSEA_ASC_TGFβ.pdf",width = 8,height = 6)
gseaplot2(ASC, "hsa04350",color = "firebrick", rel_heights=c(1, .2, .6),pvalue_table = F,title="TGF-beta signaling pathway in ASC")
dev.off()


Please Read cellchat_SS_AP.R 




