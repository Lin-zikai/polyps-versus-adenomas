library(monocle3)
library(ggplot2)
library(plyr)
library(GSVA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(KEGGREST)
library(msigdbr)
library(Seurat)
library(dplyr)
load("S:/twas/twascrc/colon_scPagwas.RData")
setwd("S:/twas/twascrc/")
seurat_object<-readRDS("merge_seurat.rds")
Idents(seurat_object)<-seurat_object$Cell_Type
seurat_object$disease_status <- ifelse(seurat_object$disease %in% "colon sessile serrated adenoma/polyp", "SS", 
                                       ifelse(seurat_object$disease %in% c("tubulovillous adenoma", "tubular adenoma"), "AP", "Other"))

setwd("S:/twas/GSVA/")

seurat_object<-subset(seurat_object,Cell_Type!="SSC")
seurat_object<-subset(seurat_object,disease_status=="AP")
AP_cds <- new_cell_data_set(as(as.matrix(seurat_object@assays$RNA@counts), 'sparseMatrix'),
                            cell_metadata = as.data.frame(seurat_object@meta.data),
                            gene_metadata = data.frame(gene_short_name = row.names(seurat_object), 
                                                       row.names = row.names(seurat_object))
)
AP_cds  <- preprocess_cds(AP_cds , num_dim = 50)
AP_cds  <- reduce_dimension(AP_cds )
AP_cds  <- cluster_cells(AP_cds )
AP_cds  <- learn_graph(AP_cds)
saveRDS(AP_cds,"AP_cds.rds")

cds<-readRDS("AP_cds.rds")
myselect <- function(cds,select.classify,my_select){
  cell_ids <- which(colData(cds)[,select.classify] == my_select)
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

AFD_lineage_cds1 <- order_cells(cds, root_pr_nodes=myselect(cds,select.classify = 'Cell_Type',my_select = "STM"))

cds_subset <- AFD_lineage_cds1
colData(cds_subset)$pseudotime <- pseudotime(cds_subset)

data_matrix <- GetAssayData(Pagwas@assays$scPagwasPaHeritability, slot = "data")

cds_exprs <- reshape2::melt(as.matrix(data_matrix))
colnames(cds_exprs) <- c("f_id", "Cell", "expression")

cds_exprs <- cds_exprs[ cds_exprs$f_id == "hsa04350",]

cds_colData <- colData(cds_subset)
cds_colData$Cell <- rownames(cds_colData)
# cds_colData <- cds_colData[cds_colData$Cell_Type == "ABS",]
cds_colData <- cds_colData[,c("pseudotime","Cell")]
cds_exprs <- merge(cds_exprs, cds_colData, by= "Cell")

cds_exprs <- as.data.frame(cds_exprs)
# saveRDS(cds_exprs,"cds_exprs.rds")
vertical_jitter = NULL
horizontal_jitter = NULL


##############################http://127.0.0.1:34113/graphics/plot_zoom_png?width=1200&height=867
disease<-data.frame(disease=seurat_object$disease_status,Cell=names(seurat_object$disease_status),Cell_Type=seurat_object$Cell_Type)
cds_exprs<-merge(cds_exprs,disease,by="Cell")
cds_exprs<-subset(cds_exprs,cds_exprs$Cell_Type!="SSC")
cds_exprs <- cds_exprs[!is.infinite(cds_exprs$pseudotime), ]
cds_exprs$expression <- scale(cds_exprs$expression)
cds_exprs$pseudotime <- scale(cds_exprs$pseudotime)

cell_type_colors <- c("ASC" = "#1f77b4", "ABS" = "#ff7f0e", "GOB" = "#2ca02c", 
                      "TAC" = "#d62728", "SSC" = "#9467bd", "CT" = "#8c564b", 
                      "STM" = "#e377c2", "TUF" = "#7f7f7f", "EE" = "#bcbd22")
disease_colors <- c("SS" = "pink")

# 绘制图形
ggplot(cds_exprs, aes(x = pseudotime, y = expression)) + 
  geom_point(aes(color = Cell_Type)) + # 散点图层，按 Cell_Type 分类
  geom_smooth(aes(x = pseudotime, y = expression, group = disease, color = disease), method = "loess", se = FALSE) + # 按 disease 分组添加平滑曲线
  scale_color_manual(values = c(cell_type_colors, disease_colors)) + # 手动指定颜色
  theme_minimal() + # 使用简洁的主题
  labs(x = "Pseudotime", y = "Expression", color = "Cell Type / Disease") # 设置坐标轴和图例的标签

##########################################################
setwd("S:/twas/twascrc/figure/正文")
pdf("Figure4 pseudotime_AP.pdf",width = 8,height = 6)
intersectgene <- c("MIR4435-2HG","PITPNC1","GDF15","TGFBR2","GUCA2A")
plot_genes_in_pseudotime(AFD_lineage_cds1[intersectgene,],
                         color_cells_by="Cell_Type",
                         min_expr=0.5)
dev.off()