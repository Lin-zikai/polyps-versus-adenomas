source("prepare.R")
setwd("S:/twas/twascrc/")

###########################################################################
# result_heng<-get_otters("S:/twas/CRC横结肠stage2")
# result_blood<-get_otters("S:/twas/CRC全血stage2")
result_heng<-readRDS("result_heng.rds")
p_heng<-0.05/nrow(result_heng)
result_heng1<-subset(result_heng,result_heng$OTTERS_p<p_heng)
twasgene<-bitr(result_heng1$Gene,"ENSEMBL","SYMBOL",OrgDb = "org.Hs.eg.db")


##########################################################################

load("S:/twas/twascrc/colon_scPagwas.RData")
varfeature<-Pagwas@assays$RNA@var.features
intersectgene<-intersect(varfeature,twasgene$SYMBOL)

setwd("S:/twas/twascrc/figure/正文")
pdf("Figure4 Venn.pdf",width = 8,height = 6)


a<-list("TWAS_Gene"=twasgene$SYMBOL,
        "scPagwas_Gene"=varfeature)
intersection_genes <- intersect(twasgene$SYMBOL, varfeature)
text_data <- data.frame(
  label = intersection_genes,
  x = rep(0.0, length(intersection_genes)), 
  y = c(-0.5,-0.3,-0.1,0.1,0.3,0.5)
)
p1=ggvenn(a,c("TWAS_Gene","scPagwas_Gene"),show_percentage = F)+
  geom_text(data = text_data, aes(x, y, label = label), size = 3.5)
p1
dev.off()


pdf("Figure4 Dotplot.pdf",width = 12,height = 8)
DotPlot(seurat_object, features = intersectgene)
dev.off()


#####################################################################################################
setwd("S:/twas/twascrc/figure/正文/")
pdf("Figure4 Barplot_intersectgene.pdf",width = 12,height =8)

intersectgene <- c("HES6","MIR4435-2HG","LIMCH1","POU2AF2","SMAD9","PITPNC1","GDF15","TGFBR2","EGFR","ERBB2","GRN","SORT1")
gene_data <- FetchData(seurat_object, vars = intersectgene)
gene_data$Cell_Type <- seurat_object$Cell_Type
gene_data$disease_status <- seurat_object$disease_status
gene_data <- subset(gene_data, gene_data$disease_status != "Other")
long_gene_data <- gene_data %>%
  pivot_longer(
    cols = intersectgene,
    names_to = "Gene",
    values_to = "Expression"
  )
test_results <- long_gene_data %>%
  group_by(Gene) %>%
  summarise(
    p.value = t.test(Expression ~ disease_status, data = .)$p.value,
    .groups = 'drop'
  ) %>%
  mutate(significance = case_when(
    p.value < 0.001 ~ '***',
    p.value < 0.01 ~ '**',
    p.value < 0.05 ~ '*',
    TRUE ~ 'ns'  # not significant
  ))

avg_expression <- long_gene_data %>%
  group_by(disease_status, Gene) %>%
  summarise(avg_expression = mean(Expression, na.rm = TRUE), .groups = "drop")
# avg_expression$avg_expression <- log(avg_expression$avg_expression + 2)
signif_positions <- avg_expression %>%
  group_by(Gene) %>%
  summarise(y.position = max(avg_expression), .groups = 'drop')
test_results <- test_results %>%
  left_join(signif_positions, by = "Gene")

p <- ggplot(avg_expression, aes(x = Gene, y = avg_expression, fill = disease_status)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7) +
  scale_fill_manual(values = c("#42B540", "#0099B3")) +
  ylab("Average Expression") +
  xlab("Gene") +
  ggtitle("") +theme_bw()+
  theme(legend.position = "bottom")
p
dev.off()

# p <- p + geom_text(data = test_results, aes(x = Gene, y = y.position, label = significance),
#                    position = position_dodge(width = 0.7), vjust = 0.5, size = 3)


#######################################################################################

group_of_cluster<-readRDS("group_cluster_marker.rds")
group_of_cluster<-group_of_cluster%>%filter(gene%in%intersectgene)
diffgene<-c("LIMA1","HES6","MIR4435-2HG","LIMCH1","POU2AF2","SMAD9","PITPNC1","GDF15","TGFBR2")
test<-data.frame()
group_of_cluster<-group_cluster_marker
for (i in diffgene){
  test1<-group_of_cluster%>%filter(gene%in%i)
  test<-rbind(test,test1)
}
heatmap_data <- reshape2::dcast(test, cluster ~ gene, value.var = "avg_log2FC")
row.names(heatmap_data) <- heatmap_data$cluster
heatmap_data$cluster <- NULL

max_value <- max(heatmap_data, na.rm = TRUE)
min_value <- min(heatmap_data, na.rm = TRUE)
breaks <- c(seq(min_value, 0, length.out = 51), seq(max_value/50, max_value, length.out = 50))


blue_palette <- colorRampPalette(c("blue", "white"))(50)
red_palette <- colorRampPalette(c("white", "red"))(50)
color_palette <- c(blue_palette, red_palette)

pheatmap(heatmap_data, color = color_palette, breaks = breaks, na_col = "grey", cluster_rows = FALSE, cluster_cols = FALSE)


###################################################################################volcano plot
setwd("S:/twas/twascrc/figure/正文/")
pdf("Figure4 volcanoplot.pdf",width = 8,height =6)
group_of_cluster<-FindMarkers(object = seurat_object,ident.1 = "SSC",ident.2 = "ASC",test.use = "wilcox", only.pos = FALSE,min.pct = 0.1)
group_of_cluster$gene<-rownames(group_of_cluster)
group_of_cluster%>%filter(gene%in%intersectgene)
volcanoplot<-group_of_cluster%>%dplyr::select(gene,P.Value=p_val_adj,logFC=avg_log2FC)
volcanoplot$change<-ifelse(volcanoplot$logFC>0,"up","down")
dataset<-volcanoplot
dataset$P.Value[dataset$P.Value==0]<-1.842251e-303
cut_off_pvalue = 0.05        #统计显著性
cut_off_logFC = 0           #差异倍数值
head(dataset)
attach(dataset)
p <- ggplot(
  # 数据、映射、颜色
  dataset, aes(x = logFC, y = -log10(P.Value), colour=change)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("#546de5","#ff4757","#d2dae2"))+
  # 辅助线
  # geom_vline(xintercept=c(-log2(1.5),log2(1.5)),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (p-value)")+
  theme_bw()+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())

dataset$label <- ifelse(dataset$gene%in%intersectgene ,as.character(dataset$gene), "")
p + geom_label_repel(data = dataset, aes(x = logFC, 
                                         y = -log10(P.Value), 
                                         label = label),
                     size = 3, box.padding = unit(0.5, "lines"),
                     point.padding = unit(0.8, "lines"), 
                     segment.color = "black", 
                     show.legend = FALSE,max.overlaps = 100000)+theme(legend.position = "bottom",legend.text = element_text(size = 12))
dev.off()

###############################
Please read pseudotime_AP.R and pseudotime_SS.R

Please read smr.R
