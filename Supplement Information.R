#FigureS1

library(CellChat)
library(patchwork)
library(cowplot)
setwd('S:/twas/cellchat/')
cellchat.SS <- readRDS("cellchat2_SS_computeCommunProb.rds")
cellchat.SS  <- filterCommunication(cellchat.SS , min.cells = 10)
cellchat.SS <- computeCommunProbPathway(cellchat.SS)
cellchat.SS <- aggregateNet(cellchat.SS)
cellchat.AP <- readRDS('cellchat2_AP_computeCommunProb.rds')
cellchat.AP  <- filterCommunication(cellchat.AP , min.cells = 10)
cellchat.AP <- computeCommunProbPathway(cellchat.AP)
cellchat.AP <- aggregateNet(cellchat.AP)
object.list <- list(SS = cellchat.SS, AP = cellchat.AP)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
setwd("S:/twas/twascrc/figure/SI")
pdf("FigureS2 Number_interactions.pdf",width = 12,height =8)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

#Figure S2
load("S:/twas/twascrc/coloncancer_scPagwas/coloncancer_scPagwas.RData")
seurat_object<-readRDS("merge_seurat.rds")
Idents(seurat_object)<-seurat_object$Cell_Type
seurat_object$disease_status <- ifelse(seurat_object$disease %in% c("colon sessile serrated adenoma/polyp"),"SS", 
                                       
                                       ifelse(seurat_object$disease %in% c("tubulovillous adenoma", "tubular adenoma"), "AP", "Other"))

intersectgene <- c("HES6","MIR4435-2HG","LIMCH1","POU2AF2","SMAD9","PITPNC1","GDF15","TGFBR2")
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
ssc_cells <- subset(seurat_object, idents = "SSC")#ssc_cells <- subset(seurat_object, idents = "ASC")
gene_data <- FetchData(ssc_cells, vars = intersectgene)
testa<-as.data.frame(Pagwas@assays$scPagwasPaHeritability@data)
hsa<-data.frame(t(testa[rownames(testa)%in%c("hsa04350","hsa04012"),]))
hsa<-subset(hsa,rownames(hsa)%in%rownames(gene_data))
gene_data<-cbind(gene_data,hsa)
intersectgene <- c("HES6","MIR4435-2HG","LIMCH1","POU2AF2","SMAD9","PITPNC1","GDF15","TGFBR2","hsa04350","hsa04012")
n <- length(intersectgene)
correlation_matrix <- matrix(nrow = n, ncol = n)
p_value_matrix <- matrix(nrow = n, ncol = n)
rownames(correlation_matrix) <- intersectgene
rownames(p_value_matrix) <- intersectgene
colnames(correlation_matrix) <- intersectgene
colnames(p_value_matrix) <- intersectgene

for (i in 1:n) {
  for (j in 1:n) {
    result <- cor.test(gene_data[[intersectgene[i]]], gene_data[[intersectgene[j]]])
    correlation_matrix[i, j] <- result$estimate
    p_value_matrix[i, j] <- result$p.value
  }
}
correlation_df <- reshape2::melt(correlation_matrix, varnames = c("Var1", "Var2"), value.name = "Correlation")
p_value_df <- reshape2::melt(p_value_matrix, varnames = c("Var1", "Var2"), value.name = "PValue")

combined_df <- merge(correlation_df, p_value_df, by = c("Var1", "Var2"))

color_scheme <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

heatmap_plot <- ggplot(combined_df, aes(x = Var2, y = Var1, fill = Correlation)) +
  geom_tile() +xlab("")+ylab("")+
  scale_fill_gradientn(colors = color_scheme, limits = c(-0.3, 0.3)) +
  geom_text(aes(label = ifelse(PValue < 0.05, "*", "")), color = "black", size = 3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(angle = 45, hjust = 1))

setwd("S:/twas/twascrc/figure/正文")
pdf("FigureS2 heatmap_SSC_correlation.pdf",width = 8,height = 6)#pdf("Figure4 heatmap_ASC_correlation.pdf",width = 8,height = 6)
print(heatmap_plot)
dev.off()


Idents(seurat_object)<-seurat_object$disease_status
pdf("FigureS2 Dotplot.pdf",width = 12,height = 8)
intersectgene <- c("HES6","MIR4435-2HG","LIMCH1","POU2AF2","SMAD9","PITPNC1","GDF15","TGFBR2")
DotPlot(seurat_object, features = intersectgene)
dev.off()

##################
smrres<-fread("S:/twas/SMR/colon_eqtl.smr")
smrres<-subset(smrres,smrres$ProbeChr==8)
smrres<-smrres%>%mutate(is_highlight=ifelse(Gene%in%c("POU5F1B","C11orf53","COLCA1","COLCA2"),"yes","no"))%>%
  mutate(is_annotate=ifelse(Gene%in%c("POU5F1B","C11orf53","COLCA1","COLCA2"),"yes","no"))
setwd("S:/twas/twascrc/figure/SI")
pdf("SMR_CHR8.pdf",width = 8,height =6)
p1 <- ggplot(smrres, aes(x=Probe_bp, y=-log10(p_SMR))) +
  geom_point(aes(color=p_SMR), alpha=0.8, size=2) +
  scale_color_gradient(low="#010019", high="#89D3ED") +
  geom_point(data=subset(smrres, is_highlight=="yes"), aes(x=Probe_bp, y=-log10(p_SMR)), color="red", size=2) +
  geom_label_repel(data=subset(smrres, is_annotate=="yes"), aes(label=Gene), size=4) + 
  theme_bw() + 
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  geom_hline(yintercept = -log10(0.05/nrow(smrres)), color = 'red', size = 1.2, linetype = "dotted") +
  xlab("Chromosome 8") 

p1
dev.off()



#Figure S3
##Please use http://vip.sangerbox.com/home.html

#FigureS4
phewasplot<-read.csv("S:/twas/phewas/colon.csv")
unique(phewasplot$id.exposure)
phewasplot<-subset(phewasplot,phewasplot$pval<0.05&phewasplot$id.exposure=="PITPNC1")
# 

p<-ggplot(phewasplot, aes(x = reorder(category, minus_log10_pval), y = minus_log10_pval, color = category)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = 'Category', y = '-log10(p-value)', title = 'Phewas of PITPNC1') +
  scale_color_hue(l = 50) 
setwd("S:/twas/twascrc/figure/SI")
pdf("Figure7 Phewas_PITPNC1.pdf",width = 12,height =8)
p
dev.off()


###Table S1

setwd("S:/twas/twascrc/")
result_heng<-get_otters("S:/twas/CRC横结肠stage2")
p_heng<-0.05/nrow(result_heng)
result_heng1<-subset(result_heng,result_heng$OTTERS_p<p_heng)
SItable<-result_heng1%>%dplyr::select(Gene:OTTERS_p)
setwd("S:/twas/twascrc/figure/SI/")
write.csv(SItable,"Table S1.csv")


##########Table S2
phewastable<-read.csv("S:/twas/phewas/colon.csv")
phewastable<-subset(phewastable,phewastable$pval<0.05&phewastable$category%in%"neoplasms")
phewastable<-dplyr::arrange(phewastable,phenostring)
phewastable<-phewastable%>%dplyr::select(Gene=id.exposure,Method=method,nsnp:pval,phenocode:category)
setwd("S:/twas/twascrc/figure/SI/")
write.csv(phewastable,"Table S2.csv",row.names = F)
