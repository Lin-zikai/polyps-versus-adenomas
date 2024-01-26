source("prepare.R")

setwd("S:/twas/RNA-seq/")
####################################################################

data<-read_xlsx("GSE76987_RightColonProcessed.xlsx")
data<-data%>%dplyr::select(Ensembl_ID,`FPKM SSA/P-1`:`FPKM SSA/P-21`,`FPKM AP-1`:`FPKM AP-10`)
datagene<-bitr(data$Ensembl_ID,"ENSEMBL","SYMBOL",OrgDb = "org.Hs.eg.db")
data<-data%>%inner_join(datagene,by = c("Ensembl_ID"="ENSEMBL"))%>%dplyr::select(SYMBOL,
`FPKM SSA/P-1`:`FPKM SSA/P-21`,`FPKM AP-1`:`FPKM AP-10`,-Ensembl_ID)%>%distinct(SYMBOL,.keep_all = T)
SSA_exp<-as.data.frame(data%>%dplyr::select(`FPKM SSA/P-1`:`FPKM SSA/P-21`))
rownames(SSA_exp)<-data$SYMBOL

AP_exp<-as.data.frame(data%>%dplyr::select(`FPKM AP-1`:`FPKM AP-10`))
rownames(AP_exp)<-data$SYMBOL
group<-rbind(data.frame(sample=colnames(SSA_exp),group="SSA"),
             data.frame(sample=colnames(AP_exp),group="AP"))
data3<-as.data.frame(data[,-1])
rownames(data3)<-as.data.frame(data)[,1]
data3<-log2(data3+1)

Group=factor(group$group,levels = c("SSA","AP"))
design=model.matrix(~0+Group)
colnames(design)<-c("SSA","AP")
design
fit<-lmFit(data3,design)
contrast.matrix<-makeContrasts(SSA-AP,levels = design)

fit2<-contrasts.fit(fit,contrast.matrix)
fit2<-eBayes(fit2)
diff=topTable(fit2,adjust.method = "fdr" ,p.value=0.05,number=30000,sort.by = "logFC")

#######################################################Boxplot

genes <- c("HES6", "MIR4435-2HG", "LIMCH1", "POU2AF2", "SMAD9", "PITPNC1")
SS_transposed <- as.data.frame(t(SSA_exp[genes, ]))
SS_transposed<-log(SS_transposed+1)
AP_transposed <- as.data.frame(t(AP_exp[genes, ]))
AP_transposed<-log(AP_transposed+1)
SS_transposed$Group<-"SS"
AP_transposed$Group<-"AP"
data_to_plot <- rbind(SS_transposed, AP_transposed)
data_to_plot<-pivot_longer(data = data_to_plot,cols = -Group,names_to = "Gene",values_to = "Expression")
p <- ggplot(data_to_plot, aes(x=Gene, y=Expression, fill=Group)) +
  geom_boxplot() +scale_fill_manual(values=c("#316291", "#EC0002"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  stat_compare_means(aes(group = Group), label = "p.signif",method = "wilcox")+theme_bw()
setwd("S:/twas/twascrc/figure/正文")
pdf("Figure7 Boxplot.pdf",width = 8,height = 6)
p
dev.off()
###############################################correlation analysis
rownames(SSA_exp)<-data$SYMBOL
SSA_exps<-SSA_exp%>%mutate(SYMBOL=rownames(SSA_exp))
rownames(AP_exp)<-data$SYMBOL
genes_of_interest <- c("TGFBR2","GDF15","HES6", "MIR4435-2HG", "LIMCH1", "POU2AF2", "SMAD9", "PITPNC1", "GRN", "SORT1")
SSA_genes <- SSA_exp %>% filter(rownames(SSA_exp) %in% genes_of_interest) %>% t()
AP_genes <- AP_exp %>% filter(rownames(AP_exp) %in% genes_of_interest) %>% t()
cor_SSA <- cor(SSA_genes, use = "complete.obs")
cor_AP <- cor(AP_genes, use = "complete.obs")

p1 <- pheatmap(cor_SSA, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
               cluster_rows = F, cluster_cols = F, main = "SSA Group", silent = TRUE, return = TRUE)

p2 <- pheatmap(cor_AP, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
               cluster_rows = F, cluster_cols = F, main = "AP Group", silent = TRUE, return = TRUE)

setwd("S:/twas/twascrc/figure/正文")
pdf("Figure7 heatmap.pdf",width = 8,height = 6)
grid.arrange(p1$gtable, p2$gtable, ncol = 2)
dev.off()

#######################survival analysis on MIR4435-2HG AND SMAD9
#Please use http://vip.sangerbox.com/home.html



############################################################Phewas scatter plot
phewasplot<-read.csv("S:/twas/phewas/colon.csv")
unique(phewasplot$id.exposure)
phewasplot<-subset(phewasplot,phewasplot$pval<0.05&phewasplot$id.exposure=="SMAD9")
# manhattan

p<-ggplot(phewasplot, aes(x = reorder(category, minus_log10_pval), y = minus_log10_pval, color = category)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = 'Category', y = '-log10(p-value)', title = 'Phewas of SMAD9') +
  scale_color_hue(l = 50) # 使用不同颜色代表不同类别
setwd("S:/twas/twascrc/figure/正文")
pdf("Figure7 Phewas_SMAD9.pdf",width = 12,height =8)
p
dev.off()

