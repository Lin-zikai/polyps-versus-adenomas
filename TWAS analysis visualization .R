source("prepare.R")
setwd("S:/twas/twascrc/")

############################################
# result_heng<-get_otters("S:/twas/CRC横结肠stage2")
# result_blood<-get_otters("S:/twas/CRC全血stage2")
result_heng<-readRDS("result_heng.rds")
p_heng<-0.05/nrow(result_heng)
result_heng1<-subset(result_heng,result_heng$OTTERS_p<p_heng)
# SItable<-result_heng1%>%dplyr::select(Gene:OTTERS_p)
# setwd("S:/twas/twascrc/figure/SI/")
# write.csv(SItable,"resulheng_table.csv")
twasgene<-bitr(result_heng1$Gene,"ENSEMBL","SYMBOL",OrgDb = "org.Hs.eg.db")


#####################################################   manhattan #########################

Gene<-bitr(result_heng$Gene,"ENSEMBL","SYMBOL",OrgDb = "org.Hs.eg.db")
manhattan_data<-merge(result_heng,Gene,by.x="Gene",by.y="ENSEMBL")
manhattan_data$BP<-as.numeric(manhattan_data$GeneEnd)-as.numeric(manhattan_data$GeneStart)+1
manhattan_data<-manhattan_data%>%dplyr::select(SNP=SYMBOL,CHR=CHROM,BP,P=OTTERS_p)%>%distinct(SNP,.keep_all = T)
data<-manhattan_data
max_bp_per_chr <- data %>% 
  group_by(CHR) %>% 
  summarise(max_bp = max(BP))
cumulative_max_bp <- cumsum(max_bp_per_chr$max_bp)

data$all <- 0
for (chr in 1:22) {
  chr_indices <- which(data$CHR == chr)
  if (chr > 1) {
    prev_sum <- cumulative_max_bp[chr - 1]
  } else {
    prev_sum <- 0
  }
  
  data$all[c(chr_indices)] <- data$BP[c(chr_indices)] + prev_sum
}
data$BPcum <- data$all
data<-data%>%mutate(is_highlight=ifelse(SNP%in%twasgene$SYMBOL,"yes","no"))%>%
  mutate(is_annotate=ifelse(SNP%in%twasgene$SYMBOL,"yes","no"))
data$CHR <- factor(data$CHR, levels = as.character(1:22))
my_colors <- c(brewer.pal(n = 9, name = "Set3"),brewer.pal(n=8,"Set2"),brewer.pal(n=5,"Set2"))
X_axis <-  data %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
X_axis<-X_axis%>%arrange(CHR)
X_axis$center[7]<-9616844
X_axis$center[9]<-11528417
X_axis$center[11]<-13823261
X_axis$center[12]<-14823261
X_axis$CHR <- factor(X_axis$CHR, levels = as.character(1:22))
setwd("S:/twas/twascrc/figure/正文")
pdf("Figure3 Manhattan.pdf",width = 16,height = 12)
p1 <- ggplot(data, aes(x=BPcum, y=-log10(P))) +
  geom_point( aes(color=data$CHR), alpha=0.8, size=1.3) +
  scale_color_manual(values = my_colors) +
  scale_x_continuous( label = X_axis$CHR, breaks= X_axis$center ) +
  scale_y_continuous(expand = c(0, 0)) +ylim(0,25)+  
  geom_point(data=subset(data, is_highlight=="yes"), color="#1D91C0", size=2) +
  geom_label_repel( data=subset(data, is_annotate=="yes"), aes(label=SNP), size=2) +
  theme_bw() + 
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  )+geom_hline(yintercept = c(-log10(p_heng)), color = c('red'), size = 1.2, linetype = c("dotted"))+xlab("Chromosome") 
p1
dev.off()