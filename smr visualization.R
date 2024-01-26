#Please read SMR method from Yang Lab(https://yanglab.westlake.edu.cn/software/smr/#Overview)


####################################
smrres<-fread("S:/twas/SMR/colon_eqtl.smr")
smrres<-subset(smrres,smrres$ProbeChr==11)
# smrres<-subset(smrres,smrres$ProbeChr==8)
smrres<-smrres%>%mutate(is_highlight=ifelse(Gene%in%c("POU5F1B","C11orf53","COLCA1","COLCA2"),"yes","no"))%>%
  mutate(is_annotate=ifelse(Gene%in%c("POU5F1B","C11orf53","COLCA1","COLCA2"),"yes","no"))
smrres$Gene[smrres$Gene=="C11orf53"]<-"POU2AF2"
setwd("S:/twas/twascrc/figure/正文")
pdf("Figure4 SMR_CHR11.pdf",width = 8,height =6)
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
  xlab("Chromosome 11") 

p1
dev.off()


