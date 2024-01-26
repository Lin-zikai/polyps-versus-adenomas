source("prepare.R")
setwd("S:/twas/cellchat/")

##############################################################################
seurat_object <- readRDS("S:/twas/twascrc/merge_seurat.rds")
seurat_object$disease_status <- ifelse(seurat_object$disease %in% "colon sessile serrated adenoma/polyp", "SS", 
                                       ifelse(seurat_object$disease %in% c("tubulovillous adenoma", "tubular adenoma"), "AP", "Other"))

seurat_object<- subset(seurat_object, subset = disease_status == "SS")#seurat_object <- subset(seurat_object, subset = disease_status == "AP")


class(seurat_object)

data.input <- seurat_object@assays$RNA@counts # normalized data matrix

meta <- seurat_object@meta.data# a dataframe with rownames containing cell mata data
unique(meta$Cell_Type)
cell.use <- rownames(meta)[meta$Cell_Type %in% c("ABS","SSC","GOB","TAC","ASC","TUF","STM","EE","CT")]
data.input <- data.input[, cell.use]
meta = meta[cell.use, ]
unique(meta$Cell_Type)
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "Cell_Type")

cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "Cell_Type") # set "labels" as default cell identity

levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
groupSize


CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)


dplyr::glimpse(CellChatDB$interaction)

CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")

cellchat@DB <- CellChatDB.use


cellchat <- subsetData(cellchat,features = NULL)

cellchat <- identifyOverExpressedGenes(cellchat)

cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = T)

saveRDS(cellchat,"cellchat_SS_computeCommunProb.rds")
saveRDS(cellchat,"cellchat_AP_computeCommunProb.rds")

###################################################################################
cellchat.SS <- readRDS("cellchat2_SS_computeCommunProb.rds")
cellchat.SS  <- filterCommunication(cellchat.SS , min.cells = 10)
cellchat.SS <- computeCommunProbPathway(cellchat.SS)
cellchat.SS <- aggregateNet(cellchat.SS)
cellchat.AP <- readRDS('cellchat2_AP_computeCommunProb.rds')
cellchat.AP  <- filterCommunication(cellchat.AP , min.cells = 10)
cellchat.AP <- computeCommunProbPathway(cellchat.AP)
cellchat.AP <- aggregateNet(cellchat.AP)


netVisual_circle(cellchat.SS@net$count, vertex.weight = groupSize, weight.scale = T,
                 label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat.AP@net$count, vertex.weight = groupSize, weight.scale = T,
                 label.edge= F, title.name = "Number of interactions")



df.net.AP <- subsetCommunication(cellchat.AP)
pathways.show.AP <- df.net.AP$pathway_name
netVisual_aggregate(cellchat.AP, signaling = pathways.show.AP[pathways.show.AP%in%"EGF"][1], layout = "circle")

df.net.SS <- subsetCommunication(cellchat.SS)
pathways.show.SS <- df.net.SS$pathway_name
netVisual_aggregate(cellchat.SS, signaling = pathways.show.SS[pathways.show.SS%in%"EGF"][1], layout = "circle")

setwd("S:/twas/twascrc/figure/正文")
pdf("Figure2 cellchatBubble_ASC.pdf",width = 8,height = 6)
netVisual_bubble(cellchat.AP, sources.use = 2, targets.use = c(1,3:8),angle.x = 45)
dev.off()

pdf("Figure2 cellchatBubble_SSC.pdf",width = 8,height = 6)
netVisual_bubble(cellchat.SS, sources.use = 5, targets.use = c(1:4,6:8),angle.x = 45)
dev.off()

###############################
p1 <- netAnalysis_contribution(cellchat.SS, signaling = pathways.show.SS[pathways.show.SS%in%"EGF"][1],
                               title =  pathways.show[pathways.show.SS%in%"EGF"][1])
p2 <- netAnalysis_contribution(cellchat.AP, signaling = pathways.show.AP[pathways.show.AP%in%"EGF"][1],
                               title =  pathways.show.AP[pathways.show.AP%in%"EGF"][1])
cowplot::plot_grid(p1, p2, align = "h",ncol=2)


#
pairLR.CXCL <- extractEnrichedLR(cellchat.SS, signaling = pathways.show.SS[pathways.show.SS%in%"EGF"][1],
                                 geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
vertex.receiver = seq(1,4) # a numeric vector
p1<-netVisual_individual(cellchat.SS, signaling = pathways.show.SS,  
                         pairLR.use = LR.show, vertex.receiver = vertex.receiver,
                         layout = 'hierarchy')

pairLR.CXCL <- extractEnrichedLR(cellchat.AP, signaling = pathways.show.AP[pathways.show.AP%in%"EGF"][1],
                                 geneLR.return = FALSE)#
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
vertex.receiver = seq(1,4) # a numeric vector
p2<-netVisual_individual(cellchat.AP, signaling = pathways.show.AP[pathways.show.AP%in%"EGF"][1],  
                         pairLR.use = LR.show, vertex.receiver = vertex.receiver,
                         layout = 'hierarchy')

cowplot::plot_grid(p1, p2 ,align = "h",ncol=2)


netVisual_individual(cellchat.AP, signaling =pathways.show.AP[pathways.show.AP%in%"EGF"][1], 
                     pairLR.use = LR.show, layout = "circle")
netVisual_individual(cellchat.SS, signaling =pathways.show.SS[pathways.show.SS%in%"EGF"][1], 
                     pairLR.use = LR.show, layout = "circle")


##########################
netVisual_bubble(cellchat.SS, sources.use = 4, 
                       targets.use = c(5:11), remove.isolate = FALSE)
netVisual_bubble(cellchat.AP, sources.use = 4, 
                       targets.use = c(5:11), remove.isolate = FALSE)


##########################################################

object.list <- list(SS = cellchat.SS, AP = cellchat.AP)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))


setwd("S:/twas/twascrc/figure/正文")
pdf("Figure2 Information_flow.pdf",width = 8,height = 6)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg2
dev.off()
























