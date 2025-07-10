title: "Morpheus"
# Preparation ================
# clear existing variables
rm(list=ls())
# ===========================
# set working directory
setwd("/DATA07/home/shlee/Mopheus/Final/")
getwd()
list.files()
# Libraries =================
library(Seurat)
library(dplyr)
library(reticulate)
library(sctransform)
library(cowplot)
library(ggplot2)
library(viridis)
library(tidyr)
library(magrittr)
library(reshape2)
library(progeny)
library(readr)
library(stringr)
library(readxl)
library(harmony)
library(EnhancedVolcano)
library(genekitr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(Seurat)
library(SeuratWrappers)
library(tradeSeq)

# Loding Data =================

# PreTotal
load('/DATA07/home/shlee/Mopheus/Final/PreTotal01.robj')
# PreTNK 
load('/DATA07/home/shlee/Mopheus/Final/PreTotal01.robj')
#PreCD8T
load('Pre_T(cm_em_rm_emra_ex)_Final01.robj')
# PreMyeloid
load('/DATA07/home/shlee/Mopheus/Final/Pre_Myeloid_Final01.robj')
# PreTcellAddGiniIndexmexta
load('/DATA07/home/shlee/Mopheus/Final/PreTcell_Giniindexaddmeta.robj')
# PreTcellSlingshot
load("/DATA07/home/shlee/Mopheus/Final/slingshot_Final_0101.robj")
#CellChat
load('/DATA07/home/shlee/Mopheus/Final/Pre_CellChatFinal_CB.robj')
load('/DATA07/home/shlee/Mopheus/Final/Pre_CellChatFinal_NCB.robj')

#################################### Marker ################################################

TNKmarkers <- c('CD3D',  'CD4', 'CD8A', 'CD8B',
                'LEF1', 'SELL',  'TCF7',  
                'GNLY', 'IFNG', 'NKG7', 'PRF1', 'GZMA', 'GZMB', 'GZMH', 'GZMK',
                'MKI67', 'TOP2A',
                'HAVCR2', 'LAG3', 'PDCD1', 'CTLA4', 'TIGIT', 'BTLA',
                'KLRC1', 'ANXA1', 'ANKRD28', 'IL7R', 'CD69', 'CD40LG', 'ZNF683', 'ITGAE',
                'FOXP3', 'IL2RA', 'IKZF2',
                'NCR1', 'NCAM1', 'TYROBP', 'KLRD1', 'KLRF1', 'KLRB1', 'CX3CR1', 'FCGR3A', 'XCL1', 'XCL2')
TNKmarkers<- rev(TNKmarkers)

TNKmarkers <- c('CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B', 'PTPRC',
                'LEF1', 'SELL', 'CCR7', 'TCF7',  
                'GNLY', 'IFNG', 'NKG7', 'PRF1', 'GZMA', 'GZMB', 'GZMH', 'GZMK',
                'MKI67', 'TOP2A',
                'HAVCR2', 'LAG3', 'PDCD1', 'CTLA4', 'TIGIT', 'BTLA',
                'KLRC1', 'ANXA1', 'ANKRD28', 'IL7R', 'CD69', 'CD40LG', 'ZNF683', 'ITGAE',
                'FOXP3', 'IL2RA', 'IKZF2',
                'NCR1', 'NCAM1', 'TYROBP', 'KLRD1', 'KLRF1', 'KLRB1', 'CX3CR1', 'FCGR3A', 'FGFBP2', 'SPON2', 'XCL1', 'XCL2')

TNKmarkers1_CD8 <- c('CCR7', 'SELL', 'TCF7', 'LEF1', 'NELL2', 'IL7R', 'ANXA1', 'FOS', 'JUN', 'ZNF683', 'TOB1',
                     'GZMK', 'GPR183', 'EOMES', 'ITM2C', 'GZMA', 'GZMB', 'GZMH', 'CX3CR1', 'FGFBP2', 'SPON2',
                     'PRF1', 'NKG7', 'IFNG', 'CCL4', 'CCL5' ,'CCL4L2', 'TIGIT', 'CTLA4', 'PDCD1', 'HAVCR2')

TNKmarkers1_CD4 <- c('SELL', 'TCF7', 'CCR7', 'LEF1', 'FHIT', 'NOSIP', 'IL7R', 'CD40LG', 'ANXA1','FOS', 'FOSB',
                     'JUN', 'GPR183', 'ZNF683', 'GZMK', 'GZMA', 'GZMH', 'CXCR3', 'XCL1', 'XCL2', 'STAT1', 
                     'STAT4', 'IFNG', 'PRF1', 'NKG7', 'CXCL13', 'CD200', 'IL21', 'CXCR5', 'BCL6', 'NMB', 
                     'ICA1', 'GNG4', 'EBI3', 'ITM2A', 'IL2RA', 'FOXP3', 'TNFRSF4', 'TNFRSF9',
                     'TNFRSF18', 'IKZF2', 'IL1R2', 'IL1R1', 'LAIR2', 'HAVCR2', 'PDCD1','CTLA4', 'TIGIT', 'IL17A', 'IL23R', 'IL22')

T_CD8 <- c('CX3CR1', 'ZNF683','IFIT1', 'IFIT2', 'IFIT3', 'TNFSF10', 'LTB', 'LTA', 'TNF', 'FASLG', 'GZMM', 'GZMK', 'GZMH', 'GZMB', 'GZMA', 'PRF1',
           'NKG7', 'GNLY', 'KLRK1', 'KLRD1', 'KLRB1', 'IFNG')

Macrophage_gene <- c('PDGFA', 'VCAN', 'MMP9', 'SPP1', 'VEGFB', 'VEGFA', 'HLA-DQB1', 'HLA-DQA1', 'HLA-DPB1', 'HLA-DQA1','HLA-DRA',
                     'HLA-C', 'HLA-B', 'HLA-A', 'CD86', 'CD80', 'GPNMB', 'CCL18', 'TREM2', 'PLD4', 'SELENOP', 'F13A1', 'STAB1',
                     'MSR1', 'CD163L1', 'AXL', 'FOLR2', 'C1QB', 'MERTK', 'MRC1', 'CD163', 'GBP1', 'IDO1', 'STAT1', 'CXCL11',
                     'CXCL10', 'CXCL9', 'S100A9', 'S100A8', 'SOCS3', 'NDRG1', 'CXCL2', 'TNF', 'CCL4', 'CCL3', 'CCL2', 'IL1B')


################# TOP Gene ######################################
TNKmaker <- c('FOXP3', 'IL2RA', 'IL1R2',
              'CXCL13', 'PTPN13', 'GNG4',
              'LEF1', 'CCR7', 'SELL',
              'LINC00892', 'TMEM173', 'TIMP1',
              'IL4I1', 'LMNA', 'GRASP',
              'EGR1', 'CATSPERB', 'FOSB',
              'TACSTD2', 'IL19', 'DGCR9',
              'DKK3', 'ENC1', 'ITM2C',
              'GZMK','GZMH', 'CCL5',
              'ZNF683', 'CAPG', 'LINC02446', 
              'CX3CR1', 'FGFBP2', 'PRSS23',
              'KIR2DL4', 'HAVCR2', 'LAG3',
              'IFIT1', 'IFIT2', 'IFIT3',
              'MKI67', 'TOP2A', 'SPC25',
              'SLC4A10', 'TRAV1-2', 'NCR3',
              'SMAP2', 'SQSTM1', 'MBD2',
              'DLL1', 'KRT81', 'XCL1',
              'SPRY2', 'TRDC', 'KLRC2',
              'PTGDS', 'CXCR2', 'MYOM2')

Myeloidmarker <- c('FCN1', 'VCAN', 'S100A9', 'S100A12',
                   'FCGR3A', 'CDKN1C', 'MTSS1', 'RHOC',
                   'OLR1', 'EREG', 'G0S2', 'IL1B',
                   'CXCL10', 'CXCL9', 'CXCL11', 'GBP5',
                   'FOLR2', 'SLC40A1', 'RNASE1', 'MAF',
                   'SPP1', 'MMP9', 'MMP19', 'SULT1C2',
                   'TREM2', 'EMP1', 'LGALS3', 'CSTB')




#################################### Figure1 ###############################################

# Figure1A
DimPlot(obj_TNK_S, group.by = 'TNK_celltype1', label = T, repel = T)

# Figure1B
dat <- table(obj_TNK_S$sample_id, obj_TNK_S$TNK_celltype1)
dat <- prop.table(table(obj_TNK_S$sample_id, obj_TNK_S$TNK_celltype1),1)
dat <- as.data.frame(dat)
colnames(dat) <- c("Batch", "T_Celltype", "Freq")

dat$ClinicalBenefit <- "none"
dat$ClinicalBenefit[dat$Batch == 'MORSCRHBOT00101'] <- 'Non-Responder'
dat$ClinicalBenefit[dat$Batch == 'MORSCRHBOT00102'] <- 'Non-Responder'
dat$ClinicalBenefit[dat$Batch == 'MORSCRHBOT00202'] <- 'Responder'
dat$ClinicalBenefit[dat$Batch == 'MORSCRHBOT00402'] <- 'Responder'
dat$ClinicalBenefit[dat$Batch == 'MORSCRHBOT00702'] <- 'Responder'
dat$ClinicalBenefit[dat$Batch == 'MORSCRHBOT00901'] <- 'Responder' 

dat$ClinicalBenefit <- factor(dat$ClinicalBenefit, levels = c("Responder", "Non-Responder"))

ggplot(dat, aes(x = T_Celltype, y = Freq, fill = ClinicalBenefit)) + 
  geom_boxplot(width = 0.5, color = "black", outlier.shape = NA, position = position_dodge(width = 0.75)) + 
  ylab("Relative frequency") + 
  xlab("") + 
  geom_signif(
    test = "t.test", 
    comparisons = list(c("Responder", "NonResponder")), 
    map_signif_level = TRUE, 
    step_increase = 0.1, 
    y_position = c(0,0.4), 
    color = "black",
    aes(group = T_Celltype) 
  ) +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black", size = 9, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(color = "black", size = 11, face = "bold"),
    legend.title = element_text(color = "black", size = 12, face = "bold"),
    legend.text = element_text(color = "black", size = 10),
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    title = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "top"
  ) +
  labs(fill = "Response")

# Figure1C
Idents(obj_CD8T) <- 'ClinicalBenefit'
markers.cluster1 <- FindMarkers(obj_CD8T, ident.1 = "CB", ident.2 = "NCB", min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
data <- markers.cluster1
EnhancedVolcano(data,
                lab = rownames(data),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pointSize = 0.3,
                labSize = 2.5,
                drawConnectors = TRUE,
                widthConnectors = 0.1,
                pCutoff = 0.05,
                FCcutoff = 0.25,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = F,
                selectLab = c('PDCD1','CTLA4','LAG3', 'JUNB', 'GNLY', 'FOS', 'IL7R', 'JUN', 'TNFAIP3',
                              'NKG7', 'PRF1', 'GZMK', 'GZMA', 'GZMB', 'ZNF683', 'HLA-DQA1', 'CCR7','DUSP8',
                              'XCL1', 'CXCL13', 'SPP1', 'IFNG', 'CD7', 'TIGIT', 'VCAM1', 'PDCL3',
                              'KLRC1', 'ANXA1', "CXCR4", 'JUND', 'DUSP4', 'NFKBIA', 'NFKBIZ', 'RPS29', 'RPS10', 'CXCL13',
                              'HAVCR2', 'TNFALP3', 'NR4A2', 'DUSP2', 'SOCS3', 'SOCS1', 'DUSP5', 'AREG', 'TNFSF9',
                              'TXNIP', 'PTGDR', 'HSPA1A', 'LPP', 'HSPA1B', 'TXK', 'S1PR5', 'CD40LG',
                              'NCAM1', 'KLRG1', 'FOSB', 'IRF4', 'LINC00861', 'RPL38', 'GIMAP7', 'FGR', 'S1PR1'))   

# Figure1D
Idents(obj_pre) <- 'ClinicalBenefit'
markers.cluster1 <- FindMarkers(obj_TNK_S, ident.1 = "CB" , ident.2 = "NCB", min.pct = 0.25, logfc.threshold = 0.1, verbose = T)
deg <- markers.cluster1
input_genes = lapply(rownames(deg), function(x) strsplit(x, '[.]')[[1]][1]) %>% unlist()
deg['ENTREZID'] = mapIds(org.Hs.eg.db, input_genes, 'ENTREZID', 'SYMBOL')
deg = deg %>% arrange(desc(avg_log2FC))
input_gsea = deg$avg_log2FC %>% as.numeric()
names(input_gsea) = deg$ENTREZID %>% as.character()

gse1 <- geneset::getGO(org = "human",ont = "bp")
gsea1 <- genGSEA(input_gsea, gse1)
plotGSEA(gsea1, plot_type = "bar", colour = c("red", "blue"), label_by = 'description')

# Figure1E
T_obj_s_Gini$ind <- "FALSE"
T_obj_s_Gini$ind[T_obj_s_Gini$TNK_celltype1 %in% c('CD8 Tcm', 'CD8 Trm', 'CD8 Tem', 'CD8 Temra', 'CD8 Tex')] <- TRUE
T_obj_s_Gini_CD8 <- subset(T_obj_s_Gini, subset = ind == TRUE)
T_obj_s_Gini_CD8@meta.data <- droplevels(T_obj_s_Gini_CD8@meta.data)

metadata <- T_obj_s_Gini_CD8@meta.data

avg_gini_index <- metadata %>%
  group_by(TNK_celltype1, sample_id, ClinicalBenefit) %>%
  summarise(mean_gini = mean(gini_index), .groups = "drop")  # gini_index의 샘플별 평균 계산

ggplot(avg_gini_index, aes(x = TNK_celltype1, y = mean_gini, fill = TNK_celltype1)) +
  geom_boxplot() +
  labs(title = "Mean Gini Index Distribution by T Cell Type",
       x = "T Cell Type",
       y = "Mean Gini Index") +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black", size = 9, face = "bold", angle = 90, hjust = 1),
    axis.text.y = element_text(color = "black", size = 11, face = "bold"),
    legend.title = element_text(color = "black", size = 12, face = "bold"),
    legend.text = element_text(color = "black", size = 10),
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    title = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "top"
  ) 

#CB_NCB
ggplot(avg_gini_index, aes(x = TNK_celltype1, y = mean_gini, fill = ClinicalBenefit)) +
  geom_boxplot() +
  geom_signif(data = avg_gini_index, test = "t.test", comparisons = list(c("CB", "NCB")), 
              map_signif_level = F, step_increase = 0.01, y_position = 0.5, color = "black") +
  labs(title = "Mean Gini Index Distribution by TNK Cell Type and Clinical Benefit",
       x = "Clinical Benefit",
       y = "Mean Gini Index") +
  theme_minimal() +
  scale_fill_manual(values = c("CB" = "red", "NCB" = "#1E90FF")) +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black", size = 9, face = "bold", angle = 90, hjust = 1),
    axis.text.y = element_text(color = "black", size = 11, face = "bold"),
    legend.title = element_text(color = "black", size = 12, face = "bold"),
    legend.text = element_text(color = "black", size = 10),
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    title = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "top"
  ) 

#################################### Figure2 ###############################################
# Libraries =================
library(monocle3)

# Figure2A 
cds <- as.cell_data_set(obj_T)


#expression_matrix <- as(as.matrix(obj_T@assays$RNA@counts), 'sparseMatrix')
#cell_metadata <- obj_T@meta.data
#gene_annotation <- data.frame(genes = row.names(expression_matrix))
#rownames(gene_annotation) <- gene_annotation$genes
#cds <- new_cell_data_set(expression_matrix,
#                         cell_metadata = cell_metadata,
#                         gene_metadata = gene_annotation)

cds <- estimate_size_factors(cds)
#cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(GENESIS_MQ[["RNA"]])
cds <- preprocess_cds(cds, num_dim = 100)
cds <- align_cds(cds, alignment_group = "sample_id")
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds, use_partition = T)
#cds <- learn_graph(cds, use_partition = T,  learn_graph_control = list(ncenter=1000))
cds <- order_cells(cds)

colData(cds)
# to gene metdata
fData(cds)
rownames(fData(cds))[1:10]
# since it misses the gene_short_name column, let's add it
fData(cds)$gene_short_name <- rownames(fData(cds))

# to get counts
counts(cds)
cds@int_colData@listData$reducedDims$UMAP <- obj_T@reductions$umap@cell.embeddings
Morpheus_trajectory <- plot_cells(cds,
                                  color_cells_by = 'pseudotime',
                                  cell_size=0.55,
                                  label_groups_by_cluster = FALSE,
                                  label_branch_points = F,
                                  label_roots = T,
                                  label_leaves = FALSE,
                                  group_label_size = 0,
                                  show_trajectory_graph = TRUE)

# Figure2B
reduced_data <- reducedDims(cds)$UMAP  
cell_meta_data <- colData(cds)

new_dats <- data.frame(
  X = reduced_data[, 1],  
  Y = reduced_data[, 2],  
  TNK_celltype1 = cell_meta_data$TNK_celltype1,  
  ClinicalBenefit = cell_meta_data$ClinicalBenefit
)
ggplot(new_dats, aes(x=X, y=Y)) + geom_point(color="grey", size=0.5) + geom_pointdensity(data=new_dats[new_dats$ClinicalBenefit =="CB",]) + scale_color_viridis()
ggplot(new_dats, aes(x=X, y=Y)) + geom_point(color="grey", size=0.5) + geom_pointdensity(data=new_dats[new_dats$ClinicalBenefit =="NCB",]) + scale_color_viridis()

# Figure2C
library(slingshot)
sc <- obj_T
sc <- FindNeighbors(sc, reduction="harmony",dims = 1:40)
sc <- FindClusters(sc, resolution = 7)
DefaultAssay(sc) = 'RNA'
Idents(sc) = 'RNA_snn_res.7'
DimPlot(sc, group.by = 'RNA_snn_res.7', label = T, repel = T)
reduction = 'umap'
sds = slingshot(Embeddings(sc, reduction)[,1:2], clusterLabels = Idents(sc), start.clus = '27')
sc@tools[['slingshot']] = SlingshotDataSet(sds)
pseudotime = slingPseudotime(sds)

curves = colnames(pseudotime)
palette = viridis(100, end = 0.95)
head(rd)
# add reduceDim
rd <- Embeddings(sc, reduction)[,1:2] %>% as.matrix()
rd <- rd[colnames(sc),]

sc$slingshot_pseudotime_curve1 = pseudotime[,1]
sc$slingshot_pseudotime_curve2 = pseudotime[,2]
sc$slingshot_pseudotime_curve3 = pseudotime[,3]
sc$slingshot_pseudotime_curve4 = pseudotime[,4]
sc$slingshot_pseudotime_curve5 = pseudotime[,5]
sc$slingshot_pseudotime_curve6 = pseudotime[,6]
sc$slingshot_pseudotime_curve7 = pseudotime[,7]
sc$slingshot_pseudotime_curve8 = pseudotime[,8]
sc$slingshot_pseudotime_curve9 = pseudotime[,9]
sc$slingshot_pseudotime_curve10 = pseudotime[,10]
sc$slingshot_pseudotime_curve11 = pseudotime[,11]
sc$slingshot_pseudotime_curve12 = pseudotime[,12]
sc$slingshot_pseudotime_curve13 = pseudotime[,13]
sc$slingshot_pseudotime_curve14 = pseudotime[,14]
sc$slingshot_pseudotime_curve15 = pseudotime[,15]
sc$slingshot_pseudotime_curve16 = pseudotime[,16]
sc$slingshot_pseudotime_curve17 = pseudotime[,17]
sc$slingshot_pseudotime_curve18 = pseudotime[,18]
sc$slingshot_pseudotime_curve19 = pseudotime[,19]
sc$slingshot_pseudotime_curve20 = pseudotime[,20]
sc$slingshot_pseudotime_curve21 = pseudotime[,21]
sc$slingshot_pseudotime_curve22 = pseudotime[,22]
sc$slingshot_pseudotime_curve23 = pseudotime[,23]

df <- sc@meta.data
psts <- slingPseudotime(sds) %>%
  as.data.frame() %>%
  mutate(cells = rownames(.),
         conditions = df$ClinicalBenefit) %>%
  pivot_longer(starts_with("Lineage"), values_to = "pseudotime", names_to = "lineages")

#TCRdata load
load("/DATA02/home/jyhwang4/mormor/mormor.rdata")
head(meta3t)
colnames(meta3t) <- c("cells", "Clonality2")

psts <- left_join(psts, meta3t, by = c("cells"))
psts1 <- subset(psts, lineages == 'Lineage1')
psts2 <- subset(psts, lineages == 'Lineage2')
psts3 <- subset(psts, lineages == 'Lineage3')
psts4 <- subset(psts, lineages == 'Lineage4')
psts5 <- subset(psts, lineages == 'Lineage5')
psts6 <- subset(psts, lineages == 'Lineage6')
psts7 <- subset(psts, lineages == 'Lineage7')
psts8 <- subset(psts, lineages == 'Lineage8')
psts9 <- subset(psts, lineages == 'Lineage9')
psts10 <- subset(psts, lineages == 'Lineage10')
psts11 <- subset(psts, lineages == 'Lineage11')
psts12 <- subset(psts, lineages == 'Lineage12')
psts13 <- subset(psts, lineages == 'Lineage13')
psts14 <- subset(psts, lineages == 'Lineage14')
psts15 <- subset(psts, lineages == 'Lineage15')
psts16 <- subset(psts, lineages == 'Lineage16')
psts17 <- subset(psts, lineages == 'Lineage17')
psts18 <- subset(psts, lineages == 'Lineage18')
psts19 <- subset(psts, lineages == 'Lineage19')
psts20 <- subset(psts, lineages == 'Lineage20')
psts21 <- subset(psts, lineages == 'Lineage21')
psts22 <- subset(psts, lineages == 'Lineage22')
psts23 <- subset(psts, lineages == 'Lineage23')

ggplot(psts1, aes(x = pseudotime,  fill = conditions)) +
  geom_density(alpha = 0.5) +
  scale_fill_brewer(type = "qual") +
  theme_bw()+
  ylim(0,0.2) +
  theme(legend.position = "bottom",
        panel.grid.major = element_blank())

FeaturePlot(sc, 
            feature = 'slingshot_pseudotime_curve3', 
            cols = c('red', 'yellow')) +
  ggtitle('') +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank()) +
  labs(colour = 'Pseudotime')

###################### Clonality ###############################
psts1_NCB <- psts3[,c("cells", "pseudotime", "Clonality2", "conditions")]
psts1_NCB <- subset(psts1_NCB, conditions %in% "NCB")
psts1_CB <- psts3[,c("cells","pseudotime", "Clonality2", "conditions")]
psts1_CB <- subset(psts1_CB, conditions %in% "CB")

colnames(psts1_NCB) <- c("cells","x", "y", "group")
colnames(psts1_CB) <- c("cells","x", "y", "group")

rownames(psts1_NCB) <- psts1_NCB$cells
rownames(psts1_CB) <- psts1_CB$cells


# Data summary function (average by x interval)
summarize_data <- function(data, bin_width = 0.5) {
  data %>%
    mutate(x_bin = floor(x / bin_width) * bin_width) %>%
    group_by(x_bin) %>%
    summarise(
      x = mean(x, na.rm = TRUE),
      weighted_mean = sum(y) / n(),
      .groups = 'drop'
    ) %>%
    na.omit()
}

# Data Summary
binned_NCB <- summarize_data(psts1_NCB, bin_width = 0.5)
binned_CB <- summarize_data(psts1_CB, bin_width = 0.5)

# Line graph function definition
plot_tcr_line <- function(data1, data2) {
  ggplot() +
    geom_line(data = data1, aes(x = x, y = weighted_mean), color = "blue", size = 1) +
    geom_line(data = data2, aes(x = x, y = weighted_mean), color = "red", size = 1) +
    ylim(0, 0.5) +  #
    xlim(0, 13)+
    xlab("Pseudotime") + ylab("Weighted Clonality (Mean)") +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12)
    )
}

# Create a line graph
p_line <- plot_tcr_line(binned_CB, binned_NCB)

# Plot Output
print(p_line)

# Figure2D 
Idents(obj_TNK_S) <- 'TNK_celltype1'

markers.cluster1 <- FindMarkers(obj_TNK_S, ident.1 = "CD8 Trm" , ident.2 = "CD8 Temra", min.pct = 0.25, logfc.threshold = 0.1, verbose = T)
deg <- markers.cluster1

input_genes = lapply(rownames(deg), function(x) strsplit(x, '[.]')[[1]][1]) %>% unlist()
deg['ENTREZID'] = mapIds(org.Hs.eg.db, input_genes, 'ENTREZID', 'SYMBOL')
deg = deg %>% arrange(desc(avg_log2FC))
input_gsea = deg$avg_log2FC %>% as.numeric()
names(input_gsea) = deg$ENTREZID %>% as.character()

gse1 <- geneset::getGO(org = "human",ont = "bp")
gsea1 <- genGSEA(input_gsea, gse1)

#gsea1 <- genGSEA(input_gsea, gse1, p_cutoff = 0.3, q_cutoff = 1)

plotGSEA(gsea1, plot_type = "bar", colour = c("red", "blue"), label_by = 'description')

# Figure2E 
Idents(obj_TNK_S) <- 'TNK_celltype1'

markers.cluster1 <- FindMarkers(obj_TNK_S, ident.1 = "CD8 Tex" , ident.2 = "CD8 Temra", min.pct = 0.25, logfc.threshold = 0.1, verbose = T)
deg <- markers.cluster1

input_genes = lapply(rownames(deg), function(x) strsplit(x, '[.]')[[1]][1]) %>% unlist()
deg['ENTREZID'] = mapIds(org.Hs.eg.db, input_genes, 'ENTREZID', 'SYMBOL')
deg = deg %>% arrange(desc(avg_log2FC))
input_gsea = deg$avg_log2FC %>% as.numeric()
names(input_gsea) = deg$ENTREZID %>% as.character()

gse1 <- geneset::getGO(org = "human",ont = "bp")
gsea1 <- genGSEA(input_gsea, gse1)

gsea1 <- genGSEA(input_gsea, gse1, p_cutoff = 0.3, q_cutoff = 1)

plotGSEA(gsea1, plot_type = "bar", colour = c("red", "blue"), label_by = 'description')

#################################### Figure3 ###############################################

# Figure3A 
DimPlot(obj_Myeloid1s_pre1, group.by ='Myeloid_celltype_final1_s', label = T)

# Figure3B
markers.cluster1 <- FindMarkers(obj_Myeloid1s_pre1_Macro, ident.1 = "CB", ident.2 = "NCB", min.pct = 0.25, logfc.threshold = 0.1, verbose = FALSE)

data <- markers.cluster1
EnhancedVolcano(data,
                lab = rownames(data),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 3.0,
                title = 'Volcano plot with EnhancedVolcano',
                legendPosition = 'right',
                legendLabSize = 14,
                drawConnectors = TRUE,
                colAlpha = 0.70,
                #gridlines.major = FALSE,
                #gridlines.minor = FALSE,
                widthConnectors = 0.1,
                colConnectors = 'grey30'
)

# Figure3C
Idents(obj_Myeloid1s_pre1) <- 'ClinicalBenefit'
DotPlot(obj_Myeloid1s_pre1, features = 'CD274')

CXCL10_Macro <- subset(obj_Myeloid1s_pre1, subset = Myeloid_celltype_final1_s == "CXCL10+ Macrophage")
CXCL10_Macro@meta.data <- droplevels(CXCL10_Macro@meta.data)
DotPlot(CXCL10_Macro, features = 'CD14')

# Figure3D
library(CellChat)
cellchat <- createCellChat(object = obj_CB, group.by = "CellType1", assay = "RNA")

#cellchat<-addMeta(cellchat, meta = "MT2")
cellchat<-setIdent(cellchat,ident.use= "CellType1") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize<-as.numeric(table(cellchat@idents)) # number of cells in each cell group
CellChatDB<-CellChatDB.human# useCellChatDB.mouseif running on mouse data
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB$interaction$annotation
#CellChatDB.use <-subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
CellChatDB.use=CellChatDB
cellchat@DB<-CellChatDB.use
cellchat<-subsetData(cellchat) # This step is necessary even if using the whole database
#future::plan("multiprocess", workers = 4) # do parallel
cellchat<-identifyOverExpressedGenes(cellchat)
cellchat<-identifyOverExpressedInteractions(cellchat)
cellchat<-smoothData(cellchat, adj = PPI.human)
cellchat<-computeCommunProb(cellchat) 
cellchat<-filterCommunication(cellchat,min.cells= 10)
cellchat<-computeCommunProbPathway(cellchat)
cellchat<-aggregateNet(cellchat)
groupSize<-as.numeric(table(cellchat@idents))

cellchat <- createCellChat(object = obj_NCB, group.by = "CellType1", assay = "RNA")

#cellchat<-addMeta(cellchat, meta = "MT2")
cellchat<-setIdent(cellchat,ident.use= "CellType1") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize<-as.numeric(table(cellchat@idents)) # number of cells in each cell group
CellChatDB<-CellChatDB.human# useCellChatDB.mouseif running on mouse data
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB$interaction$annotation
#CellChatDB.use <-subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
CellChatDB.use=CellChatDB
cellchat@DB<-CellChatDB.use
cellchat<-subsetData(cellchat) # This step is necessary even if using the whole database
#future::plan("multiprocess", workers = 4) # do parallel
cellchat<-identifyOverExpressedGenes(cellchat)
cellchat<-identifyOverExpressedInteractions(cellchat)
cellchat<-smoothData(cellchat, adj = PPI.human)
cellchat<-computeCommunProb(cellchat) 
cellchat<-filterCommunication(cellchat,min.cells= 10)
cellchat<-computeCommunProbPathway(cellchat)
cellchat<-aggregateNet(cellchat)
groupSize<-as.numeric(table(cellchat@idents))

cellchat_CB_S <- cellchat
cellchat_NCB_S <- cellchat

targetdata <- unique(obj_CD8T_S1$TNK_celltype1)
sourcedata <- unique(obj_Myeloid1s_pre1_Macro$Myeloid_celltype_final1_s)   

pathways.show <- c("CXCL") 
netVisual_aggregate(cellchat_NCB_S, signaling = pathways.show,  vertex.receiver = vertex.receiver)
netVisual_aggregate(cellchat_NCB_S, signaling = pathways.show, layout = "circle",sources.use = sourcedata,  targets.use = targetdata)

# Figure3E
netAnalysis_contribution(cellchat_CB_S, signaling = pathways.show, thresh = T)
netAnalysis_contribution(cellchat_NCB_S, signaling = pathways.show, thresh = T)

# Figure3F
Idents(obj_CD8T_S1) <- 'TNK_celltype1'
DotPlot(obj_CD8T_S1, features = "CXCR3") 

# Figure3G
DotPlot(obj_CD8T_S1, features = "CXCR3", group.by = 'TNK_celltype1', split.by = 'ClinicalBenefit') 

################################ Supplementary Figure1 ##################################################

# Supple Figure1A
DimPlot(obj_pre, group.by = "CellType", cols = use_colors)

# Supple Figure1B
mainmarkers <- c(
  'KIT', 'MS4A2', 'GATA2',                    # Mast
  'CD68', 'FCGR3A', 'MARCO', 'LYZ',           # Myeloid
  'CD79A', 'IGHG3', 'IGHA2', 'JCHAIN',        # Bcell
  'NKG7', 'GNLY', 'KLRD1','NCAM1',            # NK cells
  'TRAC', 'CD3G', 'CD3E', 'CD3D',             # Tcell
  'THY1','COL1A2','COL1A1', 'DCN',            # Fibroblast
  'RAMP2', 'FLT1', 'CLDN5', 'PECAM1',         # Endothelial
  "GPRC5A", "ALB", 'KRT19', 'EPCAM') 
mainmarkers <- rev(mainmarkers)

DotPlot(obj_pre, group.by="CellType", features = mainmarkers) + theme_bw() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = pals::brewer.orrd(10)) + 
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, face = "bold")) +
  theme(axis.text.y = element_text(face = "bold")) +
  geom_hline(yintercept =seq(1.5, 24, 1), linetype="solid", size=0.5, color="black") 

# Supple Figure1C
Idents(obj_T_S) <- 'TNK_celltype1'
DotPlot(obj_T_S, features = 'PDCD1')

# Supple Figure1D
T_CD8 <- c('CX3CR1', 'ZNF683','IFIT1', 'IFIT2', 'IFIT3', 'TNFSF10', 'LTB', 'LTA', 'TNF', 'FASLG', 'GZMM', 'GZMK', 'GZMH', 'GZMB', 'GZMA', 'PRF1',
           'NKG7', 'GNLY', 'KLRK1', 'KLRD1', 'KLRB1', 'IFNG')

obj_TNK_S$ind <- "FALSE"
obj_TNK_S$ind[obj_TNK_S$TNK_celltype1 %in% c('CD8 Tcm', 'CD8 Trm', 'CD8 Tem', 'CD8 Temra', 'CD8 Tex', 'CD8 IFN activating T')] <- TRUE
obj_TNK_S_CD8 <- subset(obj_TNK_S, subset = ind == TRUE)

obj_TNK_S_CD8@meta.data <- droplevels(obj_TNK_S_CD8@meta.data)
unique(obj_TNK_S_CD8$TNK_celltype1)

obj_TNK_S_CD8$TNK_celltype1 <- factor(obj_TNK_S_CD8$TNK_celltype1, levels = c('CD8 Tcm', 'CD8 Trm', 'CD8 Tem', 'CD8 Temra', 'CD8 Tex', 'CD8 IFN activating T'))

Idents(obj_TNK_S_CD8) <- "TNK_celltype1"
expression_lv = as.data.frame(AverageExpression(obj_TNK_S_CD8), rownames = NULL)
expression_lv = data.frame(gene = rownames(expression_lv), 'CD8_Tcm' = expression_lv$RNA.CD8.Tcm,
                           'CD8_Trm' = expression_lv$RNA.CD8.Trm, 'CD8_Tem' = expression_lv$RNA.CD8.Tem,
                           'CD8_Temra' = expression_lv$RNA.CD8.Temra, 'CD8 Tex' = expression_lv$RNA.CD8.Tex,
                           'CD8 IFN activating T' = expression_lv$RNA.CD8.IFN.activating.T)

chemo_mac = subset(expression_lv, subset = expression_lv$gene %in% T_CD8)
head(chemo_mac)

row.names(chemo_mac) = chemo_mac$gene
chemo_mac[, 1] = NULL

z_chemo_mac = apply(chemo_mac, 1, scale)
new_mat = t(z_chemo_mac)

colnames(new_mat) <- c('CD8 Tcm', 'CD8 Trm', 'CD8 Tem', 'CD8 Temra', 'CD8 Tex', 'CD8 IFN activating T')
library(pheatmap)
T_CD8 <- rev(T_CD8)
new_mat_ordered <- new_mat[T_CD8, ]

pheatmap(new_mat_ordered, main = 'TNK DEG', 
         cluster_rows = F, cluster_cols = F, border_color = 'Black', angle_col = 90,
         gaps_row = c(1, 6, 12, 17, 20, 21, 22))

# Supple Figure1E
DimPlot(obj_CD8T, group.by = 'TNK_celltype1')
################################ Supplementary Figure2 ##################################################

# Supple Figure2A
Myeloidmarker <- c('FCN1', 'VCAN', 'S100A9', 'S100A12',
                   'FCGR3A', 'CDKN1C', 'MTSS1', 'RHOC',
                   'OLR1', 'EREG', 'G0S2', 'IL1B',
                   'CXCL10', 'CXCL9', 'CXCL11', 'GBP5',
                   'FOLR2', 'SLC40A1', 'RNASE1', 'MAF',
                   'SPP1', 'MMP9', 'MMP19', 'SULT1C2',
                   'TREM2', 'EMP1', 'LGALS3', 'CSTB')

Idents(obj_Myeloid1s_pre1) <- 'Myeloid_celltype_final1_s'
expression_lv = as.data.frame(AverageExpression(obj_Myeloid1s_pre1), rownames = NULL)

expression_lv = data.frame(gene = rownames(expression_lv), 'CD14' = expression_lv$RNA.CD14..Monocyte,
                           'CD16' = expression_lv$RNA.CD16..Monocyte, 'OLR1' = expression_lv$RNA.OLR1..Monocyte,
                           'CXCL10' = expression_lv$RNA.CXCL10..Macrophage,
                           'FOLR2' = expression_lv$RNA.FOLR2..Macrophage,
                           'SPP1' = expression_lv$RNA.SPP1..Macrophage,
                           'TREM2' = expression_lv$RNA.TREM2..Macrophage)

obj_Myeloid1s_pre1$Myeloid_celltype_final1_s <- factor(obj_Myeloid1s_pre1$Myeloid_celltype_final1_s, levels = c('CD14+ Monocyte', 'CD16+ Monocyte',
                                                                                                                'OLR1+ Monocyte',  
                                                                                                                'CXCL10+ Macrophage', 'FOLR2+ Macrophage',
                                                                                                                'SPP1+ Macrophage','TREM2+ Macrophage'))

chemo_mac = subset(expression_lv, subset = expression_lv$gene %in% Myeloidmarker)
head(chemo_mac)

row.names(chemo_mac) = chemo_mac$gene
chemo_mac[, 1] = NULL

z_chemo_mac = apply(chemo_mac, 1, scale)
new_mat = t(z_chemo_mac)

colnames(new_mat) <- c('CD14+ Monocyte', 'CD16+ Monocyte',
                       'OLR1+ Monocyte', 'CXCL10+ Macrophage',
                       'FOLR2+ Macrophage','SPP1+ Macrophage','TREM2+ Macrophage')

library(pheatmap)
Myeloidmarker <- rev(Myeloidmarker)
myeloid_markers_ordered <- new_mat[Myeloidmarker, ]

pheatmap(myeloid_markers_ordered, main = 'Myeloid DEG', 
         cluster_rows = F, cluster_cols = F, border_color = 'Black', angle_col = 90)

# Supple Figure2B

obj_Myeloid1s_pre1$Myeloid_celltype_final1_s

dat <- table(obj_Myeloid1s_pre1$sample_id, obj_Myeloid1s_pre1$Myeloid_celltype_final1_s)
dat <- prop.table(table(obj_Myeloid1s_pre1$sample_id, obj_Myeloid1s_pre1$Myeloid_celltype_final1_s),1)
dat <- as.data.frame(dat)
colnames(dat) <- c("Batch", "M_Celltype", "Freq")

dat$percent <- dat$Freq * 100

dat$ClinicalBenefit <- "none"

dat$ClinicalBenefit[dat$Batch == 'MORSCRHBOT00101'] <- 'Non-Responder'
dat$ClinicalBenefit[dat$Batch == 'MORSCRHBOT00102'] <- 'Non-Responder'
dat$ClinicalBenefit[dat$Batch == 'MORSCRHBOT00202'] <- 'Responder'
dat$ClinicalBenefit[dat$Batch == 'MORSCRHBOT00402'] <- 'Responder'
dat$ClinicalBenefit[dat$Batch == 'MORSCRHBOT00702'] <- 'Responder'
dat$ClinicalBenefit[dat$Batch == 'MORSCRHBOT00901'] <- 'Responder' 

dat$ClinicalBenefit <- factor(dat$ClinicalBenefit, levels = c("Responder", "Non-Responder"))

ggplot(dat, aes(x = M_Celltype, y = Freq, fill = ClinicalBenefit)) + 
  geom_boxplot(width = 0.5, color = "black", outlier.shape = NA, position = position_dodge(width = 0.75)) + 
  ylab("Relative frequency") + 
  xlab("") + 
  ylim(0, 0.7) + 
  geom_signif(
    comparisons = list(c("Responder", "NonResponder")), 
    test = "t.test", 
    map_signif_level = TRUE, 
    step_increase = 0.1, 
    y_position = 0.5  
  ) +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black", size = 9, face = "bold", angle = 90, hjust = 1),
    axis.text.y = element_text(color = "black", size = 11, face = "bold"),
    legend.title = element_text(color = "black", size = 12, face = "bold"),
    legend.text = element_text(color = "black", size = 10),
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    title = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "top"
  ) +
  labs(fill = "Response")  

# Supple Figure 2C

CD274_enrichedScore <- AddModuleScore(obj_Myeloid1s_pre1, features = CD274_enriched, name="CD274_enriched")
CD274_scores <- CD274_enrichedScore@meta.data$CD274_enriched1
sample_id <- obj_Myeloid1s_pre1@meta.data$sample_id 

data <- data.frame(
  sample_id = sample_id,
  CD274_scores1 = CD274_scores,
  Celltype = obj_Myeloid1s_pre1@meta.data$Myeloid_celltype_final1_s
)

sample_data <- data %>%
  group_by(sample_id,Celltype) %>%
  summarize(mean_CD274 = mean(CD274_scores1, na.rm = TRUE)) %>% ungroup()

library(ggpubr)
sample_data$Celltype  <- factor(sample_data$Celltype, levels = c("CD14+ Monocyte", "CD16+ Monocyte", "OLR1+ Monocyte",
                                                                 "CXCL10+ Macrophage", "FOLR2+ Macrophage", "SPP1+ Macrophage", 
                                                                 "TREM2+ Macrophage"))

ggplot(sample_data, aes(x=Celltype, y=mean_CD274, fill = Celltype)) +  geom_boxplot() + stat_summary(fun.y=mean, geom="point", shape=23, size=1.1, fill = "red") +xlab("")+ylab("Frequency")+
  theme_bw()+
  theme(panel.background = element_blank(),
        axis.text.x =element_text(angle = 90, hjust = 1, size = 9,face="bold"),
        axis.text.y = element_text(color="black",size = 11,face="bold"),
        legend.title = element_text(color="black",size=15,face="bold"),
        legend.text = element_text(color="black",size=7),
        strip.text.x = element_text(angle = 0,face="bold",size=6),
        axis.title.x = element_text(color="black",size = 12,face="bold"),
        title = element_text(color="black",size = 12,face="bold")) +
  ggtitle("CD274 Expression by Cell Type")

# Supple Figure2C
cellchat_CB_S <- cellchat
cellchat_NCB_S <- cellchat
object.list <- list(CB = cellchat_CB_S, NCB = cellchat_NCB_S)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2






















































































