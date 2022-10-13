library(dplyr, quietly = T)
library(Seurat, quietly = T)
library(patchwork, quietly = T)
library(dplyr, quietly = T)
library(ggplot2)
# Reading in data and creating merged seurat object
sc_Ts11_C <- Mo_datain("/sc_Ts11-C/outs/filtered_feature_bc_matrix/")
sc_Ts11_C$ID <- 'sc_Ts11-C'
sc_Ts11_M <- Mo_datain("/sc_Ts11-M/outs/filtered_feature_bc_matrix/")
sc_Ts11_M$ID <- 'sc_Ts11-M'
sc_Ts11_P <- Mo_datain("/sc_Ts11-P/outs/filtered_feature_bc_matrix/")
sc_Ts11_P$ID <- 'sc_Ts11-P'
sc_Ts15_C <- Mo_datain("/sc_Ts15-C/outs/filtered_feature_bc_matrix/")
sc_Ts15_C$ID <- 'sc_Ts15-C'
sc_Ts15_P <- Mo_datain("/sc_Ts15-P/outs/filtered_feature_bc_matrix/")
sc_Ts15_P$ID <- 'sc_Ts15-P'
sc_Ts6_4_M <- Mo_datain("/sc_Ts6-4-M/outs/filtered_feature_bc_matrix/")
sc_Ts6_4_M$ID <- 'sc_Ts6-4-M'
sc_Ts6_4_P <- Mo_datain("/sc_Ts6-4-P/outs/filtered_feature_bc_matrix/")
sc_Ts6_4_P$ID <- 'sc_Ts6-4-P'
sc_Ts6_8_M <- Mo_datain("/sc_Ts6_8-M/outs/filtered_feature_bc_matrix/")
sc_Ts6_8_M$ID <- 'sc_Ts6_8-M'
sc_Ts6_8_P <- Mo_datain("/sc_Ts6-8-P/outs/filtered_feature_bc_matrix/")
sc_Ts6_8_P$ID <- 'sc_Ts6-8-P'
sc_Ts6_C <- Mo_datain("/sc_Ts6-C/outs/filtered_feature_bc_matrix/")
sc_Ts6_C$ID <- 'sc_Ts6-C'
sc_Ts6_P <- Mo_datain("/sc_Ts6-P/outs/filtered_feature_bc_matrix/")
sc_Ts6_P$ID <- 'sc_Ts6-P'
sc_Ts8_15_M1 <- Mo_datain("/sc_Ts8_15-M1/outs/filtered_feature_bc_matrix/")
sc_Ts8_15_M1$ID <- 'sc_Ts8_15-M1'
sc_Ts8_15_M2 <- Mo_datain("/sc_Ts8_15-M2/outs/filtered_feature_bc_matrix/")
sc_Ts8_15_M2$ID <- 'sc_Ts8_15-M2'
sc_Ts8_15_P <- Mo_datain("/sc_Ts8-15-P/outs/filtered_feature_bc_matrix/")
sc_Ts8_15_P$ID <- 'sc_Ts8-15-P'
sc_Ts8_C <- Mo_datain("/sc_Ts8-C/outs/filtered_feature_bc_matrix/")
sc_Ts8_C$ID <- 'sc_Ts8-C'
sc_Ts8_M <- Mo_datain("/sc_Ts8-M/outs/filtered_feature_bc_matrix/")
sc_Ts8_M$ID <- 'sc_Ts8-M'
sc_Ts8_P <- Mo_datain("/sc_Ts8-P/outs/filtered_feature_bc_matrix/")
sc_Ts8_P$ID <- 'sc_Ts8-P'
sc_WT_C <- Mo_datain("/sc_WT-C/outs/filtered_feature_bc_matrix/")
sc_WT_C$ID <- 'sc_WT-C'
sc_WT_P <- Mo_datain("/sc_WT-P/outs/filtered_feature_bc_matrix/")
sc_WT_P$ID <- 'sc_WT-P'

sce_merge <- merge(sc_Ts11_C, c(sc_Ts11_M, 
sc_Ts11_P, 
sc_Ts15_C, 
sc_Ts15_P, 
sc_Ts6_4_M, 
sc_Ts6_4_P, 
sc_Ts6_8_M, 
sc_Ts6_8_P, 
sc_Ts6_C, 
sc_Ts6_P, 
sc_Ts8_15_M1, 
sc_Ts8_15_M2, 
sc_Ts8_15_P, 
sc_Ts8_C, 
sc_Ts8_M, 
sc_Ts8_P, 
sc_WT_C, 
sc_WT_P ))
VariableFeatures(sce_merge) <- c(VariableFeatures(sc_Ts11_C),VariableFeatures(sc_Ts11_M),VariableFeatures(sc_Ts11_P),
                                VariableFeatures(sc_Ts15_C),VariableFeatures(sc_Ts15_P),VariableFeatures(sc_Ts6_4_M),VariableFeatures(sc_Ts6_4_P),
                                VariableFeatures(sc_Ts6_8_M),VariableFeatures(sc_Ts6_8_P),VariableFeatures(sc_Ts6_C),VariableFeatures(sc_Ts6_P),
                                VariableFeatures(sc_Ts8_15_M1),VariableFeatures(sc_Ts8_15_M2),VariableFeatures(sc_Ts8_15_P),VariableFeatures(sc_Ts8_C),
                                VariableFeatures(sc_Ts8_M),VariableFeatures(sc_Ts8_P),VariableFeatures(sc_WT_C),VariableFeatures(sc_WT_P))

rm(sc_Ts11_C,
   sc_Ts11_M, 
sc_Ts11_P, 
sc_Ts15_C, 
sc_Ts15_P, 
sc_Ts6_4_M, 
sc_Ts6_4_P, 
sc_Ts6_8_M, 
sc_Ts6_8_P, 
sc_Ts6_C, 
sc_Ts6_P, 
sc_Ts8_15_M1, 
sc_Ts8_15_M2, 
sc_Ts8_15_P, 
sc_Ts8_C, 
sc_Ts8_M, 
sc_Ts8_P, 
sc_WT_C, 
sc_WT_P )
Idents(sce_merge) <- sce_merge$ID
sce_merge <- set_celltype(sce_merge, new.cluster.ids = c('C','M','P','C','P','M','P','M','P','C','P','M','M','P','C','M','P','C','P'))
sce_merge$group1 <- sce_merge$celltype
Idents(sce_merge) <- sce_merge$ID
sce_merge <- set_celltype(sce_merge, new.cluster.ids = c('Ts11','Ts11','Ts11','Ts15','Ts15','Ts6','Ts6','Ts6+8','Ts6+8','Ts6','Ts6','Ts8+15',
                                                         'Ts8+15','Ts8+15','Ts8','Ts8','Ts8','WT','WT'))
sce_merge$group2 <- sce_merge$celltype
sce_merge$group3 <- paste(sce_merge$group2, sce_merge$group1, sep = '_')
sce_merge$group3 <- factor_order_change(c('WT_C','Ts6_C','Ts8_C','Ts11_C','Ts15_C','WT_P','Ts6_P','Ts8_P','Ts11_P','Ts15_P','Ts6+8_P','Ts8+15_P', 
                                          'Ts6_M','Ts8_M', 'Ts11_M', 'Ts6+8_M', 'Ts8+15_M'), sce_merge$group3)
sce_merge <- subset(sce_merge, subset = ID!= 'sc_Ts6-4-P') # Discarded because of poor quality
#  Dimensionality reduction and cell annotation
sce_merge <- RunPCA(sce_merge, verbose = FALSE)
sce_merge <- FindNeighbors(sce_merge, dims = 1:30)
sce_merge <- RunUMAP(sce_merge, dims = 1:30, min.dist = 0.2, spread = 1)
sce_merge <- FindClusters(sce_merge, verbose = FALSE)
Idents(sce_merge) <- sce_merge$seurat_clusters
options(repr.plot.width=20, repr.plot.height=12)
total_marker <- c('Ptprc','Ly6g','Csf3r','Tlr4','Siglece','Cd14','Fcgr3','Itgam','Fcgr1','Ly75','Mertk','Nkg7',
                  'Klrc1','Klrb1c','Cd3e','Cd3d','Cd9','Alas2','Gypa','Pecam1','Hes1','Col3a1',
                  'Col1a1','Sdc1','Epcam','Hsp90aa1','Sox2','Npm1','Dppa5a','Zfp42','Pou5f1','Trh','Des','Tnnc2','Myl1','Sox17','Nes','Msi1','Olig1','Gfap','Sox9','Tubb3','Dcx',
                  'Map2')
DotPlot(sce_merge, features = total_marker) + RotatedAxis()
Idents(sce_merge) <- sce_merge$seurat_clusters
sce_merge <- set_celltype(sce_merge, new.cluster.ids = c('ES_Ori','ES_Ori','Macrophage','Granulocyte','ES_NSC','ES_NC','ES_Schw','ES_FB','ES_NC','ES_NSC',
                                                         'ES_Ori','ES_FB','ES_NC','ES_NSC',
                                                         'Granulocyte','ES_Stem','ES_Stem','ES_Ori','Macrophage',
                                                         'Granulocyte','Epithelium','ES_Schw','DC','Endothelium','Erythroid','ES_NC','ES_Oligo','DC',
                                                         'Erythroid','Granulocyte','Endothelium','T&ILC','ES_Ori','ES_FB','ES_Stem','Macrophage',
                                                         'T&ILC','Fibroblast','Erythroid','ES_Musc','DC','Fibroblast','ES_Schw','Granulocyte','ES_Ori',
                                                         'Macrophage','ES_FB'))

sce_merge$celltype <- factor_order_change(c('Granulocyte',"Macrophage","DC","T&ILC","Erythroid","Endothelium","Fibroblast","Epithelium",'ES_Ori',"ES_Stem","ES_NSC",'ES_Oligo','ES_Schw',"ES_NC",'ES_FB','ES_Musc'), sce_merge$celltype)
Idents(sce_merge) <- sce_merge$celltype
# Visulizing marker genes
options(repr.plot.width=18, repr.plot.height=6)
total_marker <- c('Ptprc','Ly6g','Csf3r','Tlr4','Siglece','Cd14','Fcgr3','Itgam','Fcgr1','Ly75','Mertk','Nkg7',
                  'Klrc1','Klrb1c','Cd3e','Cd3d','Cd9','Alas2','Gypa','Pecam1','Hes1','Col3a1',
                  'Col1a1','Sdc1','Epcam','Hsp90aa1','Sox2','Npm1','Dppa5a','Zfp42','Pou5f1','Trh','Sox17','Nes','Msi1','Olig1','Gfap','Sox9','Tubb3',
                  'Dcx','Map2','Des','Tnnc2','Myl1')

DotPlot(sce_merge, features = total_marker, cols = 'RdBu') + RotatedAxis()

ggsave("/total_marker.pdf", width = 18, height = 6)
# Perpare object for inferCNV
sce_merge_temp <- subset(sce_merge, subset = celltype%in%c('ES_Ori',"ES_Stem","ES_NSC",'ES_Oligo','ES_Schw',"ES_NC",'ES_FB','ES_Musc'))
# Down sample the metadata
target_mata <- sce_merge_temp@meta.data
filter_meta <- data.frame()
celltype <- unique(target_mata$group3)

target_mata$X <- rownames(target_mata)

for (i in celltype) {
    temp <- target_mata[target_mata$group3==i,]
    if (dim(temp)[1]>1500) {
        temp1 <- sample(temp$X, 1500)
        temp2 <- temp[temp$X%in%temp1,]
        filter_meta <- rbind(filter_meta, temp2)
    } else {
        filter_meta <- rbind(filter_meta, temp)
    }
    }

target_mata <- sce_merge@meta.data
target_mata$X <- rownames(target_mata)

celltype <- c('Granulocyte',"Macrophage","DC","T&ILC","Erythroid","Endothelium","Fibroblast","Epithelium")

for (i in celltype) {
    temp <- target_mata[target_mata$celltype==i,]
    if (dim(temp)[1]>1200) {
        temp1 <- sample(temp$X, 1200)
        temp2 <- temp[temp$X%in%temp1,]
        filter_meta <- rbind(filter_meta, temp2)
    } else {
        filter_meta <- rbind(filter_meta, temp)
    }
    }

sce_merge_temp <- sce_merge[,filter_meta$X]
DefaultAssay(sce_merge_temp) <- 'SCT'
saveRDS(sce_merge_temp,'subset_for_CNV.RDS')
