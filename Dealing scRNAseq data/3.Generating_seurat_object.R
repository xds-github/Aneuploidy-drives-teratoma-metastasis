library(dplyr, quietly = T)
library(Seurat, quietly = T)
library(patchwork, quietly = T)
library(dplyr, quietly = T)
library(ggplot2)
# Reading in data and creating merged seurat object
sc_Ts11_C <- Mo_datain("/share/home/xudeshu/scanpy_dic/heter/scRNA_out/sc_Ts11-C/outs/filtered_feature_bc_matrix/")
sc_Ts11_C$ID <- 'sc_Ts11-C'
sc_Ts11_M <- Mo_datain("/share/home/xudeshu/scanpy_dic/heter/scRNA_out/sc_Ts11-M/outs/filtered_feature_bc_matrix/")
sc_Ts11_M$ID <- 'sc_Ts11-M'
sc_Ts11_P <- Mo_datain("/share/home/xudeshu/scanpy_dic/heter/scRNA_out/sc_Ts11-P/outs/filtered_feature_bc_matrix/")
sc_Ts11_P$ID <- 'sc_Ts11-P'
sc_Ts15_C <- Mo_datain("/share/home/xudeshu/scanpy_dic/heter/scRNA_out/sc_Ts15-C/outs/filtered_feature_bc_matrix/")
sc_Ts15_C$ID <- 'sc_Ts15-C'
sc_Ts15_P <- Mo_datain("/share/home/xudeshu/scanpy_dic/heter/scRNA_out/sc_Ts15-P/outs/filtered_feature_bc_matrix/")
sc_Ts15_P$ID <- 'sc_Ts15-P'
sc_Ts6_4_M <- Mo_datain("/share/home/xudeshu/scanpy_dic/heter/scRNA_out/sc_Ts6-4-M/outs/filtered_feature_bc_matrix/")
sc_Ts6_4_M$ID <- 'sc_Ts6-4-M'
sc_Ts6_4_P <- Mo_datain("/share/home/xudeshu/scanpy_dic/heter/scRNA_out/sc_Ts6-4-P/outs/filtered_feature_bc_matrix/")
sc_Ts6_4_P$ID <- 'sc_Ts6-4-P'
sc_Ts6_8_M <- Mo_datain("/share/home/xudeshu/scanpy_dic/heter/scRNA_out/sc_Ts6_8-M/outs/filtered_feature_bc_matrix/")
sc_Ts6_8_M$ID <- 'sc_Ts6_8-M'
sc_Ts6_8_P <- Mo_datain("/share/home/xudeshu/scanpy_dic/heter/scRNA_out/sc_Ts6-8-P/outs/filtered_feature_bc_matrix/")
sc_Ts6_8_P$ID <- 'sc_Ts6-8-P'
sc_Ts6_C <- Mo_datain("/share/home/xudeshu/scanpy_dic/heter/scRNA_out/sc_Ts6-C/outs/filtered_feature_bc_matrix/")
sc_Ts6_C$ID <- 'sc_Ts6-C'
sc_Ts6_P <- Mo_datain("/share/home/xudeshu/scanpy_dic/heter/scRNA_out/sc_Ts6-P/outs/filtered_feature_bc_matrix/")
sc_Ts6_P$ID <- 'sc_Ts6-P'
sc_Ts8_15_M1 <- Mo_datain("/share/home/xudeshu/scanpy_dic/heter/scRNA_out/sc_Ts8_15-M1/outs/filtered_feature_bc_matrix/")
sc_Ts8_15_M1$ID <- 'sc_Ts8_15-M1'
sc_Ts8_15_M2 <- Mo_datain("/share/home/xudeshu/scanpy_dic/heter/scRNA_out/sc_Ts8_15-M2/outs/filtered_feature_bc_matrix/")
sc_Ts8_15_M2$ID <- 'sc_Ts8_15-M2'
sc_Ts8_15_P <- Mo_datain("/share/home/xudeshu/scanpy_dic/heter/scRNA_out/sc_Ts8-15-P/outs/filtered_feature_bc_matrix/")
sc_Ts8_15_P$ID <- 'sc_Ts8-15-P'
sc_Ts8_C <- Mo_datain("/share/home/xudeshu/scanpy_dic/heter/scRNA_out/sc_Ts8-C/outs/filtered_feature_bc_matrix/")
sc_Ts8_C$ID <- 'sc_Ts8-C'
sc_Ts8_M <- Mo_datain("/share/home/xudeshu/scanpy_dic/heter/scRNA_out/sc_Ts8-M/outs/filtered_feature_bc_matrix/")
sc_Ts8_M$ID <- 'sc_Ts8-M'
sc_Ts8_P <- Mo_datain("/share/home/xudeshu/scanpy_dic/heter/scRNA_out/sc_Ts8-P/outs/filtered_feature_bc_matrix/")
sc_Ts8_P$ID <- 'sc_Ts8-P'
sc_WT_C <- Mo_datain("/share/home/xudeshu/scanpy_dic/heter/scRNA_out/sc_WT-C/outs/filtered_feature_bc_matrix/")
sc_WT_C$ID <- 'sc_WT-C'
sc_WT_P <- Mo_datain("/share/home/xudeshu/scanpy_dic/heter/scRNA_out/sc_WT-P/outs/filtered_feature_bc_matrix/")
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

