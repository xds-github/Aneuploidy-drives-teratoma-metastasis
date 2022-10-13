library(Seurat)
library(reshape2)
library(future, quietly = T)
plan("multiprocess", workers = 8)
options(future.globals.maxSize = 15000 * 1024^2)
##For mouse
Mo_datain <- function(dic = "", mito = "^mt-"){
  pbmc.counts <- Read10X(data.dir = dic)
  pbmc <- CreateSeuratObject(counts = pbmc.counts)
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = mito)
  library(scater)
  qc.lib2 <- isOutlier(pbmc$nCount_RNA, log=TRUE, type="lower")
  print("nCount_RNA")
  print(attr(qc.lib2, "thresholds"))
  qc.nexprs2 <- isOutlier(pbmc$nFeature_RNA, log=TRUE, type="lower")
  print("nFeature_RNA")
  print(attr(qc.nexprs2, "thresholds"))
  qc.mito2 <- isOutlier(pbmc$percent.mt, type="higher")
  print("percent.mt")
  print(attr(qc.mito2, "thresholds"))
  return(pbmc)
}
set_celltype <- function(sce, new.cluster.ids = c(1,2)){
  names(new.cluster.ids) <- levels(sce)
  sce <- RenameIdents(sce, new.cluster.ids)
  sce$celltype <- Idents(sce)
  return(sce)
}
factor_order_change <- function(new_order, old_factor){
  new_order <- factor(1:length(new_order),labels = new_order)
  new_factor <- factor(old_factor,levels = levels(new_order))
  return(new_factor)
}
# Initiate data
Ts8_M <- Mo_datain("/201113C_Ts8_lung/filtered_feature_bc_matrix/")
Ts8_M$ID <- 'Ts8-M'

Ts11_M <- Mo_datain("/201113B_TS11_lung/filtered_feature_bc_matrix/")
Ts11_M$ID <- 'Ts11-M'

Ts6_8_M <- Mo_datain("/201113C_Ts6_8_lung/filtered_feature_bc_matrix/")
Ts6_8_M$ID <- 'Ts6_8-M'

Ts8_15_M1 <- Mo_datain("/191522A_Ts8_lung/filtered_feature_bc_matrix/")
Ts8_15_M1$ID <- 'Ts8_15-M1'

Ts8_15_M2 <- Mo_datain("/201113A_Ts15_lung/filtered_feature_bc_matrix/")
Ts8_15_M2$ID <- 'Ts8_15-M2'

total.list <- list(Ts8_M,Ts11_M,Ts6_8_M,Ts8_15_M1,Ts8_15_M2)
samplelist <- c('Ts8_M','Ts11_M','Ts6_8_M','Ts8_15_M1','Ts8_15_M2')
rm(Ts8_M,Ts11_M,Ts6_8_M,Ts8_15_M1,Ts8_15_M2)
# Merge samples
sce_merge <- merge(x =  total.list[[1]], y = c(total.list[[2]],total.list[[3]],total.list[[4]],total.list[[5]]), add.cell.ids=c(0,1,2,3,4))
rm(total.list)
sce_merge <- subset(sce_merge, subset = percent.mt<=25)
sce_merge <- subset(sce_merge, subset = nFeature_RNA >200)
# run sctransform
sce_merge <- SCTransform(sce_merge, verbose = FALSE)
# These are now standard steps in the Seurat workflow for visualization and clustering
sce_merge <- RunPCA(sce_merge, verbose = FALSE)
sce_merge <- RunUMAP(sce_merge, dims = 1:20, verbose = FALSE)
sce_merge <- FindNeighbors(sce_merge, dims = 1:30, verbose = FALSE)
sce_merge <- FindClusters(sce_merge, verbose = FALSE, resolution = 0.5)
DimPlot(sce_merge)
VlnPlot(sce_merge, features = c('Luciferin','GFP','Mki67'))
total_marker <- c('Ptprc','Ly6g','Csf3r','Tlr4','Siglece','Cd14','Fcgr3','Itgam','Fcgr1','Fcgr2b','Ly75','Mertk','Nkg7','Klrc1','Klrb1c','Cd3e','Cd3d','Cd9','Itga2b','Alas2','Gypa','Pecam1','Hes1','Col3a1',
                  'Col1a1','Epcam','Sox2','Sdc1','Dppa5a','Zfp42','Nes','Msi1','Gfap','Sox9','Tubb3','Dcx','Map2','Hif1a','Mapkapk5','Hspb1')
DotPlot(sce_merge, features = total_marker) + RotatedAxis()
Idents(sce_merge) <- sce_merge$seurat_clusters
# Annotating clusters
sce_merge <- set_celltype(sce_merge, new.cluster.ids = c('Granulocyte','Granulocyte','DC',"Macrophage","AC_NC","AC_NSC",'AC_FB','AC',"Erythroid",'Granulocyte',"Endothelium","Erythroid","Macrophage",'AC',"T_ILC","AC_Schw","Macrophage",'AC_Epi',"Erythroid","Fibroblast",'Granulocyte','AC',"T_ILC","Epithelium","Epithelium",'AC_Epi','Granulocyte','Granulocyte',"Fibroblast",'DC',"Erythroid",'Granulocyte'))
sce_merge$celltype <- factor_order_change(c('Granulocyte',"Macrophage","DC","T_ILC","Erythroid","Endothelium",
                                            "Fibroblast","Epithelium",'AC_Epi',"AC","AC_NSC",'AC_Schw',"AC_NC","AC_FB"), sce_merge$celltype)
Idents(sce_merge) <- sce_merge$celltype
DimPlot(sce_merge, label = T)
DimPlot(sce_merge, cols = c("#023FA5","#7D87B9","#BEC1D4","#D6BCC0","#BB7784","#8E063B","#4A6FE3","#8595E1","#B5BBE3","#E6AFB9",
                            "#E07B91","#D33F6A","#11C638","#8DD593"), pt.size = 1.5)
DimPlot(sce_merge, cols = c("#4A6FE3","#8595E1","#B5BBE3","#E6AFB9",
                            "#E07B91","#D33F6A","#11C638","#8DD593","#023FA5","#7D87B9","#BEC1D4","#D6BCC0","#BB7784","#8E063B"), pt.size = 1.5)
sce_merge <- saveRDS(sce_merge,'E:/heteroploidy/single_cell/total_merge3.RDS')
# Visulizing targeted genes
FeaturePlot(sce_merge, features = 'GFP', pt.size = 4, order = T, cols = featureplot4, raster = T)
FeaturePlot(sce_merge, features = 'Luciferin', pt.size = 4, order = T, cols = featureplot4, raster = T)
FeaturePlot(sce_merge, features = 'Nanog', pt.size = 4, order = T, cols = featureplot4, raster = T)
FeaturePlot(sce_merge, features = c('Tfec','Habp2','Dkk1','Cubn','Pramel6'), pt.size = 4, order = T, cols = featureplot4, raster = T, label = T)
# Prepare for cell distribution
temp1 <- rownames(sce_merge@meta.data)
temp2 <- substring(temp1, 1,1)
temp3 <- substring(temp1, 3,20)
sce_merge$barcode <- paste(temp3, temp2, sep = '-')
write.csv(sce_merge@meta.data[,c('celltype','barcode','ID')],'meta_data.csv', row.names = F, quote = F)
sce_merge <- saveRDS(sce_merge,'/single_cell/total_merge.RDS')
