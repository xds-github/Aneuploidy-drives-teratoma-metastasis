library(infercnv)
sce_merge_temp <- readRDS("/subset_for_CNV.RDS")
raw_matrix <- as.matrix(sce_merge_temp@assays$RNA@counts)
meta_data <- sce_merge_temp@meta.data
meta_data$barcode <- row.names(meta_data)
meta_data <- meta_data[,c('barcode','celltype')]
rownames(meta_data) <- NULL
write.table(meta_data, '/share/home/xudeshu/scanpy_dic/heter/scRNA_out/cellAnnotations.txt', sep = '\t', quote = F, row.names = F, col.names = F)
