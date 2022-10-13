library(infercnv)
sce_merge_temp <- readRDS("/subset_for_CNV.RDS")
raw_matrix <- as.matrix(sce_merge_temp@assays$RNA@counts)
meta_data <- sce_merge_temp@meta.data
meta_data$barcode <- row.names(meta_data)
meta_data <- meta_data[,c('barcode','celltype')]
rownames(meta_data) <- NULL
write.table(meta_data, '/cellAnnotations.txt', sep = '\t', quote = F, row.names = F, col.names = F)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=raw_matrix,
                                    annotations_file="/cellAnnotations.txt",
                                    delim="\t",
                                    gene_order_file="/gene_order.list",
                                    ref_group_names=c('Granulocyte',"Macrophage","DC","T&ILC","Erythroid","Endothelium","Fibroblast","Epithelium"))

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="/infer_cnv_out/", 
                             cluster_by_groups=T, 
                             denoise=F,
                             num_threads = 8,
                             HMM=F) 

sce_merge_temp <- subset(sce_merge_temp, subset = celltype%in%c('ES_Ori',"ES_Stem","ES_NSC",'ES_Oligo','ES_Schw',"ES_NC",'ES_FB','ES_Musc'))
raw_matrix <- as.matrix(sce_merge_temp@assays$RNA@counts)
meta_data <- sce_merge_temp@meta.data
meta_data$barcode <- row.names(meta_data)
meta_data <- meta_data[,c('barcode','ID')]
rownames(meta_data) <- NULL

write.table(meta_data, '/scRNA_out/cellAnnotations2.txt', sep = '\t', quote = F, row.names = F, col.names = F)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=raw_matrix,
                                    annotations_file="/scRNA_out/cellAnnotations2.txt",
                                    delim="\t",
                                    gene_order_file="/share/soft/infercnv_resource/gene_order.list",
                                    ref_group_names=c('sc_WT-C','sc_WT-P'))

# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="/infer_cnv_out_AC/",
                             cluster_by_groups=T, 
                             denoise=F,
                             num_threads = 8,
                             HMM=F)
