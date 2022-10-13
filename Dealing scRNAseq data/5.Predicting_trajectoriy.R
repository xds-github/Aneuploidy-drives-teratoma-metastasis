library(Seurat)
library(reshape2)
library(monocle3)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggplot2)
sce_merge <- readRDS('/subset_for_CNV.RDS')
sce_merge <- subset(sce_merge, subset = celltype%in%c('ES_Ori',"ES_Stem","ES_NSC",'ES_Oligo','ES_Schw',"ES_NC",'ES_FB','ES_Musc'))
DefaultAssay(sce_merge) <- 'RNA'
sce_merge <- subset(sce_merge, subset = group2%in%c('Ts11','Ts15','Ts6','Ts6+8','Ts8+15','Ts8'))
data <- GetAssayData(sce_merge, assay = 'RNA', slot = 'counts')
cell_metadata <- sce_merge@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 100)
cds <- align_cds(cds, alignment_group = "ID")
cds <- reduce_dimension(cds, reduction_method = 'UMAP',
                        umap.min_dist = 0.4,
                        umap.n_neighbors = 50)
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "celltype") +
  scale_color_manual(values=c('#023FA5','#7D87B9','#E07B91','#D33F6A','#8DD593',
                              '#C6DEC7','#EAD3C6','#F0B98D'))
cds <- cluster_cells(cds)
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "partition")
cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
