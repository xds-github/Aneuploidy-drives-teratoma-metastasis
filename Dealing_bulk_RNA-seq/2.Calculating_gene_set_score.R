library(DESeq2)
library(apeglm)
library(IHW)
library(org.Hs.eg.db)
library(biomaRt)
library(curl)
library(ggplot2)
library(Seurat)
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
human_syb2mus_syb <- function(gene_list){
  gene_list <- unique(gene_list)
  refer_gene_list <- read.table('Ensembl_human_mouse_ID_change.txt', sep = '\t', header = T,check.names = F)
  refer_gene_list <- refer_gene_list[refer_gene_list$`Gene name`%in%gene_list,]
  refer_gene_list <- refer_gene_list[,c('Gene name','Mouse gene name')]
  refer_gene_list <- na.omit(refer_gene_list)
  refer_gene_list <- refer_gene_list[refer_gene_list$`Gene name`!='',]
  refer_gene_list <- refer_gene_list[refer_gene_list$`Mouse gene name`!='',]
  refer_gene_list <- refer_gene_list[!duplicated(refer_gene_list$`Mouse gene name`),]
  gene_list <- gene_list[gene_list%in%refer_gene_list$`Gene name`]
  out_list <- data.frame()
  for (i in gene_list) {
    out_list <- rbind(out_list,refer_gene_list[refer_gene_list$`Gene name`==i,])
  }
  return(out_list$`Mouse gene name`)
}
# ----------------------- Read in matrix -----------------------------------
meta_data <- read.csv('sample_info_230816.csv',check.names = F)
meta_data$group2 <- paste(meta_data$chr, meta_data$day, sep = '_')
sample_N <- meta_data$sampleName
count_matrix <- read.table(paste0('BRNA_Ts11-EB-D8-1','_ReadsPerGene.txt'),check.names = F)
count_matrix$V3 <- count_matrix$V4 <- NULL
for (i in rownames(meta_data)) {
  temp <- read.table(paste0(meta_data[i,'Dir'],meta_data[i,'fileName']),check.names = F)
  temp <- plyr::rename(temp, c('V2'=meta_data[i,'sampleName']))
  temp$V3 <- temp$V4 <- NULL
  count_matrix <- merge(count_matrix, temp, by = 'V1')  
}
rownames(count_matrix) <- count_matrix$V1
count_matrix$V1 <- count_matrix$V2 <- NULL
# ----------------------- transform ID-----------------------------------
refer_gene_list <- read.table('Ensembl_mouse_ID.txt', sep = '\t', header = T,check.names = F)
refer_gene_list <- refer_gene_list[refer_gene_list$`Gene stable ID`%in%rownames(count_matrix),]
refer_gene_list <- refer_gene_list[!duplicated(refer_gene_list$`Gene stable ID`),]
dup_list <- unique(refer_gene_list$`Gene name`[duplicated(refer_gene_list$`Gene name`)])
dup_list <- refer_gene_list[refer_gene_list$`Gene name`%in%dup_list,'Gene stable ID']
refer_gene_list$`Gene_name2` <- paste(refer_gene_list$`Gene stable ID`,refer_gene_list$`Gene name` ,sep = '_')
for (i in dup_list) {
  refer_gene_list[refer_gene_list$`Gene stable ID`==i,'Gene name'] <-  refer_gene_list[refer_gene_list$`Gene stable ID`==i,'Gene_name2']
}
refer_gene_list$`Gene_name2` <- NULL
count_matrix <- count_matrix[refer_gene_list$`Gene stable ID`,]
rownames(count_matrix) <- refer_gene_list$`Gene name`
saveRDS(count_matrix, 'total_raw_count_matrix.RDS')
# ----------------------- Generating dds object -----------------------------------
rownames(meta_data) <- meta_data$sampleName
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = meta_data, design = ~ group2)
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]
dds <- DESeq(dds)
saveRDS(dds, 'dds.RDS')
# remove batch effect and QC
vsd <- vst(dds, blind=FALSE)
# ----------------------- Generating sce object -----------------------------------
sce <- CreateSeuratObject(counts = count_matrix, meta.data = meta_data)
sce <- PercentageFeatureSet(sce, pattern = "^mt-", col.name = "percent.mt")
VlnPlot(sce, features = 'percent.mt')
sce <- SCTransform(sce)
sce$group2 <- factor_order_change(c("WT_D8","Ts6_D8","Ts8_D8","Ts11_D8","Ts15_D8"), sce$group2)
Idents(sce) <- sce$group2
# Calculating gene set score
temp <- read.csv('Targeted_GO_gene_list.csv')
gene_list <- total_list[total_list$GO=='GO:0007398','gene_list']
temp <- human_syb2mus_syb(gene_list)
gene_list <- list(temp)
sce <- AddModuleScore(sce,features = gene_list,name = "Ectoderm_development_score")
gene_list <- total_list[total_list$GO=='GO:0007498','gene_list']
temp <- human_syb2mus_syb(gene_list)
gene_list <- list(temp)
sce <- AddModuleScore(sce,features = gene_list,name = "Mesoderm_development_score")
gene_list <- total_list[total_list$GO=='GO:0007492','gene_list']
temp <- human_syb2mus_syb(gene_list)
gene_list <- list(temp)
sce <- AddModuleScore(sce,features = gene_list,name = "Endoderm_development_score")
gene_list <- total_list[total_list$GO=='GO:0019827','gene_list'] #Stem cell population maintenance
temp <- human_syb2mus_syb(gene_list)
gene_list <- list(temp)
sce <- AddModuleScore(sce,features = gene_list,name = "Stem_score")
gene_list <- total_list[total_list$GO=='GO:0036493','gene_list'] #Positive regulation of translation in response to endoplasmic reticulum stress
temp <- human_syb2mus_syb(gene_list)
gene_list <- list(temp)
sce <- AddModuleScore(sce,features = gene_list,name = "ER_stress_score")
saveRDS(sce, 'sce_with_gene_set_socre.RDS')
