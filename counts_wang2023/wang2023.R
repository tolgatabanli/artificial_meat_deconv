library(tidyverse)
library(omnideconv)
library(Seurat)
setwd("~/artifical_meat_deconv/counts_wang2023")

# read raw counts
sc <- Seurat::ReadMtx(mtx = "/nfs/data2/huji_artificial_meat/deconvolution_data/external_sc_references/Wang2023/raw_count_export_clean_metadata/matrix.mtx.gz",
                      cells = "/nfs/data2/huji_artificial_meat/deconvolution_data/external_sc_references/Wang2023/raw_count_export_clean_metadata/barcodes.tsv.gz",
                      features = "/nfs/data2/huji_artificial_meat/deconvolution_data/external_sc_references/Wang2023/raw_count_export_clean_metadata/features.tsv.gz")

# cell metadata
metadata_of_normalized <- data.table::fread("/nfs/data2/huji_artificial_meat/deconvolution_data/external_sc_references/Wang2023/GSE205347_Complete_10x_cattle_Metadata.csv.gz")

# filter to cells existing in metadata of normalized/processed data
stopifnot(length(intersect(metadata_of_normalized$V1, colnames(sc))) == nrow(metadata_of_normalized$V1))
sc <- sc[, colnames(sc) %in% metadata_of_normalized$V1]

# check rename_clusters and lables
cat("In metadata, the columns 'rename_clusters' and 'lables' differ in following number of rows: ")
metadata_of_normalized %>%
  mutate(rename_cluster = gsub("\\d+\\.", "", rename_cluster)) %>%
  distinct(rename_cluster, lables) %>%
  filter(rename_cluster != lables) %>%
  nrow()

sc <- Seurat::CreateSeuratObject(counts = sc, meta.data = metadata_of_normalized)

sc <- NormalizeData(sc)
sc <- FindVariableFeatures(sc)
sc <- ScaleData(sc)
sc <- RunPCA(sc)
sc <- RunUMAP(sc, dims = 1:30)
DimPlot(sc, group.by = "lables", label = T) + guides(color = "none")

ggsave("figures/wang2023_umap_plot.png", width = 12)

# convert IDs to ENSEMBL
# remove multimappers (symbol to multiple Ensembl IDs)
# aggregate same-mappers (multiple symbols to same Ensembl ID)
message("Mapping ids...")
library(org.Bt.eg.db)
library(AnnotationDbi)
mat <- as.matrix(GetAssayData(sc, assay = "RNA", layer = "counts"))

mixed_names <- rownames(mat)
is_ensembl <- grepl("^ENSBTAG", mixed_names)

gene_symbols <- mapIds(org.Bt.eg.db, keys = mixed_names[!is_ensembl],
                       column = "ENSEMBL", keytype = "SYMBOL", multiVals = "asNA")

mixed_names[!is_ensembl] <- ifelse(is.na(gene_symbols),
                                   mixed_names[!is_ensembl],
                                   gene_symbols)
rownames(mat) <- mixed_names

# keep only ENSEMBLs
mat <- mat[grepl("^ENSBTAG", rownames(mat)), ]

# merge counts of duplicate Ensembl IDs
aggregated_mat <- rowsum(mat, group = rownames(mat))

message(sprintf(
  "ID mapping: %d/%d rows already Ensembl IDs, %d gene symbols mapped uniquely to an Ensembl ID,
  %d symbols could not be mapped and dropped. %d ids at the end were duplicated (first not counted) and their counts are now aggregated.",
  sum(is_ensembl), length(mixed_names), sum(!is.na(gene_symbols)), sum(is.na(gene_symbols)), sum(base::duplicated(rownames(mat)))
))

# overwrite counts of the original object
sc[["RNA"]] <- CreateAssayObject(counts = aggregated_mat)
cat("Final object:\n")
sc


sc <- NormalizeData(sc)
sc <- FindVariableFeatures(sc)
sc <- ScaleData(sc)
sc <- RunPCA(sc)
sc <- RunUMAP(sc, dims = 1:30)
DimPlot(sc, group.by = "lables", label = T) + guides(color = "none")
ggsave("figures/wang2023_umap_after_id_mapping.png", width = 12)

saveRDS(sc, file = "wang_ids_mapped.RDS")

# write in mtx format for python reads
Matrix::writeMM(as(aggregated_mat, "dgCMatrix"), "wang_ids_mapped/counts.mtx")
write.csv(data.frame(gene = rownames(aggregated_mat)), "wang_ids_mapped/genes.csv", row.names = FALSE)
write.csv(sc@meta.data, "wang_ids_mapped/metadata.csv")

