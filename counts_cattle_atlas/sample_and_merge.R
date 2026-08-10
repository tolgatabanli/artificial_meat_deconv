library(Seurat)
library(parallel)
library(tidyverse)

setwd("/nfs/home/students/t.tabanli/artifical_meat_deconv/counts_cattle_atlas")

rds_files <- list.files("/nfs/data2/huji_artificial_meat/cattle_atlas/NG_scRNA", full.names = T) 
rds_files <- rds_files[endsWith(rds_files, ".rds")]
celltype_to_lineage <- read_csv("cell_type_lineage_map.csv")

# matrix of tissue x cell_type -> how many cells from each type per tissue
# results <- mclapply(rds_files, function(f) {
#   obj <- readRDS(f)
#   
#   # assign label if cell_type unknown
#   #obj$cell_type <- ifelse(obj$cell_type == "Unknown cells", obj$label, obj$cell_type)
#   ct_table <- table(obj$cell_type)
#   #ct_table <- ct_table[!names(ct_table) %in% "Unknown cells"]
#   
#   ct_occurrences <- as.vector(ct_table)
#   names(ct_occurrences) <- names(ct_table)
#   
#   tissue_to_add <- as.data.frame(as.list(ct_occurrences), check.names = FALSE)
#   if (length(unique(obj$tissue)) > 1) {
#     warning("Multiple tissues detected. Using the file name instead: ", f)
#     rownames(tissue_to_add) <- f
#   } else {
#     rownames(tissue_to_add) <- unique(obj$tissue)
#   }
#   
#   tissue_to_add
# }, mc.cores = 20)
# 
# tissue_by_ct <- as.matrix(dplyr::bind_rows(results))
# tissue_by_ct[is.na(tissue_by_ct)] <- 0
# # Normalize by column sums -> find how much each tissue contributes to each cell type
# tissue_contributions <- prop.table(tissue_by_ct, margin = 2)
# 
# readr::write_csv(tibble::rownames_to_column(as.data.frame(tissue_by_ct), "tissue"), file = "tissue_by_ct.csv", col_names = T)
# readr::write_csv(tibble::rownames_to_column(as.data.frame(tissue_contributions), "tissue"), file = "tissue_contributions.csv", col_names = T)

tissue_by_ct <- readr::read_csv("tissue_by_ct.csv", col_names = T) %>%
  column_to_rownames("tissue")
tissue_contributions <- readr::read_csv("tissue_contributions.csv", col_names = T) %>%
  column_to_rownames("tissue")

cat("Total number: ", sum(colSums(tissue_by_ct)), "\n")
# scale down by a factor of 30
sampling_matrix <- ceiling(tissue_by_ct / 30)


# sampling_matrix is a tissue x cell_type matrix, specifying
# how many cells of each type from each tissue to sample
sample_from_seurat <- function(seurat_obj, sampling_matrix) {
  tissue <- unique(seurat_obj$tissue)
  if (length(tissue) > 1) {
    warning("Multiple tissues detected. Using the first instead: ", tissue[1])
    tissue <- tissue[1]
  }
  
  partitions <- split(colnames(seurat_obj), seurat_obj$cell_type)
  #partitions[["Unknown cells"]] <- NULL # we're ignoring unknown cells to decrease analytic burden:)
  
  cols_to_select <- sapply(names(partitions),
                           function(ct_name) sample(partitions[[ct_name]],
                                                    size = sampling_matrix[tissue, ct_name]))
  cols_to_select <- unlist(cols_to_select)
  seurat_obj <- seurat_obj[, cols_to_select]
  
  return(seurat_obj)
}

sampled_list <- mclapply(rds_files, function(f) {
  set.seed(digest::digest2int(f)) # reproducibility
  
  obj <- readRDS(f)
  sampled <- sample_from_seurat(obj, sampling_matrix)
  sampled
}, mc.cores = 19)
# max file is 11 GB, min core in a node is 64, and min memory is 384 GB
# takes ca. 15 mins

final_obj <- merge(sampled_list[[1]], y = sampled_list[-1])
#final_obj$label <- trimws(final_obj$label)
final_obj$cell_type <- trimws(final_obj$cell_type)
final_obj$cell_type <- gsub("'", "", final_obj$cell_type) # Henle's loop causes issues in column parse
final_obj$cell_type <- gsub("/", "_", final_obj$cell_type) # inner/outer -> remove slash

# convert cattle IDs to ENSEMBL
# remove multimappers (symbol to multiple Ensembl IDs)
# remove same-mappers (multiple symbols to same Ensembl ID)
message("Mapping ids...")
library(org.Bt.eg.db)
library(AnnotationDbi)
mat <- as.matrix(GetAssayData(final_obj, assay = "RNA", layer = "counts"))

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
final_obj[["RNA"]] <- CreateAssayObject(counts = aggregated_mat)

# add cell lineage categories from the publication
lin <- celltype_to_lineage$cell_lineage[match(final_obj$cell_type, celltype_to_lineage$cell_type)]
final_obj$cell_lineage <- ifelse(is.na(lin), as.character(final_obj$cell_type), lin)

cat("Final object:\n")
final_obj

saveRDS(final_obj, file = "sampled_cattle_atlas.RDS")

# write in mtx format for python reads
Matrix::writeMM(as(aggregated_mat, "dgCMatrix"), "sampled_cattle_atlas/counts.mtx")
write.csv(data.frame(gene = rownames(aggregated_mat)), "sampled_cattle_atlas/genes.csv", row.names = FALSE)
write.csv(final_obj@meta.data, "sampled_cattle_atlas/metadata.csv")
