#!/usr/bin/Rscript

options(warn=-1)
setwd("~/artifical_meat_deconv")

library(conflicted)
conflict_prefer_all("dplyr")
library(omnideconv)
library(tidyverse)

counts <- data.table::fread('/nfs/home/students/t.tabanli/artifical_meat_deconv/messmer_counts/GSE211428_counts.tsv.gz')
mat <- as.matrix(counts[, -1])
rownames(mat) <- counts$V1

col_meta <- data.table::fread('/nfs/data2/huji_artificial_meat/deconvolution_data/external_sc_references/messmer2023_GSE211428/GSE211428_combined_meta.tsv.gz')
col_meta <- as.data.frame(col_meta)
col_meta$seurat_barcode <- gsub("^(.*)\\.([0-9]+)_([0-9]+)$", "Inj\\3_\\1-\\2", col_meta$Unique_barcode)
rownames(col_meta) <- col_meta$seurat_barcode
col_meta$V1 <- NULL
col_meta$seurat_barcode <- NULL

noTenocyte_index <- rownames(col_meta)[col_meta$predicted.pheno != "Tenocytes"]

# make the names path-compatible AND include timepoints
basic_cell_type_names <- gsub("/", "-", col_meta$predicted.pheno)
basic_cell_type_names[grep("Monocytes", basic_cell_type_names)] <- "Monocytes"
basic_cell_type_names <- gsub(" ", "-", basic_cell_type_names)
col_meta$predicted.pheno <- paste0(col_meta$orig.ident, "_", basic_cell_type_names)

cat("colnames(mat):", length(colnames(mat)), "\n")
cat("col_meta barcodes:", nrow(col_meta), "\n")
cat("Intersection:", length(intersect(colnames(mat), rownames(col_meta))), "\n")
cat("In mat not in col_meta:", length(setdiff(colnames(mat), rownames(col_meta))), "\n")
cat("In col_meta not in mat:", length(setdiff(rownames(col_meta), colnames(mat))), "\n")


common <- intersect(intersect(colnames(mat), rownames(col_meta)), noTenocyte_index)
col_meta <- col_meta[common, ]
col_meta <- col_meta[col_meta$orig.ident == "Injection1", ]
mat <- mat[, rownames(col_meta)]
gc()

message("starting with signature building")
omnideconv::build_model_dwls(single_cell_object = mat,
                             cell_type_annotations = col_meta$predicted.pheno,
                             dwls_method = "mast_optimized",
                             path = "sig_dwls_messmer_recreated_counts_timepoint1_noTeno/",
                             verbose = T, ncores = 12)

load('sig_dwls_messmer_recreated_counts_timepoint1_noTeno/Sig.RData')

read_tpm <- function(path) {
  bulk_matrix <- data.table::fread(path)
  
  char_cols <- names(bulk_matrix)[sapply(bulk_matrix, is.character)]
  if (length(char_cols) == 0) stop("No character column found for gene IDs.")
  
  gene_col_name <- if ("gene_id" %in% char_cols) "gene_id" else char_cols[1]
  
  if (length(char_cols) > 1) {
    warning("Multiple character columns found: ", 
            paste(char_cols, collapse = ", "),
            ". Using '", gene_col_name, "' as gene identifiers; discarding others.",
            call. = FALSE)
  }
  
  gene_col <- bulk_matrix[[gene_col_name]]
  bulk_matrix[, (char_cols) := NULL]
  
  mat <- as.matrix(bulk_matrix)
  rownames(mat) <- gene_col
  return(mat)
}


bulk_matrix_adipo <- read_tpm('/nfs/data2/huji_artificial_meat/adipogenesis/adipogenesis_tpm.tsv')
bulk_matrix_myo <- read_tpm('/nfs/data2/huji_artificial_meat/myogenesis/nfcore_out/star_salmon/salmon.merged.gene_tpm.tsv')

length(intersect(rownames(Sig), rownames(bulk_matrix_adipo)))
length(intersect(rownames(Sig), rownames(bulk_matrix_myo)))

comp_adipo <- deconvolute_dwls(bulk_gene_expression = bulk_matrix_adipo,
                               signature = Sig,
                               dwls_submethod = "DampenedWLS",
                               verbose = T)
comp_myo <- deconvolute_dwls(bulk_gene_expression = bulk_matrix_myo,
                             signature = Sig,
                             dwls_submethod = "DampenedWLS",
                             verbose = T)

saveRDS(comp_adipo, "sig_dwls_messmer_recreated_counts_timepoint1_noTeno/adipo_compositions_messmer2023_dwls_signature.Rds")
saveRDS(comp_myo, "sig_dwls_messmer_recreated_counts_timepoint1_noTeno/myo_compositions_messmer2023_dwls_signature.Rds")


