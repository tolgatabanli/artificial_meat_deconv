#!/usr/bin/Rscript
options(warn=-1)
setwd("~/artifical_meat_deconv")

#results_dir <- "messmer_counts/dwls_results/all"
results_dir <- "messmer_counts/dwls_results/no_teno"

library(conflicted)
conflict_prefer_all("dplyr")
library(omnideconv)
library(tidyverse)
library(parallel)

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
#col_meta$predicted.pheno <- paste0(col_meta$orig.ident, "_", basic_cell_type_names)
col_meta$predicted.pheno <- basic_cell_type_names

message("colnames(mat):", length(colnames(mat)), "\n")
message("col_meta barcodes:", nrow(col_meta), "\n")
message("Intersection:", length(intersect(colnames(mat), rownames(col_meta))), "\n")
message("In mat not in col_meta:", length(setdiff(colnames(mat), rownames(col_meta))), "\n")
message("In col_meta not in mat:", length(setdiff(rownames(col_meta), colnames(mat))), "\n")


common <- intersect(intersect(colnames(mat), rownames(col_meta)), noTenocyte_index)
#common <- intersect(colnames(mat), rownames(col_meta))
mat <- mat[, common]
col_meta <- col_meta[common, ]
gc()

message("starting with signature building")
omnideconv::build_model_dwls(single_cell_object = mat,
                             cell_type_annotations = col_meta$predicted.pheno,
                             dwls_method = "mast_optimized",
                             path = results_dir,
                             verbose = T, ncores = 6)

load(file.path(results_dir, 'Sig.RData'))

read_tpm <- function(path) {
  bulk_matrix <- data.table::fread(path)
  
  char_cols <- names(bulk_matrix)[sapply(bulk_matrix, is.character)]
  if (length(char_cols) == 0) stop("No character column found for gene IDs.")
  gene_col_name <- if ("gene_id" %in% char_cols) "gene_id" else char_cols[1]
  if (length(char_cols) > 1) {
    warning("Multiple character columns found: ", paste(char_cols, collapse = ", "),
            ". Using '", gene_col_name, "' as gene identifiers; discarding others.", call. = FALSE)
    }
  gene_col <- bulk_matrix[[gene_col_name]]
  bulk_matrix[, (char_cols) := NULL]
  mat <- as.matrix(bulk_matrix)
  rownames(mat) <- gene_col
  return(mat)
}


bulk_matrix_adipo <- read_tpm('/nfs/data2/huji_artificial_meat/adipogenesis/adipogenesis_tpm.tsv')
bulk_matrix_myo <- read_tpm('/nfs/data2/huji_artificial_meat/myogenesis/nfcore_out/star_salmon/salmon.merged.gene_tpm.tsv')

message("Intersection of Signature and Adipo genes: ", length(intersect(rownames(Sig), rownames(bulk_matrix_adipo))))
message("Intersection of Signature and Myo genes: ",  length(intersect(rownames(Sig), rownames(bulk_matrix_myo))))


deconvolute_dwls_parallel <- function(bulk, sig, ncores = 4, logfile = "deconv_progress.log") {
  chunks <- parallel::splitIndices(ncol(bulk), ncores)
  n_total <- ncol(bulk)
  results <- mclapply(chunks, function(idx) {
    chunk_results <- lapply(idx, function(i) {
      res <- deconvolute_dwls(bulk_gene_expression = bulk[, i, drop = FALSE],
                              signature = sig,
                              dwls_submethod = "DampenedWLS",
                              verbose = FALSE)
      cat(sprintf("Done: %s (%d/%d)\n", colnames(bulk)[i], i, n_total),
          file = logfile, append = TRUE)
      res
    })
    do.call(rbind, chunk_results)
  }, mc.cores = ncores)
  do.call(rbind, results)
}

message("___adipo___")
comp_adipo <- deconvolute_dwls_parallel(bulk_matrix_adipo, Sig, ncores = 38)
message("___myo___")
comp_myo <- deconvolute_dwls_parallel(bulk_matrix_myo, Sig, ncores = 26)

saveRDS(comp_adipo, file.path(results_dir, "adipo_compositions_messmer2023_dwls_signature.Rds"))
saveRDS(comp_myo, file.path(results_dir, "myo_compositions_messmer2023_dwls_signature.Rds"))


