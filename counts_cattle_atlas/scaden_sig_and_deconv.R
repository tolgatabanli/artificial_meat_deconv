library(Seurat)
library(omnideconv)
library(tidyverse)
library(parallel)
source("~/artifical_meat_deconv/toolkit.R")


setwd("/nfs/home/students/t.tabanli/artifical_meat_deconv/counts_cattle_atlas")

message("Reading cattle...")
cattle <- readRDS("sampled_cattle_atlas.RDS")
cattle$cell_type <- trimws(cattle$cell_type)
cattle$cell_type <- gsub("'", "", cattle$cell_type) # Henle's loop causes issues in column parse
dim(cattle)

mat <- as.matrix(GetAssayData(cattle, assay = "RNA", layer = "counts"))

# convert cattle IDs to ENSEMBL
# remove multimappers (symbol to multiple Ensembl IDs)
# remove same-mappers (multiple symbols to same Ensembl ID)
message("Mapping ids...")
library(org.Bt.eg.db)
library(AnnotationDbi)

mixed_names <- rownames(mat)
is_ensembl <- grepl("^ENSBTAG", mixed_names)

gene_symbols <- mapIds(org.Bt.eg.db, keys = mixed_names[!is_ensembl],
                       column = "ENSEMBL", keytype = "SYMBOL", multiVals = "asNA")

mixed_names[!is_ensembl] <- ifelse(is.na(gene_symbols),
                                   mixed_names[!is_ensembl],
                                   gene_symbols)
rownames(mat) <- mixed_names

# keep only deduplicated ENSEMBLs
is_ensembl_row <- grepl("^ENSBTAG", rownames(mat))
is_dup_row <- base::duplicated(rownames(mat))
mat <- mat[is_ensembl_row & !is_dup_row, ]

message(sprintf(
  "ID mapping: %d/%d rows already Ensembl IDs, %d gene symbols mapped uniquely to an Ensembl ID,
  %d symbols could not be mapped and dropped. %d ids at the end were duplicated and dropped as well.",
  sum(is_ensembl), length(mixed_names), sum(!is.na(gene_symbols)), sum(is.na(gene_symbols)), sum(is_ensembl_row & is_dup_row)
))

message("Reading tpms...")

# Read TPM matrices and keep only relevant samples.
bulk_matrix_adipo <- read_tpm('/nfs/data2/huji_artificial_meat/adipogenesis/adipogenesis_tpm.tsv')
bulk_meta_adipo <- read.csv("/nfs/data2/huji_artificial_meat/adipogenesis/sample_sheet_cleaned.csv", header = T)
bulk_matrix_adipo <- bulk_matrix_adipo[, bulk_meta_adipo$sample_id]

bulk_matrix_myo <- read_tpm('/nfs/data2/huji_artificial_meat/myogenesis/nfcore_out/star_salmon/salmon.merged.gene_tpm.tsv')
bulk_meta_myo <- read.csv("/nfs/data2/huji_artificial_meat/myogenesis/sample_sheet_minimal.csv", header = T)
# weird sample id problem
bulk_matrix_myo <- bulk_matrix_myo[, gsub(".*_", "", colnames(bulk_matrix_myo)) %in% gsub(".*_", "", bulk_meta_myo$sample_id)]

message("Starting training...")
#omnideconv::build_model_scaden(single_cell_object = mat,
#                               cell_type_annotations = cattle$cell_type,
#                               bulk_gene_expression = bulk_matrix_adipo,
#                               model_path = "sig_scaden",
#                               verbose = T)

message("Signature building finished")


comp_adipo <- omnideconv::deconvolute_scaden(signature = "sig_scaden",
                                             bulk_gene_expression = bulk_matrix_adipo,
                                             verbose = T)
comp_myo <- omnideconv::deconvolute_scaden(signature = "sig_scaden",
                                           bulk_gene_expression = bulk_matrix_myo,
                                           verbose = T)

saveRDS(comp_adipo, "sig_scaden/adipo_compositions_cattleatlas_scaden.Rds")
saveRDS(comp_myo, "sig_scaden/myo_compositions_cattleatlas_scaden.Rds")

message("Deconvolution finished")

sessionInfo()
