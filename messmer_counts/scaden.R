library(Seurat)
library(omnideconv)
library(tidyverse)
source("~/artifical_meat_deconv/toolkit.R")

setwd("/nfs/home/students/t.tabanli/artifical_meat_deconv/messmer_counts/")

# Reading reference
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

# make the names path-compatible
basic_cell_type_names <- gsub("/", "-", col_meta$predicted.pheno)
basic_cell_type_names[grep("Monocytes", basic_cell_type_names)] <- "Monocytes"
basic_cell_type_names <- gsub(" ", "-", basic_cell_type_names)
col_meta$predicted.pheno <- basic_cell_type_names

message("colnames(mat):", length(colnames(mat)), "\n")
message("col_meta barcodes:", nrow(col_meta), "\n")
message("Intersection:", length(intersect(colnames(mat), rownames(col_meta))), "\n")
message("In mat not in col_meta:", length(setdiff(colnames(mat), rownames(col_meta))), "\n")
message("In col_meta not in mat:", length(setdiff(rownames(col_meta), colnames(mat))), "\n")

common <- intersect(intersect(colnames(mat), rownames(col_meta)), noTenocyte_index)
mat <- mat[, common]
col_meta <- col_meta[common, ]
gc()

# Read bulk data
message("Reading tpms...")

# Read TPM matrices and keep only relevant samples.
bulk_matrix_adipo <- read_tpm('/nfs/data2/huji_artificial_meat/adipogenesis/adipogenesis_tpm.tsv')
bulk_meta_adipo <- read.csv("~/artifical_meat_deconv/meta_bulk_adipo.csv", header = T)
bulk_matrix_adipo <- bulk_matrix_adipo[, bulk_meta_adipo$sample_id]

bulk_matrix_myo <- read_tpm('/nfs/data2/huji_artificial_meat/myogenesis/nfcore_out/star_salmon/salmon.merged.gene_tpm.tsv')
bulk_meta_myo <- read.csv("/nfs/data2/huji_artificial_meat/myogenesis/sample_sheet_cleaned.csv", header = T)
bulk_matrix_myo <- bulk_matrix_myo[, bulk_meta_myo$sample_id]


omnideconv::build_model_scaden(single_cell_object = mat,
                               cell_type_annotations = col_meta$predicted.pheno,
                               bulk_gene_expression = bulk_matrix_adipo,
                               model_path = "sig_scaden",
                               verbose = T)
message("Signature building finished")
