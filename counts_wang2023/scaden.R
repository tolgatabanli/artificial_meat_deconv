library(Seurat)
library(omnideconv)
library(tidyverse)
library(parallel)
source("~/artifical_meat_deconv/toolkit.R")

setwd("/nfs/home/students/t.tabanli/artifical_meat_deconv/counts_wang2023")

message("Reading cattle...")
wang <- readRDS("wang_ids_mapped.RDS")
wang$lables <- trimws(wang$lables)
wang$lables <- gsub("'", "", wang$lables)
wang$lables <- gsub("/", "_", wang$lables)
wang$lables[wang$lables == "Muc-2"] <- "MuC-2"
dim(wang)

mat <- as.matrix(GetAssayData(wang, assay = "RNA", layer = "counts"))

message("Reading tpms...")

# Read TPM matrices and keep only relevant samples.
bulk_matrix_adipo <- read_tpm("/nfs/data2/huji_artificial_meat/adipogenesis/adipogenesis_tpm.tsv")
bulk_meta_adipo <- read.csv("/nfs/data2/huji_artificial_meat/adipogenesis/sample_sheet_cleaned.csv", header = T)
bulk_matrix_adipo <- bulk_matrix_adipo[, bulk_meta_adipo$sample_id]

bulk_matrix_myo <- read_tpm("/nfs/data2/huji_artificial_meat/myogenesis/nfcore_out/star_salmon/salmon.merged.gene_tpm.tsv")
colnames(bulk_matrix_myo) <- gsub(".*_", "", colnames(bulk_matrix_myo))
bulk_meta_myo <- read.csv("/nfs/data2/huji_artificial_meat/myogenesis/sample_sheet_minimal.csv", header = T) %>%
  mutate(sample_id = gsub(".*_", "", sample_id))
bulk_matrix_myo <- bulk_matrix_myo[, colnames(bulk_matrix_myo) %in% bulk_meta_myo$sample_id]

cat("Starting training...")
omnideconv::build_model_scaden(single_cell_object = mat,
                               cell_type_annotations = wang$lables,
                               bulk_gene_expression = bulk_matrix_adipo,
                               model_path = "sig_scaden",
                               verbose = T)

cat("Signature building finished")


comp_adipo <- omnideconv::deconvolute_scaden(signature = "sig_scaden",
                                             bulk_gene_expression = bulk_matrix_adipo,
                                             verbose = T)
comp_myo <- omnideconv::deconvolute_scaden(signature = "sig_scaden",
                                           bulk_gene_expression = bulk_matrix_myo,
                                           verbose = T)

saveRDS(comp_adipo, "sig_scaden/adipo_compositions_wang2023_scaden.Rds")
saveRDS(comp_myo, "sig_scaden/myo_compositions_wang2023_scaden.Rds")

cat("Deconvolution finished")

sessionInfo()
