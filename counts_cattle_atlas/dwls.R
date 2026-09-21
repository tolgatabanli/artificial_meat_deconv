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
cattle$cell_type <- gsub("/", "_", cattle$cell_type)
dim(cattle)

### patch

###

mat <- as.matrix(GetAssayData(cattle, assay = "RNA", layer = "counts"))


message("Starting signature build...")
omnideconv::build_model_dwls(single_cell_object = mat,
                             cell_type_annotations = cattle$cell_type,
                             dwls_method = "mast_optimized",
                             path = "sig_dwls_fine",
                             verbose = T, ncores = 5) # TODO: defaults: diff_cutoff = 0.5, pval_cutoff = 0.01
message("Signature building finished")


load("sig_dwls_fine/Sig.RData")

message("Reading tpms...")

# Read TPM matrices and keep only relevant samples.
bulk_matrix_adipo <- read_tpm("/nfs/data2/huji_artificial_meat/adipogenesis/adipogenesis_tpm.tsv")
bulk_meta_adipo <- read.csv("/nfs/data2/huji_artificial_meat/adipogenesis/sample_sheet_cleaned.csv", header = T)
bulk_matrix_adipo <- bulk_matrix_adipo[, bulk_meta_adipo$sample_id]

bulk_matrix_myo <- read_tpm("/nfs/data2/huji_artificial_meat/myogenesis/nfcore_out/star_salmon/salmon.merged.gene_tpm.tsv")
colnames(bulk_matrix_myo) <- gsub(".*_", "", colnames(bulk_matrix_myo))
bulk_meta_myo <- read.csv("/nfs/data2/huji_artificial_meat/myogenesis/sample_sheet_minimal.csv", header = T)
bulk_matrix_myo <- bulk_matrix_myo[, gsub(".*_", "", colnames(bulk_matrix_myo)) %in% gsub(".*_", "", bulk_meta_myo$sample_id)]

length(intersect(rownames(Sig), rownames(bulk_matrix_adipo)))
length(intersect(rownames(Sig), rownames(bulk_matrix_myo)))


message("Deconvoluting...")

# deconvolute_dwls_parallel is in toolkit.R
message("___adipo___")
comp_adipo <- deconvolute_dwls_parallel(bulk_matrix_adipo, Sig, ncores = 38)
message("___myo___")
comp_myo   <- deconvolute_dwls_parallel(bulk_matrix_myo,   Sig, ncores = 13)


saveRDS(comp_adipo, "sig_dwls_fine/adipo_compositions_cattleatlas_dwls_signature.Rds")
saveRDS(comp_myo, "sig_dwls_fine/myo_compositions_cattleatlas_dwls_signature.Rds")
