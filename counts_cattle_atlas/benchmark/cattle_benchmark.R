library(SimBu)
library(Seurat)
library(tidyverse)

setwd("/nfs/home/students/t.tabanli/artifical_meat_deconv/counts_cattle_atlas")

# Cattle Atlas
cattle <- readRDS("sampled_cattle_atlas.RDS")
cattle$cell_type <- trimws(cattle$cell_type)
cattle$cell_type <- gsub("'", "", cattle$cell_type) # Henle's loop causes issues in column parse
cattle$cell_type <- gsub("/", "_", cattle$cell_type)
dim(cattle)

# Get counts and annotation
counts <- GetAssayData(cattle, assay = "RNA", layer = "counts")
colnames(counts)
rownames(counts)
annotation <- data.frame(
  "ID" = colnames(counts),
  "cell_type" = cattle$cell_type
)

# data set specification
ds <- SimBu::dataset(
  annotation = annotation,
  count_matrix = counts,
  name = "test_dataset"
)

# cell type proportions in original data
ct_by_tissue <- read_csv("tissue_by_ct.csv")
ct_numbers <- ct_by_tissue %>%
  select(-tissue) %>%
  colSums()
names(ct_numbers) <- gsub("'", "", names(ct_numbers))
names(ct_numbers) <- gsub("/", "_", names(ct_numbers))

# keep only those cells that are in the sampled data (because of 500 threshold)
ct_numbers <- ct_numbers[names(ct_numbers) %in% annotation$cell_type]
real_props <- ct_numbers / sum(ct_numbers)
custom_prop <- data.frame(as.list(real_props), row.names = "sample1", check.names = FALSE)


# simulate bulkSeq data
simulation <- SimBu::simulate_bulk(
  data = ds,
  scenario = "custom",
  custom_scenario_data = custom_prop,
  scaling_factor = "NONE",
  ncells = 100,
  nsamples = 1,
  balance_even_mirror_scenario = 0,
  BPPARAM = BiocParallel::MulticoreParam(workers = 20),
  run_parallel = TRUE
)
