library(tidyverse)

# Cattle cell-type aggregation
if (!file.exists("converter_lineage.csv")) {
  # Read cattle, just for the cell lineage aggregation
  cattle <- readRDS("counts_cattle_atlas/sampled_cattle_atlas.RDS")
  cattle$cell_type <- trimws(cattle$cell_type)
  cattle$cell_type <- gsub("'", "", cattle$cell_type)
  cattle$cell_type <- gsub("/", "_", cattle$cell_type)
  
  # Aggregate fine cattle cell types to lineages
  converter_lineage <- unique(data.frame(
    original = c(cattle$cell_type, "Unknown"),
    new_name = c(cattle$cell_lineage, "Unknown")
  ))
  
  write.csv(converter_lineage, "converter_lineage.csv", row.names = FALSE)
} else {
  converter_lineage <- read.csv("converter_lineage.csv")
}

# Wang cell-type aggregation
wang_mapping <- c(
  "FAP-1" = "FAP", "FAP-2" = "FAP",
  "Granulocyte-1" = "Granulocyte", "Granulocyte-2" = "Granulocyte",
  "SC" = "SC",
  "MuC-1" = "MuC", "MuC-2" = "MuC",
  "LEndoC" = "EndoC", "VEndoC" = "EndoC", "Non-VEndoC" = "EndoC",
  "B cell" = "Lymphoid", "CD3 T cell" = "Lymphoid", "CD4 T cell" = "Lymphoid", "CD8 T cell" = "Lymphoid", "Proliferating T cell" = "Lymphoid", "Natural killer cell" = "Lymphoid",
  "Slow myofiber" = "Myofiber", "Fast myofiber" = "Myofiber",
  "Mo_Ma" = "Myeloid", "Tenocyte" = "Myeloid", "Neutrophil" = "Myeloid"
)
converter_wang <- unique(data.frame(
  original = c(names(wang_mapping), "Unknown"),
  new_name = c(unname(wang_mapping), "Unknown")
))

agg_to_lineage <- function(mat, converter) {
  mat <- as.matrix(mat)
  index_to_lineage <- match(colnames(mat), converter$original)
  colnames(mat) <- converter$new_name[index_to_lineage]
  t(rowsum(t(mat), group = colnames(mat)))
}

read_composition <- function(file) {
  mat <- if (endsWith(file, ".csv")) read_csv(file) else readRDS(file)
  mat <- as.data.frame(mat)
  mat[mat < 0] <- 0
  if ("index" %in% colnames(mat)) mat <- column_to_rownames(mat, "index")
  as.matrix(mat)
}

read_sample_meta <- function(file) {
  read_csv(file, col_names = TRUE) %>%
    dplyr::select(-dplyr::matches("\\.\\.\\.[0-9]+")) %>% # drops unnecessary cols (...raw.file etc.)
    mutate(day = factor(day, levels = unique(day)[order(as.numeric(gsub("[^0-9]", "", unique(day))))])) %>%
    as.data.frame() %>% janitor::clean_names()
}

# Sample metadata
raw_meta_adipo <- read_sample_meta("/nfs/data2/huji_artificial_meat/adipogenesis/sample_sheet_cleaned.csv") %>%
  dplyr::select(sample_id, day, individual, group, bio_replicate)
raw_meta_myo <- read_sample_meta("/nfs/data2/huji_artificial_meat/myogenesis/sample_sheet_cleaned.csv") %>%
  dplyr::select(sample_id, day, group, individual, bio_replicate) %>%
  filter(!grepl("Satt", sample_id)) %>%
  separate_wider_delim(group, delim = "_", names = c(NA, "medium")) %>%
  mutate(sample_id = gsub(".*_", "", sample_id))

# composition matrices: experiment -> reference -> method
comp_files <- list(
  adipo = list(
    cattle = list(
      dwls      = "counts_cattle_atlas/sig_dwls_fine/adipo_compositions_cattleatlas_dwls_signature.Rds",
      rectangle = "counts_cattle_atlas/sig_rectangle_fine/rectangle_estimations_adipo.csv",
      scaden    = "counts_cattle_atlas/sig_scaden_fine/adipo_compositions_cattleatlas_scaden.Rds"
    ),
    messmer = list(
      dwls      = "messmer_counts/sig_dwls/adipo_compositions_messmer2023_dwls_signature.Rds",
      rectangle = "messmer_counts/sig_rectangle/rectangle_estimations_adipo.csv",
      scaden    = "messmer_counts/sig_scaden/adipo_compositions_messmer2023_scaden.Rds"
    ),
    wang = list(
      dwls      = "counts_wang2023/sig_dwls/adipo_compositions_wang2023_dwls_signature.Rds",
      rectangle = "counts_wang2023/sig_rectangle/rectangle_estimations_adipo.csv",
      scaden    = "counts_wang2023/sig_scaden/adipo_compositions_wang2023_scaden.Rds"
    )
  ),
  myo = list(
    cattle = list(
      dwls      = "counts_cattle_atlas/sig_dwls_fine/myo_compositions_cattleatlas_dwls_signature.Rds",
      rectangle = "counts_cattle_atlas/sig_rectangle_fine/rectangle_estimations_myo.csv",
      scaden    = "counts_cattle_atlas/sig_scaden_fine/myo_compositions_cattleatlas_scaden.Rds"
    ),
    messmer = list(
      dwls      = "messmer_counts/sig_dwls/myo_compositions_messmer2023_dwls_signature.Rds",
      rectangle = "messmer_counts/sig_rectangle/rectangle_estimations_myo.csv",
      scaden    = "messmer_counts/sig_scaden/myo_compositions_messmer2023_scaden.Rds"
    ),
    wang = list(
      dwls      = "counts_wang2023/sig_dwls/myo_compositions_wang2023_dwls_signature.Rds",
      rectangle = "counts_wang2023/sig_rectangle/rectangle_estimations_myo.csv",
      scaden    = "counts_wang2023/sig_scaden/myo_compositions_wang2023_scaden.Rds"
    )
  )
)

comp_mats <- rapply(comp_files, read_composition, how = "replace")

# fix myo meta sample ids and keep only relevant ones
for (ref in names(comp_mats[["myo"]])) {
  for (method in c("dwls", "rectangle", "scaden")) {
    myo_sample_ids <- gsub(".*_", "", rownames(comp_mats[["myo"]][[ref]][[method]]))
    rownames(comp_mats[["myo"]][[ref]][[method]]) <- myo_sample_ids
    comp_mats[["myo"]][[ref]][[method]] <- comp_mats[["myo"]][[ref]][[method]][myo_sample_ids %in% raw_meta_myo$sample_id,]
  }
}

# make consistent names in Messmer
fix_colnames <- function(mat) {
  renames <- c("Glial-cells" = "Glial cells", "Monocytes - M?" = "Monocytes-M\u03A6", "B-T-NK-cells" = "B-T-NK cells")
  colnames(mat) <- dplyr::recode(colnames(mat), !!!renames)
  mat
}

# fix messmer cell type names
for (exp in names(comp_mats)) comp_mats[[exp]]$messmer <- lapply(comp_mats[[exp]]$messmer, fix_colnames)

# aggregate cattle cells to lineages
for (exp in names(comp_mats)) comp_mats[[exp]]$cattle <- lapply(comp_mats[[exp]]$cattle, agg_to_lineage, converter = converter_lineage)
for (exp in names(comp_mats)) comp_mats[[exp]]$wang <- lapply(comp_mats[[exp]]$wang, agg_to_lineage, converter = converter_wang)