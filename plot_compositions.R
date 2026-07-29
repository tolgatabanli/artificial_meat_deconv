library(tidyverse)
setwd("~/artifical_meat_deconv")
source("toolkit.R")

read_composition <- function(file) {
  if (endsWith(file, ".csv")) file <- read_csv(file)
  else if (endsWith(file, ".Rds")) file <- readRDS(file)
  
  mat <- as.data.frame(file)
  if ("index" %in% colnames(mat)) mat <- column_to_rownames(mat, "index")
  
  renames <- c("Glial cells" = "Glial-cells", "Monocytes" = "Monocytes - M?", "B-T-NK cells" = "B-T-NK-cells")
  renames <- renames[renames %in% colnames(mat)]
  
  return(dplyr::rename(mat, any_of(renames)))
}

# sample meta data
bulk_meta_adipo <- read_csv('/nfs/data2/huji_artificial_meat/adipogenesis/sample_sheet_cleaned.csv', col_names = T) %>%
  janitor::clean_names()
bulk_meta_myo <- read_csv('/nfs/data2/huji_artificial_meat/myogenesis/sample_sheet_cleaned.csv', col_names = T) %>%
  janitor::clean_names()

bulk_meta_adipo <- bulk_meta_adipo %>%
  dplyr::select(library_name, sample_id, title, group, bio_replicate, individual, treatment, day) %>%
  mutate(day = factor(day, levels = c("D0", "D2", "D3", "D6", "D14")))
bulk_meta_myo <- bulk_meta_myo %>%
  dplyr::select(library_name, sample_id, group, day, bio_replicate, treatment) %>%
  mutate(day = factor(day, levels = c("d0", "d3", "d7", "d21")))

comp_adipo_rectangle <- read_composition("messmer_counts/rectangle_results/rectangle_estimations_adipo.csv")
comp_adipo_dwls <- read_composition('~/artifical_meat_deconv/messmer_counts/dwls_results/no_teno/adipo_compositions_messmer2023_dwls_signature.Rds')
comp_adipo_scaden <- read_composition("~/artifical_meat_deconv/messmer_counts/sig_scaden/adipo_compositions_messmer2023_scaden.Rds")

comp_myo_rectangle <- read_composition("messmer_counts/rectangle_results/rectangle_estimations_myo.csv")
comp_myo_dwls <- read_composition('~/artifical_meat_deconv/messmer_counts/dwls_results/no_teno/myo_compositions_messmer2023_dwls_signature.Rds')
comp_myo_scaden <- read_composition("~/artifical_meat_deconv/messmer_counts/sig_scaden/myo_compositions_messmer2023_scaden.Rds")
# dwls rownames are canonical (come from actual bulk matrix)
sxx_correct <- setNames(rownames(comp_myo_dwls), gsub(".*_", "", rownames(comp_myo_dwls)))
# fix meta
bulk_meta_myo$sample_id <- sxx_correct[gsub(".*_", "", bulk_meta_myo$sample_id)]
# fix scaden (currently has meta's wrong names)
rownames(comp_myo_scaden) <- sxx_correct[gsub(".*_", "", rownames(comp_myo_scaden))]

# prepare rectangle
adipo_merged <- merge_methods(comp_adipo_rectangle, comp_adipo_dwls, comp_adipo_scaden, meta = bulk_meta_adipo, method_names = c("rectangle", "dwls", "scaden"))
colors <- pals::alphabet(n_distinct(adipo_merged$cell_type))
names(colors) <- unique(adipo_merged$cell_type)
colors["Unknown"] <- "gray70"

p <- adipo_merged %>%
  ggplot(aes(x = sample_id, y = proportion, fill = cell_type)) +
  geom_col(linewidth = 0.8) +
  coord_flip() +
  facet_grid(day + individual ~ method, scales = "free_y", space = "free_y") +
  scale_fill_manual(values = colors, name = "Cell type") +
  theme_minimal() +
  labs(y = "Proportion", x = NULL) +
  theme(
    strip.text.y.right   = element_text(angle = 0, hjust = 0.5, face = "bold"),
    strip.background.y   = element_rect(fill = "gray90", color = "black", linewidth = 0.8),
    panel.border         = element_rect(color = "gray70", fill = NA, linewidth = 0.4),
    panel.spacing.y      = unit(0.3, "lines")
  )
p

myo_merged <- merge_methods(comp_myo_rectangle, comp_myo_dwls, comp_myo_scaden, meta = bulk_meta_myo, method_names = c("rectangle", "dwls", "scaden"))

sat_samples <- myo_merged %>% filter(grepl("Sat", group))
others <- myo_merged %>% filter(!grepl("Sat", group))

p <- others %>%
  ggplot(aes(x = sample_id, y = proportion, fill = cell_type)) +
  geom_col(linewidth = 0.8) +
  coord_flip() +
  facet_grid(day ~ method, scales = "free_y", space = "free_y") +
  scale_fill_manual(values = colors, name = "Cell type") +
  theme_minimal() +
  labs(y = "Proportion", x = NULL) +
  theme(
    strip.text.y.right   = element_text(angle = 0, hjust = 0.5, face = "bold"),
    strip.background.y   = element_rect(fill = "gray90", color = "black", linewidth = 0.8),
    panel.border         = element_rect(color = "gray70", fill = NA, linewidth = 0.4),
    panel.spacing.y      = unit(0.3, "lines")
  )
p
