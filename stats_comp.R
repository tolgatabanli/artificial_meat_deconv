library(compositions)
library(vegan)
library(ggrepel)
library(tidyverse)
setwd("~/artifical_meat_deconv/")
source("toolkit.R")

# Compare DWLS and Rectangle compositions of myogenesis
# If sample_meta is given, returns a list:
#   $mat | matrix result containing the metadata
#   $metadata | data frame of the sample_meta
#   $cell_types
read_composition <- function(file, sample_meta = NULL) {
  if (endsWith(file, ".csv")) file <- read_csv(file)
  else if (endsWith(file, ".Rds")) file <- readRDS(file)
  
  mat <- as.data.frame(file)
  if ("index" %in% colnames(mat)) mat <- column_to_rownames(mat, "index")
  
  renames <- c("Glial cells" = "Glial-cells", "Monocytes" = "Monocytes - M?", "B-T-NK cells" = "B-T-NK-cells")
  renames <- renames[renames %in% colnames(mat)]
  mat <- dplyr::rename(mat, any_of(renames))
  
  if (!is.null(sample_meta)) {
    cat("Joining with metadata")
    metadata <- read_csv(sample_meta, col_names = T) %>%
      select(-any_of("...1")) %>%
      # factorize day
      mutate(day = factor(day, levels = unique(day)[order(as.numeric(gsub("[^0-9]", "", unique(day))))]))
    cell_types <- colnames(mat)
    cat("Cell types: "); cat(cell_types, sep = ", ")
    cat("columns:", colnames(metadata), sep = ", ")
    
    mat <- mat %>% rownames_to_column("sample_id") %>%
      inner_join(metadata, by = "sample_id")
      #pivot_longer(all_of(cell_types), names_to="cell_type", values_to="proportion") %>% 
      #relocate(c("sample_id", "cell_type", "proportion"))
    
    return(mat)
  }
  else return(mat)
}

comp_adipo_rectangle <- read_composition("messmer_counts/rectangle_results/rectangle_estimations_adipo.csv", "/nfs/data2/huji_artificial_meat/adipogenesis/sample_sheet_cleaned.csv")
comp_adipo_dwls <- read_composition('~/artifical_meat_deconv/messmer_counts/dwls_results/no_teno/adipo_compositions_messmer2023_dwls_signature.Rds', "/nfs/data2/huji_artificial_meat/adipogenesis/sample_sheet_cleaned.csv")

comp_myo_rectangle <- read_composition("messmer_counts/rectangle_results/rectangle_estimations_myo.csv", "/nfs/data2/huji_artificial_meat/myogenesis/sample_sheet_cleaned.csv")
comp_myo_dwls <- read_composition('~/artifical_meat_deconv/messmer_counts/dwls_results/no_teno/myo_compositions_messmer2023_dwls_signature.Rds', "/nfs/data2/huji_artificial_meat/myogenesis/sample_sheet_cleaned.csv")

# remove satellite cell experiment
comp_myo_rectangle <- comp_myo_rectangle %>% filter(!grepl("Sat", sample_id))
comp_myo_dwls <- comp_myo_dwls %>% filter(!grepl("Sat", sample_id))

# PCA of the compositions
clr <- clr(acomp(comp_myo_rectangle, parts = 2:12))
rownames(clr) <- rownames(comp_myo_rectangle)
pca          <- prcomp(clr, scale. = FALSE)
var_exp      <- round(summary(pca)$importance[2, 1:2] * 100, 1)

comp_myo_rectangle %>%
  mutate(PC1 = pca$x[, 1],
         PC2 = pca$x[, 2]) %>%
  ggplot(aes(x = PC1, y = PC2, color = day)) +
  geom_point(size = 3, alpha = 0.85) +
  labs(x = paste0("PC1 (", var_exp[1], "%)"),
       y = paste0("PC2 (", var_exp[2], "%)")) +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave("images/pca_compositions_messmer_dwls.png", width = 8)

# take common cell types
cell_types <- colnames(comp_myo_rectangle)[2:11]
merged <- rbind(comp_myo_rectangle[, cell_types], comp_myo_dwls[, cell_types])
merged_meta <- bind_rows(comp_myo_rectangle, comp_myo_dwls) %>% select(where(~all(!is.na(.)))) %>%
  mutate(method = c(rep("rectangle", nrow(comp_myo_rectangle)), rep("dwls", nrow(comp_myo_dwls))))

clr <- clr(acomp(merged))
distance_matrix <- dist(clr)
adonis2(distance_matrix ~ day + method, data = merged_meta, by = "margin")

# test each day individually
days <- droplevels(unique(merged_meta$day))

# PERMANOVA per day
pairwise_day <- lapply(days, function(d) {
  idx  <- which(droplevels(merged_meta$day) == d)
  if (length(unique(merged_meta$method[idx])) < 2) {
    warning("Skipping ", d, " due to less than 2 methods.")
    return(NULL)
  }
  
  sub_dist <- as.dist(as.matrix(distance_matrix)[idx, idx])
  sub_meta <- merged_meta[idx, ]
  res <- adonis2(sub_dist ~ method, data = sub_meta, permutations = 999)
  
  data.frame(day = d, F = res[1,"F"], R2 = res[1,"R2"], p = res[1,"Pr(>F)"])
})

do.call(rbind, Filter(Negate(is.null), pairwise_day))

# wilcox per day
clr_mat <- as.matrix(clr(acomp(merged)))

pairwise_day <- lapply(days, function(d) {
  idx      <- which(droplevels(merged_meta$day) == d)
  
  sub_clr  <- clr_mat[idx, ]
  sub_meta <- merged_meta[idx, ]
  rect_idx <- sub_meta$method == "rectangle"
  
  pvals <- apply(sub_clr, 2, function(x) {
    wilcox.test(x[rect_idx], x[!rect_idx])$p.value
  })
  
  data.frame(day = d, cell_type = cell_types, p = pvals,
             p_adj = p.adjust(pvals, method = "BH"))
})

do.call(rbind, Filter(Negate(is.null), pairwise_day))

# test with bootstrap samples, compare methods for each day

# matrix has rows as samples
bootstrap_matrix <- function(mat, size) {
  row_indices <- sample(1:nrow(clr_mat), size = size, replace = TRUE)
  return(clr_mat[row_indices,])
}

merged_meta$day <- droplevels(merged_meta$day)
df <- data.frame()
for (d in merged_meta$day) {
  for (cell_type in cell_types) {
    for (method in unique(merged_meta$method)) {
      subsample <- clr[merged_meta$method == method & merged_meta$day == d, cell_type]
      a_row <- sample(subsample, size = 100, replace = T)
      sub_df <- data.frame(day = d, cell_type = cell_type, method = method, val = a_row)
      df <- bind_rows(df, sub_df)
    }
  }
  
  return(df)
}


ggplot(df, aes(x = method, y = val)) +
  geom_violin() +
  facet_wrap(~ cell_type)

########
########
########
# Read data
meta_myo <- read.csv("artifical_meat_deconv/meta_bulk_myo.csv")

messmer_myo_timepoints <- readRDS("/nfs/home/students/t.tabanli/artifical_meat_deconv/sig_dwls_messmer_recreated_counts_noTeno/myo_compositions_messmer2023_dwls_signature.Rds")
messmer_myo_timepoints[messmer_myo_timepoints < 0] <- 0
messmer_myo_timepoints <- messmer_myo_timepoints %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(cols = -sample_id, names_to = "cell_type", values_to = "proportion") %>%
  separate_wider_delim(cell_type, delim = "_", names = c("injection", "cell_type")) %>%
  mutate(sample_id = gsub(".*_", "", sample_id))

messmer_myo_all <- messmer_myo_timepoints %>%
  group_by(sample_id, cell_type) %>%
  summarise(proportion = sum(proportion), .groups = "drop") %>%
  pivot_wider(names_from = cell_type, values_from = proportion) %>%
  column_to_rownames("sample_id") %>%
  as.matrix()

messmer_myo_tp1 <- readRDS("/nfs/home/students/t.tabanli/artifical_meat_deconv/sig_dwls_messmer_recreated_counts_timepoint1_noTeno/myo_compositions_messmer2023_dwls_signature.Rds")
messmer_myo_tp1[messmer_myo_tp1 < 0] <- 0
rownames(messmer_myo_tp1) <- gsub(".*_", "", rownames(messmer_myo_tp1))
colnames(messmer_myo_tp1) <- gsub("Injection1_", "", colnames(messmer_myo_tp1))


messmer_myo_all <- messmer_myo_all[rownames(messmer_myo_all) %in% meta_myo$sample_id, ]
messmer_myo_tp1 <- messmer_myo_tp1[rownames(messmer_myo_tp1) %in% meta_myo$sample_id, ]

# PCA for All vs Injection 1 
stacked <- rbind(messmer_myo_all, messmer_myo_tp1)
stacked[stacked < 0] <- 0
source_label <- c(rep("all", nrow(comp_all)), rep("tp1", nrow(comp_tp1)))

clr_stacked  <- as.matrix(clr(acomp(stacked)))
pca          <- prcomp(clr_stacked, scale. = FALSE)
var_exp      <- round(summary(pca)$importance[2, 1:2] * 100, 1)

pca_df <- data.frame(
  PC1    = pca$x[, 1],
  PC2    = pca$x[, 2],
  source = source_label,
  sample = c(rownames(messmer_myo_all), rownames(messmer_myo_tp1))
)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = source, label = sample)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_text_repel(size = 2.5, show.legend = FALSE, max.overlaps = 20) +
  scale_color_manual(values = c("all" = "#2C6FAC", "tp1" = "#C06010"),
                     name = "Signature") +
  labs(title   = "Aitchison PCA - Compare Compositions of All vs Timepoint 1",
       x       = paste0("PC1 (", var_exp[1], "%)"),
       y       = paste0("PC2 (", var_exp[2], "%)")) +
  theme_minimal() +
  theme(legend.position = "bottom")
p_pca

# Compare all timepoints
clr_all  <- as.matrix(clr(acomp(messmer_myo_all)))
dist_all <- dist(clr_all, method = "euclidean")

pairs_to_compare <- list(
  c("Injection1", "Injection2"),
  c("Injection3", "Injection4"),
  c("Injection3", "Injection5"),
  c("Injection4", "Injection5")
)
pairwise_results <- lapply(pairs_to_compare, function(x) {
  inj1 <- messmer_myo_timepoints %>%
    filter(injection == x[1]) %>%
    pivot_wider(names_from = cell_type, values_from = proportion) %>%
    column_to_rownames("sample_id") %>%
    select(-injection) %>%
    as.matrix()
  
  inj2 <- messmer_myo_timepoints %>%
    filter(injection == x[2]) %>%
    pivot_wider(names_from = cell_type, values_from = proportion) %>%
    column_to_rownames("sample_id") %>%
    select(-injection) %>%
    as.matrix()
  
  common_ct <- intersect(colnames(inj1), colnames(inj2))
  together  <- rbind(inj1[, common_ct], inj2[, common_ct])
  together[together < 0] <- 0
  
  meta <- rbind(
    cbind(meta_myo[match(rownames(inj1), meta_myo$sample_id), ], injection = x[1]),
    cbind(meta_myo[match(rownames(inj2), meta_myo$sample_id), ], injection = x[2])
  )
  
  clr_together  <- as.matrix(clr(acomp(together)))
  dist_together <- dist(clr_together, method = "euclidean")
  res <- adonis2(dist_together ~ injection, data = meta, permutations = 999, by = "margin")
  
  list(inj1 = x[1], inj2 = x[2],
       ct_inj1   = colnames(inj1),
       ct_inj2   = colnames(inj2),
       common_ct = common_ct,
       n_shared  = length(common_ct),
       pval   = res[1, "Pr(>F)"],
       R2     = res[1, "R2"],
       F_stat = res[1, "F"])
})

# ── Panel 1: cell type presence/absence per injection ────────────────────────
presence_df <- messmer_myo_timepoints %>%
  distinct(injection, cell_type) %>%
  mutate(status = "present") %>%
  complete(injection, cell_type, fill = list(status = "absent"))

# mark which cell types were shared in ANY comparison
shared_any <- Reduce(union, lapply(pairwise_results, `[[`, "common_ct"))

presence_df <- presence_df %>%
  mutate(shared = cell_type %in% shared_any,
         label  = case_when(
           status == "present" & shared  ~ "present & compared",
           status == "present" & !shared ~ "present, never compared",
           TRUE                          ~ "absent"
         ))

p1 <- ggplot(presence_df, aes(x = injection, y = cell_type, fill = label)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_manual(values = c(
    "present & compared"     = "#2C6FAC",
    "present, never compared"= "#A8C8E8",
    "absent"                 = "#EEEEEE"),
    name = NULL) +
  labs(title = "Cell type presence per injection", x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        panel.grid   = element_blank())
p1
# ── Panel 2: pairwise PERMANOVA summary ──────────────────────────────────────
results_df <- do.call(rbind, lapply(pairwise_results, function(r) {
  data.frame(
    pair     = paste(r$inj1, "vs", r$inj2),
    inj1     = r$inj1,
    inj2     = r$inj2,
    n_shared = r$n_shared,
    pval     = r$pval,
    R2       = round(r$R2, 3),
    F_stat   = round(r$F_stat, 2),
    sig      = cut(r$pval, c(-Inf, 0.001, 0.01, 0.05, Inf),
                   labels = c("***", "**", "*", "n.s."))
  )
}))

p2 <- ggplot(results_df, aes(x = pair, y = R2)) +
  geom_col(aes(fill = pval), width = 0.6) +
  geom_text(aes(label = sprintf("p=%.3f %s\nn=%d ct", pval, sig, n_shared)),
            vjust = -0.3, size = 2.8, lineheight = 0.9) +
  geom_text(aes(y = R2 / 2, label = sprintf("F=%.2f", F_stat)),
            color = "white", size = 2.8, fontface = "bold") +
  scale_fill_gradientn(
    colors = c("#C0392B", "#E67E22", "#BDC3C7"),
    values = scales::rescale(c(0, 0.05, 1)),
    limits = c(0, 1), name = "p-value") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(title = "Pairwise PERMANOVA (Aitchison distance)",
       x = NULL, y = "R²") +
  theme_minimal() +
  theme(axis.text.x   = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank())
p2

# PERMANOVA + PERMDISP on messmer_myo_all
meta_sub <- meta_myo[match(rownames(messmer_myo_all), meta_myo$sample_id), ]

# treatment descriptions for reference
treatment_key <- meta_sub %>%
  distinct(group, day, treatment) %>%
  arrange(day, group)

# PCA
pca    <- prcomp(clr_all, scale. = FALSE)
var_exp <- round(summary(pca)$importance[2, 1:2] * 100, 1)
pca_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], sample_id = rownames(clr_all)) %>%
  left_join(meta_sub, by = "sample_id")

ggplot(pca_df, aes(PC1, PC2, color = day)) +
  geom_point(size = 2.5) + #stat_ellipse(linetype = "dashed") +
  labs(x = paste0("PC1 (", var_exp[1], "%)"), y = paste0("PC2 (", var_exp[2], "%)")) +
  theme_minimal()

ggplot(pca_df, aes(PC1, PC2, color = group)) +
  geom_point(size = 2.5) + #stat_ellipse(linetype = "dashed") +
  labs(x = paste0("PC1 (", var_exp[1], "%)"), y = paste0("PC2 (", var_exp[2], "%)")) +
  theme_minimal()

# PERMDISP
boxplot(disp_day, main = "PERMDISP — day",   ylab = "Distance to centroid")
boxplot(disp_grp, main = "PERMDISP — group", ylab = "Distance to centroid")


