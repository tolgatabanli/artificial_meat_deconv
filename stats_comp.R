library(compositions)
library(vegan)
library(tidyverse)
setwd("~/artifical_meat_deconv/")
source("toolkit.R")

# Reading routine for composition matrices and sample metadata
source("read_compositions_and_metadata.R")


# Stats (ref x method)
library(compositions)
library(zCompositions)
library(ComplexHeatmap)

exps <- c("adipo", "myo")
refs <- c("cattle", "messmer", "wang")
methods <- c("dwls", "rectangle", "scaden")

# PERLOG inferential analysis
run_perlog <- function(exp, ref, method) {
  cat("Processing:", exp, ref, method)
  # 1. choose matrix to analyze
  mat <- new_mats[[exp]][[ref]][[method]]

  meta <- if (exp == "adipo") raw_meta_adipo else raw_meta_myo
  vars_of_interest <- if (exp == "adipo") c("day", "individual") else c("day", "medium")

  meta <- as.data.frame(meta)
  rownames(meta) <- meta$sample_id
  mat <- mat[rownames(meta), ]
  
  # 2. Heatmap of raw proportions
  ann <- meta[, vars_of_interest] %>%
    mutate(across(all_of(vars_of_interest), as.factor))
  pal_list <- lapply(vars_of_interest, function(v) {
    vals <- levels(ann[[v]])
    n <- max(length(vals), 3)
    setNames(RColorBrewer::brewer.pal(n, "Accent")[seq_along(vals)], vals)
  })
  names(pal_list) <- vars_of_interest
  pheatmap::pheatmap(mat, annotation_row = ann, annotation_colors = pal_list,
                     filename = file.path("figures_heatmapRaw", paste0(paste(exp, ref, method, sep = "_"), ".png")))
  cat("\tHeatmap of raw fracs done.\n")
  
  # 3. Part elimination and zero-replacement:
  # Keep only parts that are present in at least 80% of the samples
  # Replace remaining zeros by multiplicative replacement
  keepParts <- colMeans(mat == 0) < 0.2
  mat_processed <- mat[, keepParts]
  cat("\tParts kept: ", colnames(mat_processed), "\n")
  
  if (any(mat_processed == 0)) {
    mat_processed <- zCompositions::multRepl(mat_processed, label = 0, dl = rep(1e-5, ncol(mat_processed)))
  }
  
  mat_processed <- acomp(mat_processed)
  rownames(mat_processed) <- rownames(mat)
  
  cell_types <- colnames(mat_processed)
  
  # 6. PERLOG analysis
  # day - individual for adipogenesis, day x medium for myogenesis
  # For adipo: stratify by individual, test by day
  # For myo: stratify by medium, test by day
  groups <- apply(meta[, vars_of_interest], 1, paste, collapse = "_")
  capture.output(res <- perLog(mat_processed,
                          groups = groups,
                          p = "auto", alpha = 0.05, R = 1000,
                          posthoc.g = TRUE, posthoc.lr = TRUE, mAdj = "bonferroni"),
                        file = nullfile())
  cat("\tPERLOG done.\n")
  p_overall <- res$pvalOv
  cat("\t\tOverall p-val:", p_overall, "\n")
  
  # 7. which groups are different due to which logratio pair? -> posthocs
  group_contrast_results <- res$posthoc.groups$summary %>% as.data.frame()
  group_contrast_results[, 'significant'] <- as.factor(if_else(group_contrast_results[,"Adj p-value"] <= 0.05, 1, 0))
  logratio_significances <- res$posthoc.logratios$significance
  
  posthoc_frame <- as.data.frame(t(combn(cell_types, 2))) %>%
    cross_join(data.frame(groupContrast = rownames(group_contrast_results)))
  colnames(posthoc_frame) <- c("ct1", "ct2", "groupContrast")
  posthoc_frame$lr_significance <- 0
  
  # build the data frame for logratio-pairs x groupContrasts
  mat_num <- unclass(mat_processed)
  for (groupContrast in rownames(group_contrast_results)) {
    if (group_contrast_results[groupContrast, "Adj p-value"] > 0.05) next
    rows <- posthoc_frame$groupContrast == groupContrast
    ct1  <- posthoc_frame$ct1[rows]
    ct2  <- posthoc_frame$ct2[rows]
    
    col_idx <- colnames(mat_processed)
    i1 <- match(ct1, col_idx)
    i2 <- match(ct2, col_idx)
    
    lo <- ifelse(i1 < i2, ct1, ct2)
    hi <- ifelse(i1 < i2, ct2, ct1)
    key <- paste0(lo, "/", hi)
    flip <- ifelse(ct1 == lo, 1, -1) # names are sorted alphabetically in perlog
    
    dissimilarity <- res$disEl[groupContrast, key]
    signum <- flip * sign(res$wts[groupContrast, key])
    
    # effect size: ct1/ct2 ratio between the two groups
    group_k <- strsplit(groupContrast, " vs ")[[1]][1]
    group_l <- strsplit(groupContrast, " vs ")[[1]][2]
    
    mean_fold_change <- mapply(function(type_a, type_b) {
      lr      <- log(mat_num[, type_a] / mat_num[, type_b])
      theta_k <- mean(lr[groups == group_k])
      theta_l <- mean(lr[groups == group_l])
      exp(theta_k - theta_l)
    }, ct1, ct2)
    
    # since some composition values are replaced with a very tiny value in zero-replacement
    # some cell-types can have exactly the same values in all the samples being compared as part of a group:
    # This makes Welch's statistic undefined -> NA values
    mask_significant <- mapply(function(a, b) logratio_significances[[groupContrast]][a, b], ct1, ct2)
    posthoc_frame$lr_significance[rows] <- ifelse(!is.na(mask_significant) & mask_significant != 0, dissimilarity * signum, 0)
    posthoc_frame$log2_mean_folg_change[rows] <- ifelse(!is.na(mask_significant) & mask_significant != 0, log2(mean_fold_change), NA_real_)
  }

  # Heatmap of logratios x groupcontrast
  var1 <- vars_of_interest[1]
  var2 <- vars_of_interest[2]
  gc_anno <- group_contrast_results %>%
    rownames_to_column("groupContrast") %>%
    separate_wider_delim(groupContrast, delim = " vs ", names = c("group1", "group2")) %>%
    separate_wider_delim(group1, delim = "_", names = c(paste0(var1, "1"), paste0(var2, "1"))) %>%
    separate_wider_delim(group2, delim = "_", names = c(paste0(var1, "2"), paste0(var2, "2"))) %>%
    dplyr::select(matches(var1), matches(var2), significant) %>%
    as.data.frame()
  rownames(gc_anno) <- rownames(group_contrast_results)
  gc_anno <- filter(gc_anno, significant == 1) %>% dplyr::select(-significant)
  
  # adjust ComplexHeatmap column annotation
  vals1 <- levels(droplevels(meta[[var1]]))
  vals2 <- unique(meta[,var2])
  
  cols1 <- setNames(RColorBrewer::brewer.pal(max(length(vals1), 3), "Oranges")[seq_along(vals1)], vals1)
  cols2 <- setNames(RColorBrewer::brewer.pal(max(length(vals2), 3), "Accent")[seq_along(vals2)], vals2)
  
  col_list <- setNames(list(cols1, cols1, cols2, cols2),
                       c(paste0(var1, "1"), paste0(var1, "2"), paste0(var2, "1"), paste0(var2, "2")))
  ha <- HeatmapAnnotation(df = gc_anno, col = col_list,
                          annotation_legend_param = setNames(
                            list(list(title = var1, at = vals1),
                                 list(title = var2)),
                            c(paste0(var1, "2"),
                              paste0(var2, "1"))
                          ), # maps the name of the variable to the df and col
                          show_legend = c(F, T, T, F))
  
  # heatmap matrix
  posthoc_frame <- posthoc_frame %>%
    tidyr::unite("pair", ct1, ct2, sep = "/") %>%
    pivot_wider(id_cols = "pair", names_from = groupContrast, values_from = log2_mean_folg_change) %>%
    column_to_rownames("pair") %>%
    as.matrix()
  p_mat <- apply(posthoc_frame, 2, as.numeric)
  rownames(p_mat) <- rownames(posthoc_frame)
  
  png(file.path("figures_posthoc", paste0(paste(exp, ref, method, "p_val", signif(p_overall, 3), sep = "_"), ".png")),
      width = 720, height=720, res = 80)
  draw(
    Heatmap(p_mat[, rownames(gc_anno)], show_column_names = T, name = "log2(Mean Fold Change)",
            cluster_columns=FALSE, cluster_rows=FALSE,
            top_annotation=ha, col=circlize::colorRamp2(c(-max(abs(p_mat), na.rm=TRUE), 0, max(abs(p_mat), na.rm=TRUE)),
                                              c("darkred", "white", "darkblue")),
            show_heatmap_legend=T, rect_gp = gpar(col = "gray50", lwd = 0.5)),
    merge_legends=TRUE
  )
  dev.off()
  cat("\t\tPosthoc heatmap done.\n")
}

for (exp in exps) {
  for (ref in refs) {
    for (method in methods) {
      run_perlog(exp, ref, method)
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


