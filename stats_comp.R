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
}
cat("IT RAN THROUGH")


### backup
# 5. make geometric mean barplots
cell_types <- colnames(mat_processed)
p <- clr(mat_processed) %>% as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  inner_join(meta, by = "sample_id") %>%
  pivot_longer(cols = all_of(cell_types), names_to = "cell_type", values_to = "clr") %>%
  group_by(across(all_of(c(vars_of_interest, "cell_type")))) %>%
  summarise(clr_mean = mean(clr), .groups = "drop") %>%
  ggplot(aes(x = day, y = clr_mean)) +
  geom_col() +
  facet_grid(as.formula(paste(vars_of_interest[2], "~ cell_type")), scales = "free_x") +
  theme_bw() +
  labs(title = "", x = "Cell Type", y = "Mean CLR")
ggsave(file.path("figures_geometricMean", paste0(paste(exp, ref, method, sep = "_"), ".png")), p)
cat("\tGeometric mean barplots done.\n")



names_of_groups <- rownames(group_contrast_results)
split_pairs <- strsplit(names_of_groups, split = " vs ")

external_groups <- unique(unlist(split_pairs))
group_signif_matrix <- matrix(0, nrow = length(external_groups), ncol = length(external_groups))
rownames(group_signif_matrix) <- external_groups
colnames(group_signif_matrix) <- external_groups
for (i in 1:nrow(group_contrast_results)) {
  if (group_contrast_results[i, "Adj p-value"] > 0.05) next
  pair <- split_pairs[[i]]
  group_signif_matrix[pair[1], pair[2]] <- 1 #group_contrast_results[i, "Adj p-value"]
  group_signif_matrix[pair[2], pair[1]] <- 1
}

# red is significant
pheatmap::pheatmap(group_signif_matrix, display_numbers = F, cluster_rows = F, cluster_cols = F, legend = F)
ggsave(file.path(fig_path, paste0("heatmap_group_contrast_", fig_suffix, ".png")))






# 8. which logratios drive the differences?
# considers significant group contrasts and make a boxplot
lapply(, function(x) pheatmap::pheatmap(x))

# 
plot.princomp.acomp(princomp.acomp(mat_processed))
pc <- princomp(mat_processed)

plot(pc,type="biplot")
plot(pc,choice=c(1,3),type="biplot")
plot(pc,type="loadings")
plot(pc,type="loadings",scale.sdev=-1) # Downward
plot(pc,type="relative",scale.sdev=NA) # The directions
plot(pc,type="relative",scale.sdev=1) # one sigma Upward 
plot(pc,type="relative",scale.sdev=-1) # one sigma Downward
biplot(pc)
screeplot(pc)
loadings(pc)

# james multivariate test is unnecessary, since it requires n < p and
# is parametric.
# james(y1, y2, a = 0.05, R = 999, graph = FALSE)

# Logcontrast regression,
# but not relevant :( since this is to predict a continuous response variable
# Compositional::lc.reg(y, x, z = NULL, xnew = NULL, znew = NULL)


