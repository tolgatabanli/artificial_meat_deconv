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
refs <- c("cattle", "messmer", "wang", "messmer_full")
methods <- c("dwls", "rectangle", "scaden")

# PERLOG inferential analysis
run_perlog <- function(exp, ref, method) {
  cat("Processing:", exp, ref, method)
  # 1. choose matrix to analyze
  mat <- comp_mats[[exp]][[ref]][[method]]

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
  # For myo: compare days in scaffold P and compare scaffolds at day 0 and 21
  perlog_results <- data.frame()
  overall_results <- data.frame()
  
  if (exp == "adipo") {
    strata_props <- data.frame(filter_var = c("individual", "individual"),
                               filter_val = c("F", "M"),
                               group_var = c("day", "day"))
  } else {
    strata_props <- data.frame(filter_var = c("medium", "day", "day"),
                               filter_val = c("P", "D0", "D21"),
                               group_var = c("day", "medium", "medium"))
  }
  
  # for each individual or scaffold:
  for (i in 1:nrow(strata_props)) {
    prop <- strata_props[i, ]
    
    idx <- meta[[prop$filter_var]] == prop$filter_val
    mat_s <- mat_processed[idx, , drop = FALSE]
    meta_s <- meta[idx, ]
    groups_s <- as.character(droplevels(as.factor(meta_s[[prop$group_var]])))
    
    stratum_label <- paste0(prop$filter_var, "=", prop$filter_val)
    
    res_s <- tryCatch({
      capture.output(
        out <- perLog(mat_s, groups = groups_s,
                      p = "auto", alpha = 0.05, R = 1000,
                      posthoc.g = TRUE, posthoc.lr = TRUE, mAdj = "BH"),
        file = nullfile()
      )
      out
    }, error = function(e) {
      cat("\tPERLOG non-significant for", stratum_label, ":", conditionMessage(e), "\n")
      NULL
    })
    
    if (is.null(res_s)) {
      overall_results <- rbind(overall_results,
                               data.frame(stratum = stratum_label, test_var = prop$group_var,
                                          pvalOv = NA_real_, significant_overall = FALSE))
      next
    }
    
    overall_results <- rbind(overall_results,
                             data.frame(stratum = stratum_label, test_var = prop$group_var,
                                        pvalOv = res_s$pvalOv, significant_overall = res_s$pvalOv <= 0.05))
    
    cat("\tPERLOG done for", stratum_label, "(comparing", prop$group_var, ")\n")
    
    mat_num <- unclass(mat_s)
    
    # Save disEl + wts matrices in long format
    disEl_df <- as.data.frame(as.table(res_s$disEl))
    colnames(disEl_df) <- c("groupContrast", "pair", "disEl")
    
    wts_df <- as.data.frame(as.table(res_s$wts))
    colnames(wts_df) <- c("groupContrast", "pair", "wts")
    
    df <- merge(disEl_df, wts_df, by = c("groupContrast", "pair"))
    df$groupContrast <- as.character(df$groupContrast)
    df$pair          <- as.character(df$pair)
    
    split_pair <- strsplit(df$pair, "/")
    df$ct1 <- vapply(split_pair, `[`, character(1), 1)
    df$ct2 <- vapply(split_pair, `[`, character(1), 2)
    
    # Group contrast adjusted p-values
    group_summary  <- as.data.frame(res_s$posthoc.groups$summary)
    df$group_pval  <- group_summary[df$groupContrast, "Adj p-value"]
    
    # logratio-level significance (0/1/-1)
    df$lr_sig <- mapply(function(gc, a, b) {
      sigmat <- res_s$posthoc.logratios$significance[[gc]]
      if (is.null(sigmat)) return(NA_real_)
      sigmat[a, b]
    }, df$groupContrast, df$ct1, df$ct2)
    
    # fold change of ct1/ct2 ratio between the two groups in this contrast
    df$fold_change <- mapply(function(gc, a, b) {
      g_before <- strsplit(gc, " vs ")[[1]][1]
      g_after <- strsplit(gc, " vs ")[[1]][2]
      lr <- log(mat_num[, a] / mat_num[, b])
      exp(mean(lr[groups_s == g_after]) - mean(lr[groups_s == g_before]))
    }, df$groupContrast, df$ct1, df$ct2)
    
    df$stratum  <- stratum_label
    df$test_var <- prop$group_var
    
    perlog_results <- rbind(perlog_results, df)
  }
  
  # postprocessing 
  perlog_results <- perlog_results %>%
    mutate(
      significant = !is.na(lr_sig) & lr_sig != 0 & group_pval <= 0.05,
      fold_change = ifelse(significant, fold_change, NA_real_),
      log2fc      = log2(fold_change)
    )
  
  # Plot the comparisons
  # show also insignificant result but empty
  #all_strata <- unique(overall_results$stratum)
  #perlog_results$stratum <- factor(perlog_results$stratum, levels = all_strata)
  # insig_df <- overall_results %>%
  #   filter(!significant_overall) %>%
  #   mutate(stratum = factor(stratum, levels = all_strata))
  # following is to be able to control the panel width, by creating made-up values
  # one real pair is required to avoid an empty y-value in the common y-axis 
  # insig_anchor <- insig_df %>%
  #   mutate(groupContrast = "", pair = unique(perlog_results$pair)[1])
  
  ggplot() +
    #geom_blank(data = insig_anchor, aes(x = groupContrast, y = pair)) +
    #(if (nrow(insig_df) > 0) geom_rect(data = insig_df, aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf), fill = "gray95")) +
    geom_tile(data = perlog_results, aes(x = groupContrast, y = pair, fill = log2fc),
              color = "gray40", linewidth = 0.5) +
    facet_wrap(~ stratum, scales = "free_x", space = "free_x", drop = FALSE) +
    scale_fill_gradient2(low = "firebrick", mid = "white", high = "steelblue",
                         midpoint = 0, na.value = "grey60", name = "Log2FC") +
    theme_bw() +
    labs(x = "Group Contrast", y = "Cell-Type Pair") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  pval_summary <- paste(overall_results$stratum, signif(overall_results$pvalOv, 3),
                        sep = "=", collapse = "_")
  ggsave(filename = file.path("figures_posthoc", paste0(paste(exp, ref, method, pval_summary, sep = "_"), ".png")),
    width = 9, height = 9, units = "in", dpi = 80)
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


