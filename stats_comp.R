library(compositions)
library(tidyverse)
setwd("~/artifical_meat_deconv/")
source("toolkit.R")
suppressMessages(conflicted::conflict_prefer_all("dplyr"))

# Reading routine for composition matrices and sample metadata
source("read_compositions_and_metadata.R")
new_mats <- rapply(comp_mats, function(m) { m[m < 1e-5] <- 0; m }, how ="replace")

# Stats (ref x method)
library(compositions)
library(zCompositions)

exps <- c("adipo", "myo")
refs <- c("cattle", "messmer", "wang")
methods <- c("dwls", "rectangle", "scaden")

# PERLOG inferential analysis
all_overall_results <- data.frame()
all_group_contrast_results <- data.frame(
  exp = character(0), ref = character(0), method = character(0),
  stratum = character(0), test_var = character(0),
  groupContrast = character(0), group_pval = numeric(0),
  significant = logical(0), contrast_label = character(0)
)


run_perlog <- function(exp, ref, method) {
  cat("Processing:", exp, ref, method)
  # 1. choose matrix to analyze
  mat <- new_mats[[exp]][[ref]][[method]]

  meta <- if (exp == "adipo") raw_meta_adipo else raw_meta_myo
  vars_of_interest <- if (exp == "adipo") c("day", "individual") else c("day", "scaffold")

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
  
  # 4. PERLOG analysis
  # day - individual for adipogenesis, day x scaffold for myogenesis
  # For adipo: stratify by individual, test by day
  # For myo: compare days in scaffold P and compare scaffolds at day 0 and 21
  perlog_results <- data.frame(
    groupContrast = character(0),
    pair          = character(0),
    disEl         = numeric(0),
    wts           = numeric(0),
    ct1           = character(0),
    ct2           = character(0),
    group_pval    = numeric(0),
    lr_sig        = numeric(0),
    fold_change   = numeric(0),
    stratum       = character(0),
    test_var      = character(0)
  )
  overall_results <- data.frame()
  
  if (exp == "adipo") {
    strata_props <- data.frame(filter_var = c("individual", "individual"),
                               filter_val = c("F", "M"),
                               group_var = c("day", "day"))
  } else {
    strata_props <- data.frame(filter_var = c("scaffold", "day", "day"),
                               filter_val = c("P", "D00", "D21"),
                               group_var = c("day", "scaffold", "scaffold"))
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
    
    # Save significant result
    make_ns_row <- function() {
      data.frame(exp = exp, ref = ref, method = method,
                 stratum = stratum_label, test_var = prop$group_var,
                 groupContrast = "NS", group_pval = NA_real_,
                 significant = FALSE, contrast_label = "NS")
    }
    
    if (is.null(res_s)) {
      overall_results <- rbind(overall_results,
                               data.frame(stratum = stratum_label, test_var = prop$group_var,
                                          pvalOv = NA_real_, significant_overall = FALSE))
      all_group_contrast_results <<- rbind(all_group_contrast_results,
                                           make_ns_row())
      next
    }
    
    overall_results <- rbind(overall_results,
                             data.frame(stratum = stratum_label, test_var = prop$group_var,
                                        pvalOv = res_s$pvalOv, significant_overall = res_s$pvalOv <= 0.05))
    
    # Save to global contrast result:
    group_summary <- as.data.frame(res_s$posthoc.groups$summary)
    group_summary$groupContrast <- rownames(group_summary)
    
    sig_contrasts <- group_summary %>% filter(`Adj p-value` <= 0.05)
    
    group_contrast_results <- if (nrow(sig_contrasts) == 0) {
      make_ns_row()
    } else {
      data.frame(
        exp = exp, ref = ref, method = method,
        stratum = stratum_label, test_var = prop$group_var,
        groupContrast = sig_contrasts$groupContrast,
        group_pval = sig_contrasts$`Adj p-value`,
        significant = TRUE,
        contrast_label = paste0(sub(".*=", "", stratum_label), ": ", sig_contrasts$groupContrast)
      )
    }
    all_group_contrast_results <<- rbind(all_group_contrast_results, group_contrast_results)
    
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
    group_summary <- as.data.frame(res_s$posthoc.groups$summary)
    df$group_pval <- group_summary[df$groupContrast, "Adj p-value"]
    
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
  
  # save overall results to the global var to tabularize later
  overall_results$exp    <- exp
  overall_results$ref    <- ref
  overall_results$method <- method
  all_overall_results <<- rbind(all_overall_results, overall_results)
  
  # postprocessing 
  perlog_results <- perlog_results %>%
    mutate(
      significant = !is.na(lr_sig) & lr_sig != 0 & group_pval <= 0.05,
      fold_change = ifelse(significant, fold_change, NA_real_),
      log2fc = log2(fold_change)
    )
  
  # Plot the log-ratio vs group contrast
  num_pairs <- length(cell_types)
  if (nrow(perlog_results) > 0) {
    ggplot() +
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
           width = 2700, height = 200 * (num_pairs) * (num_pairs-1) / 2, units = "px", dpi = 300)
  }
  cat("\t\tPosthoc ggplot done.\n")
}

set.seed(2026)
for (exp in exps) {
  for (ref in refs) {
    for (method in methods) {
      run_perlog(exp, ref, method)
    }
  }
}

# p adjust for group contrasts
all_group_contrast_results <- all_group_contrast_results %>%
  group_by(exp, ref, method) %>%
  mutate(group_pval_adj = p.adjust(group_pval, method = "BH")) %>%
  ungroup()

all_group_contrast_results$significant <- all_group_contrast_results$group_pval_adj < 0.05

# Write the formatted latex table for group contrast for all ref-method combis
library(kableExtra)

format_pval <- function(p) {
  dplyr::case_when(
    is.na(p)   ~ "NS",
    p < 0.001  ~ "***",
    p < 0.01   ~ "**",
    p < 0.05   ~ "*",
    TRUE       ~ "ns"
  )
}

ref_labels <- c(messmer = "Messmer2023", wang = "Wang2023", cattle = "CattleAtlas")
method_labels <- c(dwls = "DWLS", rectangle = "Rectangle", scaden = "Scaden")

# following function is modified from an AI-generated code
build_contrast_table <- function(df, exp_name) {
  exp_label <- if (exp_name == "adipo") "Adipogenesis" else "Myogenesis"
  
  df_exp <- df %>%
    filter(exp == exp_name) %>%
    mutate(short_stratum = sub(".*=", "", stratum),
           col_key       = paste(ref_labels[ref], method_labels[method], sep = "__"))
  
  col_combos <- expand.grid(method_label = unname(method_labels),
                            ref_label    = unname(ref_labels),
                            stringsAsFactors = FALSE) %>%
    mutate(col_key = paste(ref_label, method_label, sep = "__"))
  ck <- col_combos$col_key
  
  all_strata <- unique(df_exp$short_stratum)
  
  # combos whose overall stratum test failed entirely
  ns_combos <- df_exp %>%
    filter(groupContrast == "NS") %>%
    distinct(short_stratum, col_key)
  
  # one summary row per stratum: label = stratum name, NS under failed columns
  overall_df <- lapply(all_strata, function(s) {
    failed <- ns_combos$col_key[ns_combos$short_stratum == s]
    as.data.frame(c(list(short_stratum = s, groupContrast = s),
                    setNames(as.list(ifelse(ck %in% failed, "NS", "")), ck)),
                  stringsAsFactors = FALSE)
  }) %>% bind_rows()
  
  # real contrast rows; contrast not significant for a ref/method -> "NS"
  real_df <- df_exp %>%
    filter(groupContrast != "NS") %>%
    mutate(cell_value = format_pval(group_pval)) %>%
    distinct(short_stratum, groupContrast, col_key, .keep_all = TRUE) %>%
    pivot_wider(id_cols = c(short_stratum, groupContrast),
                names_from = col_key, values_from = cell_value,
                values_fill = "NS")
  real_df[setdiff(ck, colnames(real_df))] <- "NS"
  
  # blank failed columns in real rows (NS already shown in summary row)
  for (i in seq_len(nrow(ns_combos))) {
    real_df[real_df$short_stratum == ns_combos$short_stratum[i],
            ns_combos$col_key[i]] <- ""
  }
  
  wide <- bind_rows(overall_df, real_df) %>%
    select(short_stratum, groupContrast, all_of(ck)) %>%
    arrange(short_stratum, groupContrast != short_stratum, groupContrast)
  
  group_lengths    <- rle(wide$short_stratum)$lengths
  group_start_rows <- cumsum(group_lengths) - group_lengths + 1
  hline_rows       <- head(cumsum(group_lengths), -1)
  
  header_counts <- setNames(rep(length(method_labels), length(ref_labels)),
                            unname(ref_labels))
  
  kbl(select(wide, -short_stratum),
      format = "latex", booktabs = TRUE,
      col.names = c("\\shortstack{Group\\\\Contrast}",
                    paste0("\\rotatebox[origin=c]{90}{\\makebox[1.6cm][l]{",
                           col_combos$method_label, "}}")),
      escape = FALSE,
      caption = paste0("Significant PERLOG group contrasts | ", exp_label,
                       ". * $p<0.05$, ** $p<0.01$, *** $p<0.001$. ",
                       "NS: not significant. p-values are corrected for multiple testing experiment-wise using the Benjamini-Hochberg procedure.")) %>%
    add_header_above(c(" " = 1, header_counts)) %>%
    kable_styling(latex_options = "hold_position") %>%
    row_spec(group_start_rows, bold = TRUE) %>%
    row_spec(hline_rows, hline_after = TRUE)
}

build_contrast_table(all_group_contrast_results, "adipo")
build_contrast_table(all_group_contrast_results, "myo")

# ALR with a robust denominator chosen by Procrustes
# library(easyCODA)
# mat <- new_mats[["adipo"]][["messmer"]][["dwls"]]
# keepParts <- colMeans(mat == 0) < 0.8
# mat_processed <- mat[, keepParts]
# cat("\tParts kept: ", colnames(mat_processed), "\n")
# 
# if (any(mat_processed == 0)) {
#   mat_processed <- zCompositions::multRepl(mat_processed, label = 0, dl = rep(1e-5, ncol(mat_processed)))
# }
# mat_processed <- CLOSE(mat_processed)
# res <- FINDALR(mat_processed)
# res$procrust.ref
# PLOT.LRA(LRA(mat_processed), )

# check common zero-processed matrices' cell-types per reference
# x <- rapply(comp_mats, function(m) {
#   cat("__________________\n")
#   keepParts <- colMeans(m == 0) < 0.8
#   mat_processed <- m[, keepParts]
#   cat("\tParts kept: ", colnames(mat_processed), "\n")
#   
#   if (any(mat_processed == 0)) {
#     mat_processed <- zCompositions::multRepl(mat_processed, label = 0, dl = rep(1e-5, ncol(mat_processed)))
#   }
#   
#   res <- FINDALR(mat_processed)
#   cat("\tprocrust.cor:", paste0(colnames(mat_processed), "=", round(res$procrust.cor, 3), collapse = ", "), "\n")
#   
#   mat_processed
# }, how = "replace")
# for (exp in exps) for (ref in refs) print(Reduce(intersect, sapply(x[[exp]][[ref]], colnames)))

# Further possibilities
# james multivariate test is unnecessary, since it requires n < p and
# is parametric.
# james(y1, y2, a = 0.05, R = 999, graph = FALSE)

# Logcontrast regression,
# but not relevant :( since this is to predict a continuous response variable
# Compositional::lc.reg(y, x, z = NULL, xnew = NULL, znew = NULL)


