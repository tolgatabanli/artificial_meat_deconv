library(compositions)
library(vegan)
library(tidyverse)
setwd("~/artifical_meat_deconv/")
source("toolkit.R")

# Reading routine for composition matrices and sample metadata
source("read_compositions_and_metadata.R")

new_mats <- rapply(comp_mats, function(m) { m[m < 1e-5] <- 0; m }, how ="replace")

exps <- c("adipo", "myo")
refs <- c("cattle", "wang", "messmer_full") # messmer
methods <- c("dwls", "rectangle", "scaden")
m_comb <- combn(methods, 2)

# START
process_mat <- function(mat) {
  mat <- acomp(mat)
  keepParts <- colMeans(mat == 0) < 0.2
  mat_processed <- mat#[, keepParts]
  cat("\tParts kept: ", colnames(mat_processed), "\n")
  
  if (any(mat_processed == 0)) {
    cat("\t\tMult. repl. applied.\n")
    mat_processed <- zCompositions::multRepl(mat_processed, label = 0, dl = rep(1e-5, ncol(mat_processed)))
  }
  
  mat_processed <- acomp(mat_processed)
  rownames(mat_processed) <- rownames(mat)
  
  return(mat_processed)
}
processed_mats <- rapply(new_mats, process_mat, how = "replace")

processed_intersected_mats <- processed_mats
for (exp in exps) {
  for (ref in refs) {
    mat_methods <- processed_mats[[exp]][[ref]] 
    commonTypes <- Reduce(intersect, sapply(mat_methods, colnames))
    processed_intersected_mats[[exp]][[ref]] <- lapply(mat_methods, function(m) m[,commonTypes])
  }
}

# compare DWLS and Rectangle apart from Unknowns
longer_df <- data.frame()
for (exp in exps) for (ref in refs) {
  mat_dwls <- processed_intersected_mats[[exp]][[ref]][["dwls"]]
  mat_rectangle <- processed_intersected_mats[[exp]][[ref]][["rectangle"]]
  
  meta <- if (exp == "adipo") raw_meta_adipo else raw_meta_myo
  merged <- merge_methods(mat_dwls, mat_rectangle, meta = meta, method_names = c("dwls", "rectangle"))
  merged$ref <- ref
  merged$exp <- exp
  merged <- merged[, c("exp","ref","sample_id", "method", "cell_type", "proportion", "day")]
  
  longer_df <- rbind(longer_df, merged)
}
longer_df %>%
  pivot_wider(names_from = method, values_from = proportion) %>%
  filter(dwls > 1e-5 & rectangle) %>%
  filter(cell_type != "Unknown") %>%
  ggplot(aes(x = dwls, y = rectangle)) +
    geom_point() +
    geom_smooth(method = "lm") + theme_bw() +
    facet_grid(exp ~ ref, scales = "free")
  

results_df <- data.frame()
for (exp in exps) {
  meta <- if (exp == "adipo") raw_meta_adipo else raw_meta_myo
  cat("Exp: ", exp, "\n")
  for (ref in refs) {
    cat("\tRef: ", ref, "\n")
    pairs <- list()
    r_vals <- c()
    rmse_vals <- c()
    mats <- processed_intersected_mats[[exp]][[ref]]
    # for each method pair
    for (i in 1:ncol(m_comb)) {
      m1 <- m_comb[1, i]
      m2 <- m_comb[2, i]
      cat("\t\tPair: ", m1, " vs ", m2, "\n")
      # raw proportions (> replacement value)
      mat1 <- mats[[m1]]
      mat2 <- mats[[m2]]
      merged <- merge_methods(mat1, mat2, meta = meta, method_names = c(m1, m2)) %>%
        pivot_wider(names_from = method, values_from = proportion) %>%
        filter(.data[[m1]] > 1e-5 | .data[[m2]] > 1e-5)
      
      r_val <- cor(merged[[m1]], merged[[m2]], method = "spearman")
      rmse_val <- sqrt(mean((merged[[m1]] - merged[[m2]])^2))
      r_vals <- c(r_vals, r_val)
      rmse_vals <- c(rmse_vals, rmse_val)
      pairs[[i]] <- c(m1, m2)
    }
    for (i in seq_along(pairs)) {
      results_df <- rbind(results_df, data.frame(
        exp = exp, ref = ref,
        pair = paste(pairs[[i]], collapse = "-"),
        r = r_vals[i],
        rmse = rmse_vals[i]
      ))
    }
  }
}

exp_labels <- c(adipo = "Adipogenesis", myo = "Myogenesis")
ref_labels <- c(messmer = "Messmer2023", cattle = "CattleAtlas", wang = "Wang2023")
method_labels <- c(dwls = "DWLS", rectangle = "Rectangle", scaden = "Scaden")

results_df %>% mutate(pair = str_replace_all(pair, method_labels)) %>%
ggplot(aes(x = ref, y = pair, size = r, color = rmse)) +
  geom_point() +
  geom_text(aes(label = sprintf("R=%.2f\nRMSE=%.2f", r, rmse)),
            size = 3, color = "black", hjust = 0, vjust = 0.5,
            nudge_x = 0.2, show.legend = FALSE) +
  scale_size_continuous(range = c(2, 10)) +
  scale_x_discrete(labels = ref_labels,
                   expand = expansion(mult = c(0.1, 0.35))) +
  ggh4x::facet_nested(. ~ exp, scales = "free_x", space = "free_x",
                      labeller = labeller(exp = as_labeller(exp_labels))) +
  theme_bw() +
  labs(x = "Reference Dataset", y = "Method Pairs") +
  theme_minimal() +
  theme(strip.text.x = element_text(size = 12))

