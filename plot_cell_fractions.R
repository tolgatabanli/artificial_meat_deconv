library(tidyverse)
source("~/artifical_meat_deconv/toolkit.R")
setwd("~/artifical_meat_deconv/")

# Reading routine for composition matrices and sample metadata
source("read_compositions_and_metadata.R") # -> comp_mats

# Apply detection limit of 1e-5
exps <- c("adipo", "myo")
refs <- c("cattle", "messmer", "wang", "messmer_full")
methods <- c("dwls", "rectangle", "scaden")

new_mats <- rapply(comp_mats, function(mat){
  mat[mat < 1e-5] <- 0
  mat
}, how = "replace")

# input: name of the reference data set and experiment type
# output: saves plot of barplot compositions
plot_method_comparison <- function(ref, exp, meta, facet_formula, outdir = "figures_fractions", width = 9, height = 10) {
  mats <- new_mats[[exp]][[ref]]
  
  merged <- merge_methods(mats[["rectangle"]], mats[["dwls"]], mats[["scaden"]],
                          meta = meta, method_names = c("Rectangle", "DWLS", "Scaden"))
  
  colors <- pals::alphabet(n_distinct(merged$cell_type))
  names(colors) <- unique(merged$cell_type)
  if ("Unknown" %in% names(colors)) colors["Unknown"] <- "gray90"
  
  p <- merged %>%
    ggplot(aes(x = sample_id, y = proportion, fill = cell_type)) +
    geom_col(linewidth = 0.8) +
    coord_flip() +
    facet_grid(facet_formula, scales = "free_y", space = "free_y") +
    scale_fill_manual(values = colors, name = "Cell type") +
    theme_minimal() +
    labs(y = "Proportion", x = NULL) +
    theme(axis.text.y = element_blank(),
          strip.text.y.right = element_text(angle = 0, hjust = 0.5, face = "bold"),
          strip.background.y = element_rect(fill = "gray90", color = "black", linewidth = 0.8),
          panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.4),
          panel.spacing.y = unit(0.3, "lines"))
  
  ggsave(file.path(outdir, paste0(ref, "_", exp, "_method_comparison.png")), plot = p, width = width, height = height)
  invisible(p)
}

# Messmer
plot_method_comparison("messmer", "adipo",
                       meta = raw_meta_adipo,
                       facet_formula = individual + day ~ method, width = 12, height = 5)
plot_method_comparison("messmer", "myo",
                       meta = raw_meta_myo,
                       facet_formula = scaffold + day ~ method, width = 12, height = 5)
plot_method_comparison("messmer_full", "adipo",
                       meta = raw_meta_adipo,
                       facet_formula = individual + day ~ method, width = 12, height = 5)
plot_method_comparison("messmer_full", "myo",
                       meta = raw_meta_myo,
                       facet_formula = scaffold + day ~ method, width = 12, height = 5)
# CattleAtlas
plot_method_comparison("cattle", "adipo",
                       meta = raw_meta_adipo,
                       facet_formula = individual + day ~ method, width = 12, height = 5)
plot_method_comparison("cattle", "myo",
                       meta = raw_meta_myo,
                       facet_formula = scaffold + day ~ method, width = 12, height = 5)

# Wang2023
plot_method_comparison("wang", "adipo",
                       meta = raw_meta_adipo,
                       facet_formula = individual + day ~ method, width = 12, height = 5)
plot_method_comparison("wang", "myo",
                       meta = raw_meta_myo,
                       facet_formula = scaffold + day ~ method, width = 12, height = 5)

# Cell-type-wise raw proportions of selected cell-types across days
messmerCellTypes <- c("FAPs", "SCs", "SMMCs", "Myocytes")
cattleCellTypes <- c("Stromal", "Muscle", "Epithelial")
wangCellTypes <- c("FAP-1", "FAP-2", "Fast myofiber", "SC", "Tenocyte")

plot_day_composition <- function(ref, exp, method, meta, cell_types,
                                 group_var = "day", facet_var = NULL,
                                 outdir = "figures_frac_celltype_wise") {
  day_labels <- meta %>%
    count(.data[[group_var]]) %>%
    mutate(label = paste0(.data[[group_var]], " (n=", n, ")")) %>%
    dplyr::select(all_of(group_var), label) %>%
    deframe()
  
  facet <- if (is.null(facet_var))
    facet_wrap(~ cell_type, nrow = 1)
  else if (length(facet_var) == 1)
    facet_grid(rows = vars(.data[[facet_var]]), cols = vars(cell_type))
  else
    facet_grid(reformulate("cell_type", response = paste(facet_var, collapse = "+")))
  
  p <- comp_mats[[exp]][[ref]][[method]] %>%
    as.data.frame() %>%
    rownames_to_column("sample_id") %>%
    inner_join(meta, by = "sample_id") %>%
    dplyr::select(sample_id, all_of(group_var), any_of(facet_var), all_of(cell_types)) %>%
    pivot_longer(all_of(cell_types), names_to = "cell_type", values_to = "frac") %>%
    ggplot(aes(x = .data[[group_var]], y = frac)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, height = 0, alpha = 0.6, size = 1.5) +
    facet + theme_bw()
  
  ggsave(file.path(outdir, paste0(ref, "_", exp, "_", method, "_day_composition.png")), plot = p, height = 8)
  invisible(p)
}

for (exp in c("myo", "adipo")) {
  meta <- if (exp == "myo") raw_meta_myo else raw_meta_adipo
  for (method in c("dwls", "rectangle", "scaden"))
    plot_day_composition("messmer", exp, method, meta, messmerCellTypes, facet_var = c("individual"))
}

### Check Scaden issue in cattle - Scaden estimate correlations with lineage
lineage_props <- table(cattle$cell_lineage) / sum(table(cattle$cell_lineage))
lineage_props_df <- lineage_props %>% data.frame() %>% rename(lineage = Var1, lineage_props = Freq)

bind_rows(
  lineage_props_df %>%
    inner_join(data.frame(mean_prop = colMeans(comp_mats[["adipo"]][["cattle"]][["scaden"]]),
                          lineage = names(colMeans(comp_mats[["adipo"]][["cattle"]][["scaden"]]))), by = "lineage") %>%
    mutate(exp = "Adipogenesis"),
  lineage_props_df %>%
    inner_join(data.frame(mean_prop = colMeans(comp_mats[["myo"]][["cattle"]][["scaden"]]),
                          lineage = names(colMeans(comp_mats[["myo"]][["cattle"]][["scaden"]]))), by = "lineage") %>%
    mutate(exp = "Myogenesis")) %>%
  ggplot(aes(x = lineage_props, y = mean_prop, color = exp, label = lineage)) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("Adipogenesis" = "orange3", "Myogenesis" = "darkred")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  ggrepel::geom_text_repel(max.overlaps = Inf, box.padding = 0.5, size = 4) +
  labs(x = "Proportion in the dataset", y = "Scaden Estimates", color = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(),
        panel.background = element_blank(),
        legend.position = "bottom")
ggsave("figures_fractions/cattle_scaden_scatter.png", width = 10, height = 6)


### Raw estimate distributions
exps <- c("adipo", "myo")
refs <- c("cattle", "messmer", "wang")
methods <- c("dwls", "rectangle", "scaden")

collect_method <- function(method_name) {
  do.call(rbind, lapply(exps, function(exp) {
    do.call(rbind, lapply(refs, function(ref) {
      mat <- comp_mats[[exp]][[ref]][[method_name]]
      as.data.frame(mat) %>%
        rownames_to_column("sample_id") %>%
        pivot_longer(-sample_id, names_to = "cell_type", values_to = "frac") %>%
        mutate(experiment = exp, ref = ref)
    }))
  }))
}

rectangle_df <- collect_method("rectangle")
dwls_df      <- collect_method("dwls")
scaden_df    <- collect_method("scaden")

# all methods in one plot
all_df <- bind_rows(
  rectangle_df %>% mutate(method = "Rectangle"),
  dwls_df %>% mutate(method = "DWLS"),
  scaden_df %>% mutate(method = "Scaden")
)

# draw the zeros on the very left
zero_counts <- all_df %>%
  group_by(method, experiment, ref) %>%
  summarise(n_zeros = sum(frac == 0), .groups = "drop") %>%
  mutate(vjust_val = as.numeric(factor(ref)),
         ref_label = c("messmer" = "Messmer", "cattle" = "CattleAtlas", "wang" = "Wang")[ref]) %>%
  mutate(label = paste0(ref_label, " #zeros = ", n_zeros))

x_zero_pos <- min(all_df$frac[all_df$frac > 0]) / 10

all_df %>% filter(frac > 0) %>% ggplot(aes(x = frac, fill = ref)) +
  geom_histogram(bins = 60, position = "identity", alpha = 0.7) +
  geom_text(data = zero_counts,
            aes(x = x_zero_pos, y = Inf, label = label,
                vjust = vjust_val * 1.5, hjust = 0, color = ref),
            size = 3, inherit.aes = FALSE) +
  scale_color_manual(values = RColorBrewer::brewer.pal(3, "Dark2"), guide = "none") +
  geom_vline(xintercept = 1e-05, linetype = "dashed") +
  facet_grid(method ~ experiment, scales = "free_y",
             labeller = labeller(
               experiment = c(adipo = "Adipogenesis", myo = "Myogenesis")
             )) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Dark2"),
                    labels = c("messmer" = "Messmer2023", "cattle" = "CattleAtlas", "wang" = "Wang2023"),
                    name = "Reference") +
  scale_x_log10(breaks = c(x_zero_pos, 1e-24, 1e-17, 1e-10, 1e-03, 1),
                labels = c(expression(bold("0")), expression(10^-24), expression(10^-17), expression(10^-10), expression(10^-3), expression(bold("1")))) +
  labs(y = "Number", x = "Fraction estimate") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(),
        panel.background = element_blank(),
        strip.text.x = element_text(size = 10, face = "bold"),
        strip.text.y = element_text(size = 10, face = "bold"),
        strip.background = element_blank(),
        legend.position = "bottom")

ggsave("figures_fractions/raw_estimate_distribution.png", width = 10, height = 6)


  #geom_col(data = zero_counts, aes(x = x_zero_pos, y = n_zeros, fill = ref), position = "identity", alpha = 0.7, inherit.aes = FALSE) +
  