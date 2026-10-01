# converted from Rmd and modified...

## ----library------------------------------------------------------------------
library(omnideconv)
library(tidyverse)
suppressMessages(conflicted::conflict_prefer_all("dplyr"))

library(ggpubr)
library(pals)
library(Seurat)
setwd("~/artifical_meat_deconv/messmer_counts")
source("~/artifical_meat_deconv/toolkit.R")


## ----read-counts--------------------------------------------------------------
messmer_inj1 <- readRDS("~/artifical_meat_deconv/messmer_counts/Results/Inj1.Rds")
mat <- GetAssayData(messmer_inj1, assay = "RNA", layer = "counts")
mat <- as.matrix(mat)

col_meta <- data.table::fread('/nfs/data2/huji_artificial_meat/deconvolution_data/external_sc_references/messmer2023_GSE211428/GSE211428_combined_meta.tsv.gz')
col_meta <- as.data.frame(col_meta) %>%
  filter(orig.ident == "Injection1") %>%
  mutate(Unique_barcode = gsub("\\.1_1$", "-1", Unique_barcode)) %>%
  column_to_rownames("Unique_barcode")
# make the names path-compatible, since DWLS creates dirs
col_meta$predicted.pheno <- gsub("/", "-", col_meta$predicted.pheno)

cat("colnames(mat):", length(colnames(mat)), "\n")
cat("col_meta barcodes:", nrow(col_meta), "\n")
cat("Intersection:", length(intersect(colnames(mat), rownames(col_meta))), "\n")
cat("In mat not in col_meta:", length(setdiff(colnames(mat), rownames(col_meta))), "\n")
cat("In col_meta not in mat:", length(setdiff(rownames(col_meta), colnames(mat))), "\n")

common <- intersect(colnames(mat), rownames(col_meta))
mat <- mat[, common]
col_meta <- col_meta[common, ]

gc()
#' 
## ----celltype-occurrences, fig.height=5---------------------------------------
table(cell_type = col_meta$predicted.pheno) %>%
  as.data.frame() %>%
  mutate(cell_type = gsub("Monocytes - M.$", "Monocytes-M\u03A6", cell_type)) %>%
  filter(Freq > 0) %>%
  ggplot(aes(x = cell_type, y = Freq)) +
  geom_bar(stat = "identity", fill = "white", color = "gray") +
  #scale_y_log10(labels = scales::label_log()) +
  geom_text(aes(label = Freq), vjust = -0.2) +
  labs(x = "", y = "Count", title = "Messmer2023") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 11), panel.grid.major.x = element_blank())
ggsave("Figures/messmer2023_distro_cell_types_inj1.png", width = 7, height = 6)


## ----build-signature, eval=FALSE----------------------------------------------
# omnideconv::build_model_dwls(single_cell_object = mat,
#                              cell_type_annotations = col_meta$predicted.pheno,
#                              dwls_method = "mast_optimized",
#                              path = "sig_dwls_withTeno_inj1",
#                              verbose = T, ncores = 6) # TODO: defaults: diff_cutoff = 0.5, pval_cutoff = 0.01
# gc()
# # without tenocytes
# noTenocyte_index <- col_meta$predicted.pheno != "Tenocytes"
# omnideconv::build_model_dwls(single_cell_object = mat[, noTenocyte_index],
#                              cell_type_annotations = col_meta$predicted.pheno[noTenocyte_index],
#                              dwls_method = "mast_optimized",
#                              path = "sig_dwls_inj1",
#                              verbose = T, ncores = 10)
# gc()
#' 
## ----signature-check, eval = FALSE--------------------------------------------
# for the markdown, this block is skipped and
# the signature without tenocytes is used

load('sig_dwls_withTeno_inj1/Sig.RData')

colnames(Sig) <- gsub("Monocytes - M.$", "Monocytes-M\u03A6", colnames(Sig))

sig_df <- Sig %>% as.data.frame()
sig_df <- sig_df %>%
  rownames_to_column("genes") %>%
  pivot_longer(cols = -genes, names_to = "type", values_to = "val") %>%
  pivot_wider(names_from = type, values_from = val) %>%
  drop_na()

fap_genes <- sig_df$genes[sig_df$FAPs > 0]
ten_genes <- sig_df$genes[sig_df$Tenocytes > 0]

# jaccard similarity
cat("Jaccard similarity: ")
length(intersect(fap_genes, ten_genes)) / length(union(fap_genes, ten_genes))

# jaccard similarity matrix across all cell types
cell_types <- setdiff(colnames(sig_df), "genes")
gene_sets  <- lapply(cell_types, function(ct) sig_df$genes[sig_df[[ct]] > 0])
# all nonzero weight genes for each cell type
names(gene_sets) <- cell_types

jaccard_mat <- matrix(NA_real_, nrow = length(cell_types), ncol = length(cell_types),
                      dimnames = list(cell_types, cell_types))
for (i in cell_types) {
  for (j in cell_types) {
    jaccard_mat[i, j] <- length(intersect(gene_sets[[i]], gene_sets[[j]])) /
      length(union(gene_sets[[i]], gene_sets[[j]]))
  }
}

pheatmap::pheatmap(jaccard_mat,
                   display_numbers = TRUE,
                   number_format = "%.2f",
                   cluster_rows = TRUE, cluster_cols = TRUE,
                   main = "Jaccard similarity between cell type signatures",
                   filename = "Figures/messmer2023_celltype_jaccard_heatmap_inj1.png",
                   width = 8, height = 7)

# correlation of common genes:
sig_df %>%
  filter(FAPs > 0, Tenocytes > 0) %>%
  ggplot(aes(x = FAPs, y = Tenocytes)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  scale_x_log10() +
  scale_y_log10() +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top",
           aes(label = paste(..r.label.., ..rr.label.., sep = "~`,`~")))
ggsave("Figures/messmer2023_corr_commongenes_tenocytes_vs_faps_inj1.png")

# pearson correlation matrix of common genes across all cell types
corr_mat <- matrix(NA_real_, nrow = length(cell_types), ncol = length(cell_types),
                   dimnames = list(cell_types, cell_types))
for (i in cell_types) {
  for (j in cell_types) {
    df_ij <- sig_df %>% filter(.data[[i]] > 0, .data[[j]] > 0)
    corr_mat[i, j] <- if (nrow(df_ij) >= 2) {
      cor(log10(df_ij[[i]]), log10(df_ij[[j]]), method = "pearson")
    } else {
      NA_real_
    }
  }
}

pheatmap::pheatmap(corr_mat,
                   display_numbers = TRUE,
                   number_format = "%.2f",
                   cluster_rows = TRUE, cluster_cols = TRUE,
                   color = colorRampPalette(c("blue", "white", "red"))(50),
                   breaks = seq(-1, 1, length.out = 51),
                   main = "Pearson correlation between cell type signatures",
                   filename = "Figures/messmer2023_celltype_correlation_heatmap_inj1.png",
                   width = 8, height = 7)

gc()
#' 
## ----signature-load, echo=FALSE-----------------------------------------------
# load Signature without Teno for actual deconvolution
load('sig_dwls_withTeno_inj1/Sig.RData')
colnames(Sig) <- gsub("Monocytes - M.$", "Monocytes-M\u03A6", colnames(Sig))

#' 
#' 
#' # Deconvolution


## ----read-meta----------------------------------------------------------------
bulk_matrix_adipo <- read_tpm("/nfs/data2/huji_artificial_meat/adipogenesis/adipogenesis_tpm.tsv")
bulk_meta_adipo <- read.csv("/nfs/data2/huji_artificial_meat/adipogenesis/sample_sheet_cleaned.csv", header = T) %>%
  mutate(day = factor(day, levels = unique(day)[order(as.numeric(gsub("[^0-9]", "", unique(day))))]))
bulk_matrix_adipo <- bulk_matrix_adipo[, bulk_meta_adipo$sample_id]

bulk_matrix_myo <- read_tpm("/nfs/data2/huji_artificial_meat/myogenesis/nfcore_out/star_salmon/salmon.merged.gene_tpm.tsv")
colnames(bulk_matrix_myo) <- gsub(".*_", "", colnames(bulk_matrix_myo))
bulk_meta_myo <- read.csv("/nfs/data2/huji_artificial_meat/myogenesis/sample_sheet_minimal.csv", header = T) %>%
  mutate(sample_id = gsub(".*_", "", sample_id))
bulk_matrix_myo <- bulk_matrix_myo[, gsub(".*_", "", colnames(bulk_matrix_myo)) %in% gsub(".*_", "", bulk_meta_myo$sample_id)]

## ----read-tpms, eval=FALSE----------------------------------------------------
cat("intersection of rownames of Sig and adipo bulk:")
length(intersect(rownames(Sig), rownames(bulk_matrix_adipo))) # 724
cat("intersection of rownames of Sig and myo bulk:")
length(intersect(rownames(Sig), rownames(bulk_matrix_myo)))

message("Deconvoluting...")
# deconvolute_dwls_parallel is in toolkit.R
gc()
message("___adipo___")
comp_adipo <- deconvolute_dwls_parallel(bulk_matrix_adipo, Sig, ncores = 38)
gc()
message("___myo___")
comp_myo   <- deconvolute_dwls_parallel(bulk_matrix_myo,   Sig, ncores = 13)

saveRDS(comp_adipo, "sig_dwls_withTeno_inj1/adipo_compositions_messmer2023_dwls_signature.Rds")
saveRDS(comp_myo, "sig_dwls_withTeno_inj1/myo_compositions_messmer2023_dwls_signature.Rds")



#' 
#' Lastly, we plot heatmaps with sample-clustering and annotations.
#' 
## ----heatmaps, fig.width=8, fig.height=5--------------------------------------
# adipogenesis
adipo_annot <- bulk_meta_adipo %>%
  data.frame(row.names = bulk_meta_adipo$sample_id) %>%
  dplyr::select(individual, day, group)
pheatmap::pheatmap(t(comp_adipo), cluster_rows = F, cluster_cols = T,
                   annotation_col = adipo_annot,
                   filename = "Figures/messmer2023_dwls_adipo_heatmap_inj1.png", width = 10, height = 7)


# myogenesis
myo_annot <- bulk_meta_myo %>%
  filter(sample_id %in% rownames(comp_myo)) %>%
  dplyr::select(sample_id, day, group) %>%
  column_to_rownames("sample_id")
pheatmap::pheatmap(t(comp_myo), cluster_rows = F, cluster_cols = T,
                   annotation_col = myo_annot,
                   filename = "Figures/messmer2023_dwls_myo_heatmap_inj1.png", width = 10, height = 7)
