library(Seurat)
library(omnideconv)
library(tidyverse)
library(parallel)
source("~/artifical_meat_deconv/toolkit.R")

setwd("~/artifical_meat_deconv/counts_wang2023")

message("Reading wang2023")
wang <- readRDS("wang_ids_mapped.RDS")
wang$lables <- trimws(wang$lables)
wang$lables <- gsub("'", "", wang$lables) # Henle's loop causes issues in column parse
wang$lables <- gsub("/", "_", wang$lables)
wang$lables[wang$lables == "Muc-2"] <- "MuC-2"
dim(wang)

### patch
patched_DEAnalysisMASTOptimized <- function(scdata, id, path, verbose = FALSE) {
  uniqueIds <- unique(id)
  list_lrTest.table <- as.list(rep(0, length(uniqueIds)))
  pseudo.count <- 0.1
  data.used.log2 <- log2(scdata + pseudo.count)
  colnames(data.used.log2) <- make.unique(colnames(data.used.log2))
  diff.cutoff <- 0.5
  for (i in uniqueIds) {
    mist_file <- paste0(path, "/", i, "_mist.RData")
    if (file.exists(mist_file)) {
      message("Skipping ", i, " — MIST file exists")
      load(mist_file)
      index <- which(uniqueIds == i)
      list_lrTest.table[[index]] <- cluster_lrTest.table
      next
    }
    cells.symbol.list2 <- colnames(data.used.log2)[which(id == i)]
    cells.coord.list2 <- match(cells.symbol.list2, colnames(data.used.log2))
    cells.coord.list1 <- c(1:length(id))[-cells.coord.list2]
    data.used.log2.ordered <- cbind(data.used.log2[, cells.coord.list1],
                                    data.used.log2[, cells.coord.list2])
    group.v <- c(rep(0, length(cells.coord.list1)), rep(1, length(cells.coord.list2)))
    log2.stat.result <- stat.log2.optimized(data.used.log2.ordered, group.v, pseudo.count)
    de <- log2.stat.result[log2.stat.result$log2_fc > diff.cutoff, ]
    if (verbose) dim(de)
    data.1 <- data.used.log2[, cells.coord.list1, drop = FALSE]
    data.2 <- data.used.log2[, cells.coord.list2, drop = FALSE]
    groups <- c(rep("Cluster_Other", length(cells.coord.list1)), rep(i, length(cells.coord.list2)))
    groups <- as.character(groups)
    if (dim(de)[1] > 1) {
      genes.list <- rownames(de)
      counts <- as.data.frame(cbind(data.1[genes.list, ], data.2[genes.list, ]))
      vbeta.1 <- createDataForMAST(verbose, counts, groups)
      zlm.output <- doMastZlm(verbose, vbeta.1)
      zlm.lr_pvalue <- doLrTest(verbose, zlm.output)
      cluster_lrTest.table <- createLrTestTable(i, de, zlm.lr_pvalue)
      index <- which(uniqueIds == i)
      list_lrTest.table[[index]] <- cluster_lrTest.table
      saveMastResult(path, cluster_lrTest.table, i)
    }
  }
  return(list_lrTest.table)
}
environment(patched_DEAnalysisMASTOptimized) <- asNamespace("DWLS")
assignInNamespace("DEAnalysisMASTOptimized", patched_DEAnalysisMASTOptimized, ns = "DWLS")
###

mat <- as.matrix(GetAssayData(wang, assay = "RNA", layer = "counts"))


message("Starting signature build...")
omnideconv::build_model_dwls(single_cell_object = mat,
                             cell_type_annotations = wang$lables,
                             dwls_method = "mast_optimized",
                             path = "sig_dwls",
                             verbose = T, ncores = 7)
message("Signature building finished")


load("sig_dwls/Sig.RData")


message("Reading tpms...")

# Read TPM matrices and keep only relevant samples.
# TODO: Check if meta's work!
bulk_meta_adipo <- read_csv('/nfs/data2/huji_artificial_meat/adipogenesis/sample_sheet_cleaned.csv', col_names = T) %>%
  janitor::clean_names()
bulk_meta_myo <- read_csv('/nfs/data2/huji_artificial_meat/myogenesis/sample_sheet_cleaned.csv', col_names = T) %>%
  janitor::clean_names()

bulk_meta_adipo <- bulk_meta_adipo %>%
  dplyr::select(library_name, sample_id, title, group, bio_replicate, individual, treatment, day) %>%
  mutate(day = factor(day, levels = c("D0", "D2", "D6", "D14")))
bulk_meta_myo <- bulk_meta_myo %>%
  dplyr::select(library_name, sample_id, group, day, bio_replicate, treatment) %>%
  mutate(sample_id = gsub(".*_", "", sample_id)) %>%
  mutate(day = factor(day, levels = c("d0", "d3", "d7", "d21")),
         group = gsub("(?!Sat)D\\d*", "\\", group, perl = T)) %>%
  separate_wider_delim(cols = group, delim = "_", names = c(NA, "group"), too_few = "align_end")

bulk_matrix_adipo <- read_tpm('/nfs/data2/huji_artificial_meat/adipogenesis/adipogenesis_tpm.tsv')
bulk_matrix_adipo <- bulk_matrix_adipo[, bulk_meta_adipo$sample_id]

bulk_matrix_myo <- read_tpm('/nfs/data2/huji_artificial_meat/myogenesis/nfcore_out/star_salmon/salmon.merged.gene_tpm.tsv')
colnames(bulk_matrix_myo) <- gsub(".*_", "", colnames(bulk_matrix_myo))
bulk_matrix_myo <- bulk_matrix_myo[, bulk_meta_myo$sample_id]

length(intersect(rownames(Sig), rownames(bulk_matrix_adipo)))
length(intersect(rownames(Sig), rownames(bulk_matrix_myo)))

message("Deconvoluting...")

deconvolute_dwls_parallel <- function(bulk, sig, ncores = 4, logfile = "deconv_progress.log") {
  chunks <- parallel::splitIndices(ncol(bulk), ncores)
  n_total <- ncol(bulk)
  results <- mclapply(chunks, function(idx) {
    chunk_results <- lapply(idx, function(i) {
      res <- deconvolute_dwls(bulk_gene_expression = bulk[, i, drop = FALSE],
                              signature = sig,
                              dwls_submethod = "DampenedWLS",
                              verbose = FALSE)
      cat(sprintf("Done: %s (%d/%d)\n", colnames(bulk)[i], i, n_total),
          file = logfile, append = TRUE)
      res
    })
    do.call(rbind, chunk_results)
  }, mc.cores = ncores)
  do.call(rbind, results)
}

gc()

message("___adipo___")
comp_adipo <- deconvolute_dwls_parallel(bulk_matrix_adipo, Sig, ncores = 20)
message("___myo___")
comp_myo   <- deconvolute_dwls_parallel(bulk_matrix_myo,   Sig, ncores = 20)


saveRDS(comp_adipo, "sig_dwls/adipo_compositions_wang2023_dwls_signature.Rds")
saveRDS(comp_myo, "sig_dwls/myo_compositions_wang2023_dwls_signature.Rds")

