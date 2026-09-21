read_tpm <- function(path) {
  bulk_matrix <- data.table::fread(path)
  
  char_cols <- names(bulk_matrix)[sapply(bulk_matrix, is.character)]
  if (length(char_cols) == 0) stop("No character column found for gene IDs.")
  gene_col_name <- if ("gene_id" %in% char_cols) "gene_id" else char_cols[1]
  if (length(char_cols) > 1) {
    warning("Multiple character columns found: ", paste(char_cols, collapse = ", "),
            ". Using '", gene_col_name, "' as gene identifiers; discarding others.", call. = FALSE)
  }
  gene_col <- bulk_matrix[[gene_col_name]]
  bulk_matrix[, (char_cols) := NULL]
  mat <- as.matrix(bulk_matrix)
  rownames(mat) <- gene_col
  return(mat)
}

# take arbitrary number of composition matrices from different methods and merge; rownames are samples
# mats are sample x cell_types
# meta should have a "sample_id" column matching to matrices' rownames
# output is data frame in long format with columns: sample_id, method, cell_type, proportion, and other meta columns
merge_methods <- function(..., meta, method_names = NULL) {
  mats <- list(...)
  if (is.null(method_names)) method_names <- paste0("m", seq_along(mats))
  stopifnot(length(method_names) == length(mats))
  
  cell_types <- Reduce(union, lapply(mats, colnames))
  
  mats_df <- mapply(function(mat, name) {
    df <- data.frame(sample_id = rownames(mat), mat, method = name, check.names = FALSE)
    missing_cols <- setdiff(cell_types, colnames(mat))
    df[, missing_cols] <- 0
    df[, c("sample_id", "method", cell_types)]
  }, mats, method_names, SIMPLIFY = FALSE)
  
  do.call(rbind, mats_df) %>%
    inner_join(meta, by = "sample_id") %>%
    pivot_longer(all_of(cell_types), names_to = "cell_type", values_to = "proportion") %>%
    relocate(c("cell_type", "proportion"), .after = "method")
}

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
      message("Skipping ", i, " (MIST file exists)")
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


# Samples in the bulk matrix are divided to be processed in parallel.
# Used in dwls.R scripts in reference data set deconvolutions.
deconvolute_dwls_parallel <- function(bulk, sig, ncores = 4, logfile = "dwls_deconv_progress.log") {
  chunks <- parallel::splitIndices(ncol(bulk), ncores)
  n_total <- ncol(bulk)
  results <- parallel::mclapply(chunks, function(idx) {
    chunk_results <- lapply(idx, function(i) {
      res <- deconvolute_dwls(bulk_gene_expression = bulk[, i, drop = FALSE],
                              signature = sig,
                              dwls_submethod = "DampenedWLS",
                              verbose = FALSE)
      cat(sprintf("Done: %s (%d/%d)\n", colnames(bulk)[i], i, n_total), file = logfile, append = TRUE)
      res
    })
    do.call(rbind, chunk_results)
  }, mc.cores = ncores)
  do.call(rbind, results)
}
