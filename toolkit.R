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
    df
  }, mats, method_names, SIMPLIFY = FALSE)
  
  do.call(rbind, mats_df) %>%
    inner_join(meta, by = "sample_id") %>%
    pivot_longer(all_of(cell_types), names_to = "cell_type", values_to = "proportion") %>%
    relocate(c("cell_type", "proportion"), .after = "method")
}
