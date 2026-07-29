patched_DEAnalysisMASTOptimized <- function(scdata, id, path, verbose = FALSE) {
    uniqueIds <- unique(id)
    list_lrTest.table <- as.list(rep(0, length(uniqueIds)))
    pseudo.count <- 0.1
    data.used.log2 <- log2(scdata + pseudo.count)
    colnames(data.used.log2) <- make.unique(colnames(data.used.log2))
    diff.cutoff <- 0.5
    for (i in uniqueIds) {
      mist_file <- paste0(path, "/", i, "_MIST.RData")
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
  
  assignInNamespace("DEAnalysisMASTOptimized", patched_DEAnalysisMASTOptimized, ns = "DWLS")

