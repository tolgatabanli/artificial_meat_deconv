library(Seurat)
library(parallel)

setwd("/nfs/home/students/t.tabanli/artifical_meat_deconv/counts_cattle_atlas")

rds_files <- list.files("/nfs/data2/huji_artificial_meat/cattle_atlas/NG_scRNA", full.names = T) 
rds_files <- rds_files[endsWith(rds_files, ".rds")]

sample_from_seurat <- function(seurat_obj, n = 1000) {
  sampled_cols <- sample(colnames(seurat_obj), size = n)
  seurat_obj <- seurat_obj[, sampled_cols]
  gc()
  return(seurat_obj)
}

sampled_list <- mclapply(rds_files, function(f) {
  set.seed(digest::digest2int(f))
  obj <- readRDS(f)
  sampled <- sample_from_seurat(obj)
  rm(obj); gc()
  sampled
}, mc.cores = 19)
# max file is 11 GB, min core in a node is 64, and min memory is 384 GB
# takes ca. 15 mins

final_obj <- merge(sampled_list[[1]], y = sampled_list[-1])
final_obj$label <- trimws(final_obj$label)
final_obj$cell_type <- trimws(final_obj$cell_type)

saveRDS(final_obj, file = "sampled_cattle_atlas.RDS")

# write in mtx format for python reads
mat <- GetAssayData(final_obj, layer = "counts")
Matrix::writeMM(mat, "sampled_cattle_atlas/counts.mtx")
write.csv(data.frame(gene = rownames(mat)), "sampled_cattle_atlas/genes.csv", row.names = FALSE)
write.csv(final_obj@meta.data, "sampled_cattle_atlas/metadata.csv")