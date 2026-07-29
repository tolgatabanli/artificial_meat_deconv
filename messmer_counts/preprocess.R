  ## ----setup, echo=FALSE, results="hide"-----------------------------------------------------------------------------------
  # Set working directory here:
  setwd('~/artifical_meat_deconv/messmer_counts')
  
  #Set Figure directory
  dir.create("Figures", showWarning=FALSE)
  figdir <- './Figures'
  
  #Set Results directory
  dir.create("Results", showWarning=FALSE)
  resdir <- "./Results"
  
  knitr::opts_chunk$set(error=FALSE, warning=FALSE, message=FALSE, fig.path="Figures/", tidy = TRUE, echo = TRUE)
  options(bitmapType="cairo", width=100)
  
  library('reticulate')
  library('scater')
  library('biomaRt')
  library('DropletUtils')
  library('Seurat')
  library('dplyr')
  library("GEOquery")
  library('Matrix')
  
  ## ----data.in-------------------------------------------------------------------------------------------------------------
  counts_root <- "/nfs/data2/huji_artificial_meat/deconvolution_data/external_sc_references/messmer2023_GSE211428/GSE211428_RAW"
  
  Inj1_barcodes <- read.table(file.path(counts_root, 'GSM6469858_Inj1_barcodes.tsv.gz'))
  Inj1_features <-  read.table(file.path(counts_root, 'GSM6469858_Inj1_features.tsv.gz'))
  Inj1_meta <- read.table(file.path(dirname(counts_root), 'GSE211428_Inj1_meta.tsv.gz'), header = TRUE)
  Inj1_matrix <-  readMM(file.path(counts_root, 'GSM6469858_Inj1_matrix.mtx.gz'))
  colnames(Inj1_matrix) <- Inj1_barcodes$V1
  rownames(Inj1_matrix) <- Inj1_features$V1
  
  Inj2_barcodes <- read.table(file.path(counts_root, 'GSM6469859_Inj2_barcodes.tsv.gz'))
  Inj2_features <-  read.table(file.path(counts_root, 'GSM6469859_Inj2_features.tsv.gz'))
  Inj2_meta <- read.table(file.path(dirname(counts_root), 'GSE211428_Inj2_meta.tsv.gz'), header = TRUE)
  Inj2_matrix <-  readMM(file.path(counts_root, 'GSM6469859_Inj2_matrix.mtx.gz'))
  colnames(Inj2_matrix) <- Inj2_barcodes$V1
  rownames(Inj2_matrix) <- Inj2_features$V1
  
  Inj3_barcodes <- read.table(file.path(counts_root, 'GSM6469861_Inj3_barcodes.tsv.gz'))
  Inj3_features <-  read.table(file.path(counts_root, 'GSM6469861_Inj3_features.tsv.gz'))
  Inj3_meta <- read.table(file.path(dirname(counts_root), 'GSE211428_Inj3_meta.tsv.gz'), header = TRUE)
  Inj3_matrix <-  readMM(file.path(counts_root, 'GSM6469861_Inj3_matrix.mtx.gz'))
  colnames(Inj3_matrix) <- Inj3_barcodes$V1
  rownames(Inj3_matrix) <- Inj3_features$V1
  
  Inj4_barcodes <- read.table(file.path(counts_root, 'GSM6469863_Inj4_barcodes.tsv.gz'))
  Inj4_features <-  read.table(file.path(counts_root, 'GSM6469863_Inj4_features.tsv.gz'))
  Inj4_meta <- read.table(file.path(dirname(counts_root), 'GSE211428_Inj4_meta.tsv.gz'), header = TRUE)
  Inj4_matrix <-  readMM(file.path(counts_root, 'GSM6469863_Inj4_matrix.mtx.gz'))
  colnames(Inj4_matrix) <- Inj4_barcodes$V1
  rownames(Inj4_matrix) <- Inj4_features$V1
  
  Inj5_barcodes <- read.table(file.path(counts_root, 'GSM6469865_Inj5_barcodes.tsv.gz'))
  Inj5_features <-  read.table(file.path(counts_root, 'GSM6469865_Inj5_features.tsv.gz'))
  Inj5_meta <- read.table(file.path(dirname(counts_root), 'GSE211428_Inj5_meta.tsv.gz'), header = TRUE)
  Inj5_matrix <-  readMM(file.path(counts_root, 'GSM6469865_Inj5_matrix.mtx.gz'))
  colnames(Inj5_matrix) <- Inj5_barcodes$V1
  rownames(Inj5_matrix) <- Inj5_features$V1
  
  # check rownames -> should be ENSEMBL
  ensembl <- useEnsembl(biomart = 'ensembl', dataset="btaurus_gene_ensembl", version = 116)
  
  ensemblGenes <- getBM(attributes=c('ensembl_gene_id', 'chromosome_name', 'gene_biotype',
                                     'entrezgene_id', 'description', 'transcript_length'),
                        filters="ensembl_gene_id",
                        values=rownames(Inj1_matrix), mart=ensembl)
  
  gene.meta <- ensemblGenes[match(rownames(Inj1_matrix), ensemblGenes$ensembl_gene_id),]
  rownames(gene.meta) <- ifelse(is.na(gene.meta$ensembl_gene_id),
                                rownames(Inj1_matrix),
                                gene.meta$ensembl_gene_id)
  
  all(
    identical(rownames(Inj1_matrix), rownames(Inj2_matrix)),
    identical(rownames(Inj1_matrix), rownames(Inj3_matrix)),
    identical(rownames(Inj1_matrix), rownames(Inj4_matrix)),
    identical(rownames(Inj1_matrix), rownames(Inj5_matrix))
  )
  
  ## ----make.seurat---------------------------------------------------------------------------------------------------------
  Inj1 <- CreateSeuratObject(counts = Inj1_matrix, project = "SFscRNAseq", min.cells = 3)
  Inj2 <- CreateSeuratObject(counts = Inj2_matrix, project = "SFscRNAseq", min.cells = 3)
  Inj3 <- CreateSeuratObject(counts = Inj3_matrix, project = "SFscRNAseq", min.cells = 3)
  Inj4 <- CreateSeuratObject(counts = Inj4_matrix, project = "SFscRNAseq", min.cells = 3)
  Inj5 <- CreateSeuratObject(counts = Inj5_matrix, project = "SFscRNAseq", min.cells = 3)
  
  
  ## ----seurat.meta---------------------------------------------------------------------------------------------------------
  rownames(Inj1_meta) <- Inj1_meta$Cell_barcode
  rownames(Inj2_meta) <- Inj2_meta$Cell_barcode
  rownames(Inj3_meta) <- Inj3_meta$Cell_barcode
  rownames(Inj4_meta) <- Inj4_meta$Cell_barcode
  rownames(Inj5_meta) <- Inj5_meta$Cell_barcode
  
  Inj1 <- AddMetaData(Inj1, Inj1_meta)
  Inj2 <- AddMetaData(Inj2, Inj2_meta)
  Inj3 <- AddMetaData(Inj3, Inj3_meta)
  Inj4 <- AddMetaData(Inj4, Inj4_meta)
  Inj5 <- AddMetaData(Inj5, Inj5_meta)
  
  head(Inj1@meta.data)
  
  
  ## ----gene.meta-----------------------------------------------------------------------------------------------------------
  Inj1[['RNA']] <- AddMetaData(Inj1[['RNA']], gene.meta)
  Inj2[['RNA']] <- AddMetaData(Inj2[['RNA']], gene.meta)
  Inj3[['RNA']] <- AddMetaData(Inj3[['RNA']], gene.meta)
  Inj4[['RNA']] <- AddMetaData(Inj4[['RNA']], gene.meta)
  Inj5[['RNA']] <- AddMetaData(Inj5[['RNA']], gene.meta)
  head(Inj1[['RNA']]@meta.data)
  
  
  ## ----sng-----------------------------------------------------------------------------------------------------------------
  Inj1_SNG <-  subset(Inj1, subset = Category=="SNG") 
  Inj2_SNG <-  subset(Inj2, subset = Category=="SNG") 
  Inj3_SNG <-  subset(Inj3, subset = Category=="SNG") 
  Inj4_SNG <-  subset(Inj4, subset = Category=="SNG") 
  Inj5_SNG <-  subset(Inj5, subset = Category=="SNG") 
  
  par(mfrow=c(2,3))
  barplot(table(Inj1$Category), main = 'Injection 1')
  barplot(table(Inj2$Category), main = 'Injection 2')
  barplot(table(Inj3$Category), main = 'Injection 3')
  barplot(table(Inj4$Category), main = 'Injection 4')
  barplot(table(Inj5$Category), main = 'Injection 5')
  
  
  ## ----qc, fig.width = 10--------------------------------------------------------------------------------------------------
  #mt.genes <- gene.meta[gene.meta$chromosome_name == 'MT',]
  mt.genes <- na.omit(rownames(gene.meta)[!is.na(gene.meta$chromosome_name) & gene.meta$chromosome_name == 'MT'])
  
  Inj1_SNG[["percent.mt"]] <- PercentageFeatureSet(Inj1_SNG, features = mt.genes)
  Inj2_SNG[["percent.mt"]] <- PercentageFeatureSet(Inj2_SNG, features = mt.genes)
  Inj3_SNG[["percent.mt"]] <- PercentageFeatureSet(Inj3_SNG, features = mt.genes)
  Inj4_SNG[["percent.mt"]] <- PercentageFeatureSet(Inj4_SNG, features = mt.genes)
  Inj5_SNG[["percent.mt"]] <- PercentageFeatureSet(Inj5_SNG, features = mt.genes)
  
  plot1 <- VlnPlot(Inj1_SNG, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
  plot2 <- VlnPlot(Inj2_SNG, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  plot3 <- VlnPlot(Inj3_SNG, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
  plot4 <- VlnPlot(Inj4_SNG, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  plot5 <- VlnPlot(Inj5_SNG, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  plot1 / plot2 / plot3 / plot4 / plot5  & theme(axis.text.x=element_blank(), axis.title.x=element_blank())
  
  plot1 <- FeatureScatter(Inj1_SNG, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(Inj1_SNG, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot1 + plot2 & NoLegend()
  
  
  ## ----filter--------------------------------------------------------------------------------------------------------------
  featureQC <- isOutlier(Inj1_SNG$nFeature_RNA)
  libQC <- isOutlier(Inj1_SNG$nCount_RNA)
  mtQC <- isOutlier(Inj1_SNG$percent.mt)
  
  discard <- featureQC | libQC | mtQC
  Inj1_SNG <- Inj1_SNG[,!discard]
  summary(data.frame(byFeature = featureQC, byLibrary = libQC, byMT = mtQC, Discard = discard))
  
  
  featureQC <- isOutlier(Inj2_SNG$nFeature_RNA)
  libQC <- isOutlier(Inj2_SNG$nCount_RNA)
  mtQC <- isOutlier(Inj2_SNG$percent.mt)
  
  discard <- featureQC | libQC | mtQC
  Inj2_SNG <- Inj2_SNG[,!discard]
  summary(data.frame(byFeature = featureQC, byLibrary = libQC, byMT = mtQC, Discard = discard))
  
  featureQC <- isOutlier(Inj3_SNG$nFeature_RNA)
  libQC <- isOutlier(Inj3_SNG$nCount_RNA)
  mtQC <- isOutlier(Inj3_SNG$percent.mt)
  
  discard <- featureQC | libQC | mtQC
  Inj3_SNG <- Inj3_SNG[,!discard]
  summary(data.frame(byFeature = featureQC, byLibrary = libQC, byMT = mtQC, Discard = discard))
  
  featureQC <- isOutlier(Inj4_SNG$nFeature_RNA)
  libQC <- isOutlier(Inj4_SNG$nCount_RNA)
  mtQC <- isOutlier(Inj4_SNG$percent.mt)
  
  discard <- featureQC | libQC | mtQC
  Inj4_SNG <- Inj4_SNG[,!discard]
  summary(data.frame(byFeature = featureQC, byLibrary = libQC, byMT = mtQC, Discard = discard))
  
  featureQC <- isOutlier(Inj5_SNG$nFeature_RNA)
  libQC <- isOutlier(Inj5_SNG$nCount_RNA)
  mtQC <- isOutlier(Inj5_SNG$percent.mt)
  
  discard <- featureQC | libQC | mtQC
  Inj5_SNG <- Inj5_SNG[,!discard]
  summary(data.frame(byFeature = featureQC, byLibrary = libQC, byMT = mtQC, Discard = discard))
  
  
  ## ------------------------------------------------------------------------------------------------------------------------
  saveRDS(Inj1_SNG, file = './Results/Inj1.Rds')
  saveRDS(Inj2_SNG, file = './Results/Inj2.Rds')
  saveRDS(Inj3_SNG, file = './Results/Inj3.Rds')
  saveRDS(Inj4_SNG, file = './Results/Inj4.Rds')
  saveRDS(Inj5_SNG, file = './Results/Inj5.Rds')
  
  # save merged counts
  all_five <- merge(Inj1_SNG, y = c(Inj2_SNG, Inj3_SNG, Inj4_SNG, Inj5_SNG),
                    add.cell.ids = c("Inj1", "Inj2", "Inj3", "Inj4", "Inj5"))
  
  # Export the combined counts and metadata for deposition 
  all_five <- JoinLayers(all_five)
  
  counts <- GetAssayData(all_five, layer = "counts")
  gz <- gzfile("GSE211428_counts.tsv.gz", "w")
  write.table(as.matrix(counts), file = gz, sep = "\t", quote = FALSE)
  close(gz)
  
  gz <- gzfile("GSE211428_meta.tsv.gz", "w")
  write.table(all_five@meta.data, file = gz, sep = "\t", quote = FALSE)
  close(gz)
