# Cattle Atlas

Catte atlas had originally ~114 cell types.
I used a sampling procedure to reduce the functional size of RDS objects to be able to merge the atlas.
From each cell type, 500 cells were sampled in total. Cell types with less than 500 cells were dropped out.
The sampled and merged cattle atlas includes 90 cell types over 37 tissues.

```
counts_cattle_atlas/
├── benchmark                   #> Work in progress, not complete
├── cattle_atlas.Rmd             #> General exploration: UMAP,tion plots
├── cell_type_lineage_map.csv    #> Maps fine-grained cell types to cell lineages
├── dwls.R                       #> DWLS script
├── figures                      #> Output figures
├── rectangle_signature.py       #> Rectangle script
├── sample_and_merge.R           #> Samples 500 cells per celll reference used for deconvolution
├── sampled_cattle_atlas         #> Sampled output (MTX file + metadata)
├── sampled_cattle_atlas.h5ad    #> On-disk (h5ad) format, for
├── sampled_cattle_atlas.RDS     #> RDS format, for R calls
├── scaden.R                     #> Scaden script
├── sig_dwls_fine                #> DWLS signature and fraction estimates
├── sig_rectangle_fine           #> Rectangle signature and fraction estimates
├── sig_scaden_fine              #> Scaden signature and fraction estimates
├── submit_dwls.sh               #> SLURM job for dwls.R
├── submit_merge.sh              #> SLURM job for sample_and_merge.R
├── submit_scaden.sh             #> SLURM job for scaden.R
└── tissue_by_ct.csv             #> Number of cells each tissue contributes per cell type
```
