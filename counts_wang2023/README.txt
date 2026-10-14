# Wang2023

counts_wang2023/
├── dwls.R
├── figures
├── plots.Rmd                 #> For cell-type abundance and UMAP plots
├── rectangle_signature.py
├── scaden.R
├── sig_dwls                  #> sig_* folders contain signatures and fraction estimates; available on the Zenodo repo
├── sig_rectangle
├── sig_scaden
├── submit_dwls.sh
├── submit_scaden.sh
├── wang2023.R                #> Preprocessing: reads raw counts, filters using the normalized data's metadata, maps IDs, and writes output
├── wang_ids_mapped           #> MTX files after ID mapping
├── wang_ids_mapped.h5ad      #> On-disk format of the reference dataset after ID mapping
└── wang_ids_mapped.RDS       #> RDS format of the reference dataset after ID mapping
