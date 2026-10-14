# Messmer

## counts_messmer

Scripts and folders suffixed with `_inj1` operate on the dataset containing
cells from Injection 1 only; all others use all injections. The thesis
focuses on Injection 1 only.

counts_messmer/
├── cytetyper_check.Rmd          #> Control check: compares Cytetyper's cell-type annotations against those from Messmer et al.
├── dwls_inj1.R
├── dwls.R
├── Figures
├── messmer_signature.Rmd
├── messmer_umap_reproduce.Rmd
├── preprocess.R                 #> Preprocessing script: a modified version of the original authors' script, adapted to retain Ensembl IDs
├── rectangle_signature_inj1.py
├── rectangle_signature.py
├── Results                      #> Preprocessing results: RDS files for the different injections and MTX-format files for Injection 1; available in the Zenodo repository
├── scaden.R
├── scaden_inj1.R
├── sig_dwls
├── sig_dwls_inj1
├── sig_dwls_withTeno
├── sig_dwls_withTeno_inj1
├── sig_rectangle
├── sig_rectangle_inj1
├── sig_rectangle_withTeno
├── sig_rectangle_withTeno_inj1
├── sig_scaden
├── sig_scaden_inj1
├── sig_scaden_withTeno
├── sig_scaden_withTeno_inj1
├── submit_dwls_inj1.sh
├── submit_dwls.sh
├── submit_scaden.sh
└── submit_scaden_inj1.sh
