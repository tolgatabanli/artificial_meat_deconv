# Characterising Cell Differentiation in Cultured Cattle Meat Through Cell-Type Deconvolution

This repository is for my Bachelor's Thesis on cell type deconvolution of bulk sequencing data from artificially cultured cattle meat.
Single-cell sequencing allows the characterisation of a sample on a cell level. However, the gene expression data should be filtered, and cells should be clustered and annotated with appropriate cell types,
and only then is it possible to quantify the proportion of cell types in samples.
The samples come from artificial meat cultured, i.e., (induced) stem cells from animals cultured *in vitro* to produce desired amounts of fat and muscle tissue.a
The thesis aims to characterise the amount of relevant cell types and their differentiation in adipogenesis and myogenesis experiments.

# Directory Structure
## Shared conventions across counts_* dataset directories

Each counts_cattle_atlas / counts_messmer / counts_wang2023 directory follows
the same pattern:

- Method scripts (dwls.R, rectangle_signature.py, scaden.R) run the
  corresponding deconvolution method (DWLS, Rectangle, Scaden) against the
  bulk RNA-seq samples, using that dataset's prepared single-cell reference.
- submit_*.sh — SLURM batch scripts that submit the matching method script
  as a cluster job.
- sig_* folders — hold each method's output: the computed cell-type
  signature matrix plus the resulting fraction estimates.

artificial_meat_deconv
├── converter_lineage.csv              #> Maps 90 cell types in the sampled Cattle Atlas to cell lineages
├── counts_cattle_atlas                #> Cattle Atlas processing, visualisation, and deconvolution
├── counts_messmer                     #> Messmer2023 processing, visualisation, and deconvolution
├── counts_wang2023                    #> Wang2023 processing, visualisation, and deconvolution
├── figures_agreement                  #> Method concordance: how well do estimates from different methods agree?
├── figures_fractions                  #> Cell-type fraction stacked barplots
├── figures_heatmapRaw                 #> Raw heatmaps with all cell types retained, detection limit applied
├── figures_posthoc                    #> Log-ratios of group contrasts across methods, references, and experiments
├── method_agreements.R                #> Sanity check: verifies that the cell-type gene signatures are functionally relevant
├── plot_cell_fractions.R              #> Generates the cell-type fraction stacked barplots
├── read_compositions_and_metadata.R   #> Utility script that reads fraction estimates from all references (counts_*) and all methods (DWLS, Rectangle, Scaden)
├── rectangle_environment.yml          #> Conda environment (Python packages) used for Rectangle
├── stats_comp.R                       #> Statistical analysis of compositions via log-ratios of cell types
└── toolkit.R                          #> Utility script: TPM-reading routine, cross-method composition-matrix merge function, a patched MAST function with caching, and a parallel-DWLS deconvolution function
