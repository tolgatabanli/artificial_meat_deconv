#!/bin/bash
#SBATCH --job-name=messmer_rectangle_inj1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=180G
#SBATCH --time=8:00:00
#SBATCH --partition=exbio-cpu
#SBATCH --output=job_%x_%j.out
#SBATCH --error=job_%x_%j.err
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=tolgatabanli@gmail.com

source /opt/mambaforge/etc/profile.d/conda.sh
conda activate deconv_rectangle
python3 << 'EOF'

import rectanglepy as rectangle
import pandas as pd
import anndata as ad
import pickle
import scipy.io
import os
os.chdir(os.path.expanduser("~/artifical_meat_deconv/messmer_counts"))

counts = scipy.io.mmread("Results/Inj1_counts.mtx").T.tocsr() # -> cells x genes
genes = pd.read_csv("Results/Inj1_genes.csv")["x"].values
barcodes = pd.read_csv("Results/Inj1_barcodes.csv")["x"].values
meta = pd.read_csv("Results/Inj1_meta.csv", index_col=0)

sc_adata = ad.AnnData(counts)
sc_adata.obs_names = barcodes
sc_adata.var_names = genes

# remove tenocytes from the reference
meta = meta[meta['predicted.pheno'] != 'Tenocytes']

common = sc_adata.obs_names.intersection(meta.index)
sc_adata = sc_adata[common]
sc_adata.obs["cell_type"] = meta.loc[common, "predicted.pheno"].values

# bulk matrices
bulk_adipo = pd.read_csv('/nfs/data2/huji_artificial_meat/adipogenesis/adipogenesis_tpm.tsv', index_col=0, sep='\t')
bulk_myo = pd.read_csv('/nfs/data2/huji_artificial_meat/myogenesis/nfcore_out/star_salmon/salmon.merged.gene_tpm.tsv', index_col=0, sep='\t')
bulk_myo = bulk_myo.drop(columns=['gene_name'])

# metadata
meta_adipo = pd.read_csv('/nfs/data2/huji_artificial_meat/adipogenesis/sample_sheet_cleaned.csv', index_col=0)
meta_myo = pd.read_csv('/nfs/data2/huji_artificial_meat/myogenesis/sample_sheet_cleaned.csv', index_col=0)

meta_adipo.columns = meta_adipo.columns.str.lower().str.replace(r'[^a-z0-9]', '_', regex=True).str.strip('_')
meta_myo.columns = meta_myo.columns.str.lower().str.replace(r'[^a-z0-9]', '_', regex=True).str.strip('_')

meta_adipo = meta_adipo[['library_name', 'sample_id', 'title', 'group', 'bio_replicate', 'individual', 'treatment', 'day']]
meta_adipo['day'] = pd.Categorical(meta_adipo['day'], categories=['D0', 'D2', 'D6', 'D14'], ordered=True)
meta_adipo['short_id'] = meta_adipo['library_name'].str.replace("Sample ", "S", regex=True)

meta_myo = meta_myo[['library_name', 'sample_id', 'group', 'day', 'bio_replicate', 'treatment']]
meta_myo['day'] = pd.Categorical(meta_myo['day'], categories=['d0', 'd3', 'd7', 'd21'], ordered=True)
meta_myo['group'] = meta_myo['group'].str.replace(r'(?!Sat)D\d*', '', regex=True).str.split('_').str[-1]

# adjust the sample IDs so that it's only Sxx type
bulk_myo.columns = bulk_myo.columns.str.replace(r'.*_', '', regex=True)
meta_myo['sample_id'] = meta_myo['sample_id'].str.replace(r'.*_', '', regex=True)
bulk_myo = bulk_myo.loc[:, bulk_myo.columns.isin(meta_myo['sample_id'])]

# build signatures
est_adipo, sig = rectangle.rectangle(sc_adata, bulk_adipo.T)
est_adipo, signature_result = rectangle.rectangle(sc_adata, bulk_adipo.T)


# save results as pickle
with open('sig_rectangle_inj1/rectangle_signature.pkl', 'wb') as f:
    pickle.dump(signature_result, f)


# for myo
est_myo, _ = rectangle.tl.deconvolution(signature_result, bulk_myo.T)


# also as CSV
est_adipo.reset_index().to_csv('sig_rectangle_inj1/rectangle_estimations_adipo.csv', index=False)
est_myo.reset_index().to_csv('sig_rectangle_inj1/rectangle_estimations_myo.csv', index=False)


EOF
