#!/bin/bash
#SBATCH --job-name=messmer_rectangle
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

messmer = pd.read_table('/nfs/data2/huji_artificial_meat/deconvolution_data/external_sc_references/messmer2023_GSE211428/GSE211428_counts.tsv.gz', index_col=0, sep='\t')
col_meta = pd.read_table('/nfs/data2/huji_artificial_meat/deconvolution_data/external_sc_references/messmer2023_GSE211428/GSE211428_combined_meta.tsv.gz', sep = '\t')

col_meta['seurat_barcode'] = col_meta['Unique_barcode'].str.replace(r'^(.*)\.([0-9]+)_([0-9]+)$', r'Inj\3_\1-\2', regex=True)
col_meta = col_meta.set_index('seurat_barcode').drop(columns=['V1'], errors='ignore')
col_meta['predicted.pheno'] = col_meta['predicted.pheno'].str.replace('/', '-')

print(f"mat columns: {messmer.shape[1]}")
print(f"col_meta barcodes: {len(col_meta)}")
common = messmer.columns.intersection(col_meta.index)
print(f"Intersection: {len(common)}")
print(f"In mat not in col_meta: {len(messmer.columns.difference(col_meta.index))}")
print(f"In col_meta not in mat: {len(col_meta.index.difference(messmer.columns))}")

messmer = messmer[common]
col_meta = col_meta.loc[common]

# remove tenocytes from the reference
tenocytes = col_meta[col_meta['predicted.pheno'] == 'Tenocytes']
col_meta = col_meta[~col_meta.index.isin(tenocytes.index)]
messmer = messmer.drop(columns=tenocytes.index)

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
sc_adata = ad.AnnData(messmer.T, obs=col_meta['predicted.pheno'].to_frame(name="cell_type"))
est_adipo, signature_result = rectangle.rectangle(sc_adata, bulk_adipo.T)


# save results as pickle
with open('sig_rectangle/rectangle_signature.pkl', 'wb') as f:
    pickle.dump(signature_result, f)


# for myo
est_myo, _ = rectangle.tl.deconvolution(signature_result, bulk_myo.T)


# also as CSV
est_adipo.reset_index().to_csv('sig_rectangle/rectangle_estimations_adipo.csv', index=False)
est_myo.reset_index().to_csv('sig_rectangle/rectangle_estimations_myo.csv', index=False)


EOF
