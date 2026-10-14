#!/bin/bash
#SBATCH --job-name=wang_rectangle
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=40G
#SBATCH --time=2:00:00
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
import scanpy as sc
import scipy.io
import mygene
import numpy as np


counts = scipy.io.mmread("wang_ids_mapped/counts.mtx").T.tocsr()
genes = pd.read_csv("wang_ids_mapped/genes.csv")
meta = pd.read_csv("wang_ids_mapped/metadata.csv", index_col=0)
meta['lables'] = meta['lables'].str.replace('/', '_')
meta.loc[meta['lables'] == "Muc-2", 'lables'] = "MuC-2"

#adata = ad.read_h5ad("wang_ids_mapped.h5ad")
adata = sc.AnnData(X=counts, obs=meta, var=genes.set_index("gene"))
adata.write_h5ad("wang_ids_mapped.h5ad")


# bulk matrices
bulk_adipo = pd.read_csv('/nfs/data2/huji_artificial_meat/adipogenesis/adipogenesis_tpm.tsv', index_col=0, sep='\t')
bulk_myo   = pd.read_csv('/nfs/data2/huji_artificial_meat/myogenesis/nfcore_out/star_salmon/salmon.merged.gene_tpm.tsv', index_col=0, sep='\t')
bulk_myo = bulk_myo.drop(columns=['gene_name'])


# metadata
meta_adipo = pd.read_csv('/nfs/data2/huji_artificial_meat/adipogenesis/sample_sheet_cleaned.csv', index_col=0)
meta_myo   = pd.read_csv('/nfs/data2/huji_artificial_meat/myogenesis/sample_sheet_cleaned.csv', index_col=0)

meta_adipo.columns = meta_adipo.columns.str.lower().str.replace(r'[^a-z0-9]', '_', regex=True).str.strip('_')
meta_myo.columns   = meta_myo.columns.str.lower().str.replace(r'[^a-z0-9]', '_', regex=True).str.strip('_')

meta_adipo = meta_adipo[['library_name', 'sample_id', 'title', 'group', 'bio_replicate', 'individual', 'treatment', 'day']]
meta_adipo['day'] = pd.Categorical(meta_adipo['day'], categories=['D0', 'D2', 'D6', 'D14'], ordered=True)
meta_adipo['short_id'] = meta_adipo['library_name'].str.replace("Sample ", "S", regex=True)

meta_myo = meta_myo[['library_name', 'sample_id', 'group', 'day', 'bio_replicate', 'treatment']]
meta_myo['day'] = pd.Categorical(meta_myo['day'], categories=['d0', 'd3', 'd7', 'd21'], ordered=True)
meta_myo['group'] = meta_myo['group'].str.replace(r'(?!Sat)D\d*', '', regex=True).str.split('_').str[-1]


# In[ ]:


est_adipo, signature_result = rectangle.rectangle(adata, bulk_adipo.T, cell_type_col='lables')


# In[1]:


# save results as pickle
with open('sig_rectangle/rectangle_signature.pkl', 'wb') as f:
    pickle.dump(signature_result, f)

# save adipo estimations as CSV
est_adipo.reset_index().to_csv('sig_rectangle/rectangle_estimations_adipo.csv', index=False)

# for myo deconv:
with open('sig_rectangle/rectangle_signature.pkl', 'rb') as f:
    signature_result = pickle.load(f)

est_myo, _ = rectangle.tl.deconvolution(signature_result, bulk_myo.T)
# save myo estimations as CSV
est_myo.reset_index().to_csv('sig_rectangle/rectangle_estimations_myo.csv', index=False)

EOF
