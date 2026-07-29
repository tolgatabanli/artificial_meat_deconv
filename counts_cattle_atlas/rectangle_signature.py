#!/bin/bash
#SBATCH --job-name=cattle_rectangle
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=180G
#SBATCH --time=2:00:00
#SBATCH --partition=exbio-cpu
#SBATCH --output=job_%x_%j.out
#SBATCH --error=job_%x_%j.err
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=tolgatabanli@gmail.com

source /opt/mambaforge/etc/profile.d/conda.sh
conda activate deconv_rectangle
python3 << 'EOF'

# In[1]:


import rectanglepy as rectangle
import pandas as pd
import anndata as ad
import pickle
import scanpy as sc
import scipy.io
import mygene
import numpy as np


counts = scipy.io.mmread("sampled_cattle_atlas/counts.mtx").T.tocsr()
genes = pd.read_csv("sampled_cattle_atlas/genes.csv")
meta = pd.read_csv("sampled_cattle_atlas/metadata.csv", index_col=0)

#adata = ad.read_h5ad("sampled_cattle_atlas.h5ad")
adata = sc.AnnData(X=counts, obs=meta, var=genes.set_index("gene"))
adata.write_h5ad("sampled_cattle_atlas.h5ad")

# bulk matrices
bulk_adipo = pd.read_csv('/nfs/data2/huji_artificial_meat/adipogenesis/adipogenesis_tpm.tsv', index_col=0, sep='\t')
bulk_myo   = pd.read_csv('/nfs/data2/huji_artificial_meat/myogenesis/nfcore_out/star_salmon/salmon.merged.gene_tpm.tsv', index_col=0, sep='\t')
bulk_myo = bulk_myo.drop(columns=['gene_name'])


# cell filtering:
#bulk_genes = set(bulk_adipo.index)
#sc_var_names = np.array(adata.var_names)

#min_bulk_overlap = 50
#valid = []
#for ct in adata.obs['cell_type'].unique():
#    ct_mask = adata.obs['cell_type'] == ct
#    expressed = np.array(adata[ct_mask].X.sum(axis=0)).flatten() > 0
#    sc_genes = set(sc_var_names[expressed])
#    if len(sc_genes & bulk_genes) >= min_bulk_overlap:
#        valid.append(ct)

#dropped = set(adata.obs['cell_type'].unique()) - set(valid)
#print(f"Dropping {len(dropped)} cell types: {dropped}")
#adata_filtered = adata[adata.obs['cell_type'].isin(valid)].copy()

min_cells = 50
ct_counts = adata.obs['cell_type'].value_counts()
valid_ct = ct_counts[ct_counts >= min_cells].index
adata = adata[adata.obs['cell_type'].isin(valid_ct)].copy()
print(f"Kept {len(valid_ct)}/{len(ct_counts)} cell types, {adata.n_obs} cells")

# In[3]:


mg = mygene.MyGeneInfo()

var_names = adata.var_names.tolist()
is_ensembl = [g.startswith("ENSBTAG") for g in var_names]
symbols = [g for g, e in zip(var_names, is_ensembl) if not e]

def parse_hits(results):
    mapping = {}
    for hit in results:
        if hit.get("notfound"): continue
        ensembl = hit.get("ensembl")
        # list = multiple Ensembl IDs = asNA
        gene_id = ensembl.get("gene") if isinstance(ensembl, dict) else None
        if gene_id:
            mapping[hit["query"]] = gene_id
    return mapping

symbol_to_ensembl = parse_hits(
    mg.querymany(symbols, scopes="symbol", fields="ensembl.gene", species=9913)
)

# apply mapping
new_names = [symbol_to_ensembl.get(g, g) for g in var_names]

# keep only deduplicated ENSBTAGs (first occurrence wins, mirrors R duplicated())
seen = set()
keep = []
for i, name in enumerate(new_names):
    if name.startswith("ENSBTAG") and name not in seen:
        seen.add(name)
        keep.append(i)

n_already = sum(is_ensembl)
n_mapped = len(symbol_to_ensembl)
n_unmapped = len(symbols) - n_mapped
n_dup = sum(1 for n in new_names if n.startswith("ENSBTAG")) - len(keep)

adata = adata[:, keep]

print(f"ID mapping: {n_already}/{len(var_names)} already Ensembl IDs, "
      f"{n_mapped} gene symbols mapped uniquely to an Ensembl ID, "
      f"{n_unmapped} symbols could not be mapped and dropped. "
      f"{n_dup} ids were duplicated and dropped as well.")


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


est_adipo, signature_result = rectangle.rectangle(adata, bulk_adipo.T)


# In[1]:


# save results as pickle
with open('rectangle_results/rectangle_signature.pkl', 'wb') as f:
    pickle.dump(signature_result, f)

# save adipo estimations as CSV
est_adipo.reset_index().to_csv('rectangle_results/rectangle_estimations_adipo.csv', index=False)

# for myo deconv:
with open('rectangle_results/rectangle_signature.pkl', 'rb') as f:
    signature_result = pickle.load(f)

est_myo, _ = rectangle.tl.deconvolution(signature_result, bulk_myo.T)
# save myo estimations as CSV
est_myo.reset_index().to_csv('rectangle_results/rectangle_estimations_myo.csv', index=False)

EOF
