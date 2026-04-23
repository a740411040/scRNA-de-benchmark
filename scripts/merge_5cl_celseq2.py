from pathlib import Path
import scanpy as sc
import anndata as ad
import pandas as pd
import scipy.io
import scipy.sparse as sp

IN_DIR = Path(r"D:\scRNA_benchmark\standardized_h5ad")
OUT_H5AD = Path(r"D:\scRNA_benchmark\standardized_h5ad\SCMixology_5cl_CELseq2_merged.h5ad")
BRIDGE_DIR = Path(r"D:\scRNA_benchmark\bridge")
BRIDGE_DIR.mkdir(parents=True, exist_ok=True)

files = [
    ("SCMixology_5cl_CELseq2_p1.h5ad", "p1"),
    ("SCMixology_5cl_CELseq2_p2.h5ad", "p2"),
    ("SCMixology_5cl_CELseq2_p3.h5ad", "p3"),
]

adatas = []
for fname, batch in files:
    adata = sc.read_h5ad(IN_DIR / fname)
    adata.obs["batch"] = batch
    adata.obs["dataset_id"] = "SCMixology_5cl_CELseq2_merged"
    adatas.append(adata)

adata_merged = ad.concat(adatas, join="outer", merge="same", label=None, index_unique=None)
adata_merged.layers["counts"] = adata_merged.layers["counts"].copy() if "counts" in adata_merged.layers else adata_merged.X.copy()
adata_merged.write_h5ad(OUT_H5AD)

print(adata_merged)
print("saved:", OUT_H5AD)

# 同时导出 bridge files 给 R 用
X = adata_merged.layers["counts"]
if sp.issparse(X):
    X = X.tocsr()
else:
    X = sp.csr_matrix(X)

scipy.io.mmwrite(BRIDGE_DIR / "SCMixology_5cl_CELseq2_merged_counts.mtx", X.T)

obs = adata_merged.obs.copy()
obs["cell_id"] = adata_merged.obs_names.astype(str)
obs.to_csv(BRIDGE_DIR / "SCMixology_5cl_CELseq2_merged_obs.csv", index=False)

var = adata_merged.var.copy()
var["gene_id"] = adata_merged.var_names.astype(str)
var.to_csv(BRIDGE_DIR / "SCMixology_5cl_CELseq2_merged_var.csv", index=False)

pd.DataFrame({"gene_id": adata_merged.var_names.astype(str)}).to_csv(
    BRIDGE_DIR / "SCMixology_5cl_CELseq2_merged_genes.csv", index=False
)
pd.DataFrame({"cell_id": adata_merged.obs_names.astype(str)}).to_csv(
    BRIDGE_DIR / "SCMixology_5cl_CELseq2_merged_cells.csv", index=False
)

print("bridge exported.")