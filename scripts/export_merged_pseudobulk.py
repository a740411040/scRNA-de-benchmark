from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

INPUT_H5AD = Path(r"D:\scRNA_benchmark\standardized_h5ad\SCMixology_5cl_CELseq2_merged.h5ad")
OUT_DIR = Path(r"D:\scRNA_benchmark\benchmark_runs\SCMixology_5cl_CELseq2_merged")
OUT_DIR.mkdir(parents=True, exist_ok=True)

adata = sc.read_h5ad(INPUT_H5AD)

# 确保 counts 层存在
if "counts" not in adata.layers:
    adata.layers["counts"] = adata.X.copy()

# 统一 condition
if "cell_line_demuxlet" in adata.obs.columns:
    adata.obs["condition"] = adata.obs["cell_line_demuxlet"].astype(str)
elif "cell_line" in adata.obs.columns:
    adata.obs["condition"] = adata.obs["cell_line"].astype(str)
else:
    raise ValueError("找不到 condition 列（cell_line_demuxlet / cell_line）")

# 统一 batch
if "batch" not in adata.obs.columns:
    raise ValueError("merged 数据中没有 batch 列，无法做协变量比较。")

adata.obs["batch"] = adata.obs["batch"].astype(str)

# 这里构造 pseudo-bulk sample：condition + batch
adata.obs["pb_sample"] = (
    adata.obs["condition"].astype(str) + "_" + adata.obs["batch"].astype(str)
)

counts = adata.layers["counts"]
if sp.issparse(counts):
    counts = counts.tocsr()

sample_ids = adata.obs["pb_sample"].values
unique_samples = pd.unique(sample_ids)

pb_mat = []
pb_obs = []

for sid in unique_samples:
    idx = np.where(sample_ids == sid)[0]
    sub_counts = counts[idx]
    summed = np.asarray(sub_counts.sum(axis=0)).ravel()
    pb_mat.append(summed)

    pb_obs.append({
        "sample": sid,
        "condition": str(adata.obs.iloc[idx[0]]["condition"]),
        "batch": str(adata.obs.iloc[idx[0]]["batch"]),
        "n_cells": len(idx),
    })

pb_mat = np.vstack(pb_mat)
pb_obs = pd.DataFrame(pb_obs)
pb_counts = pd.DataFrame(pb_mat, index=pb_obs["sample"], columns=adata.var_names)

pb_counts.to_csv(OUT_DIR / "pseudobulk_counts.csv")
pb_obs.to_csv(OUT_DIR / "pseudobulk_metadata.csv", index=False)

print("saved:")
print(OUT_DIR / "pseudobulk_counts.csv")
print(OUT_DIR / "pseudobulk_metadata.csv")
print("\nmetadata:")
print(pb_obs)