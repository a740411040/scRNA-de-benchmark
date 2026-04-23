from pathlib import Path
import numpy as np
import scanpy as sc

INPUT_H5AD = Path(r"D:\scRNA_benchmark\standardized_h5ad\SCMixology_5cl_10x.h5ad")
OUT_DIR = Path(r"D:\scRNA_benchmark\pipeline_output\SCMixology_5cl_10x")
OUT_DIR.mkdir(parents=True, exist_ok=True)

adata = sc.read_h5ad(INPUT_H5AD)

# 保存 raw counts
if "counts" not in adata.layers:
    adata.layers["counts"] = adata.X.copy()

# 基础 QC 指标
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

# 基础过滤
if "n_genes_by_counts" in adata.obs.columns:
    adata = adata[adata.obs["n_genes_by_counts"] >= 200].copy()
if "pct_counts_mt" in adata.obs.columns:
    adata = adata[adata.obs["pct_counts_mt"] <= 20].copy()

sc.pp.filter_genes(adata, min_cells=3)

# 下采样到 5000 cells
target_n_cells = 5000
if adata.n_obs > target_n_cells:
    np.random.seed(1234)
    keep_idx = np.random.choice(adata.n_obs, target_n_cells, replace=False)
    adata = adata[keep_idx].copy()

# Normalize + log1p
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers["lognorm"] = adata.X.copy()

# HVG
sc.pp.highly_variable_genes(adata, flavor="cell_ranger", n_top_genes=3000)
adata = adata[:, adata.var["highly_variable"]].copy()

# PCA / neighbors / UMAP / Leiden
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)

# 保存
adata.write_h5ad(OUT_DIR / "processed.h5ad")

# 画图
sc.pl.umap(adata, color=["condition", "celltype", "leiden"], save=False, show=False)