from pathlib import Path
from itertools import combinations

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

IN_DIR = Path(r"D:\scRNA_benchmark\standardized_h5ad")
OUT_ROOT = Path(r"D:\scRNA_benchmark\benchmark_runs")
OUT_ROOT.mkdir(parents=True, exist_ok=True)

# 第一批建议先跑这些
TARGET_FILES = [
    "SCMixology_5cl_CELseq2_p1.h5ad",
    "SCMixology_5cl_CELseq2_p2.h5ad",
    "SCMixology_5cl_CELseq2_p3.h5ad",
    "SCMixology_RNAmix_sce2.h5ad",
    "SCMixology_RNAmix_sce8.h5ad",
]

SUMMARY = []


def choose_condition_col(adata):
    # 优先更可信的细胞系标签
    for c in ["cell_line_demuxlet", "cell_line", "condition"]:
        if c in adata.obs.columns:
            return c
    return "condition"


def make_sample_col(adata):
    obs = adata.obs.copy()

    if "sample" in obs.columns and obs["sample"].nunique() > 1:
        return adata

    # 尽量构造一个最小可用 sample
    if "batch" in obs.columns and obs["batch"].nunique() > 1:
        obs["sample"] = obs["platform"].astype(str) + "_" + obs["batch"].astype(str)
    elif "cell_number" in obs.columns:
        obs["sample"] = obs["platform"].astype(str) + "_" + obs["cell_number"].astype(str)
    else:
        obs["sample"] = obs["platform"].astype(str) + "_sample1"

    adata.obs = obs
    return adata


def preprocess_adata(adata, target_n_cells=5000):
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    # mt 基因
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    # 基础过滤
    if "n_genes_by_counts" in adata.obs.columns:
        adata = adata[adata.obs["n_genes_by_counts"] >= 200].copy()
    if "pct_counts_mt" in adata.obs.columns:
        adata = adata[adata.obs["pct_counts_mt"] <= 20].copy()

    sc.pp.filter_genes(adata, min_cells=3)

    # 下采样
    if adata.n_obs > target_n_cells:
        np.random.seed(1234)
        keep_idx = np.random.choice(adata.n_obs, target_n_cells, replace=False)
        adata = adata[keep_idx].copy()

    # normalize + log1p
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["lognorm"] = adata.X.copy()

    # HVG
    sc.pp.highly_variable_genes(adata, flavor="cell_ranger", n_top_genes=3000)
    adata = adata[:, adata.var["highly_variable"]].copy()

    # PCA / UMAP / Leiden
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.5)

    return adata


def run_scanpy_pairwise(adata, cond_col, group1, group2, outdir: Path):
    tmp = adata[adata.obs[cond_col].isin([group1, group2])].copy()
    tmp.obs[cond_col] = pd.Categorical(
        tmp.obs[cond_col].astype(str),
        categories=[group1, group2]
    )
    tmp.X = tmp.layers["lognorm"].copy()

    comparison = f"{group1}_vs_{group2}"

    # Wilcoxon
    sc.tl.rank_genes_groups(
        tmp,
        groupby=cond_col,
        groups=[group1],
        reference=group2,
        method="wilcoxon",
        corr_method="benjamini-hochberg",
        tie_correct=True,
        pts=True,
    )
    res_wil = sc.get.rank_genes_groups_df(tmp, group=group1)
    res_wil["gene"] = res_wil["names"].astype(str)
    res_wil["dataset"] = str(adata.obs["dataset_id"].iloc[0])
    res_wil["method"] = "scanpy_wilcoxon"
    res_wil["comparison"] = comparison
    res_wil.to_csv(outdir / f"scanpy_wilcoxon_{comparison}.csv", index=False)

    # t-test
    sc.tl.rank_genes_groups(
        tmp,
        groupby=cond_col,
        groups=[group1],
        reference=group2,
        method="t-test",
        corr_method="benjamini-hochberg",
        pts=True,
    )
    res_t = sc.get.rank_genes_groups_df(tmp, group=group1)
    res_t["gene"] = res_t["names"].astype(str)
    res_t["dataset"] = str(adata.obs["dataset_id"].iloc[0])
    res_t["method"] = "scanpy_ttest"
    res_t["comparison"] = comparison
    res_t.to_csv(outdir / f"scanpy_ttest_{comparison}.csv", index=False)

    return res_wil, res_t


def run_scanpy_all_pairwise(adata, cond_col, outdir: Path):
    levels = sorted(adata.obs[cond_col].astype(str).unique())
    if len(levels) < 2:
        return None

    pair_list = list(combinations(levels, 2))
    all_rows = []

    for group1, group2 in pair_list:
        print(f"Scanpy pairwise: {group1} vs {group2}")
        res_wil, res_t = run_scanpy_pairwise(adata, cond_col, group1, group2, outdir)

        all_rows.append({
            "comparison": f"{group1}_vs_{group2}",
            "wilcoxon_n": res_wil.shape[0],
            "ttest_n": res_t.shape[0],
        })

    pair_df = pd.DataFrame(all_rows)
    pair_df.to_csv(outdir / "scanpy_pairwise_summary.csv", index=False)
    return pair_df


def export_pseudobulk(adata, outdir: Path):
    obs = adata.obs.copy()

    # 为不同数据集定制更细的 pseudo-bulk sample
    dataset_name = str(obs["dataset_id"].iloc[0]) if "dataset_id" in obs.columns else "unknown"

    # 默认保证 condition / celltype 存在
    if "condition" not in obs.columns:
        obs["condition"] = "unknown"
    if "celltype" not in obs.columns:
        obs["celltype"] = obs["condition"].astype(str)

    # 1) CELseq2 五细胞系：用 celltype + cell_number
    if dataset_name.startswith("SCMixology_5cl_CELseq2") and "cell_number" in obs.columns:
        obs["pb_sample"] = (
            obs["celltype"].astype(str) + "_" + obs["cell_number"].astype(str)
        )

    # 2) RNAmix：用板位 p_row + p_column
    elif dataset_name.startswith("SCMixology_RNAmix") and all(c in obs.columns for c in ["p_row", "p_column"]):
        obs["pb_sample"] = (
            obs["p_row"].astype(str) + "_" + obs["p_column"].astype(str)
        )

    # 3) 如果 batch 存在且 batch 不是单一值，用 condition + batch
    elif "batch" in obs.columns and obs["batch"].astype(str).nunique() > 1:
        obs["pb_sample"] = (
            obs["condition"].astype(str) + "_" + obs["batch"].astype(str)
        )

    # 4) 如果 sample 本身就有多个值，用 condition + sample
    elif "sample" in obs.columns and obs["sample"].astype(str).nunique() > 1:
        obs["pb_sample"] = (
            obs["condition"].astype(str) + "_" + obs["sample"].astype(str)
        )

    # 5) 兜底：无法构造有效 pseudo-bulk
    else:
        obs["pb_sample"] = obs["condition"].astype(str) + "_sample1"

    adata.obs = obs

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
            "celltype": str(adata.obs.iloc[idx[0]]["celltype"]),
            "batch": str(adata.obs.iloc[idx[0]].get("batch", "unknown")),
            "platform": str(adata.obs.iloc[idx[0]].get("platform", "unknown")),
            "n_cells": len(idx),
        })

    pb_mat = np.vstack(pb_mat)
    pb_obs = pd.DataFrame(pb_obs)
    pb_counts = pd.DataFrame(pb_mat, index=pb_obs["sample"], columns=adata.var_names)

    pb_counts.to_csv(outdir / "pseudobulk_counts.csv")
    pb_obs.to_csv(outdir / "pseudobulk_metadata.csv", index=False)

    return pb_counts, pb_obs


def main():
    for fname in TARGET_FILES:
        fpath = IN_DIR / fname
        if not fpath.exists():
            print(f"跳过，不存在: {fpath}")
            continue

        print("=" * 90)
        print("处理:", fname)
        adata = sc.read_h5ad(fpath)

        dataset_name = fpath.stem
        outdir = OUT_ROOT / dataset_name
        outdir.mkdir(parents=True, exist_ok=True)

        cond_col = choose_condition_col(adata)

        # 把标准 condition 同步成更可信列
        adata.obs["condition"] = adata.obs[cond_col].astype(str)

        # 如果没 celltype，则先用 condition 占位
        if "celltype" not in adata.obs.columns:
            adata.obs["celltype"] = adata.obs["condition"].astype(str)

        adata = preprocess_adata(adata)
        adata = make_sample_col(adata)

        # condition 至少要有两组
        if adata.obs["condition"].astype(str).nunique() < 2:
            print("分组不足，跳过 DE:", fname)
            continue

        # 保存处理后对象
        adata.write_h5ad(outdir / "processed.h5ad")

        # Scanpy pairwise DE
        pair_df = run_scanpy_all_pairwise(adata, "condition", outdir)
        if pair_df is None or pair_df.empty:
            print("没有可运行的 pairwise 比较:", fname)
            continue

        print(pair_df)

        # pseudo-bulk
        pb_counts, pb_obs = export_pseudobulk(adata, outdir)

        SUMMARY.append({
            "dataset": dataset_name,
            "n_cells_after_qc": adata.n_obs,
            "n_genes_after_hvg": adata.n_vars,
            "condition_col_used": cond_col,
            "condition_n": adata.obs["condition"].nunique(),
            "sample_n": adata.obs["sample"].nunique(),
            "celltype_n": adata.obs["celltype"].nunique(),
            "pseudobulk_n": pb_obs.shape[0],
            "scanpy_pairwise_n": pair_df.shape[0],
        })

    summary_df = pd.DataFrame(SUMMARY)
    summary_df.to_csv(OUT_ROOT / "benchmark_batch_summary.csv", index=False)
    print("\n已保存:", OUT_ROOT / "benchmark_batch_summary.csv")
    print(summary_df)


if __name__ == "__main__":
    main()