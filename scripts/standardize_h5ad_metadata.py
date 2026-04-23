from pathlib import Path
import re
import numpy as np
import pandas as pd
import scanpy as sc

IN_DIR = Path(r"D:\scRNA_benchmark\converted_h5ad")
OUT_DIR = Path(r"D:\scRNA_benchmark\standardized_h5ad")
OUT_DIR.mkdir(parents=True, exist_ok=True)

SUMMARY_ROWS = []


def infer_platform(dataset_id: str) -> str:
    x = dataset_id.lower()
    if "10x" in x:
        return "10x"
    if "celseq2" in x or "celseq" in x:
        return "CELseq2"
    if "dropseq" in x or "drop-seq" in x:
        return "Dropseq"
    if "rnamix" in x:
        return "plate_mix"
    if "9cellmix" in x:
        return "plate_mix"
    return "unknown"


def first_existing_col(df: pd.DataFrame, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None


def build_standard_fields(adata):
    obs = adata.obs.copy()
    dataset_id = str(obs["dataset_id"].iloc[0]) if "dataset_id" in obs.columns else "unknown_dataset"

    # -------------------------
    # 1. cell_id
    # -------------------------
    if "cell_id" not in obs.columns:
        obs["cell_id"] = obs.index.astype(str)

    # -------------------------
    # 2. platform
    # -------------------------
    if "platform" not in obs.columns:
        obs["platform"] = infer_platform(dataset_id)

    # -------------------------
    # 3. batch
    # -------------------------
    batch_col = first_existing_col(obs, ["batch", "Batch", "plate", "pool"])
    if batch_col is not None:
        obs["batch_std"] = obs[batch_col].astype(str)
    else:
        obs["batch_std"] = "unknown"

    # -------------------------
    # 4. celltype / cell identity
    #    优先选最“像标签”的列
    # -------------------------
    celltype_col = first_existing_col(
        obs,
        [
            "cell_line",
            "cell_line_demuxlet",
            "demuxlet_cls",
            "cell_type",
            "celltype",
            "CellType",
            "traj",
            "mix",
        ],
    )

    if celltype_col is not None:
        obs["celltype_std"] = obs[celltype_col].astype(str)
    else:
        # 对 9cellmix 这类数据，H1975/H2228/HCC827 三列可能是 membership/概率
        trio_cols = [c for c in ["H1975", "H2228", "HCC827"] if c in obs.columns]
        if len(trio_cols) == 3:
            tmp = obs[trio_cols].copy()

            def choose_celltype(row):
                # 若是布尔/0-1型，直接选最大者
                vals = pd.to_numeric(row, errors="coerce")
                if vals.notna().sum() > 0:
                    return vals.idxmax()
                return "unknown"

            obs["celltype_std"] = tmp.apply(choose_celltype, axis=1).astype(str)
        else:
            obs["celltype_std"] = "unknown"

    # -------------------------
    # 5. condition
    #    benchmark 时最重要，按数据集类型尽量自动推
    # -------------------------
    condition_col = first_existing_col(
        obs,
        [
            "condition",
            "Condition",
            "group",
            "Group",
            "mix",
            "mRNA_amount",
            "demuxlet_cls",
        ],
    )

    if condition_col is not None:
        obs["condition_std"] = obs[condition_col].astype(str)
    else:
        # 对 5cl / 10x 这种数据，先用 cell line 作为 condition 候选
        if "cell_line" in obs.columns:
            obs["condition_std"] = obs["cell_line"].astype(str)
        elif "cell_line_demuxlet" in obs.columns:
            obs["condition_std"] = obs["cell_line_demuxlet"].astype(str)
        elif all(c in obs.columns for c in ["H1975", "H2228", "HCC827"]):
            trio_cols = ["H1975", "H2228", "HCC827"]

            def choose_condition(row):
                vals = pd.to_numeric(row[trio_cols], errors="coerce")
                if vals.notna().sum() > 0:
                    return vals.idxmax()
                return "unknown"

            obs["condition_std"] = obs.apply(choose_condition, axis=1).astype(str)
        else:
            obs["condition_std"] = "unknown"

    # -------------------------
    # 6. sample
    #    没有 donor/sample 时，构造一个“最小可用 sample”
    # -------------------------
    sample_col = first_existing_col(
        obs,
        [
            "sample",
            "Sample",
            "sample_id",
            "batch",
            "cell_number",
            "p_column",
            "p_row",
        ],
    )

    if sample_col is not None:
        obs["sample_std"] = obs[sample_col].astype(str)
    else:
        # 按数据集场景构造 sample
        if "batch_std" in obs.columns and obs["batch_std"].nunique() > 1:
            obs["sample_std"] = (
                obs["platform"].astype(str) + "_" + obs["batch_std"].astype(str)
            )
        elif all(c in obs.columns for c in ["p_row", "p_column"]):
            obs["sample_std"] = (
                obs["platform"].astype(str)
                + "_"
                + obs["p_row"].astype(str)
                + "_"
                + obs["p_column"].astype(str)
            )
        else:
            obs["sample_std"] = obs["platform"].astype(str) + "_sample1"

    # -------------------------
    # 7. 清理一些值
    # -------------------------
    for c in ["celltype_std", "condition_std", "sample_std", "batch_std", "platform"]:
        obs[c] = obs[c].astype(str).replace({"nan": "unknown", "None": "unknown"})

    adata.obs = obs

    # 统一再补几个标准列名，后续 pipeline 直接用
    adata.obs["celltype"] = adata.obs["celltype_std"]
    adata.obs["condition"] = adata.obs["condition_std"]
    adata.obs["sample"] = adata.obs["sample_std"]
    adata.obs["batch"] = adata.obs["batch_std"]

    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    if "gene_id" not in adata.var.columns:
        adata.var["gene_id"] = adata.var_names.astype(str)

    if "gene_symbol" not in adata.var.columns:
        adata.var["gene_symbol"] = adata.var_names.astype(str)

    return adata


def summarize_adata(adata, out_name: str):
    obs = adata.obs
    row = {
        "file": out_name,
        "dataset_id": obs["dataset_id"].iloc[0] if "dataset_id" in obs.columns else out_name.replace(".h5ad", ""),
        "n_cells": adata.n_obs,
        "n_genes": adata.n_vars,
        "condition_n": obs["condition"].nunique() if "condition" in obs.columns else np.nan,
        "sample_n": obs["sample"].nunique() if "sample" in obs.columns else np.nan,
        "celltype_n": obs["celltype"].nunique() if "celltype" in obs.columns else np.nan,
        "batch_n": obs["batch"].nunique() if "batch" in obs.columns else np.nan,
        "platform_n": obs["platform"].nunique() if "platform" in obs.columns else np.nan,
        "condition_levels": ";".join(map(str, sorted(obs["condition"].unique()[:20]))) if "condition" in obs.columns else "",
        "celltype_levels": ";".join(map(str, sorted(obs["celltype"].unique()[:20]))) if "celltype" in obs.columns else "",
    }
    SUMMARY_ROWS.append(row)


def main():
    files = sorted(IN_DIR.glob("*.h5ad"))
    print(f"检测到 {len(files)} 个 h5ad 文件。")

    for f in files:
        print("=" * 80)
        print("读取:", f.name)
        adata = sc.read_h5ad(f)

        adata = build_standard_fields(adata)

        out_path = OUT_DIR / f.name
        adata.write_h5ad(out_path)

        summarize_adata(adata, f.name)

        print(adata)
        print("condition:", adata.obs["condition"].value_counts().head().to_dict())
        print("sample:", adata.obs["sample"].value_counts().head().to_dict())
        print("celltype:", adata.obs["celltype"].value_counts().head().to_dict())
        print("已保存:", out_path)

    summary_df = pd.DataFrame(SUMMARY_ROWS)
    summary_path = OUT_DIR / "dataset_summary_standardized.csv"
    summary_df.to_csv(summary_path, index=False)
    print("\n标准化摘要已保存到：", summary_path)
    print(summary_df)


if __name__ == "__main__":
    main()