from pathlib import Path
import pandas as pd
import re

ROOT = Path(r"D:\scRNA_benchmark\benchmark_runs")
OUT = Path(r"D:\scRNA_benchmark\benchmark_results")
OUT.mkdir(parents=True, exist_ok=True)

rows = []

# --------------------------------------------------
# 工具函数：统一 gene 列
# --------------------------------------------------
def ensure_gene_column(df: pd.DataFrame):
    out = df.copy()

    if "gene" in out.columns:
        out["gene"] = out["gene"].astype(str)
        return out

    if "names" in out.columns:
        out["gene"] = out["names"].astype(str)
        return out

    return None


# --------------------------------------------------
# 1. 收集 Scanpy pairwise 结果
# 文件名格式：
#   scanpy_ttest_A549_vs_H1975.csv
#   scanpy_wilcoxon_A549_vs_H1975.csv
# --------------------------------------------------
for dataset_dir in ROOT.iterdir():
    if not dataset_dir.is_dir():
        continue

    dataset = dataset_dir.name

    for f in dataset_dir.glob("scanpy_*.csv"):
        if f.name == "scanpy_pairwise_summary.csv":
            continue

        m = re.match(r"^(scanpy_ttest|scanpy_wilcoxon)_(.+)\.csv$", f.name)
        if not m:
            continue

        method = m.group(1)
        comparison = m.group(2)

        df = pd.read_csv(f)
        df = ensure_gene_column(df)
        if df is None:
            print("Scanpy 文件没有 gene/names 列，跳过:", f)
            continue

        sub = df.copy()
        sub["dataset"] = dataset
        sub["method"] = method
        sub["comparison"] = comparison
        rows.append(sub)

        print(f"[Scanpy] 已收集: {dataset} | {method} | {comparison}")


# --------------------------------------------------
# 2. 收集 R pseudo-bulk pairwise 结果
# 目录：
#   r_pseudobulk_de_pairwise
# 文件名格式：
#   edgeR_A549_vs_H1975.csv
#   DESeq2_A549_vs_H1975.csv
#   limma_voom_A549_vs_H1975.csv
# --------------------------------------------------
for dataset_dir in ROOT.iterdir():
    if not dataset_dir.is_dir():
        continue

    dataset = dataset_dir.name
    rdir = dataset_dir / "r_pseudobulk_de_pairwise"
    if not rdir.exists():
        continue

    print(f"\n检查 pseudo-bulk 结果目录: {rdir}")

    for f in rdir.glob("*.csv"):
        if f.name in ["pairwise_summary.csv", "pseudobulk_metadata_filtered.csv"]:
            continue

        m = re.match(r"^(edgeR|DESeq2|limma_voom)_(.+)\.csv$", f.name)
        if not m:
            print("  未匹配到 pseudo-bulk 结果文件模式，跳过:", f.name)
            continue

        method = m.group(1)
        comparison = m.group(2)

        df = pd.read_csv(f)
        df = ensure_gene_column(df)
        if df is None:
            print("  没有 gene/names 列，跳过:", f.name)
            continue

        sub = df.copy()
        sub["dataset"] = dataset
        sub["method"] = method
        sub["comparison"] = comparison
        rows.append(sub)

        print(f"  已收集: {dataset} | {method} | {comparison}")


# --------------------------------------------------
# 3. 收集 R Seurat / MAST pairwise 结果
# 目录：
#   r_seurat_mast_pairwise
# 文件名格式：
#   seurat_wilcoxon_A549_vs_H1975.csv
#   MAST_A549_vs_H1975.csv
# --------------------------------------------------
for dataset_dir in ROOT.iterdir():
    if not dataset_dir.is_dir():
        continue

    dataset = dataset_dir.name
    rdir = dataset_dir / "r_seurat_mast_pairwise"
    if not rdir.exists():
        continue

    print(f"\n检查 Seurat/MAST 结果目录: {rdir}")

    for f in rdir.glob("*.csv"):
        if f.name == "seurat_mast_pairwise_summary.csv":
            continue

        m = re.match(r"^(seurat_wilcoxon|MAST)_(.+)\.csv$", f.name)
        if not m:
            print("  未匹配到 Seurat/MAST 结果文件模式，跳过:", f.name)
            continue

        method = m.group(1)
        comparison = m.group(2)

        df = pd.read_csv(f)
        df = ensure_gene_column(df)
        if df is None:
            print("  没有 gene/names 列，跳过:", f.name)
            continue

        sub = df.copy()
        sub["dataset"] = dataset
        sub["method"] = method
        sub["comparison"] = comparison
        rows.append(sub)

        print(f"  已收集: {dataset} | {method} | {comparison}")


# --------------------------------------------------
# 4. 保存汇总结果
# --------------------------------------------------
if len(rows) == 0:
    print("\n没有收集到任何结果。")
else:
    all_df = pd.concat(rows, ignore_index=True)

    out_file = OUT / "all_de_results_long.csv"
    all_df.to_csv(out_file, index=False)

    print("\n已保存:", out_file)

    print("\nmethod 计数：")
    print(all_df["method"].value_counts())

    print("\n去重后的 dataset / method / comparison：")
    print(
        all_df[["dataset", "method", "comparison"]]
        .drop_duplicates()
        .sort_values(["dataset", "method", "comparison"])
        .head(200)
    )