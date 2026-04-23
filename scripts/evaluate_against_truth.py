from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score

RESULTS_PATH = Path(r"D:\scRNA_benchmark\benchmark_results\all_de_results_long.csv")
TRUTH_DIR = Path(r"D:\scRNA_benchmark\truth_tables")
OUT_DIR = Path(r"D:\scRNA_benchmark\benchmark_results")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# -------------------------
# 1. 读取方法结果
# -------------------------
res = pd.read_csv(RESULTS_PATH, low_memory=False)
print("方法结果 shape:", res.shape)
print(res[["dataset", "method", "comparison"]].drop_duplicates().head(20))

# -------------------------
# 2. 读取 truth tables
# -------------------------
truth_tables = {}
for f in TRUTH_DIR.glob("*.csv"):
    df = pd.read_csv(f)
    truth_tables[f.stem] = df

print("\ntruth tables:", list(truth_tables.keys()))

# -------------------------
# 3. 定义 truth positive
# -------------------------
def build_truth_set(truth_df, fdr_cutoff=0.05, logfc_cutoff=0.05):
    df = truth_df.copy()
    df["FDR"] = pd.to_numeric(df["FDR"], errors="coerce")
    df["logFC"] = pd.to_numeric(df["logFC"], errors="coerce")
    df["gene"] = df["gene"].astype(str)

    pos = df[
        (df["FDR"] < fdr_cutoff) &
        (df["logFC"].abs() > logfc_cutoff)
    ].copy()

    return set(pos["gene"])

truth_sets = {
    name: build_truth_set(df, fdr_cutoff=0.05, logfc_cutoff=0.05)
    for name, df in truth_tables.items()
}

for k, v in truth_sets.items():
    print(f"{k}: {len(v)} truth positive genes")

# -------------------------
# 4. comparison -> truth table 映射
# 只评估包含 H1975 / H2228 / HCC827 的比较
# scanpy 的 single_pair 先自动跳过
# -------------------------
def choose_truth_table(comparison: str):
    comparison = str(comparison)
    if "H1975" in comparison:
        return "H1975_DEtable"
    if "H2228" in comparison:
        return "H2228_DEtable"
    if "HCC827" in comparison:
        return "HCC827_DEtable"
    return None

# -------------------------
# 5. 统一生成 score
# score 越大，越支持 DE
# -------------------------
def first_valid_numeric_series(df, candidates):
    """
    按顺序寻找“存在且不是全 NaN”的数值列
    """
    for c in candidates:
        if c in df.columns:
            s = pd.to_numeric(df[c], errors="coerce")
            if s.notna().sum() > 0:
                return c, s
    return None, pd.Series(np.nan, index=df.index)


def first_valid_numeric_series(df, candidates):
    for c in candidates:
        if c in df.columns:
            s = pd.to_numeric(df[c], errors="coerce")
            if s.notna().sum() > 0:
                return c, s
    return None, pd.Series(np.nan, index=df.index)


def make_score(df, dataset=None, method=None, comparison=None):
    # 0. 针对 Seurat / MAST：优先用显著性 * 效应大小
    if method in ["seurat_wilcoxon", "MAST"]:
        p_col, p_s = first_valid_numeric_series(df, ["p_val_adj", "p_val", "padj", "pvalue"])
        fc_col, fc_s = first_valid_numeric_series(df, ["avg_log2FC", "avg_logFC", "logFC", "log2FoldChange"])

        if p_col is not None and fc_col is not None:
            score = -np.log10(p_s + 1e-300) * fc_s.abs()
            print(f"[score] {dataset} | {method} | {comparison} -> using -log10({p_col}) * abs({fc_col})")
            return score

        if p_col is not None:
            print(f"[score] {dataset} | {method} | {comparison} -> using -log10({p_col})")
            return -np.log10(p_s + 1e-300)

        if fc_col is not None:
            print(f"[score] {dataset} | {method} | {comparison} -> using abs({fc_col})")
            return fc_s.abs()

    # 1. 统计量类：优先
    col, s = first_valid_numeric_series(df, ["scores", "stat", "t", "LR", "F"])
    if col is not None:
        print(f"[score] {dataset} | {method} | {comparison} -> using {col}")
        return s

    # 2. fold change 类
    col, s = first_valid_numeric_series(
        df,
        ["avg_log2FC", "avg_logFC", "log2FoldChange", "logFC", "logfoldchanges"]
    )
    if col is not None:
        print(f"[score] {dataset} | {method} | {comparison} -> using abs({col})")
        return s.abs()

    # 3. p / FDR 类
    col, s = first_valid_numeric_series(
        df,
        ["p_val_adj", "pvals_adj", "padj", "FDR", "adj.P.Val",
         "p_val", "PValue", "P.Value", "pvalue", "pvals"]
    )
    if col is not None:
        print(f"[score] {dataset} | {method} | {comparison} -> using -log10({col})")
        return -np.log10(s + 1e-300)

    print(f"[score] {dataset} | {method} | {comparison} -> no valid score column")
    return pd.Series(np.nan, index=df.index)

# -------------------------
# 6. 统一 gene 列
# -------------------------
def ensure_gene_column(df):
    out = df.copy()

    if "gene" in out.columns:
        out["gene"] = out["gene"].astype(str)
        return out

    if "names" in out.columns:
        out["gene"] = out["names"].astype(str)
        return out

    return None

# -------------------------
# 7. 计算 metrics
# -------------------------
metrics_rows = []

for (dataset, method, comparison), df in res.groupby(["dataset", "method", "comparison"]):
    truth_name = choose_truth_table(comparison)
    if truth_name is None:
        continue
    if truth_name not in truth_sets:
        continue

    sub = ensure_gene_column(df)
    if sub is None:
        continue

    sub["gene"] = sub["gene"].astype(str)
    sub["score"] = make_score(sub, dataset=dataset, method=method, comparison=comparison)
    sub = sub.dropna(subset=["score"]).copy()

    truth_set = truth_sets[truth_name]
    sub["label"] = sub["gene"].isin(truth_set).astype(int)
    

    if sub["label"].nunique() < 2:
        continue

    threshold = sub["score"].quantile(0.95)
    sub["pred"] = (sub["score"] >= threshold).astype(int)

    tp = ((sub["pred"] == 1) & (sub["label"] == 1)).sum()
    fp = ((sub["pred"] == 1) & (sub["label"] == 0)).sum()
    fn = ((sub["pred"] == 0) & (sub["label"] == 1)).sum()
    tn = ((sub["pred"] == 0) & (sub["label"] == 0)).sum()

    tpr = tp / (tp + fn) if (tp + fn) > 0 else np.nan
    fpr = fp / (fp + tn) if (fp + tn) > 0 else np.nan

    auroc = roc_auc_score(sub["label"], sub["score"])
    auprc = average_precision_score(sub["label"], sub["score"])

    metrics_rows.append({
        "dataset": dataset,
        "method": method,
        "comparison": comparison,
        "truth_table": truth_name,
        "n_genes": len(sub),
        "n_truth_genes": int(sub["label"].sum()),
        "AUROC": auroc,
        "AUPRC": auprc,
        "TPR": tpr,
        "FPR": fpr
    })

metrics_df = pd.DataFrame(metrics_rows)

if len(metrics_df) == 0:
    print("\n没有生成任何可评估结果。")
    metrics_df.to_csv(OUT_DIR / "benchmark_metrics.csv", index=False)
    print("已保存空文件:", OUT_DIR / "benchmark_metrics.csv")
else:
    metrics_df["benchmark_score"] = (
        0.35 * metrics_df["AUROC"] +
        0.35 * metrics_df["AUPRC"] +
        0.20 * metrics_df["TPR"] +
        0.10 * (1 - metrics_df["FPR"])
    )

    metrics_df.to_csv(OUT_DIR / "benchmark_metrics.csv", index=False)

    print("\n已保存:", OUT_DIR / "benchmark_metrics.csv")
    print(
        metrics_df.sort_values(
            ["dataset", "comparison", "benchmark_score"],
            ascending=[True, True, False]
        ).head(50)
    )