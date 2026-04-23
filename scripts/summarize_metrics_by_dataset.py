from pathlib import Path
import pandas as pd

METRICS = Path(r"D:\scRNA_benchmark\benchmark_results\benchmark_metrics.csv")
OUT = Path(r"D:\scRNA_benchmark\benchmark_results")
OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(METRICS)

# 1. 方法总体均值
method_summary = (
    df.groupby("method")[["AUROC", "AUPRC", "TPR", "FPR", "benchmark_score"]]
    .mean()
    .reset_index()
    .sort_values("benchmark_score", ascending=False)
)
method_summary.to_csv(OUT / "summary_by_method.csv", index=False)

# 2. 数据集 × 方法均值
dataset_method_summary = (
    df.groupby(["dataset", "method"])[["AUROC", "AUPRC", "TPR", "FPR", "benchmark_score"]]
    .mean()
    .reset_index()
    .sort_values(["dataset", "benchmark_score"], ascending=[True, False])
)
dataset_method_summary.to_csv(OUT / "summary_by_dataset_method.csv", index=False)

# 3. comparison × 方法均值
comparison_method_summary = (
    df.groupby(["comparison", "method"])[["AUROC", "AUPRC", "TPR", "FPR", "benchmark_score"]]
    .mean()
    .reset_index()
    .sort_values(["comparison", "benchmark_score"], ascending=[True, False])
)
comparison_method_summary.to_csv(OUT / "summary_by_comparison_method.csv", index=False)

print("saved:")
print(OUT / "summary_by_method.csv")
print(OUT / "summary_by_dataset_method.csv")
print(OUT / "summary_by_comparison_method.csv")
print("\nTop methods overall:")
print(method_summary)