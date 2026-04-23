from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

NO_COV = Path(r"D:\scRNA_benchmark\benchmark_results\benchmark_metrics_no_covariate.csv")
WITH_COV = Path(r"D:\scRNA_benchmark\benchmark_results\benchmark_metrics_with_covariate.csv")
OUT = Path(r"D:\scRNA_benchmark\benchmark_results\plots")
OUT.mkdir(parents=True, exist_ok=True)

df0 = pd.read_csv(NO_COV)
df1 = pd.read_csv(WITH_COV)

df0["model"] = "no_covariate"
df1["model"] = "with_covariate"

df = pd.concat([df0, df1], ignore_index=True)

# 1. 方法平均表
mean_df = (
    df.groupby(["model", "method"])[["AUROC", "AUPRC", "TPR", "FPR", "benchmark_score"]]
    .mean()
    .reset_index()
)
mean_df.to_csv(OUT / "covariate_method_mean_summary.csv", index=False)
print(mean_df)

# 2. benchmark_score 对比图
pivot_score = mean_df.pivot(index="method", columns="model", values="benchmark_score")
pivot_score.plot(kind="bar", figsize=(8, 5))
plt.ylabel("Benchmark score")
plt.title("Covariate effect on benchmark score")
plt.xticks(rotation=30, ha="right")
plt.tight_layout()
plt.savefig(OUT / "covariate_benchmark_score_barplot.png", dpi=300)
plt.close()

# 3. AUROC 对比图
pivot_auroc = mean_df.pivot(index="method", columns="model", values="AUROC")
pivot_auroc.plot(kind="bar", figsize=(8, 5))
plt.ylabel("AUROC")
plt.title("Covariate effect on AUROC")
plt.xticks(rotation=30, ha="right")
plt.tight_layout()
plt.savefig(OUT / "covariate_AUROC_barplot.png", dpi=300)
plt.close()

# 4. FPR 对比图
pivot_fpr = mean_df.pivot(index="method", columns="model", values="FPR")
pivot_fpr.plot(kind="bar", figsize=(8, 5))
plt.ylabel("FPR")
plt.title("Covariate effect on FPR")
plt.xticks(rotation=30, ha="right")
plt.tight_layout()
plt.savefig(OUT / "covariate_FPR_barplot.png", dpi=300)
plt.close()

# 5. delta 表
merge_df = pd.merge(
    df0,
    df1,
    on=["dataset", "method", "comparison", "truth_table"],
    suffixes=("_no_cov", "_with_cov")
)

merge_df["delta_benchmark_score"] = merge_df["benchmark_score_with_cov"] - merge_df["benchmark_score_no_cov"]
merge_df["delta_AUROC"] = merge_df["AUROC_with_cov"] - merge_df["AUROC_no_cov"]
merge_df["delta_AUPRC"] = merge_df["AUPRC_with_cov"] - merge_df["AUPRC_no_cov"]
merge_df["delta_TPR"] = merge_df["TPR_with_cov"] - merge_df["TPR_no_cov"]
merge_df["delta_FPR"] = merge_df["FPR_with_cov"] - merge_df["FPR_no_cov"]

merge_df.to_csv(OUT / "covariate_delta_table.csv", index=False)

delta_mean = (
    merge_df.groupby("method")[["delta_benchmark_score", "delta_AUROC", "delta_AUPRC", "delta_TPR", "delta_FPR"]]
    .mean()
    .reset_index()
)
delta_mean.to_csv(OUT / "covariate_delta_mean_summary.csv", index=False)
print(delta_mean)

# 6. Δbenchmark_score 图
plt.figure(figsize=(8, 5))
plt.bar(delta_mean["method"], delta_mean["delta_benchmark_score"])
plt.ylabel("Δ benchmark score (with - no)")
plt.title("Mean covariate gain by method")
plt.xticks(rotation=30, ha="right")
plt.tight_layout()
plt.savefig(OUT / "covariate_delta_benchmark_score_barplot.png", dpi=300)
plt.close()

# 7. ΔFPR 图
plt.figure(figsize=(8, 5))
plt.bar(delta_mean["method"], delta_mean["delta_FPR"])
plt.ylabel("Δ FPR (with - no)")
plt.title("Mean covariate effect on FPR")
plt.xticks(rotation=30, ha="right")
plt.tight_layout()
plt.savefig(OUT / "covariate_delta_FPR_barplot.png", dpi=300)
plt.close()

print("saved to:", OUT)