from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

METRICS_PATH = Path(r"D:\scRNA_benchmark\benchmark_results\benchmark_metrics.csv")
OUT_DIR = Path(r"D:\scRNA_benchmark\benchmark_results\plots")
OUT_DIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(METRICS_PATH)

# 方法均值
mean_df = (
    df.groupby("method")[["AUROC", "AUPRC", "TPR", "FPR", "benchmark_score"]]
    .mean()
    .reset_index()
    .sort_values("benchmark_score", ascending=False)
)

print(mean_df)

# 1. 综合分条形图
plt.figure(figsize=(8, 5))
plt.bar(mean_df["method"], mean_df["benchmark_score"])
plt.ylabel("Mean benchmark score")
plt.title("Overall benchmark score by method")
plt.xticks(rotation=30, ha="right")
plt.tight_layout()
plt.savefig(OUT_DIR / "mean_benchmark_score_barplot.png", dpi=300)
plt.close()

# 2. AUROC 箱线图
plt.figure(figsize=(8, 5))
df.boxplot(column="AUROC", by="method")
plt.title("AUROC distribution by method")
plt.suptitle("")
plt.ylabel("AUROC")
plt.xticks(rotation=30)
plt.tight_layout()
plt.savefig(OUT_DIR / "AUROC_boxplot.png", dpi=300)
plt.close()

# 3. AUPRC 箱线图
plt.figure(figsize=(8, 5))
df.boxplot(column="AUPRC", by="method")
plt.title("AUPRC distribution by method")
plt.suptitle("")
plt.ylabel("AUPRC")
plt.xticks(rotation=30)
plt.tight_layout()
plt.savefig(OUT_DIR / "AUPRC_boxplot.png", dpi=300)
plt.close()

# 4. TPR 箱线图
plt.figure(figsize=(8, 5))
df.boxplot(column="TPR", by="method")
plt.title("TPR distribution by method")
plt.suptitle("")
plt.ylabel("TPR")
plt.xticks(rotation=30)
plt.tight_layout()
plt.savefig(OUT_DIR / "TPR_boxplot.png", dpi=300)
plt.close()

# 5. FPR 箱线图
plt.figure(figsize=(8, 5))
df.boxplot(column="FPR", by="method")
plt.title("FPR distribution by method")
plt.suptitle("")
plt.ylabel("FPR")
plt.xticks(rotation=30)
plt.tight_layout()
plt.savefig(OUT_DIR / "FPR_boxplot.png", dpi=300)
plt.close()

# 6. comparison x method 热图
heat_df = df.pivot_table(
    index="comparison",
    columns="method",
    values="benchmark_score",
    aggfunc="mean"
)

plt.figure(figsize=(9, max(5, 0.35 * len(heat_df))))
plt.imshow(heat_df.values, aspect="auto")
plt.xticks(range(len(heat_df.columns)), heat_df.columns, rotation=45, ha="right")
plt.yticks(range(len(heat_df.index)), heat_df.index)
plt.colorbar(label="Benchmark score")
plt.title("Benchmark score heatmap")
plt.tight_layout()
plt.savefig(OUT_DIR / "benchmark_score_heatmap.png", dpi=300)
plt.close()

print("plots saved to:", OUT_DIR)