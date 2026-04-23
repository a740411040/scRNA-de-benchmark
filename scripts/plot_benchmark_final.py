from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

METRICS = Path(r"D:\scRNA_benchmark\benchmark_results\benchmark_metrics.csv")
OUT = Path(r"D:\scRNA_benchmark\benchmark_results\plots")
OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(METRICS)

# 1. 方法总体均值
mean_df = (
    df.groupby("method")[["AUROC", "AUPRC", "TPR", "FPR", "benchmark_score"]]
    .mean()
    .reset_index()
    .sort_values("benchmark_score", ascending=False)
)

# 图1：总体条形图
plt.figure(figsize=(8, 5))
plt.bar(mean_df["method"], mean_df["benchmark_score"])
plt.ylabel("Mean benchmark score")
plt.title("Overall benchmark score by method")
plt.xticks(rotation=30, ha="right")
plt.tight_layout()
plt.savefig(OUT / "fig1_overall_benchmark_score.png", dpi=300)
plt.close()

# 图2：AUROC 箱线图
plt.figure(figsize=(8, 5))
df.boxplot(column="AUROC", by="method")
plt.title("AUROC by method")
plt.suptitle("")
plt.ylabel("AUROC")
plt.xticks(rotation=30)
plt.tight_layout()
plt.savefig(OUT / "fig2_AUROC_boxplot.png", dpi=300)
plt.close()

# 图3：AUPRC 箱线图
plt.figure(figsize=(8, 5))
df.boxplot(column="AUPRC", by="method")
plt.title("AUPRC by method")
plt.suptitle("")
plt.ylabel("AUPRC")
plt.xticks(rotation=30)
plt.tight_layout()
plt.savefig(OUT / "fig3_AUPRC_boxplot.png", dpi=300)
plt.close()

# 图4：FPR 箱线图
plt.figure(figsize=(8, 5))
df.boxplot(column="FPR", by="method")
plt.title("FPR by method")
plt.suptitle("")
plt.ylabel("FPR")
plt.xticks(rotation=30)
plt.tight_layout()
plt.savefig(OUT / "fig4_FPR_boxplot.png", dpi=300)
plt.close()

# 图5：dataset × method 热图
heat_df = df.pivot_table(
    index="dataset",
    columns="method",
    values="benchmark_score",
    aggfunc="mean"
)

plt.figure(figsize=(9, max(4, 0.6 * len(heat_df))))
plt.imshow(heat_df.values, aspect="auto")
plt.xticks(range(len(heat_df.columns)), heat_df.columns, rotation=45, ha="right")
plt.yticks(range(len(heat_df.index)), heat_df.index)
plt.colorbar(label="Mean benchmark score")
plt.title("Dataset × Method benchmark score")
plt.tight_layout()
plt.savefig(OUT / "fig5_dataset_method_heatmap.png", dpi=300)
plt.close()

print("plots saved to:", OUT)