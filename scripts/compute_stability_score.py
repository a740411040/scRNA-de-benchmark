from pathlib import Path
import pandas as pd

METRICS = Path(r"D:\scRNA_benchmark\benchmark_results\benchmark_metrics.csv")
OUT = Path(r"D:\scRNA_benchmark\benchmark_results")

df = pd.read_csv(METRICS)

stab = (
    df.groupby("method")[["benchmark_score", "AUROC", "FPR"]]
    .std()
    .reset_index()
    .rename(columns={
        "benchmark_score": "sd_benchmark",
        "AUROC": "sd_AUROC",
        "FPR": "sd_FPR"
    })
)

# 反向归一化：波动越小，稳定性越高
for col in ["sd_benchmark", "sd_AUROC", "sd_FPR"]:
    maxv = stab[col].max()
    if maxv == 0:
        stab["stability_" + col] = 1.0
    else:
        stab["stability_" + col] = 1 - stab[col] / maxv

stab["stability_score"] = (
    0.5 * stab["stability_sd_benchmark"] +
    0.3 * stab["stability_sd_AUROC"] +
    0.2 * stab["stability_sd_FPR"]
)

mean_perf = (
    df.groupby("method")[["benchmark_score", "AUROC", "AUPRC", "TPR", "FPR"]]
    .mean()
    .reset_index()
)

final = pd.merge(mean_perf, stab[["method", "stability_score"]], on="method")
final["robust_overall_score"] = 0.7 * final["benchmark_score"] + 0.3 * final["stability_score"]
final = final.sort_values("robust_overall_score", ascending=False)

final.to_csv(OUT / "method_stability_summary.csv", index=False)
print(final)