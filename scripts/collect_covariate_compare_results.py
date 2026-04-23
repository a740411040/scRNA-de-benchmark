from pathlib import Path
import pandas as pd

ROOT = Path(r"D:\scRNA_benchmark\benchmark_runs\SCMixology_5cl_CELseq2_merged\covariate_compare")
OUT = Path(r"D:\scRNA_benchmark\benchmark_results")
OUT.mkdir(parents=True, exist_ok=True)

rows_no = []
rows_with = []

def add_result(df, method, comparison, dataset, bucket):
    if "gene" not in df.columns:
        return
    sub = df.copy()
    sub["dataset"] = dataset
    sub["method"] = method
    sub["comparison"] = comparison
    bucket.append(sub)

# 这里 comparison 先固定为同一个 pair
# 因为你当前 covariate compare 脚本是只对一个 pair 做的
comparison_name = "A549_vs_H1975"
dataset_name = "SCMixology_5cl_CELseq2_merged"

file_map = {
    "edgeR_no_cov.csv": ("edgeR", "no"),
    "edgeR_with_cov.csv": ("edgeR", "with"),
    "DESeq2_no_cov.csv": ("DESeq2", "no"),
    "DESeq2_with_cov.csv": ("DESeq2", "with"),
    "limma_voom_no_cov.csv": ("limma_voom", "no"),
    "limma_voom_with_cov.csv": ("limma_voom", "with"),
    "MAST_no_cov.csv": ("MAST", "no"),
    "MAST_with_cov.csv": ("MAST", "with"),
}

for fname, (method, cov_type) in file_map.items():
    f = ROOT / fname
    if not f.exists():
        print("missing:", f)
        continue

    df = pd.read_csv(f)
    if cov_type == "no":
        add_result(df, method, comparison_name, dataset_name, rows_no)
    else:
        add_result(df, method, comparison_name, dataset_name, rows_with)

if len(rows_no) > 0:
    no_df = pd.concat(rows_no, ignore_index=True)
    no_path = OUT / "all_de_results_no_covariate.csv"
    no_df.to_csv(no_path, index=False)
    print("saved:", no_path)
    print(no_df["method"].value_counts())

if len(rows_with) > 0:
    with_df = pd.concat(rows_with, ignore_index=True)
    with_path = OUT / "all_de_results_with_covariate.csv"
    with_df.to_csv(with_path, index=False)
    print("saved:", with_path)
    print(with_df["method"].value_counts())