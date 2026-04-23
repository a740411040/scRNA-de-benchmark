from pathlib import Path
import pandas as pd
import scipy.io
import scipy.sparse as sp
import anndata as ad

BRIDGE_DIR = Path(r"D:\scRNA_benchmark\bridge")
OUT_DIR = Path(r"D:\scRNA_benchmark\converted_h5ad")
OUT_DIR.mkdir(parents=True, exist_ok=True)

def build_h5ad_from_bridge(prefix: str):
    print(f"正在组装: {prefix}")

    counts_path = BRIDGE_DIR / f"{prefix}_counts.mtx"
    obs_path = BRIDGE_DIR / f"{prefix}_obs.csv"
    var_path = BRIDGE_DIR / f"{prefix}_var.csv"
    genes_path = BRIDGE_DIR / f"{prefix}_genes.csv"
    cells_path = BRIDGE_DIR / f"{prefix}_cells.csv"

    X = scipy.io.mmread(counts_path).tocsr()

    obs = pd.read_csv(obs_path)
    var = pd.read_csv(var_path)
    genes = pd.read_csv(genes_path)
    cells = pd.read_csv(cells_path)

    # Matrix Market from R writeMM(counts) => 行是gene，列是cell
    # AnnData需要 cell x gene，因此要转置
    X = X.T.tocsr()

    obs.index = cells["cell_id"].astype(str).values
    var.index = genes["gene_id"].astype(str).values

    obs = obs.loc[cells["cell_id"].astype(str).values].copy()
    var = var.loc[genes["gene_id"].astype(str).values].copy()

    adata = ad.AnnData(X=X, obs=obs, var=var)

    # 保存 raw counts
    adata.layers["counts"] = adata.X.copy()

    # 补统一字段
    adata.obs["dataset_id"] = prefix
    adata.obs["cell_id"] = adata.obs_names.astype(str)

    if "gene_id" not in adata.var.columns:
        adata.var["gene_id"] = adata.var_names.astype(str)

    # 尝试补 gene_symbol
    if "Symbol" in adata.var.columns:
        adata.var["gene_symbol"] = adata.var["Symbol"].astype(str)
    elif "symbol" in adata.var.columns:
        adata.var["gene_symbol"] = adata.var["symbol"].astype(str)
    else:
        adata.var["gene_symbol"] = adata.var_names.astype(str)

    out_path = OUT_DIR / f"{prefix}.h5ad"
    adata.write_h5ad(out_path)

    print(f"已保存: {out_path}")
    print(adata)
    print()

def main():
    # 自动扫描所有 bridge 前缀
    count_files = sorted(BRIDGE_DIR.glob("*_counts.mtx"))
    prefixes = [f.name.replace("_counts.mtx", "") for f in count_files]

    print("检测到以下 prefix：")
    for p in prefixes:
        print(" -", p)
    print()

    for prefix in prefixes:
        build_h5ad_from_bridge(prefix)

if __name__ == "__main__":
    main()