from __future__ import annotations

"""
scRNA-seq 差异表达（DE）benchmark：Python / Scanpy 侧流程
-------------------------------------------------------
这个脚本专门负责：
1) 读取与预处理 AnnData
2) 构建原生单细胞检验所需输入（log-normalized）
3) 构建跨语言桥接文件（给 R / Seurat / Bioconductor）
4) 人为噪声注入（dropout / gaussian / mislabel）
5) 运行 Scanpy 默认 t-test / Wilcoxon
6) 运行 MEMENTO（给出尽量稳定的包装；API 版本变化时可小改）
7) 评估 AUROC / AUPRC / TPR / FPR
8) 绘图

设计原则
--------
- 所有“原始计数”统一存储在 adata.layers['counts']
- Scanpy 原生 DE 使用 log1p 归一化矩阵
- 伪 bulk 与 R 侧方法一律使用原始 count
- 噪声注入尽量在 counts / lognorm 两个层面都保留，便于不同方法公平比较

注意
----
- 对大数据集建议 QC 后再随机下采样到 ~5000 cells
- 对高斯噪声，建议在 HVG 或有限基因集合上使用，以控制内存
- MEMENTO 的公开 API 在不同小版本间可能有轻微变动；本脚本给出的是 v0.1.x 思路
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Sequence
import json
import subprocess
import tempfile
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.io import mmwrite
from sklearn.metrics import average_precision_score, roc_auc_score
import matplotlib.pyplot as plt


# -----------------------------------------------------------------------------
# 基础工具函数
# -----------------------------------------------------------------------------

def _to_csr(x):
    """确保矩阵为 CSR，方便按行切片与稀疏写出。"""
    if sp.issparse(x):
        return x.tocsr()
    return sp.csr_matrix(x)


def _dense_from_any(x) -> np.ndarray:
    """小规模数据转 dense。建议在 5000 cells / 2000~3000 genes 内使用。"""
    if sp.issparse(x):
        return x.toarray()
    return np.asarray(x)


def _get_layer_matrix(adata: sc.AnnData, layer: str):
    if layer == "X":
        return adata.X
    if layer not in adata.layers:
        raise KeyError(f"层 {layer!r} 不存在。当前层：{list(adata.layers.keys())}")
    return adata.layers[layer]


def _ensure_counts_layer(adata: sc.AnnData, counts_layer: str = "counts") -> None:
    """如果没有 counts 层，则默认把当前 X 视为原始 count。"""
    if counts_layer not in adata.layers:
        adata.layers[counts_layer] = _to_csr(adata.X.copy())


def load_benchmark_adata(adata_path: str | Path) -> sc.AnnData:
    """
    统一读取 benchmark 输入。

    支持：
    - .h5ad：直接读取
    - .RData：调用 R 将对象临时转换成 h5ad 后再读取
    """
    adata_path = Path(adata_path)
    suffix = adata_path.suffix.lower()

    if suffix == ".h5ad":
        return sc.read_h5ad(adata_path)

    if suffix == ".rdata":
        with tempfile.TemporaryDirectory(prefix="scmixology_rdata_") as tmpdir:
            tmp_h5ad = Path(tmpdir) / "converted_from_rdata.h5ad"
            r_code = """
args <- commandArgs(trailingOnly = TRUE)
input_path <- args[1]
output_path <- args[2]
obj_env <- new.env(parent = emptyenv())
loaded_names <- load(input_path, envir = obj_env)
loaded_objects <- mget(loaded_names, envir = obj_env)

pick_first <- function(class_name) {
  idx <- vapply(loaded_objects, function(x) inherits(x, class_name), logical(1))
  if (any(idx)) loaded_objects[[which(idx)[1]]] else NULL
}

obj <- pick_first("SingleCellExperiment")
if (is.null(obj)) {
  se_obj <- pick_first("SummarizedExperiment")
  if (!is.null(se_obj)) {
    obj <- SingleCellExperiment::as(se_obj, "SingleCellExperiment")
  }
}
if (is.null(obj)) {
  seu_obj <- pick_first("Seurat")
  if (!is.null(seu_obj)) {
    obj <- Seurat::as.SingleCellExperiment(seu_obj)
  }
}
if (is.null(obj)) {
  stop(
    paste0(
      "RData 文件中没有可转换的 SingleCellExperiment / SummarizedExperiment / Seurat 对象。已加载对象：",
      paste(loaded_names, collapse = ", ")
    )
  )
}
zellkonverter::writeH5AD(obj, output_path)
"""
            try:
                subprocess.run(
                    ["Rscript", "-e", r_code, str(adata_path), str(tmp_h5ad)],
                    check=True,
                    capture_output=True,
                    text=True,
                )
            except FileNotFoundError as exc:
                raise RuntimeError(
                    "检测到 .RData 输入，但当前环境没有可用的 Rscript，无法自动转换为 h5ad。"
                ) from exc
            except subprocess.CalledProcessError as exc:
                stderr = exc.stderr.strip() if exc.stderr else str(exc)
                raise RuntimeError(
                    "读取 .RData 失败。请确认 R 已安装，并且可用包包括 zellkonverter、"
                    "SingleCellExperiment，以及在需要时的 Seurat。原始错误："
                    f"{stderr}"
                ) from exc

            return sc.read_h5ad(tmp_h5ad)

    raise ValueError(f"目前只支持 .h5ad 或 .RData 输入：{adata_path}")


# -----------------------------------------------------------------------------
# 预处理
# -----------------------------------------------------------------------------

def qc_and_downsample(
    adata: sc.AnnData,
    *,
    counts_layer: str = "counts",
    min_genes_per_cell: int = 200,
    min_cells_per_gene: int = 3,
    mt_prefix: str = "MT-",
    max_mt_pct: float = 20.0,
    target_n_cells: Optional[int] = 5000,
    seed: int = 1234,
) -> sc.AnnData:
    """
    标准 QC：
    - 细胞至少检测到 min_genes_per_cell 个基因
    - 基因至少出现在 min_cells_per_gene 个细胞中
    - 过滤线粒体比例过高的细胞
    - 如细胞过多，随机下采样到 target_n_cells
    """
    adata = adata.copy()
    _ensure_counts_layer(adata, counts_layer)

    # 用 counts 做 QC 更稳妥
    adata.X = _get_layer_matrix(adata, counts_layer).copy()
    adata.var["mt"] = adata.var_names.str.upper().str.startswith(mt_prefix.upper())
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    adata = adata[adata.obs["n_genes_by_counts"] >= min_genes_per_cell].copy()
    adata = adata[adata.obs["pct_counts_mt"] <= max_mt_pct].copy()
    sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)

    if target_n_cells is not None and adata.n_obs > target_n_cells:
        rng = np.random.default_rng(seed)
        keep = rng.choice(adata.n_obs, size=target_n_cells, replace=False)
        adata = adata[keep].copy()

    # 恢复 X 为 counts，避免后续误用
    adata.X = _get_layer_matrix(adata, counts_layer).copy()
    return adata


def make_native_and_bridge_layers(
    adata: sc.AnnData,
    *,
    counts_layer: str = "counts",
    target_sum: float = 1e4,
    n_top_hvgs: Optional[int] = None,
) -> sc.AnnData:
    """
    为 benchmark 统一准备两个输入视图：
    1) counts: 给 pseudo-bulk / edgeR / DESeq2 / limma-voom / nebula
    2) lognorm: 给 Scanpy / Seurat 默认单细胞检验

    可选：只保留 HVGs，适合大样本噪声 stress test。
    """
    adata = adata.copy()
    _ensure_counts_layer(adata, counts_layer)

    # 先基于原始 count 构造 lognorm
    adata.X = _get_layer_matrix(adata, counts_layer).copy()
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    adata.layers["lognorm"] = _to_csr(adata.X.copy())

    if n_top_hvgs is not None:
        # 用原始 count 做 HVG 选择更符合单细胞常规流程
        tmp = adata.copy()
        tmp.X = _get_layer_matrix(tmp, counts_layer).copy()
        sc.pp.highly_variable_genes(tmp, flavor="seurat_v3", n_top_genes=n_top_hvgs, layer=None)
        keep = tmp.var["highly_variable"].to_numpy()
        adata = adata[:, keep].copy()

    # 默认把 X 指向 lognorm，便于 Scanpy 直接运行
    adata.X = _get_layer_matrix(adata, "lognorm").copy()
    return adata


def export_bridge_files(
    adata: sc.AnnData,
    outdir: str | Path,
    *,
    counts_layer: str = "counts",
    write_h5ad: bool = True,
    prefix: str = "benchmark_input",
) -> None:
    """
    导出给 R 侧使用的桥接文件：
    - h5ad：R 可用 zellkonverter::readH5AD 读取
    - 10x 风格 mtx/tsv：R / Python 都通用
    - obs / var 元数据 csv
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    _ensure_counts_layer(adata, counts_layer)

    if write_h5ad:
        adata.write_h5ad(outdir / f"{prefix}.h5ad")

    counts = _to_csr(_get_layer_matrix(adata, counts_layer))
    mmwrite(outdir / f"{prefix}.matrix.mtx", counts.T)  # 10x 风格：gene x cell

    pd.DataFrame({"gene": adata.var_names}).to_csv(
        outdir / f"{prefix}.features.tsv", sep="\t", index=False, header=False
    )
    pd.DataFrame({"cell": adata.obs_names}).to_csv(
        outdir / f"{prefix}.barcodes.tsv", sep="\t", index=False, header=False
    )
    adata.obs.to_csv(outdir / f"{prefix}.obs.csv")
    adata.var.to_csv(outdir / f"{prefix}.var.csv")


# -----------------------------------------------------------------------------
# 噪声注入
# -----------------------------------------------------------------------------

def inject_dropout_noise(
    adata: sc.AnnData,
    *,
    counts_layer: str = "counts",
    output_layer: str = "counts_dropout",
    dropout_rate: float = 0.2,
    low_expr_quantile: float = 0.3,
    seed: int = 1234,
) -> sc.AnnData:
    """
    人为增加 dropout：
    只对“低表达的非零值”施加置零操作，更接近真实单细胞稀疏性。

    参数
    ----
    dropout_rate:
        在候选低表达非零条目中，被强制置 0 的概率。
    low_expr_quantile:
        仅考虑非零表达值中最底部 low_expr_quantile 分位数那部分。
    """
    adata = adata.copy()
    _ensure_counts_layer(adata, counts_layer)
    mat = _to_csr(_get_layer_matrix(adata, counts_layer)).tocoo(copy=True)

    if mat.nnz == 0:
        adata.layers[output_layer] = mat.tocsr()
        return adata

    rng = np.random.default_rng(seed)
    threshold = np.quantile(mat.data, low_expr_quantile)
    candidate_idx = np.where(mat.data <= threshold)[0]

    if candidate_idx.size == 0:
        adata.layers[output_layer] = mat.tocsr()
        return adata

    to_drop = candidate_idx[rng.random(candidate_idx.size) < dropout_rate]
    keep_mask = np.ones(mat.data.shape[0], dtype=bool)
    keep_mask[to_drop] = False

    new_mat = sp.csr_matrix(
        (mat.data[keep_mask], (mat.row[keep_mask], mat.col[keep_mask])),
        shape=mat.shape,
    )
    adata.layers[output_layer] = new_mat
    return adata


def inject_gaussian_noise(
    adata: sc.AnnData,
    *,
    input_layer: str = "lognorm",
    output_layer: str = "lognorm_gaussian",
    output_count_layer: str = "counts_gaussian",
    cell_fraction: float = 0.1,
    sigma: float = 0.25,
    seed: int = 1234,
) -> sc.AnnData:
    """
    随机选取一部分细胞，在 log-normalized 表达矩阵上叠加高斯噪声。

    为什么在 lognorm 上加噪？
    --------------------
    因为高斯背景噪声更自然地对应“连续表达空间”的扰动；
    随后再通过 expm1 + round 回写出近似 count 层，供 count-based 方法使用。

    注意
    ----
    - 建议在 QC + HVG 之后再用此函数，避免超大矩阵 densify 占内存。
    """
    adata = adata.copy()
    x = _dense_from_any(_get_layer_matrix(adata, input_layer))
    rng = np.random.default_rng(seed)

    n_cells = adata.n_obs
    n_pick = max(1, int(np.floor(n_cells * cell_fraction)))
    picked = rng.choice(n_cells, size=n_pick, replace=False)

    x_noisy = x.copy()
    x_noisy[picked, :] += rng.normal(loc=0.0, scale=sigma, size=(n_pick, x.shape[1]))
    x_noisy = np.clip(x_noisy, a_min=0.0, a_max=None)

    adata.layers[output_layer] = _to_csr(x_noisy)
    approx_counts = np.rint(np.expm1(x_noisy)).astype(np.int64)
    approx_counts[approx_counts < 0] = 0
    adata.layers[output_count_layer] = _to_csr(approx_counts)
    adata.obs[f"{output_layer}_flag"] = False
    adata.obs.iloc[picked, adata.obs.columns.get_loc(f"{output_layer}_flag")] = True
    return adata


def inject_mislabeling(
    adata: sc.AnnData,
    *,
    label_col: str,
    mislabel_rate: float = 0.1,
    output_col: str = "celltype_mislabeled",
    seed: int = 1234,
) -> sc.AnnData:
    """
    将 X% 的细胞标签随机错配为其它类别。
    这是评估算法对 annotation 错误敏感性的关键操作。
    """
    adata = adata.copy()
    if label_col not in adata.obs:
        raise KeyError(f"obs 中不存在标签列 {label_col!r}")

    labels = adata.obs[label_col].astype(str).copy()
    uniq = sorted(labels.unique().tolist())
    if len(uniq) < 2:
        raise ValueError("标签类别不足 2，无法进行错配。")

    rng = np.random.default_rng(seed)
    n_flip = max(1, int(np.floor(adata.n_obs * mislabel_rate)))
    picked = rng.choice(adata.n_obs, size=n_flip, replace=False)

    new_labels = labels.to_numpy().copy()
    for idx in picked:
        current = new_labels[idx]
        candidates = [x for x in uniq if x != current]
        new_labels[idx] = rng.choice(candidates)

    adata.obs[output_col] = new_labels
    adata.obs[f"{output_col}_flag"] = False
    adata.obs.iloc[picked, adata.obs.columns.get_loc(f"{output_col}_flag")] = True
    return adata


# -----------------------------------------------------------------------------
# DE 方法封装：Scanpy / MEMENTO
# -----------------------------------------------------------------------------

def run_scanpy_rank_genes_groups(
    adata: sc.AnnData,
    *,
    group_col: str,
    case_label: str,
    control_label: str,
    method: str = "wilcoxon",
    use_layer: str = "lognorm",
    corr_method: str = "benjamini-hochberg",
) -> pd.DataFrame:
    """
    Scanpy 原生 DE：支持 t-test / wilcoxon。
    这里显式做二组比较，避免 'rest' 参考组引入额外复杂性。
    """
    if method not in {"wilcoxon", "t-test", "t-test_overestim_var"}:
        raise ValueError("method 必须是 'wilcoxon' / 't-test' / 't-test_overestim_var'")

    tmp = adata[adata.obs[group_col].isin([case_label, control_label])].copy()
    tmp.obs[group_col] = pd.Categorical(tmp.obs[group_col], categories=[case_label, control_label])
    tmp.X = _get_layer_matrix(tmp, use_layer).copy()

    sc.tl.rank_genes_groups(
        tmp,
        groupby=group_col,
        groups=[case_label],
        reference=control_label,
        method=method,
        corr_method=corr_method,
        tie_correct=True,
        pts=True,
    )

    df = sc.get.rank_genes_groups_df(tmp, group=case_label)
    df = df.rename(
        columns={
            "names": "gene",
            "scores": "score",
            "logfoldchanges": "logFC",
            "pvals": "pval",
            "pvals_adj": "padj",
        }
    )
    df["method"] = f"scanpy_{method}"
    df["case"] = case_label
    df["control"] = control_label
    return df


def run_memento_example(
    adata: sc.AnnData,
    *,
    group_col: str,
    case_label: str,
    control_label: str,
    replicate_col: Optional[str] = None,
    counts_layer: str = "counts",
    capture_rate: float = 0.07,
) -> pd.DataFrame:
    """
    MEMENTO 示例包装。

    说明
    ----
    - MEMENTO 官方实现与 Scanpy 兼容，但不同小版本 API 细节可能变化。
    - 此处给出 v0.1.x 常见思路：setup -> 分组 -> 计算 moments -> 单基因检验。
    - 若你本地安装版本函数名略有不同，请按照 memento 文档将下面 2~3 行函数名做微调。
    """
    try:
        import memento  # type: ignore
    except Exception as exc:
        raise ImportError(
            "未能导入 memento。请先安装：pip install git+https://github.com/yelabucsf/scrna-parameter-estimation.git"
        ) from exc

    tmp = adata[adata.obs[group_col].isin([case_label, control_label])].copy()
    _ensure_counts_layer(tmp, counts_layer)
    tmp.X = _get_layer_matrix(tmp, counts_layer).copy()
    tmp.obs[group_col] = pd.Categorical(tmp.obs[group_col], categories=[control_label, case_label])
    tmp.obs["capture_rate"] = capture_rate

    # ---- 下面几步是 memento 的典型分析框架 ----
    # 注：不同版本可能是 setup_memento / create_groups / compute_1d_moments / ht_1d_moments
    # 如果你的版本函数签名不同，通常只需改动参数名。
    memento.setup_memento(tmp, q_column="capture_rate")

    if replicate_col is None:
        # 没有 donor / sample 时，退化为按条件分组；
        # 更推荐在多样本设计中使用 replicate_col。
        tmp.obs["memento_group"] = tmp.obs[group_col].astype(str)
        memento.create_groups(tmp, label_columns=["memento_group"])
    else:
        memento.create_groups(tmp, label_columns=[group_col, replicate_col])

    memento.compute_1d_moments(tmp, min_perc_group=0.7)

    # 常见做法：对条件项做检验
    # 若 API 有变动，请把 binary_test_1d / get_1d_ht_result 对应替换成你版本的函数。
    memento.binary_test_1d(
        tmp,
        treatment_col=group_col,
        num_boot=500,
        verbose=0,
        treatment=case_label,
        control=control_label,
    )
    res = memento.get_1d_ht_result(tmp)

    # 统一输出字段名
    col_map = {
        "gene": "gene",
        "de_coef": "logFC",
        "de_pval": "pval",
        "de_fdr": "padj",
    }
    rename_cols = {k: v for k, v in col_map.items() if k in res.columns}
    res = res.rename(columns=rename_cols).copy()
    if "gene" not in res.columns:
        # 某些版本 gene 作为 index
        res = res.reset_index().rename(columns={res.columns[0]: "gene"})
    if "logFC" not in res.columns:
        # 若只有效应量列名不同，保留一个近似统一名
        effect_cols = [c for c in res.columns if "coef" in c.lower() or "effect" in c.lower()]
        if effect_cols:
            res = res.rename(columns={effect_cols[0]: "logFC"})
    if "padj" not in res.columns and "pval" in res.columns:
        res["padj"] = pd.Series(res["pval"]).rank(method="average") / len(res)  # 占位；建议本地改为 BH

    res["method"] = "memento"
    res["case"] = case_label
    res["control"] = control_label
    return res


# -----------------------------------------------------------------------------
# 评估指标
# -----------------------------------------------------------------------------

def evaluate_de_result(
    de_df: pd.DataFrame,
    *,
    gold_genes: Sequence[str] | set[str],
    gene_col: str = "gene",
    p_col: str = "padj",
    score_col: Optional[str] = None,
    alpha: float = 0.05,
) -> dict:
    """
    在已知 gold standard DE gene 列表时，计算：
    - AUROC
    - AUPRC
    - TPR (Sensitivity)
    - FPR

    默认打分逻辑
    ------------
    score = -log10(adjusted_pvalue + 1e-300)
    即：显著性越强，分数越高。
    """
    gold_genes = set(gold_genes)
    x = de_df.copy()
    x = x.drop_duplicates(subset=[gene_col]).copy()
    x[gene_col] = x[gene_col].astype(str)
    x["is_gold"] = x[gene_col].isin(gold_genes).astype(int)

    if score_col is None:
        x["_score_"] = -np.log10(pd.to_numeric(x[p_col], errors="coerce").fillna(1.0).clip(lower=1e-300))
        score_col = "_score_"

    y_true = x["is_gold"].to_numpy()
    y_score = pd.to_numeric(x[score_col], errors="coerce").fillna(0.0).to_numpy()

    # 当 gold / non-gold 只有一类时，AUROC / AUPRC 不可定义
    if len(np.unique(y_true)) < 2:
        auroc = np.nan
        auprc = np.nan
    else:
        auroc = roc_auc_score(y_true, y_score)
        auprc = average_precision_score(y_true, y_score)

    called = pd.to_numeric(x[p_col], errors="coerce").fillna(1.0).to_numpy() < alpha
    tp = int(((called == 1) & (y_true == 1)).sum())
    fp = int(((called == 1) & (y_true == 0)).sum())
    fn = int(((called == 0) & (y_true == 1)).sum())
    tn = int(((called == 0) & (y_true == 0)).sum())

    tpr = tp / (tp + fn) if (tp + fn) > 0 else np.nan
    fpr = fp / (fp + tn) if (fp + tn) > 0 else np.nan

    return {
        "AUROC": auroc,
        "AUPRC": auprc,
        "TPR": tpr,
        "FPR": fpr,
        "TP": tp,
        "FP": fp,
        "FN": fn,
        "TN": tn,
        "n_genes": int(x.shape[0]),
    }


def rank_methods(metrics_df: pd.DataFrame) -> pd.DataFrame:
    """
    推荐的综合打分逻辑：
        total_score = 0.35*AUROC + 0.35*AUPRC + 0.20*TPR + 0.10*(1-FPR)

    为什么这样定权重？
    ------------------
    - AUROC: 全局排序能力
    - AUPRC: 稀疏正类（DE gene 通常只占少数）更敏感
    - TPR: 找得到多少真阳性
    - 1-FPR: 惩罚乱报
    """
    out = metrics_df.copy()
    for c in ["AUROC", "AUPRC", "TPR", "FPR"]:
        out[c] = pd.to_numeric(out[c], errors="coerce")
    out["composite_score"] = (
        0.35 * out["AUROC"]
        + 0.35 * out["AUPRC"]
        + 0.20 * out["TPR"]
        + 0.10 * (1.0 - out["FPR"])
    )
    return out.sort_values("composite_score", ascending=False)


# -----------------------------------------------------------------------------
# 绘图
# -----------------------------------------------------------------------------

def plot_metric_vs_noise(
    metrics_df: pd.DataFrame,
    *,
    x_col: str,
    y_col: str,
    method_col: str = "method",
    title: Optional[str] = None,
    out_png: Optional[str | Path] = None,
) -> None:
    """不同噪声强度下的 AUROC / AUPRC 折线图。"""
    fig, ax = plt.subplots(figsize=(8, 5))
    for method, sub in metrics_df.groupby(method_col):
        sub = sub.sort_values(x_col)
        ax.plot(sub[x_col], sub[y_col], marker="o", label=method)
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.set_title(title or f"{y_col} vs {x_col}")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    if out_png is not None:
        fig.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_fpr_boxplot(
    metrics_df: pd.DataFrame,
    *,
    method_col: str = "method",
    fpr_col: str = "FPR",
    out_png: Optional[str | Path] = None,
    title: str = "False Positive Rate across scenarios",
) -> None:
    """不同方法在多场景下 FPR 的箱线图。"""
    methods = list(metrics_df[method_col].dropna().unique())
    data = [metrics_df.loc[metrics_df[method_col] == m, fpr_col].dropna().to_numpy() for m in methods]

    fig, ax = plt.subplots(figsize=(max(8, len(methods) * 0.7), 5))
    ax.boxplot(data, labels=methods, patch_artist=False)
    ax.set_ylabel("FPR")
    ax.set_title(title)
    ax.tick_params(axis="x", rotation=45)
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    if out_png is not None:
        fig.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close(fig)


# -----------------------------------------------------------------------------
# 一个最小可运行示例
# -----------------------------------------------------------------------------

def minimal_example(
    adata_path: str | Path,
    outdir: str | Path,
    *,
    group_col: str,
    case_label: str,
    control_label: str,
    label_col_for_noise: Optional[str] = None,
):
    """
    给你一个可直接照改的最小主流程：
    1) 读取 h5ad
    2) QC + 下采样
    3) 构建 counts/lognorm
    4) 注入 dropout / gaussian / mislabel 三种噪声
    5) 跑 scanpy_wilcoxon 与 scanpy_ttest
    6) 导出桥接文件给 R
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    adata = load_benchmark_adata(adata_path)
    adata = qc_and_downsample(adata, target_n_cells=5000)
    adata = make_native_and_bridge_layers(adata, n_top_hvgs=3000)

    # 噪声版本
    d1 = inject_dropout_noise(adata, dropout_rate=0.2, low_expr_quantile=0.3)
    d2 = inject_gaussian_noise(adata, cell_fraction=0.2, sigma=0.25)
    if label_col_for_noise is not None:
        d3 = inject_mislabeling(adata, label_col=label_col_for_noise, mislabel_rate=0.1)
    else:
        d3 = adata.copy()

    for name, x in {
        "clean": adata,
        "dropout20": d1,
        "gaussian20": d2,
        "mislabel10": d3,
    }.items():
        export_bridge_files(x, outdir / name, prefix=name)

        # Scanpy 默认 t-test / Wilcoxon
        for method in ["wilcoxon", "t-test"]:
            res = run_scanpy_rank_genes_groups(
                x,
                group_col=group_col,
                case_label=case_label,
                control_label=control_label,
                method=method,
            )
            res.to_csv(outdir / f"{name}.{method}.csv", index=False)


if __name__ == "__main__":
    # 这里不强行做复杂 CLI，避免你改起来麻烦。
    # 推荐直接在 notebook / python 脚本中 import 本文件中的函数。
    example_input = r"E:\data\SCMixology_5cellline_extension.RData"
    print(
        "benchmark_scanpy_pipeline.py 已加载。\n"
        "现在 adata_path 同时支持 .h5ad 和 .RData。\n"
        f"示例输入文件：{example_input}\n"
        "建议用法：import 本脚本中的函数，在 notebook 中按数据集逐个调用。"
    )
