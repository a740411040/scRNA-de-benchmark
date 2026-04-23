# scRNA-seq Differential Expression Benchmark

A benchmarking project for **single-cell RNA-seq differential expression (DE) analysis methods**, covering both **native single-cell approaches** and **pseudo-bulk approaches**, with an additional comparison of **batch covariate modeling**.

## Overview

This project systematically evaluates the performance of multiple scRNA-seq DE methods under a unified benchmarking framework. The benchmark focuses on three key questions:

1. Which method provides the best overall balance between sensitivity and false-positive control?
2. How do native single-cell methods compare with pseudo-bulk methods?
3. Does adding a batch covariate consistently improve performance?

The benchmark includes:

- **Native single-cell methods**: Scanpy t-test, Scanpy Wilcoxon, Seurat Wilcoxon, MAST  
- **Pseudo-bulk methods**: edgeR, DESeq2, limma-voom  
- **Covariate comparison**: no covariate vs. batch-adjusted models  
- **Stability analysis**: robustness across datasets and comparisons  

---

## Key Findings

- **Best overall method:** **edgeR**  
  It achieved the best balance across ranking quality, recall, false-positive control, and stability.

- **Best default method in Python/Scanpy workflows:** **Scanpy Wilcoxon**  
  It was consistently more stable and reliable than Scanpy t-test.

- **High-recall but more aggressive methods:** **Seurat Wilcoxon** and **MAST**  
  These methods recovered more true positives, but also produced higher false-positive rates.

- **Most batch-sensitive method:** **limma-voom**  
  It showed the largest performance gain after adding batch as a covariate.

---

## Benchmark Design

### Datasets

The main benchmark is based on the **SCMixology 5-cell-line CELseq2** datasets:

- `SCMixology_5cl_CELseq2_p1`
- `SCMixology_5cl_CELseq2_p2`
- `SCMixology_5cl_CELseq2_p3`

These datasets were selected because they:

- have clear cell line labels,
- support pairwise DE comparison,
- can be matched to a gold-standard truth table.

For batch covariate analysis, the three CELseq2 batches were merged into:

- `SCMixology_5cl_CELseq2_merged`

Gold-standard DE truth was derived from:

- `SCMixology_90cell_DEtruth`

Additional datasets such as `SCMixology_RNAmix_*` and `SCMixology_10x*` were used for supplementary exploration, but not as the main pseudo-bulk scoring datasets.

---

## Methods Evaluated

| Category | Method | Input | Notes |
| :-- | :-- | :-- | :-- |
| Native single-cell | Scanpy t-test | log-normalized matrix | Pairwise DE in Scanpy |
| Native single-cell | Scanpy Wilcoxon | log-normalized matrix | More robust Scanpy default |
| Native single-cell | Seurat Wilcoxon | Seurat log-normalized matrix | `FindMarkers(test.use = "wilcox")` |
| Native single-cell / regression | MAST | Seurat log-normalized matrix | Supports batch covariates |
| Pseudo-bulk | edgeR | aggregated raw counts | QLF test |
| Pseudo-bulk | DESeq2 | aggregated raw counts | standard generalized model |
| Pseudo-bulk | limma-voom | aggregated raw counts | batch-sensitive in this benchmark |

---

## Preprocessing

To ensure fair comparison, all datasets were processed using a unified pipeline:

- preserve raw count matrix in `adata.layers["counts"]`
- filter cells with `n_genes_by_counts >= 200`
- filter cells with `pct_counts_mt <= 20`
- keep genes expressed in at least 3 cells
- downsample to 5000 cells if needed
- normalize with `normalize_total(target_sum=1e4)` + `log1p`
- identify 3000 highly variable genes
- perform PCA / neighbors / UMAP / Leiden for structure inspection

For benchmarking:

- **native single-cell methods** used the log-normalized matrix
- **pseudo-bulk methods** used aggregated raw counts

---

## Evaluation Metrics

The benchmark uses the following metrics:

- **AUROC**: global ranking performance
- **AUPRC**: enrichment of true positives under class imbalance
- **TPR**: true positive recovery rate
- **FPR**: false positive rate

Predicted positives were defined as the **top 5% genes** ranked by each method-specific score.

### Benchmark Score

The composite benchmark score is defined as:

```text
0.35 × AUROC + 0.35 × AUPRC + 0.20 × TPR + 0.10 × (1 − FPR)
```

This weighting emphasizes ranking quality while still accounting for recall and false-positive control.

---

## Main Benchmark Results

### Overall ranking

| Rank | Method | Benchmark Score | Interpretation |
| :-- | :-- | :--: | :-- |
| 1 | **edgeR** | **0.431** | Best overall balance |
| 2 | **Seurat Wilcoxon** | 0.365 | High recall, more aggressive |
| 3 | **MAST** | 0.363 | Similar to Seurat, supports regression |
| 4 | **Scanpy Wilcoxon** | 0.357 | Best Scanpy-native default |
| 5 | **Scanpy t-test** | 0.351 | Moderate baseline |
| 6 | **DESeq2** | 0.257 | Weak recovery under current setting |
| 7 | **limma-voom** | 0.231 | Underperforms without covariate adjustment |

### Interpretation

- **edgeR** showed the strongest overall performance, with the best AUPRC and strong control of false positives.
- **Seurat Wilcoxon** and **MAST** achieved higher recall but were more likely to report false positives.
- **Scanpy Wilcoxon** consistently outperformed Scanpy t-test and is the preferred default within Python-first workflows.
- **DESeq2** and **limma-voom** were less competitive in the main benchmark under the current setup.

---

## Batch Covariate Analysis

To assess the effect of explicit batch modeling, the three CELseq2 batches were merged and analyzed with and without `batch` as a covariate.

### Score changes after adding batch

| Method | No Covariate | With Batch Covariate | ΔScore | Interpretation |
| :-- | :--: | :--: | :--: | :-- |
| edgeR | 0.497 | 0.490 | -0.007 | already robust |
| DESeq2 | 0.396 | 0.389 | -0.007 | little benefit |
| **limma-voom** | 0.385 | **0.457** | **+0.071** | strongest gain |
| MAST | 0.390 | 0.393 | +0.003 | slight improvement |

### Interpretation

- **limma-voom** benefited the most from batch correction, indicating that its performance was strongly affected by technical noise when batch was not modeled.
- **MAST** showed a small improvement.
- **edgeR** and **DESeq2** did not gain additional advantage from batch modeling in this benchmark.

This suggests that **covariate modeling should be treated as a method-dependent option rather than a universally beneficial default**.

---

## Stability Analysis

To quantify robustness across datasets and comparisons, a `stability_score` was calculated based on variation in:

- benchmark score
- AUROC
- FPR

A combined metric was then defined as:

```text
robust_overall_score = 0.7 × benchmark_score + 0.3 × stability_score
```

### Robust overall ranking

| Method | Benchmark Score | Stability Score | Robust Overall Score |
| :-- | :--: | :--: | :--: |
| **edgeR** | 0.431 | 0.423 | **0.429** |
| MAST | 0.363 | 0.459 | 0.392 |
| Seurat Wilcoxon | 0.365 | 0.437 | 0.387 |
| Scanpy t-test | 0.351 | 0.426 | 0.374 |
| Scanpy Wilcoxon | 0.357 | 0.395 | 0.368 |
| DESeq2 | 0.257 | 0.408 | 0.302 |
| limma-voom | 0.231 | 0.199 | 0.222 |

**Conclusion:** edgeR remained the top method even after jointly considering performance and stability.

---

## Recommended Usage

| Scenario | Recommended Method | Rationale |
| :-- | :-- | :-- |
| Replicates available, want best overall balance | **edgeR** | strongest combined performance |
| Strong concern about false positives | **edgeR** / **Scanpy Wilcoxon** | more conservative |
| High-recall candidate screening | **Seurat Wilcoxon** / **MAST** | more aggressive recovery |
| Python-first workflow without reliable replicates | **Scanpy Wilcoxon** | best native Scanpy choice |
| Strong suspected batch effect | compare **edgeR / limma-voom / MAST** covariate models | benefit is method-dependent |

---

## Main Output Files

Key outputs generated in this project include:

- `benchmark_metrics.csv`  
  Main benchmark metrics across methods

- `benchmark_metrics_no_covariate.csv`  
  Metrics without batch covariate

- `benchmark_metrics_with_covariate.csv`  
  Metrics with batch covariate

- `method_stability_summary.csv`  
  Stability and robustness summary

- `fig1_overall_benchmark_score.png` to `fig5_dataset_method_heatmap.png`  
  Main benchmark figures

- `covariate_*`  
  Covariate comparison figures and delta plots

---

## Take-Home Messages

- **edgeR** is the best default pseudo-bulk method in this benchmark.
- **Scanpy Wilcoxon** is the best default native method for Python workflows.
- **Seurat Wilcoxon** and **MAST** are useful when recall matters more than strict false-positive control.
- **Batch correction is not universally beneficial**; its effect depends strongly on the method.
- **limma-voom** is the method that benefits most from explicit batch modeling in this benchmark.

---

## Citation

If you use this benchmark framework, results, or evaluation design in your own work, please cite the accompanying report or repository as appropriate.

---

## Notes

This README is a GitHub-oriented summary of the full benchmark report. For detailed experimental design, truth mapping, pairwise comparisons, and figure interpretation, please refer to the full report.
