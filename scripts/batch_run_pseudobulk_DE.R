# ============================================
# 批量运行 pseudo-bulk DE：pairwise edgeR / DESeq2 / limma-voom
# 适合多 condition 数据集（例如 5 个 cell line）
# ============================================

suppressPackageStartupMessages({
  library(edgeR)
  library(DESeq2)
  library(limma)
})

ROOT_DIR <- "D:/scRNA_benchmark/benchmark_runs"

TARGET_DATASETS <- c(
  "SCMixology_5cl_CELseq2_p1",
  "SCMixology_5cl_CELseq2_p2",
  "SCMixology_5cl_CELseq2_p3",
  "SCMixology_RNAmix_sce2",
  "SCMixology_RNAmix_sce8"
)

# ----------------------------
# 读取 pseudo-bulk 输入
# ----------------------------
read_pb_input <- function(dataset_dir) {
  counts_path <- file.path(dataset_dir, "pseudobulk_counts.csv")
  meta_path   <- file.path(dataset_dir, "pseudobulk_metadata.csv")

  if (!file.exists(counts_path) || !file.exists(meta_path)) {
    stop("缺少 pseudo-bulk 输入文件: ", dataset_dir)
  }

  counts <- read.csv(counts_path, row.names = 1, check.names = FALSE)
  meta <- read.csv(meta_path, stringsAsFactors = FALSE)

  rownames(meta) <- meta$sample
  counts <- counts[meta$sample, , drop = FALSE]

  list(counts = counts, meta = meta)
}

# ----------------------------
# 过滤 pseudo-bulk 样本
# 这里把 min_cells 改小，否则第二个 replicate 会被过滤掉
# ----------------------------
filter_pb_samples <- function(counts, meta, min_cells = 3, min_reps_per_condition = 2) {
  meta2 <- meta[meta$n_cells >= min_cells, , drop = FALSE]
  if (nrow(meta2) == 0) return(NULL)

  counts2 <- counts[meta2$sample, , drop = FALSE]

  if (length(unique(meta2$condition)) < 2) return(NULL)

  cond_tab <- table(meta2$condition)
  keep_conditions <- names(cond_tab)[cond_tab >= min_reps_per_condition]

  meta2 <- meta2[meta2$condition %in% keep_conditions, , drop = FALSE]
  if (nrow(meta2) == 0) return(NULL)

  counts2 <- counts2[meta2$sample, , drop = FALSE]

  if (length(unique(meta2$condition)) < 2) return(NULL)

  list(counts = counts2, meta = meta2)
}

# ----------------------------
# 单次 pairwise edgeR
# ----------------------------
run_edger_pair <- function(counts, meta) {
  y <- DGEList(counts = t(counts))
  keep <- filterByExpr(y, group = meta$condition)
  y <- y[keep, , keep.lib.sizes = FALSE]

  y <- calcNormFactors(y)
  design <- model.matrix(~ condition, data = meta)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = 2)

  res <- topTags(qlf, n = Inf)$table
  res$gene <- rownames(res)
  res <- res[, c("gene", setdiff(colnames(res), "gene"))]
  res
}

# ----------------------------
# 单次 pairwise DESeq2
# ----------------------------
run_deseq2_pair <- function(counts, meta, cond1, cond2) {
  dds <- DESeqDataSetFromMatrix(
    countData = t(counts),
    colData = meta,
    design = ~ condition
  )

  dds <- dds[rowSums(counts(dds)) > 1, ]
  dds <- DESeq(dds)

  res <- results(dds, contrast = c("condition", cond2, cond1))
  res <- as.data.frame(res)
  res$gene <- rownames(res)
  res <- res[, c("gene", setdiff(colnames(res), "gene"))]
  res
}

# ----------------------------
# 单次 pairwise limma-voom
# ----------------------------
run_limma_pair <- function(counts, meta) {
  y <- DGEList(counts = t(counts))
  keep <- filterByExpr(y, group = meta$condition)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)

  design <- model.matrix(~ condition, data = meta)
  v <- voom(y, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)

  res <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
  res$gene <- rownames(res)
  res <- res[, c("gene", setdiff(colnames(res), "gene"))]
  res
}

# ----------------------------
# 主流程
# ----------------------------
summary_rows <- list()

for (dataset_name in TARGET_DATASETS) {
  cat("\n", strrep("=", 80), "\n", sep = "")
  cat("处理数据集:", dataset_name, "\n")

  dataset_dir <- file.path(ROOT_DIR, dataset_name)
  if (!dir.exists(dataset_dir)) {
    cat("目录不存在，跳过:", dataset_dir, "\n")
    next
  }

  outdir <- file.path(dataset_dir, "r_pseudobulk_de_pairwise")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  obj <- read_pb_input(dataset_dir)
  counts <- obj$counts
  meta <- obj$meta

  filtered <- filter_pb_samples(counts, meta, min_cells = 3, min_reps_per_condition = 2)
  if (is.null(filtered)) {
    cat("过滤后没有足够的 condition/replicates，跳过。\n")
    next
  }

  counts <- filtered$counts
  meta <- filtered$meta
  meta$condition <- as.factor(meta$condition)
  meta$platform  <- as.factor(meta$platform)

  cat("过滤后样本数:", nrow(meta), "\n")
  cat("condition 分布:\n")
  print(table(meta$condition))

  # 保存过滤后的 metadata
  write.csv(meta, file.path(outdir, "pseudobulk_metadata_filtered.csv"), row.names = FALSE)

  cond_levels <- sort(unique(as.character(meta$condition)))
  pair_list <- combn(cond_levels, 2, simplify = FALSE)

  pair_summary <- list()

  for (pair in pair_list) {
    cond1 <- pair[1]
    cond2 <- pair[2]

    cat("\n--- pairwise:", cond1, "vs", cond2, "---\n")

    sub_meta <- meta[meta$condition %in% c(cond1, cond2), , drop = FALSE]
    sub_counts <- counts[sub_meta$sample, , drop = FALSE]

    sub_meta$condition <- factor(sub_meta$condition, levels = c(cond1, cond2))

    # 每组至少 2 个 replicate
    cond_tab <- table(sub_meta$condition)
    if (any(cond_tab < 2)) {
      cat("该 pair replicates 不足，跳过。\n")
      next
    }

    pair_prefix <- paste0(cond1, "_vs_", cond2)

    # edgeR
    edger_res <- tryCatch(
      run_edger_pair(sub_counts, sub_meta),
      error = function(e) {
        cat("edgeR 失败:", conditionMessage(e), "\n")
        return(NULL)
      }
    )
    if (!is.null(edger_res)) {
      write.csv(edger_res,
                file.path(outdir, paste0("edgeR_", pair_prefix, ".csv")),
                row.names = FALSE)
    }

    # DESeq2
    deseq2_res <- tryCatch(
      run_deseq2_pair(sub_counts, sub_meta, cond1, cond2),
      error = function(e) {
        cat("DESeq2 失败:", conditionMessage(e), "\n")
        return(NULL)
      }
    )
    if (!is.null(deseq2_res)) {
      write.csv(deseq2_res,
                file.path(outdir, paste0("DESeq2_", pair_prefix, ".csv")),
                row.names = FALSE)
    }

    # limma
    limma_res <- tryCatch(
      run_limma_pair(sub_counts, sub_meta),
      error = function(e) {
        cat("limma-voom 失败:", conditionMessage(e), "\n")
        return(NULL)
      }
    )
    if (!is.null(limma_res)) {
      write.csv(limma_res,
                file.path(outdir, paste0("limma_voom_", pair_prefix, ".csv")),
                row.names = FALSE)
    }

    pair_summary[[pair_prefix]] <- data.frame(
      dataset = dataset_name,
      comparison = pair_prefix,
      n_samples = nrow(sub_meta),
      cond1_n = unname(cond_tab[1]),
      cond2_n = unname(cond_tab[2]),
      stringsAsFactors = FALSE
    )
  }

  if (length(pair_summary) > 0) {
    pair_df <- do.call(rbind, pair_summary)
    write.csv(pair_df, file.path(outdir, "pairwise_summary.csv"), row.names = FALSE)
    summary_rows[[dataset_name]] <- pair_df
  }
}

if (length(summary_rows) > 0) {
  summary_df <- do.call(rbind, summary_rows)
  write.csv(summary_df, file.path(ROOT_DIR, "r_pseudobulk_pairwise_summary.csv"), row.names = FALSE)
  cat("\n已保存总汇总文件:\n")
  print(summary_df)
} else {
  cat("\n没有成功处理的 pairwise 结果。\n")
}