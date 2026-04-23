suppressPackageStartupMessages({
  library(edgeR)
  library(DESeq2)
  library(limma)
})

ROOT <- "D:/scRNA_benchmark/benchmark_runs/SCMixology_5cl_CELseq2_merged"
OUT_DIR <- file.path(ROOT, "covariate_compare")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

counts_path <- file.path(ROOT, "pseudobulk_counts.csv")
meta_path   <- file.path(ROOT, "pseudobulk_metadata.csv")

if (!file.exists(counts_path) || !file.exists(meta_path)) {
  stop("缺少 pseudobulk 输入文件，请先生成 merged dataset 的 pseudobulk_counts.csv 和 pseudobulk_metadata.csv")
}

counts <- read.csv(counts_path, row.names = 1, check.names = FALSE)
meta <- read.csv(meta_path, stringsAsFactors = FALSE)

rownames(meta) <- meta$sample
counts <- counts[meta$sample, , drop = FALSE]

# 基础清理
meta$condition <- as.factor(meta$condition)
meta$batch <- as.factor(meta$batch)

cat("样本数:", nrow(meta), "\n")
cat("condition:\n")
print(table(meta$condition))
cat("batch:\n")
print(table(meta$batch))

# 确保 batch 至少两个水平
has_batch_cov <- nlevels(meta$batch) > 1

# 只保留每个 condition 至少 2 个 pseudo-bulk 样本
cond_tab <- table(meta$condition)
keep_cond <- names(cond_tab)[cond_tab >= 2]
meta <- meta[meta$condition %in% keep_cond, , drop = FALSE]
counts <- counts[meta$sample, , drop = FALSE]

meta$condition <- droplevels(meta$condition)
meta$batch <- droplevels(meta$batch)

if (nlevels(meta$condition) < 2) {
  stop("过滤后 condition 不足两组，无法比较。")
}

cat("\n过滤后样本数:", nrow(meta), "\n")
cat("过滤后 condition:\n")
print(table(meta$condition))
cat("过滤后 batch:\n")
print(table(meta$batch))

# 这里为了演示，先只取前两个 condition 做 no_cov / with_cov 对比
lev <- levels(meta$condition)
g1 <- lev[1]
g2 <- lev[2]

sub_meta <- meta[meta$condition %in% c(g1, g2), , drop = FALSE]
sub_counts <- counts[sub_meta$sample, , drop = FALSE]

sub_meta$condition <- factor(sub_meta$condition, levels = c(g1, g2))
sub_meta$batch <- droplevels(factor(sub_meta$batch))

cat("\n用于协变量比较的 pair:", g1, "vs", g2, "\n")
print(table(sub_meta$condition, sub_meta$batch))

# ---------- edgeR ----------
run_edger <- function(counts, meta, with_cov = FALSE) {
  y <- DGEList(counts = t(counts))
  design <- if (with_cov && nlevels(meta$batch) > 1) {
    model.matrix(~ batch + condition, data = meta)
  } else {
    model.matrix(~ condition, data = meta)
  }

  keep <- filterByExpr(y, group = meta$condition)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)

  coef_id <- ncol(design)
  qlf <- glmQLFTest(fit, coef = coef_id)
  res <- topTags(qlf, n = Inf)$table
  res$gene <- rownames(res)
  res
}

# ---------- DESeq2 ----------
run_deseq2 <- function(counts, meta, with_cov = FALSE) {
  dds <- DESeqDataSetFromMatrix(
    countData = t(counts),
    colData = meta,
    design = if (with_cov && nlevels(meta$batch) > 1) ~ batch + condition else ~ condition
  )

  dds <- dds[rowSums(counts(dds)) > 1, ]
  dds <- DESeq(dds)

  lev <- levels(meta$condition)
  res <- results(dds, contrast = c("condition", lev[2], lev[1]))
  res <- as.data.frame(res)
  res$gene <- rownames(res)
  res
}

# ---------- limma-voom ----------
run_limma <- function(counts, meta, with_cov = FALSE) {
  y <- DGEList(counts = t(counts))
  design <- if (with_cov && nlevels(meta$batch) > 1) {
    model.matrix(~ batch + condition, data = meta)
  } else {
    model.matrix(~ condition, data = meta)
  }

  keep <- filterByExpr(y, group = meta$condition)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  v <- voom(y, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)

  coef_id <- ncol(design)
  res <- topTable(fit, coef = coef_id, number = Inf, sort.by = "P")
  res$gene <- rownames(res)
  res
}

# 保存 metadata，确保 no_cov 和 with_cov 用的是同一批样本
write.csv(sub_meta, file.path(OUT_DIR, "covariate_compare_metadata.csv"), row.names = FALSE)

# no covariate
write.csv(run_edger(sub_counts, sub_meta, FALSE), file.path(OUT_DIR, "edgeR_no_cov.csv"), row.names = FALSE)
write.csv(run_deseq2(sub_counts, sub_meta, FALSE), file.path(OUT_DIR, "DESeq2_no_cov.csv"), row.names = FALSE)
write.csv(run_limma(sub_counts, sub_meta, FALSE), file.path(OUT_DIR, "limma_voom_no_cov.csv"), row.names = FALSE)

# with covariate
if (has_batch_cov && nlevels(sub_meta$batch) > 1) {
  write.csv(run_edger(sub_counts, sub_meta, TRUE), file.path(OUT_DIR, "edgeR_with_cov.csv"), row.names = FALSE)
  write.csv(run_deseq2(sub_counts, sub_meta, TRUE), file.path(OUT_DIR, "DESeq2_with_cov.csv"), row.names = FALSE)
  write.csv(run_limma(sub_counts, sub_meta, TRUE), file.path(OUT_DIR, "limma_voom_with_cov.csv"), row.names = FALSE)
} else {
  cat("\n当前 pair 没有足够 batch 水平，跳过 with_cov 版本。\n")
}