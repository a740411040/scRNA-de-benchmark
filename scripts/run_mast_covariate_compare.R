suppressPackageStartupMessages({
  library(Seurat)
  library(MAST)
  library(Matrix)
})

ROOT_BRIDGE <- "D:/scRNA_benchmark/bridge"
OUT_DIR <- "D:/scRNA_benchmark/benchmark_runs/SCMixology_5cl_CELseq2_merged/covariate_compare"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

counts <- readMM(file.path(ROOT_BRIDGE, "SCMixology_5cl_CELseq2_merged_counts.mtx"))
obs <- read.csv(file.path(ROOT_BRIDGE, "SCMixology_5cl_CELseq2_merged_obs.csv"), stringsAsFactors = FALSE)
genes <- read.csv(file.path(ROOT_BRIDGE, "SCMixology_5cl_CELseq2_merged_genes.csv"), stringsAsFactors = FALSE)
cells <- read.csv(file.path(ROOT_BRIDGE, "SCMixology_5cl_CELseq2_merged_cells.csv"), stringsAsFactors = FALSE)

# ---------- 关键修复：不要直接用重复 cell_id 做 rownames ----------
raw_cell_ids <- as.character(cells$cell_id)

# 用顺序号强制生成唯一列名，避免重复
unique_cell_ids <- paste0(raw_cell_ids, "__", seq_along(raw_cell_ids))

# 保证维度一致
stopifnot(nrow(obs) == length(unique_cell_ids))
stopifnot(ncol(counts) == length(unique_cell_ids))
stopifnot(nrow(counts) == nrow(genes))

# 同步到 obs / counts
obs$cell_id_raw <- raw_cell_ids
obs$cell_id_unique <- unique_cell_ids
rownames(obs) <- unique_cell_ids

rownames(counts) <- as.character(genes$gene_id)
colnames(counts) <- unique_cell_ids

# 创建对象
seu <- CreateSeuratObject(counts = counts, meta.data = obs)
seu <- NormalizeData(seu, verbose = FALSE)

# 统一 condition
if ("cell_line_demuxlet" %in% colnames(seu@meta.data)) {
  seu$condition <- seu$cell_line_demuxlet
} else if ("cell_line" %in% colnames(seu@meta.data)) {
  seu$condition <- seu$cell_line
} else {
  stop("找不到 condition 列（cell_line_demuxlet / cell_line）")
}

# 如果没有 batch，就报错，因为你现在就是要做协变量比较
if (!"batch" %in% colnames(seu@meta.data)) {
  stop("obs 中没有 batch 列，无法做协变量比较。")
}

Idents(seu) <- "condition"

run_mast_one <- function(seu, with_cov = FALSE) {
  lev <- sort(unique(as.character(seu$condition)))
  g1 <- lev[1]
  g2 <- lev[2]

  sub <- subset(seu, idents = c(g1, g2))
  Idents(sub) <- "condition"

  res <- FindMarkers(
    sub,
    ident.1 = g1,
    ident.2 = g2,
    test.use = "MAST",
    logfc.threshold = 0,
    min.pct = 0,
    latent.vars = if (with_cov) "batch" else NULL
  )
  res$gene <- rownames(res)
  res
}

write.csv(run_mast_one(seu, FALSE), file.path(OUT_DIR, "MAST_no_cov.csv"), row.names = FALSE)
write.csv(run_mast_one(seu, TRUE),  file.path(OUT_DIR, "MAST_with_cov.csv"), row.names = FALSE)