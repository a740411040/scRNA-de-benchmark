suppressPackageStartupMessages({
  library(Seurat)
  library(MAST)
  library(Matrix)
})

ROOT_BRIDGE <- "D:/scRNA_benchmark/bridge"
ROOT_OUT <- "D:/scRNA_benchmark/benchmark_runs"

TARGET_DATASETS <- c(
  "SCMixology_5cl_CELseq2_p1",
  "SCMixology_5cl_CELseq2_p2",
  "SCMixology_5cl_CELseq2_p3",
  "SCMixology_RNAmix_sce2",
  "SCMixology_RNAmix_sce8"
)

read_bridge_to_seurat <- function(prefix) {
  counts_path <- file.path(ROOT_BRIDGE, paste0(prefix, "_counts.mtx"))
  obs_path <- file.path(ROOT_BRIDGE, paste0(prefix, "_obs.csv"))
  genes_path <- file.path(ROOT_BRIDGE, paste0(prefix, "_genes.csv"))
  cells_path <- file.path(ROOT_BRIDGE, paste0(prefix, "_cells.csv"))

  counts <- readMM(counts_path)
  obs <- read.csv(obs_path, stringsAsFactors = FALSE)
  genes <- read.csv(genes_path, stringsAsFactors = FALSE)
  cells <- read.csv(cells_path, stringsAsFactors = FALSE)

  rownames(obs) <- cells$cell_id
  rownames(counts) <- genes$gene_id
  colnames(counts) <- cells$cell_id

  seu <- CreateSeuratObject(counts = counts, meta.data = obs)

  # 统一 condition
  if ("cell_line_demuxlet" %in% colnames(seu@meta.data)) {
    seu$condition <- seu$cell_line_demuxlet
  } else if ("cell_line" %in% colnames(seu@meta.data)) {
    seu$condition <- seu$cell_line
  } else if ("mix" %in% colnames(seu@meta.data)) {
    seu$condition <- seu$mix
  } else {
    seu$condition <- "unknown"
  }

  # 统一 celltype
  if ("cell_line" %in% colnames(seu@meta.data)) {
    seu$celltype <- seu$cell_line
  } else if ("cell_line_demuxlet" %in% colnames(seu@meta.data)) {
    seu$celltype <- seu$cell_line_demuxlet
  } else {
    seu$celltype <- seu$condition
  }

  # 统一 batch
  if (!"batch" %in% colnames(seu@meta.data)) {
    seu$batch <- "unknown"
  }

  seu
}

run_one_pair <- function(sub, g1, g2, dataset_name, outdir) {
  comp <- paste0(g1, "_vs_", g2)

  # Seurat Wilcoxon
  res_wil <- tryCatch(
    FindMarkers(
      object = sub,
      ident.1 = g1,
      ident.2 = g2,
      test.use = "wilcox",
      logfc.threshold = 0,
      min.pct = 0,
      min.cells.group = 3
    ),
    error = function(e) {
      message("[Seurat Wilcoxon failed] ", dataset_name, " ", comp, " : ", conditionMessage(e))
      return(NULL)
    }
  )

  if (!is.null(res_wil) && nrow(res_wil) > 0) {
    res_wil$gene <- rownames(res_wil)
    res_wil$dataset <- dataset_name
    res_wil$method <- "seurat_wilcoxon"
    res_wil$comparison <- comp
    write.csv(
      res_wil,
      file.path(outdir, paste0("seurat_wilcoxon_", comp, ".csv")),
      row.names = FALSE
    )
  }

  # MAST 需要先 NormalizeData
  res_mast <- tryCatch(
    FindMarkers(
      object = sub,
      ident.1 = g1,
      ident.2 = g2,
      test.use = "MAST",
      logfc.threshold = 0,
      min.pct = 0,
      min.cells.group = 3,
      latent.vars = if ("batch" %in% colnames(sub@meta.data) && length(unique(sub$batch)) > 1) "batch" else NULL
    ),
    error = function(e) {
      message("[MAST failed] ", dataset_name, " ", comp, " : ", conditionMessage(e))
      return(NULL)
    }
  )

  if (!is.null(res_mast) && nrow(res_mast) > 0) {
    res_mast$gene <- rownames(res_mast)
    res_mast$dataset <- dataset_name
    res_mast$method <- "MAST"
    res_mast$comparison <- comp
    write.csv(
      res_mast,
      file.path(outdir, paste0("MAST_", comp, ".csv")),
      row.names = FALSE
    )
  }

  data.frame(
    dataset = dataset_name,
    comparison = comp,
    stringsAsFactors = FALSE
  )
}

run_pairwise_seurat_mast <- function(seu, dataset_name, outdir) {
  # 预处理给 Seurat/MAST 用
  seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000, verbose = FALSE)

  Idents(seu) <- "condition"
  cond_levels <- sort(unique(as.character(seu$condition)))
  cond_levels <- cond_levels[cond_levels != "unknown"]

  if (length(cond_levels) < 2) {
    message("条件不足，跳过: ", dataset_name)
    return(NULL)
  }

  pair_list <- combn(cond_levels, 2, simplify = FALSE)
  summary_list <- list()

  for (pair in pair_list) {
    g1 <- pair[1]
    g2 <- pair[2]

    sub <- subset(seu, idents = c(g1, g2))
    Idents(sub) <- "condition"

    # 每组至少 3 个细胞
    tab <- table(Idents(sub))
    if (length(tab) < 2 || any(tab < 3)) {
      message("细胞数不足，跳过: ", dataset_name, " ", g1, " vs ", g2)
      next
    }

    summary_list[[paste0(g1, "_vs_", g2)]] <- run_one_pair(sub, g1, g2, dataset_name, outdir)
  }

  if (length(summary_list) > 0) {
    pair_df <- do.call(rbind, summary_list)
    write.csv(pair_df, file.path(outdir, "seurat_mast_pairwise_summary.csv"), row.names = FALSE)
    return(pair_df)
  }

  NULL
}

for (dataset_name in TARGET_DATASETS) {
  cat("\n", strrep("=", 80), "\n", sep = "")
  cat("处理数据集:", dataset_name, "\n")

  outdir <- file.path(ROOT_OUT, dataset_name, "r_seurat_mast_pairwise")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  seu <- read_bridge_to_seurat(dataset_name)
  pair_df <- run_pairwise_seurat_mast(seu, dataset_name, outdir)

  if (!is.null(pair_df)) {
    print(pair_df)
  } else {
    cat("没有成功的 Seurat/MAST pairwise 结果:", dataset_name, "\n")
  }
}