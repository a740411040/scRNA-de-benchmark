# ==============================
# 批量导出 SCMixology bridge files
# ==============================

setwd("D:/scRNA_benchmark")

library(SingleCellExperiment)
library(Matrix)

dir.create("D:/scRNA_benchmark/bridge", showWarnings = FALSE, recursive = TRUE)
dir.create("D:/scRNA_benchmark/truth_tables", showWarnings = FALSE, recursive = TRUE)

# --------------------------------
# 通用函数：导出 SingleCellExperiment 为 bridge files
# --------------------------------
export_sce_bridge <- function(sce, prefix, outdir = "D:/scRNA_benchmark/bridge") {
  cat("正在导出:", prefix, "\n")

  cnt <- counts(sce)
  if (!inherits(cnt, "dgCMatrix")) {
    cnt <- as(cnt, "dgCMatrix")
  }

  # counts matrix
  writeMM(cnt, file = file.path(outdir, paste0(prefix, "_counts.mtx")))

  # obs metadata
  obs <- as.data.frame(colData(sce))
  obs$cell_id <- colnames(sce)
  write.csv(
    obs,
    file = file.path(outdir, paste0(prefix, "_obs.csv")),
    row.names = FALSE
  )

  # var metadata
  var <- as.data.frame(rowData(sce))
  var$gene_id <- rownames(sce)
  write.csv(
    var,
    file = file.path(outdir, paste0(prefix, "_var.csv")),
    row.names = FALSE
  )

  # gene names
  write.csv(
    data.frame(gene_id = rownames(sce)),
    file = file.path(outdir, paste0(prefix, "_genes.csv")),
    row.names = FALSE
  )

  # cell names
  write.csv(
    data.frame(cell_id = colnames(sce)),
    file = file.path(outdir, paste0(prefix, "_cells.csv")),
    row.names = FALSE
  )

  cat("导出完成:", prefix, "\n\n")
}

# --------------------------------
# 1. SCMixology_5cellline_extension.RData
# --------------------------------
rm(list = ls()[!ls() %in% c("export_sce_bridge")])
load("D:/scRNA_benchmark/raw_data/SCMixology_5cellline_extension.RData")

export_sce_bridge(sc_Celseq2_5cl_p1, "SCMixology_5cl_CELseq2_p1")
export_sce_bridge(sc_Celseq2_5cl_p2, "SCMixology_5cl_CELseq2_p2")
export_sce_bridge(sc_Celseq2_5cl_p3, "SCMixology_5cl_CELseq2_p3")
export_sce_bridge(sce_sc_10x_5cl_qc, "SCMixology_5cl_10x")

# --------------------------------
# 2. SCMixology_9cellmix.RData
# --------------------------------
rm(list = ls()[!ls() %in% c("export_sce_bridge")])
load("D:/scRNA_benchmark/raw_data/SCMixology_9cellmix.RData")

export_sce_bridge(sce_9cells_qc,  "SCMixology_9cellmix_main")
export_sce_bridge(sce_POP_sel_qc, "SCMixology_9cellmix_POP_sel")
export_sce_bridge(sce_SC1_qc,     "SCMixology_9cellmix_SC1")
export_sce_bridge(sce_SC2_qc,     "SCMixology_9cellmix_SC2")
export_sce_bridge(sce_SC3_qc,     "SCMixology_9cellmix_SC3")
export_sce_bridge(sce_SC4_qc,     "SCMixology_9cellmix_SC4")

# --------------------------------
# 3. SCMixology_RNAmix.RData
# --------------------------------
rm(list = ls()[!ls() %in% c("export_sce_bridge")])
load("D:/scRNA_benchmark/raw_data/SCMixology_RNAmix.RData")

export_sce_bridge(sce2_qc, "SCMixology_RNAmix_sce2")
export_sce_bridge(sce8_qc, "SCMixology_RNAmix_sce8")

# --------------------------------
# 4. SCMixology_90cell_DEtruth.RData
# --------------------------------
rm(list = ls()[!ls() %in% c("export_sce_bridge")])
load("D:/scRNA_benchmark/raw_data/SCMixology_90cell_DEtruth.RData")

write.csv(H1975_DEtable,
          file = "D:/scRNA_benchmark/truth_tables/H1975_DEtable.csv",
          row.names = FALSE)

write.csv(H2228_DEtable,
          file = "D:/scRNA_benchmark/truth_tables/H2228_DEtable.csv",
          row.names = FALSE)

write.csv(HCC827_DEtable,
          file = "D:/scRNA_benchmark/truth_tables/HCC827_DEtable.csv",
          row.names = FALSE)

cat("90cell_DEtruth 已导出为 csv。\n\n")

# --------------------------------
# 5. 打印 bridge 文件概况
# --------------------------------
bridge_files <- list.files("D:/scRNA_benchmark/bridge", full.names = TRUE)
truth_files <- list.files("D:/scRNA_benchmark/truth_tables", full.names = TRUE)

cat("bridge 文件数量:", length(bridge_files), "\n")
cat("truth 文件数量:", length(truth_files), "\n\n")

print(head(bridge_files, 20))
print(truth_files)