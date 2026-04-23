args <- commandArgs(trailingOnly = TRUE)

cat("收到参数:\n")
print(args)

if (length(args) < 2) {
  stop("用法: Rscript convert_rdata_to_h5ad.R input.RData output.h5ad")
}

input_file <- normalizePath(args[1], mustWork = FALSE)
output_file <- args[2]

cat("输入文件:", input_file, "\n")
cat("输出文件:", output_file, "\n")
cat("当前工作目录:", getwd(), "\n")

if (!file.exists(input_file)) {
  stop(paste("输入文件不存在:", input_file))
}

out_dir <- dirname(output_file)
if (!dir.exists(out_dir)) {
  cat("输出目录不存在，正在创建:", out_dir, "\n")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

cat("开始 load RData...\n")
load(input_file)

objs <- ls()
cat("RData 中对象列表:\n")
print(objs)

target_obj <- NULL
target_name <- NULL

for (x in objs) {
  obj <- get(x)
  cls <- class(obj)
  cat("检查对象:", x, " 类别:", paste(cls, collapse = ", "), "\n")

  if ("Seurat" %in% cls) {
    target_obj <- obj
    target_name <- x
    break
  }

  if ("SingleCellExperiment" %in% cls) {
    target_obj <- obj
    target_name <- x
    break
  }
}

if (is.null(target_obj)) {
  stop("未找到 Seurat 或 SingleCellExperiment 对象。")
}

cat("找到目标对象:", target_name, "\n")
cat("对象类型:", class(target_obj), "\n")

if ("Seurat" %in% class(target_obj)) {
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(SingleCellExperiment))
  suppressPackageStartupMessages(library(zellkonverter))
  sce_obj <- as.SingleCellExperiment(target_obj)
  cat("开始 writeH5AD...\n")
  writeH5AD(sce_obj, output_file)
} else if ("SingleCellExperiment" %in% class(target_obj)) {
  suppressPackageStartupMessages(library(zellkonverter))
  cat("开始 writeH5AD...\n")
  writeH5AD(target_obj, output_file)
} else {
  stop("目标对象既不是 Seurat 也不是 SingleCellExperiment。")
}

cat("writeH5AD 执行完毕。\n")
cat("输出文件是否存在:", file.exists(output_file), "\n")
if (file.exists(output_file)) {
  cat("输出文件大小(bytes):", file.info(output_file)$size, "\n")
}