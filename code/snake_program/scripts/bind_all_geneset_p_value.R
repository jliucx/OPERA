geneset_files <- snakemake@input[["geneset_files"]]
geneset_output <- snakemake@output[["geneset_merged"]]

# 排序函数
extract_index <- function(paths) {
  as.integer(gsub("\\D+", "", basename(paths)))
}

geneset_files_sorted <- geneset_files[order(extract_index(geneset_files))]

# 合并数据
geneset_df <- do.call(rbind, lapply(geneset_files_sorted, readRDS))

# 保存合并结果
saveRDS(geneset_df, geneset_output)
