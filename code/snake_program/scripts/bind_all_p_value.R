geneset_files <- snakemake@input[["geneset_files"]]
genewise_files <- snakemake@input[["genewise_files"]]
geneset_output <- snakemake@output[["geneset_merged"]]
genewise_output <- snakemake@output[["genewise_merged"]]

# 排序函数
extract_index <- function(paths) {
  as.integer(gsub("\\D+", "", basename(paths)))
}

geneset_files_sorted <- geneset_files[order(extract_index(geneset_files))]
genewise_files_sorted <- genewise_files[order(extract_index(genewise_files))]

# 合并数据
geneset_df <- do.call(rbind, lapply(geneset_files_sorted, readRDS))
genewise_df <- do.call(rbind, lapply(genewise_files_sorted, readRDS))

# 保存合并结果
saveRDS(geneset_df, geneset_output)
saveRDS(genewise_df, genewise_output)
