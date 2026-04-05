suppressPackageStartupMessages({
  library(presto)
  library(sceptre)
  library(sceptredata)
  library(readr)
  library(dplyr)
  library(Matrix)
  library(parallel)
  library(BH)
  library(Rcpp)
  library(peakRAM)
  library(sgt)
  library(doParallel)
  library(foreach)
  library(speedglm)
  library(matrixStats)
  library(bigmemory)
})


suppressMessages(sourceCpp("/project/xuanyao/jiaming/Getting_started/code/crt_sample.cpp"))
suppressMessages(source("/project/xuanyao/jiaming/Getting_started/code/functions_optimized.R"))


combined_test_tensor <- readRDS(snakemake@input[["combined_test_result"]])
gene_module_index_matrix <- readRDS(snakemake@input[["gene_module_index_matrix"]])


max_statistics <- (-(min_p_val_fast_low_ram(-test_normalizer(combined_test_tensor), gene_module_index_matrix)))


print(snakemake@output[["max_combined_statistics"]])
dir.create(dirname(snakemake@output[[1]]), recursive = TRUE, showWarnings = FALSE)

saveRDS(max_statistics, snakemake@output[[1]])

print("finished!")

