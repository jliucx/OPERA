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

normalized_test_tensor <- readRDS(snakemake@input[["normalized_test_result"]])
gene_module_index_matrix <- readRDS(snakemake@input[["gene_module_index_matrix"]])
advanced <- snakemake@params[["advanced"]]


  if(advanced){
    combined_test_tensor <- combine_test_statistics_tensor_advanced(normalized_test_tensor, gene_module_index_matrix)
  } else {
    combined_test_tensor <- combine_test_statistics_tensor_low_RAM(normalized_test_tensor, gene_module_index_matrix)
  }


print(snakemake@output[["combined_test_result"]])
dir.create(dirname(snakemake@output[[1]]), recursive = TRUE, showWarnings = FALSE)

saveRDS(combined_test_tensor, snakemake@output[[1]])

print("finished!")

