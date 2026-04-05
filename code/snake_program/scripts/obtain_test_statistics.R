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




grna_assignment_matrix_combined <- readRDS(snakemake@input[["chunk"]])
grna_assignment_matrix <- readRDS(snakemake@input[["guide"]])
response_matrix <- readRDS(snakemake@input[["response_matrix"]])
n_permute <- snakemake@params[["n_permute"]]
test_type <- snakemake@params[["test_type"]]


if (test_type == "t_test") {
  if (moi == "high") {

      test_result <- t_test_pooled_variance(response_matrix, grna_assignment_matrix_combined)

  } else {

      test_result <- t_test_sparse_lowmoi(response_matrix, grna_assignment_matrix_combined, negative_vector = 0)

  }


} else if (test_type == "wilcoxon_test") {

    test_result <- t_test_wilcoxon(response_matrix, grna_assignment_matrix_combined)

}


  normalized_test_tensor <- test_normalizer_rcpp(
    convert_matrix_to_tensor(test_result,
    n_permute,
    grna_assignment_matrix)
  )


dir.create(dirname(snakemake@output[[1]]), recursive = TRUE, showWarnings = FALSE)

saveRDS(normalized_test_tensor, snakemake@output[[1]])

print("finished!")
