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

# Read inputs
chunk <- readRDS(snakemake@input[["chunk"]])
covariate <- readRDS(snakemake@input[["covariate"]])

# Read params
n_permute <- snakemake@params[["n_permute"]]
use_resample <- snakemake@params[["use_resample"]]

grna_assignment_matrix_combined <- permutation_or_resampling(chunk, n_permute, covariate, use_resample)

# Save output
saveRDS(grna_assignment_matrix_combined, snakemake@output[[1]])

print("finished!")
