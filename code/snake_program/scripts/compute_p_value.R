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


suppressMessages(sourceCpp("/project/xuanyao/jiaming/Getting_started/code/fit_skew_normal.cpp"))
suppressMessages(sourceCpp("/project/xuanyao/jiaming/Getting_started/code/crt_sample.cpp"))
suppressMessages(source("/project/xuanyao/jiaming/Getting_started/code/functions_optimized.R"))


combined_test_tensor <- readRDS(snakemake@input[["combined_test_statistics"]])
normalized_test_tensor <- readRDS(snakemake@input[["normalized_test_tensor"]])
gene_module_index_matrix <- readRDS(snakemake@input[["gene_module_index_matrix"]])
guide<-readRDS(snakemake@input[["guide"]])


n_permute <-snakemake@params[["n_permute"]]
use_approximation <- snakemake@params[["use_approximation"]]
use_approximation_gene <- snakemake@params[["use_approximation_gene"]]
approximation_type <- snakemake@params[["approximation_type"]]
approximation_threshold <- snakemake@params[["approximation_threshold"]]

max_statistics <- (-(min_p_val_fast_low_ram(-test_normalizer(combined_test_tensor), gene_module_index_matrix)))

precompute_p <- compute_final_p_val(max_statistics, n_permute,1)
genewise_p<-compute_final_p_val(normalized_test_tensor,n_permute,3)


if (use_approximation) {
  threshold <- approximation_threshold

    indices <- which(!is.na(precompute_p) & precompute_p < threshold, arr.ind = TRUE)

    if (length(indices) > 0) {
      for (k in seq_len(nrow(indices))) {
        i <- indices[k, 2]
        j <- indices[k, 1]

        p_value_1 <- max_statistics[-1, i, j]
        p_value_2 <- max_statistics[1, i, j]
        p_value_1 <- p_value_1[is.finite(p_value_1)]


        if (length(p_value_1) <= 1 || length(unique(p_value_1)) <= 1) {
          result <- 1
        } else {

          if (approximation_type=="normal"){
            result <- fit_and_evaluate_skew_normal(p_value_1, p_value_2, 1, TRUE)
          } else if (approximation_type=="t"){
            result <- fit_and_evaluate_skew_t(p_value_1,p_value_2,1)
          }
        }
        precompute_p[j, i] <- result
      }
    } else {
      print("No indices meet the condition precompute_p < threshold")
    }
}

if (use_approximation_gene) {

  threshold <- approximation_threshold

    indices <- which(!is.na(genewise_p) & genewise_p < threshold, arr.ind = TRUE)

    if (length(indices) > 0) {
      for (k in seq_len(nrow(indices))) {
        i <- indices[k, 2]
        j <- indices[k, 1]

        p_value_1 <- normalized_test_tensor[-1, i, j]
        p_value_2 <- normalized_test_tensor[1, i, j]
        p_value_1 <- p_value_1[is.finite(p_value_1)]


        if (length(p_value_1) <= 1 || length(unique(p_value_1)) <= 1) {
          result <- 1
        } else {

          if (approximation_type=="normal"){
            result <- fit_and_evaluate_skew_normal(p_value_1, p_value_2, 3, TRUE)
          } else if (approximation_type=="t"){
            result <- fit_and_evaluate_skew_t(p_value_1,p_value_2,3)
          }
        }
        genewise_p[j, i] <- result

      }
    } else {
      print("No indices meet the condition precompute_p < threshold")
    }

}

rownames(genewise_p)<- rownames(guide)
rownames(precompute_p)<- rownames(guide)
colnames(genewise_p) <- rownames(gene_module_index_matrix)

print(snakemake@output[["geneset_p_value"]])
dir.create(dirname(snakemake@output[["geneset_p_value"]]), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(snakemake@output[["genewise_p_value"]]), recursive = TRUE, showWarnings = FALSE)

saveRDS(precompute_p, snakemake@output[["geneset_p_value"]])
saveRDS(genewise_p, snakemake@output[["genewise_p_value"]])

print("finished!")

