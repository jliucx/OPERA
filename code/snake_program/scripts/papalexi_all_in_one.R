suppressPackageStartupMessages({
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
  library(heavytailcombtest)

})


suppressMessages(sourceCpp("/project/xuanyao/jiaming/Getting_started/code/crt_sample.cpp"))
suppressMessages(source("/project/xuanyao/jiaming/Getting_started/code/functions_optimized.R"))


# 读取输入文件
i <- as.integer(snakemake@wildcards[["chunk"]])


grna_target_data_frame_omit <- readRDS(snakemake@input[["grna_target_dataframe"]])
gene_module_index_matrix <- readRDS(snakemake@input[["gene_module_index_matrix"]])
discovery_pairs_cis<- readRDS(snakemake@input[["discovery_pairs_cis"]])
response_matrix_NTC <- readRDS(snakemake@input[["response_matrix_NTC"]])
NTC_guides <- readRDS(snakemake@input[["non_targeting_guide"]])
grna_names <- readRDS(snakemake@input[["grna_names"]])
covariate_NTC <- readRDS(snakemake@input[["covariate_NTC"]])


# Read params
n_permute <- snakemake@params[["n_permute"]]
use_resample <- snakemake@params[["use_resample"]]
trt_cells_cutoff<-snakemake@params[["trt_cells_cutoff"]]
control_cells_cutoff<-snakemake@params[["control_cells_cutoff"]]
use_approximation <- snakemake@params[["use_approximation"]]
use_approximation_gene <- snakemake@params[["use_approximation_gene"]]
approximation_type <- snakemake@params[["approximation_type"]]
approximation_threshold <- snakemake@params[["approximation_threshold"]]
test_type <- snakemake@params[["test_type"]]
use_ACAT <- snakemake@params[["use_ACAT"]]
guide_or_target <- snakemake@params[["guide_or_target"]]
grna_name <- grna_names[i]

########################### STEP 1 Generate guide index matrix #########################
start_step1 <- Sys.time()

if (i %in% NTC_guides) {
  response_matrix <- response_matrix_NTC
  covariate <- covariate_NTC
  NTC_guides_length <- readRDS(snakemake@input[["ntc_guide_length"]])

  # 找到 i 在 NTC_guides 中的位置
  i_pos <- which(NTC_guides == i)
  start <- if (i_pos == 1) 1 else sum(NTC_guides_length[1:(i_pos - 1)]) + 1
  end <- start + NTC_guides_length[i_pos] - 1

  # 直接构造 chunk
  chunk <- sparseMatrix(i = rep(1L, end - start + 1),
                        j = start:end,
                        x = rep(1L, end - start + 1),
                        dims = c(1, ncol(response_matrix)))
} else {
  response_matrix_main <- readRDS(snakemake@input[["response_matrix"]])
  covariate_main <- readRDS(snakemake@input[["covariate"]])
  response_matrix <- cbind(response_matrix_main,response_matrix_NTC)
  covariate <- rbind(covariate_main,covariate_NTC)

  n_main_cols <- ncol(response_matrix_main)

  chunk <- sparseMatrix(i = rep(1L, n_main_cols),
                        j = 1:n_main_cols,
                        x = rep(1L, n_main_cols),
                        dims = c(1, ncol(response_matrix)))
}


############################################################################################
#############################################################################################
####################################### Papalexi only ###########################################


covariate<-covariate[,c(1,2,3,5)]
covariate[,1]<-log(covariate[,1])
covariate[,2]<-log(covariate[,2])


####################################### repologle only ###########################################
# covariate[,"response_n_umis"] <- log(covariate[,"response_n_umis"]+1)
# covariate[,"response_n_nonzero"] <- log(covariate[,"response_n_nonzero"]+1)
#
# if(use_resample){
#   covariate<-covariate[,c("batch","response_n_umis","response_n_nonzero")]
# } else{
#   covariate<-covariate[,"batch"]
# }


############################################ Nadig only ##########################################

# covariate[,1] <- as.factor(covariate[,1])
# covariate[,2]<-log(covariate[,2]+1)
# covariate[,3]<-log(covariate[,3]+1)
# 
# if (!use_resample){
#   covariate <- covariate[,1]
# }

######################################################################################################
###################################################################################################
####################################################################################################



cat("dimension of the gene expression matrix: ",dim(response_matrix),"\n")
cat("dimension of the grna matrix: ",dim(chunk),"\n")
cat("class of the grna matrix: ",class(chunk),"\n")
cat("number of non-0 element: ",sum(chunk),"\n")


rownames(chunk) <- grna_name
end_step1 <- Sys.time()
cat("Time for STEP 1 (generate guide index matrix):", difftime(end_step1, start_step1, units = "secs"), "\n")

################################ STEP 2 Perturbation Gene qc ##################################

start_step2 <- Sys.time()

if(guide_or_target=="guide"){
  chunk_rownames <- rownames(chunk)
  matched_targets <- grna_target_data_frame_omit[grna_target_data_frame_omit$grna_id %in% chunk_rownames, ]

  unique_targets <- unique(matched_targets$grna_target)
} else {
  matched_targets <- rownames(chunk)


  unique_targets <- unique(matched_targets)
  if (length(unique_targets) != 1) {
    stop("The rows of the chunk_matrix do not refer to exactly one unique target.")
  }

  unique_target <- unique_targets[1]

}


if (length(unique_targets) != 1) {
  unique_target <- NA_character_
  cis_genes <- character(0)
} else {
  unique_target <- unique_targets[1]
  cis_genes <- discovery_pairs_cis[discovery_pairs_cis$grna_target == unique_target, ]$response_id
}

rows_to_zero <- rownames(gene_module_index_matrix) %in% cis_genes



cells_with_grna <- which(chunk[1,] != 0)

all_cells <- seq_len(ncol(chunk))
control_cells_idx <- setdiff(all_cells, cells_with_grna)


trt_cells <-  response_matrix[, cells_with_grna, drop = FALSE]
ctrl_cells <- response_matrix[, control_cells_idx, drop = FALSE]

gene_do_not_pass <- rownames(response_matrix)[
  which(Matrix::rowSums(trt_cells != 0) < trt_cells_cutoff | Matrix::rowSums(ctrl_cells != 0) < control_cells_cutoff)
]

rows_to_qc <- rownames(gene_module_index_matrix) %in% gene_do_not_pass

gene_module_index_matrix[rows_to_zero, ] <- 0
gene_module_index_matrix[rows_to_qc, ] <- 0
end_step2 <- Sys.time()
cat("Time for STEP 2 (perturbation gene QC):", difftime(end_step2, start_step2, units = "secs"), "\n")
print(paste0("cis genes: ",sum(rows_to_zero), ", qc genes: ", sum(rows_to_qc),"\n"))


####################################### STEP 3 Permutation / Resample ###########################
start_step3 <- Sys.time()

grna_assignment_matrix_combined <- permutation_or_resampling(chunk, n_permute, covariate, use_resample)

end_step3 <- Sys.time()
cat("Time for STEP 3 (permutation/resample):", difftime(end_step3, start_step3, units = "secs"), "\n")
########################################## STEP 4 Obtain Test Statistics #################################
start_step4 <- Sys.time()


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
                           chunk)
)


combined_test_tensor <- combine_test_statistics_tensor_low_RAM(normalized_test_tensor, gene_module_index_matrix)


end_step4 <- Sys.time()
cat("Time for STEP 4 (test statistics):", difftime(end_step4, start_step4, units = "secs"), "\n")

#################################################### STEP 5 Obtain P Val ###################################################
start_step5 <- Sys.time()

if (!use_ACAT){
max_statistics <- (-(min_p_val_fast_low_ram(-test_normalizer(combined_test_tensor), gene_module_index_matrix)))
precompute_p <- compute_final_p_val(max_statistics, n_permute,1)
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

rownames(precompute_p)<- grna_name
} else{

  precompute_p <- compute_final_p_val(combined_test_tensor, n_permute,1)


  if (use_approximation) {
    threshold <- approximation_threshold

    indices <- which(!is.na(precompute_p) & precompute_p < threshold, arr.ind = TRUE)

    if (length(indices) > 0) {
      for (k in seq_len(nrow(indices))) {
        i <- indices[k, 2]
        j <- indices[k, 1]

        p_value_1 <- combined_test_tensor[-1, i, j]
        p_value_2 <- combined_test_tensor[1, i, j]
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


  rownames(precompute_p)<- grna_name
  precompute_p <- combine_p_using_ACAT(precompute_p)


}

genewise_p<-compute_final_p_val(normalized_test_tensor,n_permute,3)
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

rownames(genewise_p)<- grna_name
colnames(genewise_p) <- rownames(gene_module_index_matrix)


dir.create(dirname(snakemake@output[["geneset_p_value"]]), recursive = TRUE, showWarnings = FALSE)
saveRDS(precompute_p, snakemake@output[["geneset_p_value"]])
dir.create(dirname(snakemake@output[["genewise_p_value"]]), recursive = TRUE, showWarnings = FALSE)
saveRDS(genewise_p, snakemake@output[["genewise_p_value"]])

print("finish")

end_step5 <- Sys.time()
cat("Time for STEP 5 (compute p-values):", difftime(end_step5, start_step5, units = "secs"), "\n")

