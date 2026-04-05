setwd("/project/xuanyao/jiaming/Getting_started")
args <- commandArgs(trailingOnly = TRUE)


suppressPackageStartupMessages({
  library(sceptre)
  library(readr)
  library(dplyr)
  library(Matrix)
  library(ondisc)
})





exp_mat <- readRDS("/project/xuanyao/jiaming/Getting_started/data/Nadig/full_response_matrix.rds")

covariate_full <-readRDS("/project/xuanyao/jiaming/Getting_started/data/Nadig/covariate.rds")

grna_assignment_matrix <-  readRDS("/project/xuanyao/jiaming/Getting_started/data/Nadig/grna_assignment_matrix.rds")

Repologle_associated_genes <- readRDS("/project/xuanyao/jiaming/Getting_started/data/GWAS/asthma/GCST010043/flames+200kb/Nadig_associated_perturbation_flames_200kb.rds")

grna_names <- rownames(grna_assignment_matrix)

MCH_guide_idx_ordered <- match(Repologle_associated_genes, grna_names)


keep <- !is.na(MCH_guide_idx_ordered)
MCH_guide_idx <- MCH_guide_idx_ordered[keep]
Repologle_associated_genes <- Repologle_associated_genes[keep]


NTC_guide_idx <- grep("non-targeting", grna_names)

idx <- MCH_guide_idx
idx_2 <-  NTC_guide_idx

response_train_list <- vector("list", length(idx)+length(idx_2))
response_test_list  <- vector("list", length(idx))
cov_train_list      <- vector("list", length(idx)+length(idx_2))





for (k in seq_along(idx)) {
  
  indx <- idx[k]
  
  
  
  infected_cell <- which(grna_assignment_matrix[indx,]!=0)
  
  cell_to_keep <- infected_cell
  
  n <- length(infected_cell)
  message("guide=", grna_names[indx], " infected cells=", n)
  
  infected_cell <- sample(infected_cell, n)
  n_train <- floor(n / 2)
  
  train_cells <- infected_cell[seq_len(n_train)]
  test_cells  <- infected_cell[(n_train + 1):n]
  
  # 切 response matrix（genes x cells）
  response_train_list[[k]] <- exp_mat[, train_cells, drop = FALSE]
  response_test_list[[k]]  <- exp_mat[, test_cells,  drop = FALSE]
  
  # 切 covariate dataframe（cells x covariates）
  cov_train_list[[k]] <- covariate_full[train_cells, , drop = FALSE]
  
  
}

for (k in seq_along(idx_2)) {
  
  indx <- idx_2[k]
  
  pos <- length(idx)+k
  
  infected_cell <- which(grna_assignment_matrix[indx,]!=0)
  
  cell_to_keep <- infected_cell
  
  
  response_train_list[[pos]] <- exp_mat[, cell_to_keep, drop = FALSE]
  cov_train_list[[pos]] <- covariate_full[cell_to_keep, , drop = FALSE]
  
  
}


guide_names <- grna_names[c(idx,idx_2)]

response_OPERA <- do.call(cbind, response_train_list)
covariate_OPERA <- do.call(rbind,cov_train_list)


ncols_each <- sapply(response_train_list, ncol)
guide_per_col <- rep(guide_names, times = ncols_each)


guide_levels <- guide_names
i <- match(guide_per_col, guide_levels)
j <- seq_along(guide_per_col)

guide_by_col <- sparseMatrix(
  i = i, j = j, x = 1L,
  dims = c(length(guide_levels), length(guide_per_col)),
  dimnames = list(guide_levels, colnames(response_OPERA))
)

response_topic <- do.call(cbind, response_test_list)

rownames(response_topic) <-  rownames(exp_mat)
rownames(response_OPERA)<-  rownames(exp_mat)

colnames(covariate_OPERA)<-c("gem_group","response_n_umis","response_n_nonzero")

saveRDS(response_OPERA,"/project/xuanyao/jiaming/Getting_started/data/Nadig/standard_input_format_asthma/flames+200kb/response_matrix.rds")
saveRDS(covariate_OPERA,"/project/xuanyao/jiaming/Getting_started/data/Nadig/standard_input_format_asthma/flames+200kb/covariate.rds")
saveRDS(guide_by_col,"/project/xuanyao/jiaming/Getting_started/data/Nadig/standard_input_format_asthma/flames+200kb/grna_assignment_matrix.rds")
saveRDS(response_topic,"/project/xuanyao/jiaming/Getting_started/data/Nadig/topic_model/flames+200kb/asthma_perturbed_expression_matrix.rds")




print("finished!")

