setwd("/project/xuanyao/jiaming/Getting_started")
args <- commandArgs(trailingOnly = TRUE)


suppressPackageStartupMessages({
  library(sceptre)
  library(readr)
  library(dplyr)
  library(Matrix)
  library(ondisc)
})




sceptre_obj <- read_ondisc_backed_sceptre_object(
  sceptre_object_fp = "/project/xuanyao/jiaming/Getting_started/data/repologle/sceptre_object.rds",
  response_odm_file_fp = "/project/xuanyao/jiaming/Getting_started/data/repologle/gene.odm",
  grna_odm_file_fp = "/project/xuanyao/jiaming/Getting_started/data/repologle/grna.odm")

print("construct sceptre object success!")


exp_list  <- vector("list", length(rownames(sceptre_obj@response_matrix[[1]])))

for (k in seq_along(rownames(sceptre_obj@response_matrix[[1]]))) {
  exp_list[[k]] <- sceptre_obj@response_matrix[[1]][k,]
  print(k)
}

exp_mat <- do.call(rbind, exp_list)
rownames(exp_mat) <- rownames(sceptre_obj@response_matrix[[1]])
covariate_full <- sceptre_obj@covariate_data_frame



Repologle_associated_genes <- readRDS("/project/xuanyao/jiaming/Getting_started/data/GWAS/MCH/flames+200kb_mapping/Repologle_associated_perturbations.rds")

grna_names <- rownames(sceptre_obj@grna_matrix[[1]])

MCH_guide_idx_ordered <- match(Repologle_associated_genes, grna_names)


keep <- !is.na(MCH_guide_idx_ordered)
MCH_guide_idx <- MCH_guide_idx_ordered[keep]
Repologle_associated_genes <- Repologle_associated_genes[keep]


NTC_guide_idx <- grep("non[-_]targeting", grna_names)

idx <- MCH_guide_idx
idx_2 <-  NTC_guide_idx
  
response_train_list <- vector("list", length(idx)+length(idx_2))
response_test_list  <- vector("list", length(idx))
cov_train_list      <- vector("list", length(idx)+length(idx_2))





for (k in seq_along(idx)) {
  
  indx <- idx[k]


  
  infected_cell <- which(sceptre_obj@grna_matrix[[1]][indx,]!=0)
  
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
  
  infected_cell <- which(sceptre_obj@grna_matrix[[1]][indx,]!=0)
  
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
  
rownames(response_topic) <-  rownames(sceptre_obj@response_matrix[[1]])
rownames(response_OPERA)<-  rownames(sceptre_obj@response_matrix[[1]])

saveRDS(response_OPERA,"/project/xuanyao/jiaming/Getting_started/data/repologle/standard_input_format_MCH/flames+200kb/response_matrix.rds")
saveRDS(covariate_OPERA,"/project/xuanyao/jiaming/Getting_started/data/repologle/standard_input_format_MCH/flames+200kb/covariate.rds")
saveRDS(guide_by_col,"/project/xuanyao/jiaming/Getting_started/data/repologle/standard_input_format_MCH/flames+200kb/grna_assignment_matrix.rds")
saveRDS(response_topic,"/project/xuanyao/jiaming/Getting_started/data/repologle/topic_model/flames+200kb/MCH_perturbed_expression_matrix.rds")

  


print("finished!")

