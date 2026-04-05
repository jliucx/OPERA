library(Matrix)


# 读取输入文件
chunk <- readRDS(snakemake@input[["chunk"]])
grna_target_data_frame_omit <- readRDS(snakemake@input[["grna_target_dataframe"]])
gene_module_index_matrix <- readRDS(snakemake@input[["gene_module_index_matrix"]])
discovery_pairs_cis<- readRDS(snakemake@input[["discovery_pairs_cis"]])
response_matrix<-readRDS(snakemake@input[["response_matrix"]])

trt_cells_cutoff<-snakemake@params[["trt_cells_cutoff"]]
control_cells_cutoff<-snakemake@params[["control_cells_cutoff"]]
guide_or_target<-snakemake@params[["guide_or_target"]]

#############################################

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




#############################################
if (length(unique_targets) != 1) {
  unique_target <- NA_character_
  cis_genes <- character(0)
} else {
  unique_target <- unique_targets[1]
  cis_genes <- discovery_pairs_cis[discovery_pairs_cis$grna_target == unique_target, ]$response_id
}

print(unique_target)
print(head(cis_genes))


rows_to_zero <- rownames(gene_module_index_matrix) %in% cis_genes

if (ncol(chunk) != ncol(response_matrix)) {
  stop(sprintf(
    "Dimension mismatch: ncol(chunk)=%d but ncol(response_matrix)=%d. They must match (cells).",
    ncol(chunk), ncol(response_matrix)
  ))
}

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

print(paste0("cis genes: ",sum(rows_to_zero), ", qc genes: ", sum(rows_to_qc)))


saveRDS(gene_module_index_matrix, snakemake@output[["modified"]])

print("finished!")

