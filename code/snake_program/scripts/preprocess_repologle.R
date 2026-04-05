inputs <- snakemake@input
outs <- snakemake@output
params <- snakemake@params

i <- as.integer(snakemake@wildcards[["chunk"]])

setwd("/project/xuanyao/jiaming/Getting_started")


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
  library(ondisc)
})

suppressMessages(source("/project/xuanyao/jiaming/Getting_started/code/functions_optimized.R"))

lapply(outs, function(path) dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE))


# Load input files


sceptre_obj <- read_ondisc_backed_sceptre_object(
sceptre_object_fp = inputs[["sceptre_rds"]],
response_odm_file_fp = inputs[["gene_odm"]],
grna_odm_file_fp = inputs[["grna_odm"]])

print("construct sceptre object success!")

control_cell = readRDS(inputs[["control_cell"]])
gene_to_analyze = readRDS(inputs[["gene_to_analyze"]])


cell_filter <- params[["cell_filter"]]
rank <- params[["rank"]]
filter_type= params[["filter_type"]]
filter_umi_cutoff= params[["filter_umi_cutoff"]]
filter_percentage= params[["filter_percentage"]]
grna_to_analyze = params[["grna_to_analyze"]]


covariate_full <- sceptre_obj@covariate_data_frame

if(cell_filter){

  if (filter_type=="umi"){
    cell_to_filter <- which(covariate_full$response_n_umis < filter_umi_cutoff)
  } else {
    lower1 <- quantile(covariate_full[,"response_n_umis"], filter_percentage)
    upper1 <- quantile(covariate_full[,"response_n_umis"], 1-filter_percentage)

    lower2 <- quantile(covariate_full[,"response_n_nonzero"], filter_percentage)
    upper2 <- quantile(covariate_full[,"response_n_nonzero"], 1-filter_percentage)

    cell_to_filter <- which(covariate_full[,"response_n_umis"] < lower1 | covariate_full[,"response_n_umis"] > upper1 |
                              covariate_full[,"response_n_nonzero"] < lower2 | covariate_full[,"response_n_nonzero"] > upper2)
  }

} else {
  cell_to_filter <- integer(0)

}




  infected_cell <- which(sceptre_obj@grna_matrix[[1]][i,]!=0)

  # cell_to_keep <- union(infected_cell,control_cell)
  cell_to_keep <- infected_cell

  cell_to_keep <- setdiff(cell_to_keep,cell_to_filter)

  ###############################################

  print(paste0("infected cells: ",length(infected_cell)," number of cells to use: ",length(cell_to_keep)))
  covariate_temp <- covariate_full[cell_to_keep, ]

  expression_matrix <- matrix(
    nrow = length(gene_to_analyze),
    ncol = length(cell_to_keep)
  )
  print(paste0("dim of expression matrix: ", dim(expression_matrix)))
  rownames(expression_matrix) <- gene_to_analyze


  for (g in seq_along(gene_to_analyze)) {
    gene <- gene_to_analyze[g]
   expression <- sceptre_obj@response_matrix[[1]][gene,][cell_to_keep]
    expression_matrix[g, ] <- expression
  }




  response_matrix <- as(t(t(expression_matrix) / covariate_temp[, "response_n_umis"]),"dgCMatrix")

  if (rank){
    response_matrix <- rank_matrix_v2(response_matrix)
  }

  saveRDS(response_matrix, outs[["response_chunk"]])






print("finished!")

