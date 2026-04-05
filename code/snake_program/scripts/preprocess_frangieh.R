inputs <- snakemake@input
outs <- snakemake@output
params <- snakemake@params


library(Matrix)

setwd("/project/xuanyao/jiaming/Getting_started")
print(snakemake@params)

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



# Load input files
covariate <- readRDS(inputs[["covariate"]])
response_matrix_raw <- readRDS(inputs[["response"]])
grna_assignment_matrix <- readRDS(inputs[["grna"]])

cell_filter <- params[["cell_filter"]]
rank <- params[["rank"]]
filter_type= params[["filter_type"]]
filter_umi_cutoff= params[["filter_umi_cutoff"]]
filter_percentage= params[["filter_percentage"]]
grna_to_analyze = params[["grna_to_analyze"]]

grna_assignment_matrix <- grna_assignment_matrix[grna_to_analyze,,drop=F]

if(cell_filter){
  if (filter_type=="umi"){
    cell_to_filter <- which(covariate$response_n_umis < filter_umi_cutoff)
  } else {
    lower1 <- quantile(covariate[,"n_umis"], filter_percentage)
    upper1 <- quantile(covariate[,"n_umis"], 1-filter_percentage)

    lower2 <- quantile(covariate[,"n_nonzero"], filter_percentage)
    upper2 <- quantile(covariate[,"n_nonzero"], 1-filter_percentage)

    cell_to_filter <- which(covariate[,"n_umis"] < lower1 | covariate[,"n_umis"] > upper1 |
                              covariate[,"n_nonzero"] < lower2 | covariate[,"n_nonzero"] > upper2)
  }


  response_matrix_raw <- response_matrix_raw[, -cell_to_filter]
  covariate <- covariate[-cell_to_filter, ]
  grna_assignment_matrix <- grna_assignment_matrix[, -cell_to_filter, drop = FALSE]
}
response_matrix <- t(t(response_matrix_raw) / covariate[, "n_umis"])
# response_matrix <- t(t(response_matrix_raw) / covariate_full[, 1])

if (rank){
  response_matrix <- rank_matrix_v2(response_matrix)
}

covariate <- as.data.frame(covariate)

#####################################################################################################################################################
covariate[,"n_umis"]<-log(covariate[,"n_umis"]+1)
covariate[,"n_nonzero"]<-log(covariate[,"n_nonzero"]+1)

covariate <- covariate[,c(1,2,6,7,9)]
#########################################################################################################################################################
lapply(outs, function(path) dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE))

# Save outputs
saveRDS(response_matrix, outs[["response"]])
saveRDS(grna_assignment_matrix, outs[["grna"]])
saveRDS(covariate, outs[["covariate"]])
write(nrow(grna_assignment_matrix), file = outs[["grna_nrow"]])
cat("Running on node:", Sys.info()[["nodename"]], "\n")
print("finished!")

