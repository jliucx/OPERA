library(sceptre)
library(sceptredata)
library(readr)
library(dplyr)
library(peakRAM)
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

set.seed(0)
setwd("/project/xuanyao/jiaming/Getting_started/data/frangieh/ifn_gamma/")



grna_assignment_matrix <- readRDS("hallmark_set/union_grna_assignment_matrix.rds")
response_matrix <- readRDS("complete_response_matrix.rds")
grna_target_data_frame <-readRDS("target_target_data_frame.rds")
gene_within_hallmark_geneset <- readRDS("hallmark_set/gene_within_hallmark_geneset.rds")

extra_covariate <-  readRDS("cell_covariate.rds")

sceptre_object <- import_data(
  response_matrix = response_matrix,
  grna_matrix = grna_assignment_matrix,
  grna_target_data_frame = grna_target_data_frame,
  moi = "high",
  extra_covariates = extra_covariate

)

print("finish initiallize")
discovery_pairs_trans <- construct_trans_pairs(
  sceptre_object = sceptre_object
)
discovery_pairs_trans <- discovery_pairs_trans[discovery_pairs_trans$response_id%in%gene_within_hallmark_geneset,]

peak_ram <- peakRAM({
sceptre_object <- sceptre_object |>
  set_analysis_parameters(
    discovery_pairs = discovery_pairs_trans,
    side = "both",
    multiple_testing_alpha=0.05,
    control_group = "complement",
    resampling_mechanism = "crt",
    formula=formula(~ log(n_nonzero) + log(n_umis)+s_score+g2m_score+batch),
    grna_integration_strategy = "singleton"
  ) |>
  assign_grnas(method = "thresholding",
               threshold = 1) |>
  run_qc() |>
    run_calibration_check(parallel = T) |>
    run_discovery_analysis(parallel = T)

})

setwd("/project/xuanyao/jiaming/Getting_started")

print(peak_ram)
# write the results to dis
write_outputs_to_directory(sceptre_object, "/project/xuanyao/jiaming/paper/output/sceptre/frangieh_gamma/union_crt_hallmark_set_genes_match_snakemake/")
