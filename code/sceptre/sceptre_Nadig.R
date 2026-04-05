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
args <- commandArgs(trailingOnly = TRUE)

i <- as.integer(args[1])


setwd("/project/xuanyao/jiaming/Getting_started/data/Nadig/")



grna_assignment_matrix <- readRDS("grna_assignment_matrix.rds")
response_matrix <- readRDS("response_matrix_cNMF70_hallmark.rds")
grna_target_data_frame <-readRDS("grna_target_dataframe.rds")

extra_covariate <-  readRDS("covariate.rds")

sceptre_object <- import_data(
  response_matrix = response_matrix,
  grna_matrix = grna_assignment_matrix,
  grna_target_data_frame = grna_target_data_frame,
  moi = "low",
  extra_covariates = extra_covariate

)

print("finish initiallize")

discovery_pairs_trans_split100 <- readRDS("/project/xuanyao/jiaming/Getting_started/data/Nadig/discovery_pairs_trans_split100.rds")

discovery_pairs_trans<-discovery_pairs_trans_split100[[i]]


peak_ram <- peakRAM({
  sceptre_object <- sceptre_object |>
    set_analysis_parameters(
      discovery_pairs = discovery_pairs_trans,
      side = "both",
      grna_integration_strategy="singleton",
      multiple_testing_alpha=0.05,
      control_group = "nt_cells",
      resampling_mechanism = "permutations",
      formula=formula(~ log(UMI_count) + log(UMI_nonzero)+as.factor(gem_group))
    ) |>
    assign_grnas(method = "thresholding",
                 threshold = 1) |>
    run_qc() |>
    run_discovery_analysis(parallel=T)

})

setwd("/project/xuanyao/jiaming/Getting_started")

print(peak_ram)
# write the results to dis
write_outputs_to_directory(sceptre_object, paste0("output/Nadig/hallmark_cNMF70/set_",i))
