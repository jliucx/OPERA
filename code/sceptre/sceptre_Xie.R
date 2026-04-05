library(sceptre)
library(readr)
library(dplyr)
library(peakRAM)
setwd("/project/xuanyao/jiaming/Getting_started")
set.seed(0)

grna_assignment_matrix <- readRDS("/project/xuanyao/jiaming/Getting_started/data/Xie/program_data/grna_assignment_matrix.rds")
response_matrix <- readRDS("/project/xuanyao/jiaming/Getting_started/data/Xie/program_data/response_matrix.rds")
grna_target_data_frame <-readRDS("/project/xuanyao/jiaming/Getting_started/data/Xie/program_data/target_target_data_frame.rds"
)

log_covariate <- readRDS("/project/xuanyao/jiaming/Getting_started/data/Xie/program_data/log_covariate.rds")
colnames(log_covariate) <- c("a","b","c","d","e","f")
sceptre_object <- import_data(
  response_matrix = response_matrix,
  grna_matrix = grna_assignment_matrix,
  grna_target_data_frame = grna_target_data_frame,
  moi = "high",
  extra_covariates = log_covariate
  
)




discovery_pairs_trans <- construct_trans_pairs(
  sceptre_object = sceptre_object
)


peak_ram <- peakRAM({
  sceptre_object <- sceptre_object |>
    set_analysis_parameters(
      discovery_pairs = discovery_pairs_trans,
      side = "both",
      grna_integration_strategy="singleton",
      multiple_testing_alpha=0.05,
      formula=formula(~ a+b+c+d+e+as.factor(f))
    ) |>
    assign_grnas(method = "thresholding",
                 threshold = 1) |>
    # assign_grnas() |>
    run_qc(n_nonzero_trt_thresh=50,n_nonzero_cntrl_thresh=50, response_n_nonzero_range=c(0.01,0.99), response_n_umis_range=c(0.01,0.99),p_mito_threshold =1) |>
    run_calibration_check(parallel = T) |>
    run_discovery_analysis(parallel = T)
})

print(peak_ram)

# write the results to dis
write_outputs_to_directory(sceptre_object, "/project/xuanyao/jiaming/paper/output/sceptre/Xie/union_qc50/")