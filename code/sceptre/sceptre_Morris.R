library(sceptre)
library(sceptredata)
library(readr)
library(dplyr)
library(peakRAM)
set.seed(0)
setwd("/project/xuanyao/jiaming/Getting_started")

cell_to_remove <- readRDS("/project/xuanyao/jiaming/Getting_started/data/Morris_adjusted/GDO_threshold/cell_to_remove.rds")

grna_assignment_matrix <- readRDS("/project/xuanyao/jiaming/Getting_started/data/Morris_adjusted/program_data/union_grna_assignment_matrix.rds")
response_matrix <- readRDS("/project/xuanyao/jiaming/Getting_started/data/Morris_adjusted/program_data/response_matrix.rds")
grna_target_data_frame <-readRDS("/project/xuanyao/jiaming/Getting_started/data/Morris_adjusted/program_data/target_target_data_frame.rds")

log_covariate <- readRDS("/project/xuanyao/jiaming/Getting_started/data/Morris_adjusted/program_data/log_covariate.rds")
colnames(log_covariate) <- c("a","b","c","d")
sceptre_object <- import_data(
  response_matrix = response_matrix,
  grna_matrix = grna_assignment_matrix,
  grna_target_data_frame = grna_target_data_frame,
  moi = "high",
  extra_covariates = log_covariate
  
)


positive_control_pairs <- construct_positive_control_pairs(sceptre_object)

discovery_pairs_trans <- construct_trans_pairs(
  sceptre_object = sceptre_object,
  pairs_to_exclude = "none"
)
peak_ram <- peakRAM({
  sceptre_object <- sceptre_object |>
    set_analysis_parameters(
      discovery_pairs = discovery_pairs_trans,
      positive_control_pairs = positive_control_pairs,
      side = "both",
      grna_integration_strategy="singleton",
      multiple_testing_alpha=0.05,
      formula=formula(~ a+b+c+d)
    ) |>
    assign_grnas(method = "thresholding",
                 threshold = 1) |>
    # assign_grnas() |>
    run_qc(n_nonzero_trt_thresh=7,n_nonzero_cntrl_thresh=7, response_n_nonzero_range=c(0,1), response_n_umis_range=c(0,1),additional_cells_to_remove=cell_to_remove,p_mito_threshold =1) |>
    run_calibration_check(parallel = F) |>
    run_power_check(parallel = F) |>
    run_discovery_analysis(parallel = F)
})

print(peak_ram)

# write the results to dis
write_outputs_to_directory(sceptre_object, "/project/xuanyao/jiaming/paper/output/sceptre/Morris/union/")
