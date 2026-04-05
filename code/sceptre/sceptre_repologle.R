library(sceptre)
library(sceptredata)
library(readr)
library(dplyr)
library(peakRAM)
set.seed(0)
setwd("/project/xuanyao/jiaming/Getting_started")

grna_assignment_matrix <- readRDS("/project/xuanyao/jiaming/Getting_started/data/repologle/standard_input_format_MCH/flames_200kb/grna_assignment_matrix.rds")
response_matrix <- as.matrix(readRDS("/project/xuanyao/jiaming/Getting_started/data/repologle/standard_input_format_MCH/flames_200kb/response_matrix.rds"))
grna_target_data_frame <-readRDS("/project/xuanyao/jiaming/Getting_started/data/repologle/standard_input_format_MCH/flames_200kb/grna_target_dataframe.rds")
grna_target_data_frame<-grna_target_data_frame[grna_target_data_frame$grna_id %in% rownames(grna_assignment_matrix),]

covariate <- readRDS("/project/xuanyao/jiaming/Getting_started/data/repologle/standard_input_format_MCH/flames_200kb/covariate.rds")[,c(1,4,5)]
colnames(covariate) <- c("a","b","c")
sceptre_object <- import_data(
  response_matrix = response_matrix,
  grna_matrix = grna_assignment_matrix,
  grna_target_data_frame = grna_target_data_frame,
  moi = "low",
  extra_covariates = covariate
  
)



discovery_pairs_trans <- construct_trans_pairs(
  sceptre_object = sceptre_object,
  pairs_to_exclude = "none"
  
)

de_70_0.01_2_0.5 <- readRDS("/project/xuanyao/jiaming/Getting_started/data/repologle/standard_input_format_MCH/flames_200kb/de_70_0.01_2_0.5.rds")

discovery_pairs_trans<- discovery_pairs_trans[discovery_pairs_trans$response_id %in% unique(unlist(de_70_0.01_2_0.5)),]

peak_ram <- peakRAM({
  sceptre_object <- sceptre_object |>
    set_analysis_parameters(
      discovery_pairs = discovery_pairs_trans,
      side = "both",
      grna_integration_strategy="singleton",
      multiple_testing_alpha=0.05,
      formula=formula(~ as.factor(a)+log(b+1)+log(c+1)),
      resampling_mechanism="crt"
    ) |>
    assign_grnas(method = "thresholding",
                 threshold = 1) |>
    # assign_grnas() |>
    run_qc(n_nonzero_trt_thresh=7,n_nonzero_cntrl_thresh=7,p_mito_threshold =1) |>
    run_calibration_check(parallel = F) |>
    run_power_check(parallel = F) |>
    run_discovery_analysis(parallel = F)
})

print(peak_ram)

# write the results to dis
write_outputs_to_directory(sceptre_object, "/project/xuanyao/jiaming/paper/output/sceptre/repologle_K70")
