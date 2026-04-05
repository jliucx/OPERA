############################################
## Libraries + C++ sources
############################################




suppressPackageStartupMessages({
  library(splatter)
  library(SingleCellExperiment)
  library(Matrix)
  library(Rcpp)
  library(sceptre)
  library(readr)
  library(dplyr)
  library(parallel)
  library(BH)
  library(peakRAM)
  library(speedglm)
  library(matrixStats)
  library(heavytailcombtest)
  library(ACAT)
  library(limma)
  library(edgeR)
  
  # library(sgt)   # <- 如果没装会报错；需要的话再打开
})

source("/project/xuanyao/jiaming/OPERA/scripts/functions_optimized.R")
sourceCpp("/project/xuanyao/jiaming/OPERA/scripts/fit_skew_normal.cpp")
sourceCpp("/project/xuanyao/jiaming/OPERA/scripts/crt_sample.cpp")


args <- commandArgs(trailingOnly = TRUE)
v <- as.integer(args[1])  # effect_fc index
u <- as.integer(args[2])  # down_gene_frac index

simulate_perturbseq_splatter <- function(
    n_genes = 100,
    n_cells = 10000,
    infected_frac = 0.20,
    down_gene_frac = 0.70,
    effect_fc = 2,
    effect_sdlog = 0.10
) {
  stopifnot(infected_frac > 0 && infected_frac < 1)
  stopifnot(down_gene_frac >= 0 && down_gene_frac <= 1)
  stopifnot(effect_fc >= 1)
  
  params <- newSplatParams(nGenes = n_genes, batchCells = n_cells)
  
  sce <- splatSimulateGroups(
    params,
    group.prob  = c(1 - infected_frac, infected_frac),
    de.prob     = c(0, down_gene_frac),
    de.downProb = c(0, 0.5),              
    de.facLoc   = log(effect_fc),       
    de.facScale = effect_sdlog,
    verbose = FALSE
  )
  
  infected <- (colData(sce)$Group == "Group2")
  Y <- counts(sce)
  rownames(Y) <- paste0("gene", seq_len(n_genes))
  colnames(Y) <- paste0("cell", seq_len(n_cells))
  grna_assignment_matrix <- Matrix(
    as.integer(infected),
    nrow = 1, ncol = n_cells, sparse = TRUE
  )
  rownames(grna_assignment_matrix) <- "gRNA_sim"
  colnames(grna_assignment_matrix) <- colnames(Y)
  
  list(
    sce = sce,
    Y = Y,
    infected = infected,
    grna_assignment_matrix = grna_assignment_matrix
  )
}


############################################
## OPERA wrapper: returns numeric p-value
############################################

OPERA <- function(response_matrix, covariate, grna_assignment_matrix, n_permute) {
  
  grna_assignment_matrix_combined <- permutation_or_resampling(
    grna_assignment_matrix, n_permute, covariate, use_resample = TRUE
  )
  
  test_result <- t_test_wilcoxon(response_matrix, grna_assignment_matrix_combined)
  
  normalized_test_tensor <- test_normalizer_rcpp(
    convert_matrix_to_tensor(test_result, n_permute, grna_assignment_matrix)
  )
  
  gene_module_index_matrix <- matrix(1, nrow = n_genes, ncol = 1)
  colnames(gene_module_index_matrix) <- "GENESET"
  
  combined_test_tensor <- combine_test_statistics_tensor_low_RAM(
    normalized_test_tensor, gene_module_index_matrix
  )
  
  precompute_p <- compute_final_p_val(combined_test_tensor, n_permute, 1)
  
  threshold <- 0.05
  indices <- which(!is.na(precompute_p) & precompute_p < threshold, arr.ind = TRUE)
  
  if (length(indices) > 0) {
    for (k in seq_len(nrow(indices))) {
      i <- indices[k, 2]
      j <- indices[k, 1]
      
      p_value_1 <- combined_test_tensor[-1, i, j]
      p_value_2 <- combined_test_tensor[ 1, i, j]
      p_value_1 <- p_value_1[is.finite(p_value_1)]
      
      if (length(p_value_1) <= 1 || length(unique(p_value_1)) <= 1) {
        result <- 1
      } else {
        result <- fit_and_evaluate_skew_t(p_value_1, p_value_2, 1)
      }
      
      precompute_p[j, i] <- result
    }
  } else {
    message("No indices meet the condition precompute_p < threshold")
  }
  
  rownames(precompute_p) <- rownames(grna_assignment_matrix)
  
  ACAT_p <- combine_p_using_ACAT(precompute_p)
  as.numeric(ACAT_p)[1]
}

limma_acat <- function(Y_counts, infected, covariate_df = NULL) {
  # Y_counts: genes x cells counts
  # infected: logical length n_cells
  # covariate_df: data.frame with columns response_n_nonzero, response_n_umis (optional)
  
  y <- as.matrix(Y_counts)  # n_genes=100 OK
  group <- factor(ifelse(infected, "infected", "ctrl"), levels = c("ctrl","infected"))
  
  if (is.null(covariate_df)) {
    design <- model.matrix(~ group)
  } else {
    # 用你 sceptre 的 covariates
    design <- model.matrix(~ group + log(response_n_nonzero) + log(response_n_umis), data = covariate_df)
  }
  
  dge <- edgeR::DGEList(counts = y)
  dge <- edgeR::calcNormFactors(dge)
  
  v <- limma::voom(dge, design, plot = FALSE)
  fit <- limma::lmFit(v, design)
  fit <- limma::eBayes(fit)
  
  # 取 infected 的系数（groupinfected）
  coef_name <- "groupinfected"
  if (!coef_name %in% colnames(fit$coefficients)) {
    stop("limma: cannot find coef 'groupinfected'. Check design matrix colnames.")
  }
  
  p <- fit$p.value[, coef_name]
  
  # ACAT 需要 (0,1) 内的 p
  p[p >= 1] <- 0.99999
  p[p <= 0] <- 1e-250
  

as.numeric(ACAT::ACAT(p))[1]
  
}

############################################
## Output paths + data
############################################
master_seed <- 0L
seed_i <- function(i, master = master_seed) {
  as.integer((master + 104729L * i) %% (2^31 - 1L)) + 1L
}
n_sim <- 1000
n_genes <- 100
n_cells <- 10000
infected_frac <- 0.01

down_gene_frac_vec <-c(0.05,0.1,0.2,0.5)
effect_fc_vec <- c(1.02,1.03,1.04)

down_gene_frac <- down_gene_frac_vec[u]
effect_fc <-effect_fc_vec[v]
# effect_fc <-1

effect_sdlog <- 0
n_permute_opera <- 10000


out_log  <- paste0("/project/xuanyao/jiaming/paper/output/simulation_with_limma/",100*down_gene_frac,"_density_",effect_fc,"_fc")
dir.create(out_log, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(out_log, paste0("10000_resample_sim_log_",master_seed,".tsv"))

grna_target_data_frame <- readRDS("/project/xuanyao/jiaming/paper/simulation_data/grna_target_dataframe.rds")
discovery_pairs <- readRDS("/project/xuanyao/jiaming/paper/simulation_data/discovery_pairs.rds")


############################################
## Stable seed per iteration (independent of n_sim)
############################################




############################################
## Write header once
############################################

if (!file.exists(out_file)) {
  header <- data.frame(
    sim=integer(), seed=integer(),
    n_genes=integer(), n_cells=integer(),
    infected_frac=double(), de_gene_frac=double(),
    effect_fc=double(), effect_sdlog=double(), n_permute=integer(),
    n_infected=integer(),
    ACAT_p=double(), OPERA_p=double(),Sceptre_bonferoni=double(),Limma_ACAT_p=double(),
    runtime_sec=double(),
    ok=integer(), err=character(),
    stringsAsFactors = FALSE
  )
  write.table(header, out_file, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
}




############################################
## Resume logic
############################################

if (file.exists(out_file)) {
  done <- read.table(out_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  done_i <- done$sim[done$ok == 1]
} else {
  done_i <- integer(0)
}


############################################
## Main loop
############################################

for (i in setdiff(seq_len(n_sim), done_i)) {

  
  t0 <- Sys.time()
  
  row <- tryCatch({
    
    seeds <- seed_i(i)
    set.seed(seeds)
    
    sim <- simulate_perturbseq_splatter(
      n_genes = n_genes, n_cells = n_cells,
      infected_frac = infected_frac, down_gene_frac = down_gene_frac,
      effect_fc = effect_fc, effect_sdlog = effect_sdlog
    )
    
    n_infected <- sum(sim$infected)
    
    print(sum(sim$Y))
    # ---- sceptre ----
    sceptre_object <- import_data(
      response_matrix = sim$Y,
      grna_matrix = sim$grna_assignment_matrix,
      grna_target_data_frame = grna_target_data_frame,
      moi = "high"
    )
    
    sceptre_object <- sceptre_object |>
      set_analysis_parameters(
        discovery_pairs = discovery_pairs,
        side = "both",
        multiple_testing_alpha = 0.05,
        control_group = "complement",
        resampling_mechanism = "crt",
        formula = formula(~ log(response_n_nonzero) + log(response_n_umis)),
        grna_integration_strategy = "singleton"
      ) |>
      assign_grnas(method = "thresholding", threshold = 1) |>
      run_qc(n_nonzero_trt_thresh = 0, n_nonzero_cntrl_thresh = 0) |>
      run_discovery_analysis(parallel = TRUE)
    
    discovery_result <- get_result(
      sceptre_object = sceptre_object,
      analysis = "run_discovery_analysis"
    )
    
    # Avoid exact 1's for ACAT
    p <- discovery_result$p_value
    p[p >= 1] <- 0.99999
    p[p <= 0] <- 1e-250
    
    ACAT_p <- as.numeric(ACAT(p))[1]

    BH_p <- min(1, length(p) * min(p))
    
    
    # ---- OPERA ----
    response_matrix <- t(t(sim$Y) / sceptre_object@covariate_data_frame[, "response_n_umis"])
    covariate <- as.data.frame(sceptre_object@covariate_matrix[, c(2, 3)])
    
    OPERA_p <- OPERA(response_matrix, covariate, sim$grna_assignment_matrix, n_permute_opera)
    Limma_ACAT_p <- limma_acat(sim$Y, sim$infected, sceptre_object@covariate_data_frame[, c("response_n_nonzero","response_n_umis")])
    
    
    runtime_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    
    data.frame(
      sim=i, seed=seeds,
      n_genes=n_genes, n_cells=n_cells,
      infected_frac=infected_frac, de_gene_frac=down_gene_frac,
      effect_fc=effect_fc, effect_sdlog=effect_sdlog, n_permute=n_permute_opera,
      n_infected=n_infected,
      ACAT_p=ACAT_p, OPERA_p=OPERA_p,Sceptre_bonferoni=BH_p,Limma_ACAT_p=Limma_ACAT_p,
      runtime_sec=runtime_sec,
      ok=1L, err="",
      stringsAsFactors = FALSE
    )
    
  }, error = function(e) {
    
    runtime_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    
    data.frame(
      sim=i, seed=seed_i(i),
      n_genes=n_genes, n_cells=n_cells,
      infected_frac=infected_frac, de_gene_frac=down_gene_frac,
      effect_fc=effect_fc, effect_sdlog=effect_sdlog, n_permute=n_permute_opera,
      n_infected=NA_integer_,
      ACAT_p=NA_real_, OPERA_p=NA_real_,Sceptre_bonferoni=NA_real_,Limma_ACAT_p=NA_real_,
      runtime_sec=runtime_sec,
      ok=0L, err=conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })
  
  write.table(row, out_file, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)
  
  message("sim ", i, " done; ok=", row$ok,
          " ACAT=", signif(row$ACAT_p, 3),
          " OPERA=", signif(row$OPERA_p, 3))
  
  print(paste0("finish simulation ", i))
}
