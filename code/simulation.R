suppressPackageStartupMessages({
  
library(splatter)
library(SingleCellExperiment)
library(Matrix)
library(Rcpp)
library(sceptre)
library(readr)
library(dplyr)
library(Matrix)
library(parallel)
library(BH)
library(Rcpp)
library(peakRAM)
library(sgt)
library(speedglm)
library(matrixStats)
library(heavytailcombtest)
library(ACAT)
})
source("/project/xuanyao/jiaming/OPERA/scripts/functions_optimized.R")
sourceCpp("/project/xuanyao/jiaming/OPERA/scripts/fit_skew_normal.cpp")
sourceCpp("/project/xuanyao/jiaming/OPERA/scripts/crt_sample.cpp")



simulate_perturbseq_splatter_restrict_DE <- function(
    n_genes = 10000,
    n_cells = 10000,
    infected_frac = 0.20,
    down_gene_frac = 0.70,      # 这里更准确叫 de_gene_frac（DE比例），但我保持你原名
    effect_fc = 2,
    effect_sdlog = 0.10,
    de_pool_n = 100,
    up_prob = 0.5               # DE gene 上调概率（其余为下调）
) {
  stopifnot(de_pool_n <= n_genes, de_pool_n > 0)
  stopifnot(infected_frac > 0 && infected_frac < 1)
  stopifnot(down_gene_frac >= 0 && down_gene_frac <= 1)
  stopifnot(effect_fc >= 1)
  stopifnot(up_prob >= 0 && up_prob <= 1)
  
  params <- newSplatParams(nGenes = n_genes, batchCells = n_cells)
  
  sce <- splatSimulateGroups(
    params,
    group.prob  = c(1 - infected_frac, infected_frac),
    de.prob     = c(0, 0),
    verbose = FALSE
  )
  
  infected <- (colData(sce)$Group == "Group2")
  colData(sce)$infected <- infected
  colData(sce)$gRNA <- ifelse(infected, "gRNA_sim", "NTC")
  
  # 2) only sample DE genes from 1:de_pool_n
  n_de <- round(de_pool_n * down_gene_frac)
  de_idx <- sort(sample.int(de_pool_n, n_de))
  de_genes <- rep(FALSE, n_genes)
  de_genes[de_idx] <- TRUE
  
  # 3) assign random direction + effect size
  #    fc_vals >= 1 roughly centered around effect_fc
  de_factor <- rep(1, n_genes)     # multiplicative factor on mean
  de_dir <- rep(0L, n_genes)       # +1 up, -1 down, 0 non-DE
  
  if (n_de > 0) {
    fc_vals <- rlnorm(n_de, meanlog = log(effect_fc), sdlog = effect_sdlog)
    is_up <- rbinom(n_de, size = 1, prob = up_prob) == 1
    
    # up: *fc, down: /fc
    signed_fac <- ifelse(is_up, fc_vals, 1 / fc_vals)
    
    de_factor[de_idx] <- signed_fac
    de_dir[de_idx] <- ifelse(is_up, 1L, -1L)
  }
  
  # 4) rebuild mu and resample counts
  Y0 <- counts(sce)
  lib <- Matrix::colSums(Y0); lib <- lib / mean(lib)
  
  ctrl <- !infected
  gene_mu <- rowMeans(Y0[, ctrl, drop = FALSE])
  
  mu <- outer(gene_mu, lib)
  
  if (any(infected) && n_de > 0) {
    mu[de_idx, infected] <- mu[de_idx, infected] * de_factor[de_idx]
  }
  
  nb_size <- 10
  Y <- matrix(
    rnbinom(n_genes * n_cells, mu = as.vector(mu), size = nb_size),
    nrow = n_genes, ncol = n_cells
  )
  Y <- Matrix(Y, sparse = TRUE)
  rownames(Y) <- paste0("gene", seq_len(n_genes))
  colnames(Y) <- paste0("cell", seq_len(n_cells))
  
  grna_assignment_matrix <- Matrix(as.integer(infected), nrow = 1, ncol = n_cells, sparse = TRUE)
  rownames(grna_assignment_matrix) <- "gRNA_sim"
  colnames(grna_assignment_matrix) <- colnames(Y)
  
  list(
    Y = Y,
    infected = infected,
    de_genes = de_genes,
    de_idx = de_idx,
    de_factor = de_factor,   # >1 up, <1 down
    de_dir = de_dir,         # +1 up, -1 down
    grna_assignment_matrix = grna_assignment_matrix
  )
}


OPERA<-function(response_matrix,covariate,grna_assignment_matrix,n_permute){
  
  grna_assignment_matrix_combined <- permutation_or_resampling(grna_assignment_matrix, n_permute, covariate, use_resample=T)
  test_result <- t_test_wilcoxon(response_matrix, grna_assignment_matrix_combined)
  normalized_test_tensor <- test_normalizer_rcpp(
    convert_matrix_to_tensor(test_result,
                             n_permute,
                             grna_assignment_matrix)
  )
  gene_module_index_matrix <- matrix(1,nrow=100,ncol=1)
  colnames(gene_module_index_matrix) <- "GENESET"
  combined_test_tensor <- combine_test_statistics_tensor_low_RAM(normalized_test_tensor, gene_module_index_matrix)
  
  
  precompute_p <- compute_final_p_val(combined_test_tensor, n_permute,1)
  
  
  
  threshold <- 0.05
  
  indices <- which(!is.na(precompute_p) & precompute_p < threshold, arr.ind = TRUE)
  
  if (length(indices) > 0) {
    for (k in seq_len(nrow(indices))) {
      i <- indices[k, 2]
      j <- indices[k, 1]
      
      p_value_1 <- combined_test_tensor[-1, i, j]
      p_value_2 <- combined_test_tensor[1, i, j]
      p_value_1 <- p_value_1[is.finite(p_value_1)]
      
      
      if (length(p_value_1) <= 1 || length(unique(p_value_1)) <= 1) {
        result <- 1
      } else {
        
        
        result <- fit_and_evaluate_skew_t(p_value_1,p_value_2,1)
        
      }
      precompute_p[j, i] <- result
    }
  } else {
    print("No indices meet the condition precompute_p < threshold")
  }
  
  
  
  
  rownames(precompute_p)<- rownames(grna_assignment_matrix)
  
  
  
  ACAT_p <- combine_p_using_ACAT(precompute_p)
  
  ACAT_p_num <- as.numeric(ACAT_p)[1]
  return(ACAT_p_num)
  
}

out_log <- "/project/xuanyao/jiaming/paper/output/simulation/20_density_1_fc"
dir.create(out_log, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(out_log, "sim_log.tsv")

grna_target_data_frame <- readRDS("/project/xuanyao/jiaming/paper/simulation_data/grna_target_dataframe.rds")
discovery_pairs <- readRDS("/project/xuanyao/jiaming/paper/simulation_data/discovery_pairs.rds")

master_seed <- 0L
seed_i <- function(i, master = master_seed) {
  # 生成一个稳定的 32-bit 正整数 seed（避免 0）
  as.integer((master + 104729L * i) %% (2^31 - 1L)) + 1L
}



# 表头（如果文件不存在就写一次）
if (!file.exists(out_file)) {
  header <- data.frame(
    sim=integer(), seed=integer(),
    n_genes=integer(), n_cells=integer(),
    infected_frac=double(), de_pool_n=integer(), de_gene_frac=double(),
    effect_fc=double(), effect_sdlog=double(), up_prob=double(), n_permute=integer(),
    n_infected=integer(),
    n_de=integer(), n_up=integer(), n_down=integer(),
    ACAT_p=double(), OPERA_p=double(),
    runtime_sec=double(),
    ok=integer(), err=character(),
    stringsAsFactors = FALSE
  )
  write.table(header, out_file, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
}

n_sim = 20
n_genes <- 5000
n_cells=10000
infected_frac=0.01
down_gene_frac=0.2
effect_fc=1


if (file.exists(out_file)) {
  done <- read.table(out_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  done_i <- done$sim[done$ok == 1]
} else {
  done_i <- integer(0)
}


for (i in setdiff(seq_len(n_sim), done_i)) {
  
  t0 <- Sys.time()
  
  row <- tryCatch({
    seeds <- seed_i(i)
    set.seed(seeds)
    
    sim <- simulate_perturbseq_splatter_restrict_DE(
      n_genes = n_genes, n_cells = n_cells,
      infected_frac = infected_frac, down_gene_frac = down_gene_frac,
      effect_fc = effect_fc, effect_sdlog = 0,
      de_pool_n = 100, up_prob = 0.5
    )
    
    # 记录模拟真实情况
    n_infected <- sum(sim$infected)
    n_de <- length(sim$de_idx)
    n_up <- sum(sim$de_dir[sim$de_idx] ==  1L)
    n_down <- sum(sim$de_dir[sim$de_idx] == -1L)
    
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
      run_qc(n_nonzero_trt_thresh=0, n_nonzero_cntrl_thresh=0) |>
      run_discovery_analysis(parallel = TRUE)
    
    discovery_result <- get_result(
      sceptre_object = sceptre_object,
      analysis = "run_discovery_analysis"
    )
    discovery_result$p_value[discovery_result$p_value==1] <- 0.99999
    ACAT_p <- ACAT(discovery_result$p_value)
    
    # ---- OPERA ----
    response_matrix <- t(t(sim$Y) / sceptre_object@covariate_data_frame[, "response_n_umis"])[1:100, ]
    covariate <- as.data.frame(sceptre_object@covariate_matrix[, c(2,3)])
    
    OPERA_p <- OPERA(response_matrix, covariate, sim$grna_assignment_matrix, 3000)
    
    runtime_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    
    data.frame(
      sim=i, seed=seeds,
      n_genes=n_genes, n_cells=n_cells,
      infected_frac=infected_frac, de_pool_n=100, de_gene_frac=down_gene_frac,
      effect_fc=effect_fc, effect_sdlog=0, up_prob=0.5, n_permute=3000,
      n_infected=n_infected,
      n_de=n_de, n_up=n_up, n_down=n_down,
      ACAT_p=ACAT_p, OPERA_p=OPERA_p,
      runtime_sec=runtime_sec,
      ok=1L, err="",
      stringsAsFactors = FALSE
    )
    
  }, error = function(e) {
    
    runtime_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    
    data.frame(
      sim=i, seed=seeds,
      n_genes=n_genes, n_cells=n_cells,
      infected_frac=infected_frac, de_pool_n=100, de_gene_frac=down_gene_frac,
      effect_fc=effect_fc, effect_sdlog=0, up_prob=0.5, n_permute=3000,
      n_infected=NA_integer_,
      n_de=NA_integer_, n_up=NA_integer_, n_down=NA_integer_,
      ACAT_p=NA_real_, OPERA_p=NA_real_,
      runtime_sec=runtime_sec,
      ok=0L, err=conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })
  
  # 追加写入（不重复写表头）
  write.table(row, out_file, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)
  
  message("sim ", i, " done; ok=", row$ok, " ACAT=", signif(row$ACAT_p,3), " OPERA=", signif(row$OPERA_p,3))
  
  print(paste0("finish simulation ",i))
}

