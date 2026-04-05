library(sceptre)
library(sceptredata)
library(readr)
library(dplyr)
library(peakRAM)
library(future)
set.seed(0)

# ── 配置：要测试的核数 ──────────────────────────────────────────────
max_cores <- parallel::detectCores(logical = FALSE)
cores_to_test <- unique(c(1, 2, 4, 8, 16, 32, max_cores))
cores_to_test <- sort(cores_to_test[cores_to_test <= max_cores])

cat(sprintf("检测到物理核数: %d\n", max_cores))
cat(sprintf("将测试核数: %s\n\n", paste(cores_to_test, collapse = ", ")))

# ── 读取数据（只读一次） ──────────────────────────────────────────────
data_dir <- "/project/xuanyao/jiaming/Getting_started/data/repologle/standard_input_format_MCH/flames_200kb"

cat("正在读取数据...\n")
grna_assignment_matrix <- readRDS(file.path(data_dir, "grna_assignment_matrix.rds"))
response_matrix        <- as.matrix(readRDS(file.path(data_dir, "response_matrix.rds")))
grna_target_data_frame <- readRDS(file.path(data_dir, "grna_target_dataframe.rds"))
grna_target_data_frame <- grna_target_data_frame[
  grna_target_data_frame$grna_id %in% rownames(grna_assignment_matrix), ]
covariate <- readRDS(file.path(data_dir, "covariate.rds"))[, c(1, 4, 5)]
colnames(covariate) <- c("a", "b", "c")
de_70 <- readRDS(file.path(data_dir, "de_70_0.01_2_0.5.rds"))
cat("数据读取完成。\n\n")

# ── 单次运行函数 ──────────────────────────────────────────────────────
run_sceptre <- function(n_cores) {
  if (n_cores == 1) {
    future::plan(future::sequential)
    use_parallel <- FALSE
  } else {
    future::plan(future::multisession, workers = n_cores)
    use_parallel <- TRUE
  }
  on.exit(future::plan(future::sequential))

  # 构建 sceptre_object（每次独立构建，避免状态残留）
  sceptre_obj <- import_data(
    response_matrix        = response_matrix,
    grna_matrix            = grna_assignment_matrix,
    grna_target_data_frame = grna_target_data_frame,
    moi                    = "low",
    extra_covariates       = covariate
  )

  discovery_pairs_trans <- construct_trans_pairs(
    sceptre_object  = sceptre_obj,
    pairs_to_exclude = "none"
  )
  discovery_pairs_trans <- discovery_pairs_trans[
    discovery_pairs_trans$response_id %in% unique(unlist(de_70)), ]

  # 计时 + 内存
  wall_start <- proc.time()[["elapsed"]]
  ram_result <- peakRAM({
    sceptre_obj <- sceptre_obj |>
      set_analysis_parameters(
        discovery_pairs          = discovery_pairs_trans,
        side                     = "both",
        grna_integration_strategy = "singleton",
        multiple_testing_alpha   = 0.05,
        formula                  = formula(~ as.factor(a) + log(b + 1) + log(c + 1)),
        resampling_mechanism     = "crt"
      ) |>
      assign_grnas(method = "thresholding", threshold = 1) |>
      run_qc(n_nonzero_trt_thresh = 7, n_nonzero_cntrl_thresh = 7, p_mito_threshold = 1) |>
      run_calibration_check(parallel = use_parallel) |>
      run_power_check(parallel = use_parallel) |>
      run_discovery_analysis(parallel = use_parallel)
  })
  wall_time <- proc.time()[["elapsed"]] - wall_start

  list(
    n_cores       = n_cores,
    wall_time_sec = wall_time,
    peak_ram_mb   = ram_result$Peak_RAM_Used_MiB
  )
}

# ── 依次测试各核数 ────────────────────────────────────────────────────
results <- vector("list", length(cores_to_test))

for (i in seq_along(cores_to_test)) {
  nc <- cores_to_test[i]
  cat(sprintf("[%d/%d] 测试 %d 核...\n", i, length(cores_to_test), nc))
  t0 <- Sys.time()
  res <- run_sceptre(nc)
  t1 <- Sys.time()
  results[[i]] <- res
  cat(sprintf("  挂钟时间: %.1f 秒  峰值内存: %.1f MiB\n\n",
              res$wall_time_sec, res$peak_ram_mb))
}

# ── 汇总与指标计算 ────────────────────────────────────────────────────
df <- do.call(rbind, lapply(results, as.data.frame))

baseline_time <- df$wall_time_sec[df$n_cores == 1]  # 单核作为基准
df$speedup          <- baseline_time / df$wall_time_sec          # 加速比
df$parallel_efficiency <- df$speedup / df$n_cores                # 并行效率
df$core_hours       <- df$wall_time_sec * df$n_cores / 3600      # 核时（小时）

cat("══════════════════════════════════════════════════════\n")
cat("Benchmark 结果汇总\n")
cat("══════════════════════════════════════════════════════\n")
print(df, digits = 3, row.names = FALSE)

# ── 保存结果 ──────────────────────────────────────────────────────────
out_dir <- "/project/xuanyao/jiaming/paper/output/sceptre/benchmark"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

result_file <- file.path(out_dir, paste0("benchmark_repologle_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"))
write.csv(df, result_file, row.names = FALSE)
cat(sprintf("\n结果已保存至: %s\n", result_file))

# ── 打印最优核数建议 ──────────────────────────────────────────────────
best_efficiency <- df[which.max(df$parallel_efficiency), ]
best_time       <- df[which.min(df$wall_time_sec), ]
best_corehours  <- df[which.min(df$core_hours), ]

cat("\n── 建议 ──────────────────────────────────────────────\n")
cat(sprintf("最快方案:       %d 核，%.1f 秒\n",
            best_time$n_cores, best_time$wall_time_sec))
cat(sprintf("最省核时方案:   %d 核，%.4f 核时\n",
            best_corehours$n_cores, best_corehours$core_hours))
cat(sprintf("并行效率最高:   %d 核，效率 %.1f%%\n",
            best_efficiency$n_cores, best_efficiency$parallel_efficiency * 100))
