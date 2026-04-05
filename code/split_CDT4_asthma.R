setwd("/project/xuanyao/jiaming/Getting_started")
args <- commandArgs(trailingOnly = TRUE)

suppressPackageStartupMessages({
  library(sceptre)
  library(readr)
  library(dplyr)
  library(Matrix)
  library(ondisc)
})

sceptre_obj <- read_ondisc_backed_sceptre_object(
  sceptre_object_fp = "/project/xuanyao/zining/data/perturb_seq/cd4_perturb_seq/Stim8hr_merged_sgRNA/sceptre_object_donor_only.rds",
  response_odm_file_fp = "/project/xuanyao/zining/data/perturb_seq/cd4_perturb_seq/Stim8hr_merged_sgRNA/gene.odm",
  grna_odm_file_fp = "/project/xuanyao/zining/data/perturb_seq/cd4_perturb_seq/Stim8hr_merged_sgRNA/grna.odm")

print("construct sceptre object success!")

# ── 设置checkpoint目录 ────────────────────────────────────────
ckpt_dir <- "/project/xuanyao/jiaming/Getting_started/data/CDT4/checkpoint/flames+200kb"
dir.create(ckpt_dir, recursive = TRUE, showWarnings = FALSE)


# ── 1. 读取 & 匹配 ─────────────────────────────────────────────
Repologle_associated_genes <- readRDS(
  "/project/xuanyao/jiaming/Getting_started/data/GWAS/asthma/GCST010043/flames+200kb/CDT4_associated_perturbation.rds")

grna_names     <- rownames(sceptre_obj@grna_matrix[[1]])
covariate_full <- sceptre_obj@covariate_data_frame
grna_base      <- sub("-\\d+$", "", grna_names)

MCH_guide_idx     <- which(grna_base %in% Repologle_associated_genes)
MCH_guide_targets <- grna_base[MCH_guide_idx]

NTC_guide_idx <- grep("^NTC-", grna_names, ignore.case = TRUE)
NTC_names     <- grna_names[NTC_guide_idx]

print(paste("Found", length(MCH_guide_idx), "MCH guides covering",
            length(unique(MCH_guide_targets)), "targets"))
print(paste("Found", length(NTC_guide_idx), "NTC guides"))


# ── 2. 找cells ────────────────────────────────────────────────
unique_targets <- unique(MCH_guide_targets)

target_cells_list <- vector("list", length(unique_targets))
names(target_cells_list) <- unique_targets

for (tgt in unique_targets) {
  guide_idxs <- MCH_guide_idx[MCH_guide_targets == tgt]
  cells <- c()
  for (gi in guide_idxs) {
    cells <- union(cells, which(sceptre_obj@grna_matrix[[1]][gi, ] != 0))
  }
  target_cells_list[[tgt]] <- cells
  message("target=", tgt, " cells=", length(cells))
}

# ★ 过滤掉cells < 20的target
cells_count <- sapply(target_cells_list, length)
excluded    <- names(cells_count[cells_count < 20])
if (length(excluded) > 0) {
  message("Excluding targets with < 20 cells: ", paste(excluded, collapse = ", "))
}
keep              <- cells_count >= 20
target_cells_list <- target_cells_list[keep]
unique_targets    <- unique_targets[keep]
print(paste("Targets remaining after filtering:", length(unique_targets)))

# NTC guides
NTC_cells_list <- vector("list", length(NTC_guide_idx))
names(NTC_cells_list) <- NTC_names

for (k in seq_along(NTC_guide_idx)) {
  gi <- NTC_guide_idx[k]
  NTC_cells_list[[k]] <- which(sceptre_obj@grna_matrix[[1]][gi, ] != 0)
  message("NTC=", NTC_names[k], " cells=", length(NTC_cells_list[[k]]))
}

all_needed_cells <- sort(unique(c(unlist(target_cells_list), unlist(NTC_cells_list))))
print(paste("Total cells to extract:", length(all_needed_cells)))


# ── 3. 读取expression（断点续跑）─────────────────────────────
all_gene_names <- rownames(sceptre_obj@response_matrix[[1]])
n_genes        <- length(all_gene_names)
exp_mat_ckpt   <- file.path(ckpt_dir, "exp_mat.rds")

if (file.exists(exp_mat_ckpt)) {
  message(">>> Loading exp_mat from checkpoint...")
  exp_mat <- readRDS(exp_mat_ckpt)
} else {
  message(">>> Building exp_mat from scratch...")
  exp_list <- vector("list", n_genes)
  for (k in seq_len(n_genes)) {
    row_full      <- sceptre_obj@response_matrix[[1]][k, ]
    exp_list[[k]] <- row_full[all_needed_cells]
    if (k %% 500 == 0) message("gene ", k, "/", n_genes)
  }
  exp_mat <- do.call(rbind, exp_list)
  rownames(exp_mat) <- all_gene_names
  colnames(exp_mat) <- as.character(all_needed_cells)
  saveRDS(exp_mat, exp_mat_ckpt)
  message(">>> exp_mat saved to checkpoint")
}

print(paste("exp_mat dim:", nrow(exp_mat), "x", ncol(exp_mat)))


# ── 4. Train/test split（断点续跑）───────────────────────────
n_targets <- length(unique_targets)
n_NTC     <- length(NTC_cells_list)

response_train_list <- vector("list", n_targets + n_NTC)
response_test_list  <- vector("list", n_targets)
cov_train_list      <- vector("list", n_targets + n_NTC)
guide_names_out     <- character(n_targets + n_NTC)

for (k in seq_along(unique_targets)) {
  tgt      <- unique_targets[k]
  tgt_ckpt <- file.path(ckpt_dir, paste0("target_", tgt, ".rds"))
  
  if (file.exists(tgt_ckpt)) {
    message(">>> Skipping ", tgt, " (checkpoint found)")
    saved <- readRDS(tgt_ckpt)
    response_train_list[[k]] <- saved$train
    response_test_list[[k]]  <- saved$test
    cov_train_list[[k]]      <- saved$cov
    guide_names_out[k]       <- tgt
    
  } else {
    cells <- target_cells_list[[tgt]]
    
    # 确保cells都在exp_mat中
    cells <- cells[as.character(cells) %in% colnames(exp_mat)]
    
    n       <- length(cells)
    cells_s <- sample(cells, n)
    n_train <- floor(n / 2)
    
    train_cells <- cells_s[seq_len(n_train)]
    test_cells  <- cells_s[(n_train + 1):n]
    
    train_cols <- as.character(train_cells)
    test_cols  <- as.character(test_cells)
    
    r_train <- exp_mat[, train_cols, drop = FALSE]
    r_test  <- exp_mat[, test_cols,  drop = FALSE]
    cov     <- covariate_full[train_cells, , drop = FALSE]
    
    saveRDS(list(train = r_train, test = r_test, cov = cov), tgt_ckpt)
    
    response_train_list[[k]] <- r_train
    response_test_list[[k]]  <- r_test
    cov_train_list[[k]]      <- cov
    guide_names_out[k]       <- tgt
    message("target=", tgt, " train=", n_train, " test=", n - n_train)
  }
}

# NTC（断点续跑）
for (k in seq_along(NTC_cells_list)) {
  pos      <- n_targets + k
  ntc_name <- NTC_names[k]
  ntc_ckpt <- file.path(ckpt_dir, paste0("NTC_", ntc_name, ".rds"))
  
  if (file.exists(ntc_ckpt)) {
    message(">>> Skipping ", ntc_name, " (checkpoint found)")
    saved <- readRDS(ntc_ckpt)
    response_train_list[[pos]] <- saved$train
    cov_train_list[[pos]]      <- saved$cov
    guide_names_out[pos]       <- ntc_name
    
  } else {
    cells     <- NTC_cells_list[[k]]
    cell_cols <- as.character(cells)
    r_train   <- exp_mat[, cell_cols, drop = FALSE]
    cov       <- covariate_full[cells, , drop = FALSE]
    
    saveRDS(list(train = r_train, cov = cov), ntc_ckpt)
    
    response_train_list[[pos]] <- r_train
    cov_train_list[[pos]]      <- cov
    guide_names_out[pos]       <- ntc_name
    message("NTC=", ntc_name, " cells=", length(cells))
  }
}


# ── 5. 合并 & 构建guide assignment ────────────────────────────
response_OPERA  <- do.call(cbind, response_train_list)
covariate_OPERA <- do.call(rbind, cov_train_list)

ncols_each    <- sapply(response_train_list, ncol)
guide_per_col <- rep(guide_names_out, times = ncols_each)

i <- match(guide_per_col, guide_names_out)
j <- seq_along(guide_per_col)

guide_by_col <- sparseMatrix(
  i = i, j = j, x = 1L,
  dims = c(length(guide_names_out), length(guide_per_col)),
  dimnames = list(guide_names_out, colnames(response_OPERA))
)

response_topic <- do.call(cbind, response_test_list)
rownames(response_topic) <- all_gene_names
rownames(response_OPERA)  <- all_gene_names


# ── 6. 保存 ───────────────────────────────────────────────────
saveRDS(response_OPERA,
        "/project/xuanyao/jiaming/Getting_started/data/CDT4/standard_input_format_asthma/flames+200kb/response_matrix.rds")
saveRDS(covariate_OPERA,
        "/project/xuanyao/jiaming/Getting_started/data/CDT4/standard_input_format_asthma/flames+200kb/covariate.rds")
saveRDS(guide_by_col,
        "/project/xuanyao/jiaming/Getting_started/data/CDT4/standard_input_format_asthma/flames+200kb/grna_assignment_matrix.rds")
saveRDS(response_topic,
        "/project/xuanyao/jiaming/Getting_started/data/CDT4/topic_model/flames+200kb/asthma_perturbed_expression_matrix.rds")

print("finished!")