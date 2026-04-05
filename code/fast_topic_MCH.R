#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Matrix)
  library(fastTopics)
})

library(fastTopics)

mat <- readRDS("/project/xuanyao/jiaming/Getting_started/data/repologle/topic_model/flames_200kb/MCH_perturbed_expression_matrix.rds")

# fastTopics 需要 cell x gene
mat_t <- t(mat)

fit <- fit_topic_model(
  mat_t,
  k              = 50,
  numiter.main   = 200,
  numiter.refine = 200,
  control.main   = list(numiter = 4, nc = 8),    # 8核
  control.refine = list(numiter = 4, extrapolate = TRUE, nc = 8),
  verbose = "progressbar"
)

saveRDS(fit, "/project/xuanyao/jiaming/Getting_started/data/repologle/topic_model/flames_200kb/k50_iter200.rds")

de <- de_analysis(
  fit  = fit,        # 上一步的 topic model fit
  X    = mat_t,      # cell x gene 原始count矩阵
  shrink.method = "ash",
  lfc.stat      = "le",   # least extreme LFC
  verbose       = TRUE
)

saveRDS(de, "/project/xuanyao/jiaming/Getting_started/data/repologle/topic_model/flames_200kb/de_50.rds")

