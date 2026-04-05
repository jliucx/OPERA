#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Matrix)
  library(fastTopics)
})

library(fastTopics)
set.seed(0)

# mat <- readRDS("/project/xuanyao/jiaming/Getting_started/data/Nadig/topic_model/flames+200kb/asthma_perturbed_expression_matrix.rds")
# 
# # fastTopics 需要 cell x gene
# mat_t <- t(mat)
# 
# 
# 
# fit <- fit_topic_model(
#   mat_t,
#   k              = 15,
#   numiter.main   = 1000,
#   numiter.refine = 1000,
#   control.main   = list(numiter = 4, nc = 16),    # 8核
#   control.refine = list(numiter = 4, extrapolate = TRUE, nc = 16),
#   verbose = "progressbar"
# )
# 
# saveRDS(fit, "/project/xuanyao/jiaming/Getting_started/data/Nadig/topic_model/flames+200kb/k15_iter1000.rds")
# 
# 
# de <- de_analysis(
#   fit           = fit,
#   X             = mat_t,
#   shrink.method = "ash",
#   lfc.stat      = "le",
#   control       = list(nc = 16),   # ← 并行核数放在 control 里
#   verbose       = TRUE
# 
# )
# 
# saveRDS(de, "/project/xuanyao/jiaming/Getting_started/data/Nadig/topic_model/flames+200kb/de_15.rds")
# 
# 
# 
# 
# 
# fit <- fit_topic_model(
#   mat_t,
#   k              = 10,
#   numiter.main   = 1000,
#   numiter.refine = 1000,
#   control.main   = list(numiter = 4, nc = 16),    # 8核
#   control.refine = list(numiter = 4, extrapolate = TRUE, nc = 16),
#   verbose = "progressbar"
# )
# 
# saveRDS(fit, "/project/xuanyao/jiaming/Getting_started/data/Nadig/topic_model/flames+200kb/k10_iter1000.rds")
# 
# de <- de_analysis(
#   fit           = fit,
#   X             = mat_t,
#   shrink.method = "ash",
#   lfc.stat      = "le",
#   control       = list(nc = 16),   # ← 并行核数放在 control 里
#   verbose       = TRUE
#   
# )
# 
# saveRDS(de, "/project/xuanyao/jiaming/Getting_started/data/Nadig/topic_model/flames+200kb/de_10.rds")
# 
# 
# 
# 
# fit <- fit_topic_model(
#   mat_t,
#   k              = 20,
#   numiter.main   = 1000,
#   numiter.refine = 1000,
#   control.main   = list(numiter = 4, nc = 16),    # 8核
#   control.refine = list(numiter = 4, extrapolate = TRUE, nc = 16),
#   verbose = "progressbar"
# )
# 
# saveRDS(fit, "/project/xuanyao/jiaming/Getting_started/data/Nadig/topic_model/flames+200kb/k20_iter1000.rds")
# 
# de <- de_analysis(
#   fit           = fit,
#   X             = mat_t,
#   shrink.method = "ash",
#   lfc.stat      = "le",
#   control       = list(nc = 16),   # ← 并行核数放在 control 里
#   verbose       = TRUE
#   
# )
# 
# saveRDS(de, "/project/xuanyao/jiaming/Getting_started/data/Nadig/topic_model/flames+200kb/de_20.rds")
# 
# 
# 


########################################  Repologle #####################################################


# mat <- readRDS("/project/xuanyao/jiaming/Getting_started/data/repologle/topic_model/flames+200kb/MCH_perturbed_expression_matrix.rds")
# mat_t <- t(mat)
# 
# fit <- fit_topic_model(
#   mat_t,
#   k              = 50,
#   numiter.main   = 200,
#   numiter.refine = 200,
#   control.main   = list(numiter = 4, nc = 32),    # 8核
#   control.refine = list(numiter = 4, extrapolate = TRUE, nc = 32),
#   verbose = "progressbar"
# )
# 
# saveRDS(fit, "/project/xuanyao/jiaming/Getting_started/data/repologle/topic_model/flames+200kb/k50_iter200.rds")
# 
# de <- de_analysis(
#   fit           = fit,
#   X             = mat_t,
#   shrink.method = "ash",
#   lfc.stat      = "le",
#   control       = list(nc = 32),   # ← 并行核数放在 control 里
#   verbose       = TRUE
# 
# )
# 
# saveRDS(de, "/project/xuanyao/jiaming/Getting_started/data/repologle/topic_model/flames+200kb/de_50.rds")

###############################################  CDT4 #############################################################




mat <-  readRDS("/project/xuanyao/jiaming/Getting_started/data/CDT4/topic_model/flames+200kb/asthma_perturbed_expression_matrix.rds")
mat_t <- t(mat)



# 
# fit <- fit_topic_model(
#   mat_t,
#   k              = 30,
#   numiter.main   = 200,
#   numiter.refine = 200,
#   control.main   = list(numiter = 4, nc = 32),    # 8核
#   control.refine = list(numiter = 4, extrapolate = TRUE, nc = 32),
#   verbose = "progressbar"
# )
# 
# saveRDS(fit, "/project/xuanyao/jiaming/Getting_started/data/CDT4/topic_model/flames+200kb/k30_iter200.rds")


fit <- readRDS("/project/xuanyao/jiaming/Getting_started/data/CDT4/topic_model/flames+200kb/k30_iter200.rds")

de <- de_analysis(
  fit           = fit,
  X             = mat_t,
  shrink.method = "ash",
  lfc.stat      = "le",
  control       = list(nc = 32),   # ← 并行核数放在 control 里
  verbose       = TRUE
  
)

saveRDS(de, "/project/xuanyao/jiaming/Getting_started/data/CDT4/topic_model/flames+200kb/de_30.rds")



fit <- fit_topic_model(
  mat_t,
  k              = 50,
  numiter.main   = 200,
  numiter.refine = 200,
  control.main   = list(numiter = 4, nc = 32),    # 8核
  control.refine = list(numiter = 4, extrapolate = TRUE, nc = 32),
  verbose = "progressbar"
)

saveRDS(fit, "/project/xuanyao/jiaming/Getting_started/data/CDT4/topic_model/flames+200kb/k50_iter200.rds")

de <- de_analysis(
  fit           = fit,
  X             = mat_t,
  shrink.method = "ash",
  lfc.stat      = "le",
  control       = list(nc = 32),   # ← 并行核数放在 control 里
  verbose       = TRUE
  
)

saveRDS(de, "/project/xuanyao/jiaming/Getting_started/data/CDT4/topic_model/flames+200kb/de_50.rds")
