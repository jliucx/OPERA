library(Matrix)

grna_matrix <- readRDS(snakemake@input[["grna"]])
description <- snakemake@wildcards[["description"]]

output_dir <- dirname(snakemake@output[[1]])
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 拆分保存
for (i in seq_len(nrow(grna_matrix))) {
  out_path <- snakemake@output[[i]]
  saveRDS(grna_matrix[i, , drop = FALSE], file = out_path)
}

print("finished!")
