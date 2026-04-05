#!/usr/bin/env bash
set -euo pipefail

SNAKEMAKE=~/miniconda3/envs/snakemake7/bin/snakemake
JOBS=50
CLUSTER_CONFIG="config/cluster_config_2.json"

"${SNAKEMAKE}" \
  -s Snakefile_2 \
  --cluster-config "${CLUSTER_CONFIG}" \
  --cluster "module load R/4.1.0 && sbatch \
    -p {cluster.partition} \
    --account={cluster.account} \
    -t {cluster.time} \
    -N {cluster.nodes} \
    --mem={cluster.mem} \
    -o {cluster.output} \
    -e {cluster.error}" \
  --jobs "${JOBS}" \
  --rerun-incomplete \
  --keep-going
