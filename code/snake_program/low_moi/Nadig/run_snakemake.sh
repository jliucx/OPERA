#!/usr/bin/env bash
set -euo pipefail

SNAKEMAKE=~/miniconda3/envs/snakemake7/bin/snakemake
JOBS=100
CLUSTER_CONFIG="config/cluster_config.json"

"${SNAKEMAKE}" \
  -s Snakefile \
  --cluster-config "${CLUSTER_CONFIG}" \
  --cluster "module load R/4.1.0 && sbatch \
    -p {cluster.partition} \
    --qos {cluster.qos} \
    --account={cluster.account} \
    -t {cluster.time} \
    -N {cluster.nodes} \
    --mem={cluster.mem} \
    -o {cluster.output} \
    -e {cluster.error}" \
  --jobs "${JOBS}" \
  --rerun-incomplete \
  --keep-going
