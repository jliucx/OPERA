#!/bin/bash
set -euo pipefail

n_i=3   
n_j=4   

for j in $(seq 1 $n_j); do
  for i in $(seq 1 $n_i); do
    echo "Rscript --vanilla simulation_2.R $i $j" > ender
    file="line_${i}_${j}.sh"
    cat parallel_header.sh ender > "$file"
    chmod u+x "$file"
    sbatch "$file"
    rm -f "$file" ender
  done
done
