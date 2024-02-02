#!/usr/bin/env bash

module load conda_R
for seed in {1..100}; do
  if ! test -f "flex_CIs_${seed}.rda"; then
    sbatch --export="seed=${seed}" -J "flexA${seed}" run_flex.sh
  fi
done
for seed in {1..100}; do
  if ! test -f "sflex_CIs_${seed}.rda"; then
    sbatch --export="seed=${seed}" -J "sflexA${seed}" run_flex_short.sh
  fi
done
for seed in {1..100}; do
  if ! test -f "S_CIs_${seed}.rda"; then
    sbatch --export="seed=${seed}" -J "SA${seed}" run_S.sh
  fi
done
for seed in {1..100}; do
  if ! test -f "sS_CIs_${seed}.rda"; then
    sbatch --export="seed=${seed}" -J "sSA${seed}" run_S_short.sh
  fi
done
for seed in {1..100}; do
  if ! test -f "para_CIs_${seed}.rda"; then
    sbatch --export="seed=${seed}" -J "paraA${seed}" run_para.sh
  fi
done
