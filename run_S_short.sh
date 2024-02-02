#!/usr/bin/env bash

module load conda_R
Rscript "main_S_short.R" $seed
