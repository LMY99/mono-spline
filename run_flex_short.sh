#!/usr/bin/env bash

module load conda_R
Rscript "main_flex_short.R" $seed
