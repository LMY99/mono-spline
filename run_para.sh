#!/usr/bin/env bash

module load conda_R
Rscript "main_para.R" $seed
