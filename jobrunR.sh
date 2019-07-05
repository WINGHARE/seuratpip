#!/bin/bash
#SBATCH --job-name="run_seurat"
#SBATCH --output="logs/run_seurat.%j.%N.out"
#SBATCH -N 1
#SBATCH -p LM
#SBATCH --ntasks-per-node 28
#SBATCH -t 24:00:00
#SBATCH --mem 128GB
# echo commands to stdout 
#./nohu*sh
dates=$(date +%s)
Rscript --vanilla seurat_monocle.v1.R WT1_raw_feature_bc_matrix.h5 $dates > logs/$(date +%s).WT1_raw_feature_bc_matrix.out 2>&1 &
Rscript --vanilla seurat_monocle.v1.R WT3_raw_feature_bc_matrix.h5 $dates > logs/$(date +%s).WT3_raw_feature_bc_matrix.out 2>&1 &
Rscript --vanilla seurat_monocle.v1.R WT4_raw_feature_bc_matrix.h5 $dates > logs/$(date +%s).WT4_raw_feature_bc_matrix.out 2>&1 &
Rscript --vanilla seurat_monocle.v1.R KO1_filtered_feature_bc_matrix.h5 $dates > logs/$(date +%s).KO1_filtered_feature_bc_matrix.out 2>&1 &
Rscript --vanilla seurat_monocle.v1.R KO2_filtered_feature_bc_matrix.h5 $dates > logs/$(date +%s).KO2_filtered_feature_bc_matrix.out 2>&1 &
Rscript --vanilla seurat_monocle.v1.R KO3_filtered_feature_bc_matrix.h5 $dates > logs/$(date +%s).KO3_filtered_feature_bc_matrix.out 2>&1 &
Rscript --vanilla seurat_monocle.v1.R KO4_filtered_feature_bc_matrix.h5 $dates > logs/$(date +%s).KO4_filtered_feature_bc_matrix.out 2>&1 &
echo "submit to background, date:"
echo $dates
#Rscript seurat_monocle.v1.R
# move to your appropriate pylon5 directory
# this job assumes:
#  - all input data is stored in this directory
#  - all output should be stored in this directory
#This job runs with 2 nodes, 24 cores per node for a total of 48 cores.
#ibrun in verbose mode will give binding detai
