#!/bin/bash
# My first script
dates=$(date +%s)
#nohup Rscript --vanilla seurat_monocle.v1.R WT1_raw_feature_bc_matrix.h5 $dates > logs/$(date +%s).WT1_raw_feature_bc_matrix.out 2>&1 &
#nohup Rscript --vanilla seurat_monocle.v1.R WT3_raw_feature_bc_matrix.h5 $dates > logs/$(date +%s).WT3_raw_feature_bc_matrix.out 2>&1 &
#nohup Rscript --vanilla seurat_monocle.v1.R WT4_raw_feature_bc_matrix.h5 $dates > logs/$(date +%s).WT4_raw_feature_bc_matrix.out 2>&1 &
#nohup Rscript --vanilla seurat_monocle.v1.R KO1_filtered_feature_bc_matrix.h5 $dates > logs/$(date +%s).KO1_filtered_feature_bc_matrix.out 2>&1 &
#nohup Rscript --vanilla seurat_monocle.v1.R KO2_filtered_feature_bc_matrix.h5 $dates > logs/$(date +%s).KO2_filtered_feature_bc_matrix.out 2>&1 &
#nohup Rscript --vanilla seurat_monocle.v1.R KO3_filtered_feature_bc_matrix.h5 $dates > logs/$(date +%s).KO3_filtered_feature_bc_matrix.out 2>&1 &
nohup Rscript --vanilla seurat_monocle.v1.R data/GSE69405.csv $dates > logs/$(date +%s).GSE69405.out 2>&1 &
echo "submit to background, date:"
echo $dates
