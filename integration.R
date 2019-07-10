# set working directory, you may change the directory first.
#setwd("")
# loading required packege



args <- commandArgs(trailingOnly = TRUE)

# test if there is at least one argument: if not, load file 1 as imput
if (length(args)==0) {
    print("No input file name input, use ko1 filename")
  
    args <-
      c(
    "data/KO1_raw_feature_bc_matrix.h5,data/KO2_raw_feature_bc_matrix.h5,data/KO3_raw_feature_bc_matrix.h5,data/KO4_raw_feature_bc_matrix.h5,data/WT1_raw_feature_bc_matrix.h5,data/WT3_raw_feature_bc_matrix.h5,data/WT4_raw_feature_bc_matrix.h5"
    )
    #filename <- "data/WT1_raw_feature_bc_matrix.h5"
    
}
    
print(args)
filename = args[1]
    #sfname = gsub(".h5","",filename)
    #sfname = gsub(".csv","",sfname)
    #sfname = gsub("data/","",sfname)
filenames <- unlist(strsplit(filename, ","))
num.files <- length(filenames)
sfnames <- unlist(strsplit(gsub("[data/]|[.h5]","",filename),","))

    
if(num.files<=1){
  print("Not enough files to intergate")
}
sfname <- paste(gsub("[.h5]|[.csv]|[data/]|[.txt]","",filenames[[1]]),"integrated",sep="_")


if(!require(Seurat)) {
  install.packages("Seurat")
} else{
  library(Seurat)
}
library(dplyr)
################
# input hdf5 ###
################

if(length(grep(".h5",filename)>0)){
    
    my.raw.data <- list()
    for(i in 1:num.files){
      my.raw.data[[i]]<-Read10X_h5(filenames[[i]])
    }
}else if(length(grep("[.csv][.txt]",filename)>0)){
    my.raw.data<- read.csv(filename,sep="\t",header=TRUE, row.names = 1)
}else{
  my.raw.data <- Read10X(data.dir = filename)
  
}
#my.raw.data<-Read10X_h5(filename) # check `?Read10X_h5` for help

my.object <- list()

for(i in 1:num.files){
  my.object[[i]]<- CreateSeuratObject(my.raw.data[[i]])
  my.object[[i]]<- NormalizeData(my.object[[i]],verbose = FALSE)
  my.object[[i]] <- subset(my.object, subset = nFeature_RNA > 200 & nFeature_RNA < 3500)
  my.object[[i]] <- FindVariableFeatures(my.object[[i]],selection.method = "vst"
                                         ,nfeatures = 500, verbose =  FALSE)
  
  
}
names(my.object) <- sfnames
data.all.anchors <- FindIntegrationAnchors(my.object,dims=1:30)
#################
#preprocessing###
#################
# filtering out data with high count, low count, and high MT- counts. 
#my.object[["percent.mt"]] <- PercentageFeatureSet(my.object, pattern = "^MT-")

pdf(file=paste("figs/",sfname,"_vlnplot.pdf",sep = ""))
VlnPlot(my.object, features = c("nFeature_RNA", "nCount_RNA"#,"percent.mt"
                                ), ncol = 2)
dev.off()

# if(sfname!="GSE69405"){
# my.object <- subset(my.object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 #& percent.mt < 5
# )
# }

switch(sfname, 

       GSE108394={
         my.object <- subset(my.object, subset = nFeature_RNA > 2000 & nFeature_RNA < 8000)
       } 
)

#################
my.object<-NormalizeData(my.object,normalization.method = "LogNormalize",scale.factor = 10000)
##########################
# find high variable gene#
##########################
my.object<-FindVariableFeatures(my.object,selection.method = "vst",nfeatures = 5000)
# before PCA, scale data to eliminate extreme value affect.
top10 <- head(VariableFeatures(my.object), 10)
all.gene<-rownames(my.object)
my.object<-ScaleData(my.object,features = all.gene)
# after scaling, perform PCA
my.object<-RunPCA(my.object,rev.pca = F,features = VariableFeatures(object = my.object))
DimPlot(my.object,reduction = "pca")
print(my.object[["pca"]], dims = 1:5, nfeatures = 5)

pdf(file=paste("figs/",sfname,"_elbowplot.pdf",sep = ""))
ElbowPlot(my.object)# check the dims range for the following analysis
dev.off()

pdf(file=paste("figs/",sfname,"_dimheatmap.pdf",sep = ""))
DimHeatmap(my.object, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

###########################################
# CORE part: Run TSNE and UMAP if no result provided#############
###########################################

if (file.exists(paste("output/",sfname,".rds", sep = ""))){
  my.object<-readRDS(file = paste("output/",sfname,".rds", sep = ""))
}else{
  my.object<-RunTSNE(my.object,dims = 1:30,perplexity=10,dim.embed = 3)
  # run umap to get high dimension scatter plot at 2 dimensional coordinate system.
  my.object<-RunUMAP(object = my.object,dims = 1:30)
  
  #clustering by using Seurat KNN.
  # clustering by using KNN, this is seurat cluster algorithm, this part only for cell categorization
  # here has one problems: whether should we define the clustering number?
  my.object<-FindNeighbors(my.object,k.param = 6,dims = 1:30) # k.parm determine k number of clusters
  # find clustering, there will be changing the default cell type, if want to use customized cell type. 
  # use Idents() function.
  my.object<-FindClusters(my.object,resolution = 0.5)
  # Save RDS files
  saveRDS(my.object,file = paste("output/",sfname,".rds", sep = ""))
}
  

pdf(file=paste("figs/",sfname,"_dimplot.pdf",sep = ""))
DimPlot(my.object,reduction = "tsne") # if want to show t-sne, replace "tsne" to "umap"
dev.off()

# Find markers for all culsters and each cluster
my.object.markers.all<-FindAllMarkers(my.object,only.pos = T)
my.object.markers.each <- list()
for(c in unique(my.object$seurat_clusters))
{
  c.marker <- FindMarkers(my.object, ident.1 = c, min.pct = 0.25)
  my.object.markers.each[[paste("my.object.markers.cluster.",c,sep = "")]] <- c.marker

}
#my.object.markers.group.top2<-my.object.markers.all %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

#my.object.markers.group.top2<-my.object.markers.all %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

# Find top markers in a list
find_topk_marker<-function(marker,k=1){
  results<-rownames(head(marker,n=k)[1])
  return(results)
}
find_topk_marker_each<-function(markers.each,k=1){
  
  results <- c()
  for(i in 1:length(markers.each)){
    results <- c(results,(find_topk_marker(markers.each[[i]],k)))
  }
  return(results)
}


pdf(file=paste("figs/",sfname,"_markers_vlnplot.pdf",sep = ""))
VlnPlot(my.object, features = find_topk_marker(my.object.markers.all,k=9), slot = "counts", log = TRUE,pt.size=0)
dev.off()

pdf(file=paste("figs/",sfname,"_markers_vlnplotvs.pdf",sep = ""))
VlnPlot(my.object, features = find_topk_marker_each(my.object.markers.each,1), slot = "counts", log = TRUE,pt.size=0)
dev.off()
