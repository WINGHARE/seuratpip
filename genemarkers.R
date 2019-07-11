# set working directory, you may change the directory first.
#setwd("")
# loading required packege

if(!require(Seurat)) {
  install.packages("Seurat")
} else{
  library(Seurat)
}
library(dplyr)
library(optparse)

# args <- commandArgs(trailingOnly = TRUE)

option_list = list(
  make_option(c("-f", "--file"), type="character", default="data/WT1_raw_feature_bc_matrix.h5", 
              help="dataset file names, seperated by commas", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-lb", "--filterl"), type ="integer",default=500, action="store_true",
              help="The lower bound of the filter nFerature RNA [default= %default]"),
  make_option(c("-rb", "--filterr"), type ="integer", default=5000, action="store_true",
              help="The upper bound of the filter nFerature RNA", metavar="integer")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# test if there is at least one argument: if not, load file 1 as imput

filename = opt$file
    #sfname = gsub(".h5","",filename)
    #sfname = gsub(".csv","",sfname)
    #sfname = gsub("data/","",sfname)
sfname = gsub("[.h5]|[.csv]|[data/]|[.txt]","",filename)





################
# input hdf5 ###
################

if(length(grep(".h5",filename)>0)){
    my.raw.data<-Read10X_h5(filename) # check `?Read10X_h5` for help
}else if(length(grep("[.csv][.txt]",filename)>0)){
    my.raw.data<- read.csv(filename,sep="\t",header=TRUE, row.names = 1)

    if(sfname == "GSE69405"){
	    gname.69405<-my.raw.data$gene_name
	    my.raw.data <- my.raw.data[,54:143]# ncol(my.raw.data)]
	    #rownames(my.raw.data) <- gname.69405
    }
}else{
  my.raw.data <- Read10X(data.dir = filename)
  
}
#my.raw.data<-Read10X_h5(filename) # check `?Read10X_h5` for help
my.object<-CreateSeuratObject(my.raw.data) 
#################
#preprocessing###
#################
# filtering out data with high count, low count, and high MT- counts. 
#my.object[["percent.mt"]] <- PercentageFeatureSet(my.object, pattern = "^MT-")

pdf(file=paste("figs/",sfname,"_vlnplot.pdf",sep = ""))
VlnPlot(my.object, features = c("nFeature_RNA", "nCount_RNA"#,"percent.mt"
                                ), ncol = 2)
dev.off()

if(sfname!="GSE69405"){
my.object <- subset(my.object, subset = nFeature_RNA > opt$filterl & nFeature_RNA < opt$filterr)
}

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
  my.object<-RunTSNE(my.object,dims = 1:30,perplexity=10,dim.embed = 2)
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
