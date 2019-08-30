# set working directory, you may change the directory first.
#setwd("")
# loading required packege
library(Seurat)
library(dplyr)
library(cowplot)
library(optparse)
args <- commandArgs(trailingOnly = TRUE)
print(args)
names.default <-c("data/KO1_raw_feature_bc_matrix.h5,data/KO2_raw_feature_bc_matrix.h5,data/KO3_raw_feature_bc_matrix.h5,data/KO4_raw_feature_bc_matrix.h5,data/WT1_raw_feature_bc_matrix.h5,data/WT3_raw_feature_bc_matrix.h5,data/WT4_raw_feature_bc_matrix.h5")

option_list = list(
  make_option(c("-f", "--file"), type="character", default=names.default, 
              help="dataset file names, seperated by commas", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character"),
  make_option("--filterl", type ="integer", default=as.integer(500),
              help="The lower bound of the filter nFerature RNA [default= %default]", metavar="number"),
  make_option("--filterr", type ="integer", default=as.integer(5000), 
              help="The upper bound of the filter nFerature RNA", metavar="number")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser,args=args);

#args <- commandArgs(trailingOnly = TRUE)

##################  
# Argument tests #
##################

# Test if there is at least one argument: if not, use the default arguments 
# if (length(args)==0) {
#     print("No input file name input, use default filename")
#     args <-c("data/KO1_raw_feature_bc_matrix.h5,data/KO2_raw_feature_bc_matrix.h5,data/KO3_raw_feature_bc_matrix.h5,data/KO4_raw_feature_bc_matrix.h5,data/WT1_raw_feature_bc_matrix.h5,data/WT3_raw_feature_bc_matrix.h5,data/WT4_raw_feature_bc_matrix.h5")
# }
# print(args)
# filename = args[1]

filename <- opt$file
filenames <- unlist(strsplit(filename, ","))
num.files <- length(filenames)
sfnames <- unlist(strsplit(gsub("[data/]|[.h5]","",filename),","))

    
if(num.files<=1){
  print("Not enough files to intergate")
}
sfname <- paste(gsub("[.h5]|[.csv]|[data/]|[.txt]","",filenames[[1]]),"integrated",sep="_")


if (file.exists(paste("output/",sfname,".rds", sep = ""))){
  data.all.integrated<-readRDS(file = paste("output/",sfname,".rds", sep = ""))
}else{ # Start else
# Load data files base on different type of files
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
}
my.object <- list()

# Preprocess the raw files
for(i in 1:num.files){
  my.object[[i]]<- CreateSeuratObject(my.raw.data[[i]])
  my.object[[i]] <- RenameCells(my.object[[i]], add.cell.id = i)
  #my.object[[i]] <- NormalizeData(my.object[[i]],verbose = FALSE)
  #my.object[[i]]<-NormalizeData(my.object[[i]],normalization.method = "LogNormalize",
                                #scale.factor = 10000)
  
  #my.object[[i]] <- subset(my.object[[i]], subset = nFeature_RNA > opt$filterl & nFeature_RNA < opt$filterr)
  #my.object[[i]] <- FindVariableFeatures(my.object[[i]],selection.method = "vst"
                                         #,nfeatures = 500, verbose =  FALSE)
  
}

data.all.integrated<- merge(my.object[[1]], y = c(pbmc4k, pbmc8k)
                            , add.cell.ids = c("3K", "4K", "8K")
                            , project = "Tregko")

saveRDS(data.all.integrated, file = paste("output/",sfname,"combined.rds",sep = ""))