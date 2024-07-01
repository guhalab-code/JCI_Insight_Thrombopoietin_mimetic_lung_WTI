#Importing libraries
library(Seurat)
library(dplyr)
library(patchwork)
library(SeuratWrappers)
library(tidyverse)
library(ggplot2)
library(cowplot)
#Setting the theme
theme_set(theme_cowplot())

#Setting up my working directory
setwd('C:/Users/bemal/OneDrive/Pulpit/BM-scRNAseq')

#Importing scRNA-seq data from each sample
LM_naive_M.data <- Read10X(data.dir = "1/sample_filtered_feature_bc_matrix")
LM_naive_F.data <- Read10X(data.dir = "2/sample_filtered_feature_bc_matrix")

KO_naive_M.data <- Read10X(data.dir = "3/sample_filtered_feature_bc_matrix")
KO_naive_F.data <- Read10X(data.dir = "4/sample_filtered_feature_bc_matrix")

LM_TPOm_M.data <- Read10X(data.dir = "5/sample_filtered_feature_bc_matrix")
LM_TPOm_F.data <- Read10X(data.dir = "6/sample_filtered_feature_bc_matrix")

KO_TPOm_F1.data <- Read10X(data.dir = "7/sample_filtered_feature_bc_matrix")
KO_TPOm_F2.data <- Read10X(data.dir = "8/sample_filtered_feature_bc_matrix")

LM_naive_IR_M.data <- Read10X(data.dir = "9/sample_filtered_feature_bc_matrix")
LM_naive_IR_F.data <- Read10X(data.dir = "10/sample_filtered_feature_bc_matrix")

KO_naive_IR_M.data <- Read10X(data.dir = "11/sample_filtered_feature_bc_matrix")
KO_naive_IR_F.data <- Read10X(data.dir = "12/sample_filtered_feature_bc_matrix")

LM_TPOm_IR_M.data <- Read10X(data.dir = "13/sample_filtered_feature_bc_matrix")
LM_TPOm_IR_F.data <- Read10X(data.dir = "14/sample_filtered_feature_bc_matrix")

KO_TPOm_IR_M.data <- Read10X(data.dir = "15/sample_filtered_feature_bc_matrix")
KO_TPOm_IR_F.data <- Read10X(data.dir = "16/sample_filtered_feature_bc_matrix")


#Initialize the Seurat object for sample 
LM_naive_M <- CreateSeuratObject(counts = LM_naive_M.data[[1]], project = "LM_naive_M", min.cells = 3, min.features =  200)
LM_naive_M
LM_naive_F <- CreateSeuratObject(counts = LM_naive_F.data[[1]], project = "LM_naive_F", min.cells = 3, min.features =  200)
LM_naive_F

KO_naive_M <- CreateSeuratObject(counts = KO_naive_M.data[[1]], project = "KO_naive_M", min.cells = 3, min.features =  200)
KO_naive_M
KO_naive_F <- CreateSeuratObject(counts = KO_naive_F.data[[1]], project = "KO_naive_F", min.cells = 3, min.features =  200)
KO_naive_F

LM_TPOm_M <- CreateSeuratObject(counts = LM_TPOm_M.data[[1]], project = "LM_TPOm_M", min.cells = 3, min.features =  200)
LM_TPOm_M
LM_TPOm_F <- CreateSeuratObject(counts = LM_TPOm_F.data[[1]], project = "LM_TPOm_F", min.cells = 3, min.features =  200)
LM_TPOm_F

KO_TPOm_F1 <- CreateSeuratObject(counts = KO_TPOm_F1.data[[1]], project = "KO_TPOm_F1", min.cells = 3, min.features =  200)
KO_TPOm_F1
KO_TPOm_F2 <- CreateSeuratObject(counts = KO_TPOm_F2.data[[1]], project = "KO_TPOm_F2", min.cells = 3, min.features =  200)
KO_TPOm_F2

LM_naive_IR_M <- CreateSeuratObject(counts = LM_naive_IR_M.data[[1]], project = "LM_naive_IR_M", min.cells = 3, min.features =  200)
LM_naive_IR_M
LM_naive_IR_F <- CreateSeuratObject(counts = LM_naive_IR_F.data[[1]], project = "LM_naive_IR_F", min.cells = 3, min.features =  200)
LM_naive_IR_F

KO_naive_IR_M <- CreateSeuratObject(counts = KO_naive_IR_M.data[[1]], project = "KO_naive_IR_M", min.cells = 3, min.features =  200)
KO_naive_IR_M
KO_naive_IR_F <- CreateSeuratObject(counts = KO_naive_IR_F.data[[1]], project = "KO_naive_IR_F", min.cells = 3, min.features =  200)
KO_naive_IR_F

LM_TPOm_IR_M <- CreateSeuratObject(counts = LM_TPOm_IR_M.data[[1]], project = "LM_TPOm_IR_M", min.cells = 3, min.features =  200)
LM_TPOm_IR_M
LM_TPOm_IR_F <- CreateSeuratObject(counts = LM_TPOm_IR_F.data[[1]], project = "LM_TPOm_IR_F", min.cells = 3, min.features =  200)
LM_TPOm_IR_F

KO_TPOm_IR_M <- CreateSeuratObject(counts = KO_TPOm_IR_M.data[[1]], project = "KO_TPOm_IR_M", min.cells = 3, min.features =  200)
KO_TPOm_IR_M
KO_TPOm_IR_F <- CreateSeuratObject(counts = KO_TPOm_IR_F.data[[1]], project = "KO_TPOm_IR_F", min.cells = 3, min.features =  200)
KO_TPOm_IR_F



#Merging Seurat objects without deleting normalization

obj.merged <- merge(LM_naive_M, y = c(LM_naive_F, 
                                      KO_naive_M, 
                                      KO_naive_F, 
                                      LM_TPOm_M, 
                                      LM_TPOm_F, 
                                      KO_TPOm_F1, 
                                      KO_TPOm_F2, 
                                      LM_naive_IR_M, 
                                      LM_naive_IR_F,
                                      KO_naive_IR_M,
                                      KO_naive_IR_F,
                                      LM_TPOm_IR_M,
                                      LM_TPOm_IR_F,
                                      KO_TPOm_IR_M,
                                      KO_TPOm_IR_F), 
                    add.cell.ids = c("LM_naive_M", 
                                     "LM_naive_F", 
                                     "KO_naive_M", 
                                     "KO_naive_F", 
                                     "LM_TPOm_M", 
                                     "LM_TPOm_F", 
                                     "KO_TPOm_F1", 
                                     "KO_TPOm_F2", 
                                     "LM_naive_IR_M", 
                                     "LM_naive_IR_F",
                                     "KO_naive_IR_M",
                                     "KO_naive_IR_F",
                                     "LM_TPOm_IR_M",
                                     "LM_TPOm_IR_F",
                                     "KO_TPOm_IR_M",
                                     "KO_TPOm_IR_F"), 
                    project = "BM", merge.data = TRUE)


obj.merged

head(colnames(obj.merged))
tail(colnames(obj.merged))

table(obj.merged$orig.ident)


#Saving the Seurat object
saveRDS(obj.merged, file = "Seurat_analysis/obj.merged.rds")



