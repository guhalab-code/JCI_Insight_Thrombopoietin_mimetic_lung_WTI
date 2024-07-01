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
setwd('C:/Users/bemal/OneDrive/Pulpit/Jeb-scRNAseq')


#Importing scRNA-seq data from each sample
na.data <- Read10X(data.dir = "na_filtered_feature_bc_matrix")
veh.data <- Read10X(data.dir = "veh_filtered_feature_bc_matrix")
tpo.data <- Read10X(data.dir = "tpo_filtered_feature_bc_matrix")


#Initialize the Seurat object for sample 
naive_lung <- CreateSeuratObject(counts = na.data, project = "na.data", min.cells = 3, min.features =  200)
veh_lung <- CreateSeuratObject(counts = veh.data, project = "veh.data", min.cells = 3, min.features =  200)
tpo_lung <- CreateSeuratObject(counts = tpo.data, project = "tpo.data", min.cells = 3, min.features =  200)


#Merging Seurat objects without deleting normalization

obj.merged <- merge(naive_lung, y = c(veh_lung, 
                                      tpo_lung), 
                    add.cell.ids = c("naive_lung", 
                                     "veh_lung", 
                                     "tpo_lung"),
                    project = "WTI", merge.data = TRUE)
              
obj.merged

head(colnames(obj.merged))
tail(colnames(obj.merged))

table(obj.merged$orig.ident)


#Saving the Seurat object
saveRDS(obj.merged, file = "obj.merged.rds")



