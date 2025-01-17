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

#Load .rds file saved after last analysis
x <- readRDS(file = "obj.merged.rds")


#Adding percentage of mice mitochondrial genes [NOTE: change if humans to "^MT-"]
x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "mt-")

#Adding percentage of mice ribosomal genes [NOTE: change if humans to "^MT-"]
x[["percent.ribosomal"]] <- PercentageFeatureSet(x, pattern = "^Rp[sl][[:digit:]]|^Rplp[[:digit:]]|^Rpsa")

# Show QC metrics for the first 5 cells
head(x@meta.data, 100)
table(obj.merged$orig.ident)


# Visualize QC metrics as a violin plot and saving the plots in 'Figures' folder [NOTE: add pt.size=0 if you want without dots]
png(filename = "Figures/Vlnplots-before.png", width=1600, height=700)
VlnPlot(x, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0.001, raster=FALSE)
dev.off()

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
png(filename = "Figures/FeatureScatter-before.png", width=1600, height=1400)
plot1 <- FeatureScatter(x, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size=0.001, raster=FALSE)
plot2 <- FeatureScatter(x, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size=0.001, raster=FALSE)
plot1 + plot2
dev.off()
rm(plot1, plot2)


#Filtering the cells 
x <- subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
table(x$orig.ident)

#Plotting after data filtering
png(filename = "Figures/FeatureScatter-after_filtering.png", width=1600, height=700)
plot1_2 <- FeatureScatter(x, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size=0.001, raster=FALSE)
plot2_2 <- FeatureScatter(x, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size=0.001, raster=FALSE)
plot1_2 + plot2_2
dev.off()
rm(plot1_2, plot2_2)

# Visualize QC metrics as a violin plot and saving the plots in 'Figures' folder [NOTE: add pt.size=0 if you want without dots]
png(filename = "Figures/Vlnplots-after.png", width=1600, height=700)
VlnPlot(x, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0.001, raster=FALSE)
dev.off()

#Normalization 
x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000, verbose = TRUE)

#Find variable features
x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top20 <- head(VariableFeatures(x), 20)

# plot variable features with and without labels
png(filename = "Figures/VariableFeaturePlot.png", width=1200, height=700)
plot1_3 <- VariableFeaturePlot(x)
plot2_3 <- LabelPoints(plot = plot1_3, points = top20, repel = TRUE)
plot1_3 + plot2_3
dev.off()
rm(plot1_3, plot2_3)

#Scaling
x <- ScaleData(x)

#Perform linear dimensional reduction
x <- RunPCA(x, features = VariableFeatures(object = x))

# Examine and visualize PCA results a few different ways
print(x[["pca"]], dims = 1:5, nfeatures = 5)

#Loadings to PCA
png(filename = "Figures/PCA_loadings.png", width=900, height=900)
VizDimLoadings(x, dims = 1:5, reduction = "pca")
dev.off()

#PCA plot
png(filename = "Figures/PCA_plot.png", width=600, height=600)
DimPlot(x, reduction = "pca", raster=FALSE)
dev.off()

#Heatmap
png(filename = "Figures/Heatmap.png", width=600, height=600)
DimHeatmap(x, dims = 1, cells = 500)
dev.off()

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
x <- JackStraw(x, num.replicate = 100)
x <- ScoreJackStraw(x, dims = 1:12)

#JackStraw plot
png(filename = "Figures/JackStraw_plot.png", width=600, height=300)
JackStrawPlot(x, dims = 1:12)
dev.off()

#Elbow plot
png(filename = "Figures/Elbow_plot.png", width=600, height=300)
ElbowPlot(x)
dev.off()


# Integrate layers
x <- IntegrateLayers(object = x, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)

# rejoin layers
x[["RNA"]] <- JoinLayers(x[["RNA"]])

#Cluster the cells [NOTE: can change resolution to get more/less clusters]
x <- FindNeighbors(x, dims = 1:7)
x <- FindClusters(x, resolution =0.2)

# Look at cluster IDs of the first 5 cells
head(Idents(x), 5)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
x <- RunUMAP(x, dims = 1:10)

# note that you can set `label = TRUE` or use the Label Clusters function to help label
# individual clusters

#Umap plot
png(filename = "Figures/UMP_plot_02_dim7_Integrated.png", width=600, height=600)
DimPlot(x, reduction = "umap", raster=FALSE)
dev.off()

#Umap plot with samples name
png(filename = "Figures/UMP_plot_with_origID.png", width=600, height=600)
DimPlot(x, reduction = "umap", group.by = 'orig.ident')
dev.off()

Idents(object = x) <- "seurat_clusters"
levels(x)

# find markers for every cluster compared to all remaining cells, report only the positive ones
x.markers <- FindAllMarkers(x, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)

#Renaming identitites
Idents(object = x) <- "orig.ident"
levels(x)

sample.names <- c("Naive", 
                  "WTI_Vehicle", 
                  "WTI_TPOm")

names(sample.names) <- levels(x)
x2 <- RenameIdents(x, sample.names)
levels(x2)

# Stash cell identity classes
x[["samples"]] <- Idents(object = x2)

table(x@meta.data$samples, x@meta.data$orig.ident)

#Removing temporary version of Seurat object
rm(x2)

#Umap plot with samples name
png(filename = "Figures/UMP_plot_with_samples.png", width=600, height=600)
DimPlot(x, reduction = "umap", group.by = 'samples', raster=FALSE)
dev.off()

png(filename = "Figures/UMP_plot_with_clusters.png", width=600, height=400)
DimPlot(x, reduction = "umap", group.by = 'seurat_clusters', label = TRUE, raster=FALSE)
dev.off()

#Saving markers for each cluster
selected_markers <- print(x.markers %>% group_by(cluster) %>% slice_max(n = 15, order_by = avg_log2FC), n = 300)

write.csv(selected_markers, "selected_lung_markers.csv")

#Saving the Seurat object
saveRDS(x, file = "x.rds")

################################################################################################################





