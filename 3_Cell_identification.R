#Setting up my working directory
setwd('C:/Users/bemal/OneDrive/Pulpit/Jeb-scRNAseq')


#Importing libraries
library(Seurat)
library(dplyr)
library(patchwork)
library(SeuratWrappers)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(dplyr)
library(ggrepel)
#Setting the theme
theme_set(theme_cowplot())

#Load .rds file saved after last analysis
x <- readRDS(file = "x.rds")

Idents(object = x2) <- "seurat_clusters"
levels(x2)

#Violin plot
png(filename = "Figures/Feature_plot_Mpl.png", width=400, height=400)
FeaturePlot(x, features = c("Mpl"), pt.size = 0, raster=FALSE)
dev.off()

# endothelial cells marking
png(filename = "Figures/Feature_plot_Kdr.png", width=400, height=400)
FeaturePlot(x, features = c("Kdr"), pt.size = 0, raster=FALSE)
dev.off()

png(filename = "Figures/Feature_plot_Vwf.png", width=400, height=400)
FeaturePlot(x, features = c("Vwf"), pt.size = 0, raster=FALSE)
dev.off()

png(filename = "Figures/Feature_plot_Vegfc.png", width=400, height=400)
FeaturePlot(x, features = c("Vegfc"), pt.size = 0, raster=FALSE)
dev.off()

png(filename = "Figures/Feature_plot_Car4.png", width=400, height=400)
FeaturePlot(x, features = c("Car4"), pt.size = 0, raster=FALSE)
dev.off()

png(filename = "Violin_plot_Axl.png", width=400, height=400)
VlnPlot(x, features = c("Axl"), pt.size = 0, raster=FALSE)
dev.off()

png(filename = "Figures/Feature_plot_Cd34.png", width=400, height=400)
FeaturePlot(x, features = c("Cd34"), pt.size = 0, raster=FALSE)
dev.off()

#Mesenchymal cell markers
png(filename = "Figures/Feature_plot_Pdgfra.png", width=400, height=400)
FeaturePlot(x, features = c("Pdgfra"), pt.size = 0, raster=FALSE)
dev.off()

png(filename = "Figures/Feature_plot_Col13a1.png", width=400, height=400)
FeaturePlot(x, features = c("Col13a1"), pt.size = 0, raster=FALSE)
dev.off()

png(filename = "Figures/Feature_plot_Col14a1.png", width=400, height=400)
FeaturePlot(x, features = c("Col14a1"), pt.size = 0, raster=FALSE)
dev.off()

png(filename = "Figures/Feature_plot_Pdgfrb.png", width=400, height=400)
FeaturePlot(x, features = c("Pdgfrb"), pt.size = 0, raster=FALSE)
dev.off()

png(filename = "Figures/Feature_plot_Lgr6.png", width=400, height=400)
FeaturePlot(x, features = c("Lgr6"), pt.size = 0, raster=FALSE)
dev.off()

png(filename = "Figures/Feature_plot_Lgr5.png", width=400, height=400)
FeaturePlot(x, features = c("Lgr5"), pt.size = 0, raster=FALSE)
dev.off()


# Epithelial / AT markers
png(filename = "Figures/Feature_plot_Aqp5.png", width=400, height=400)
FeaturePlot(x, features = c("Aqp5"), pt.size = 0, raster=FALSE)
dev.off()

png(filename = "Figures/Feature_plot_Ager.png", width=400, height=400)
FeaturePlot(x, features = c("Ager"), pt.size = 0, raster=FALSE)
dev.off()

png(filename = "Figures/Feature_plot_Krt8.png", width=400, height=400)
FeaturePlot(x, features = c("Krt8"), pt.size = 0, raster=FALSE)
dev.off()

png(filename = "Figures/Feature_plot_Sftpc.png", width=400, height=400)
FeaturePlot(x, features = c("Sftpc"), pt.size = 0, raster=FALSE)
dev.off()

png(filename = "Figures/Feature_plot_Epcam.png", width=400, height=400)
FeaturePlot(x, features = c("Epcam"), pt.size = 0, raster=FALSE)
dev.off()

png(filename = "Figures/Feature_plot_Krt19.png", width=400, height=400)
FeaturePlot(x, features = c("Krt19"), pt.size = 0, raster=FALSE)
dev.off()

png(filename = "Figures/Feature_plot_Lrrn4.png", width=400, height=400)
FeaturePlot(x, features = c("Lrrn4"), pt.size = 0, raster=FALSE)
dev.off()

png(filename = "Figures/Feature_plot_Rbp1.png", width=400, height=400)
FeaturePlot(x, features = c("Rbp1"), pt.size = 0, raster=FALSE)
dev.off()


#Changing id to cluster numbers
Idents(object = x) <- "seurat_clusters"

x <- ScaleData(x, features = rownames(x))
#Heatmap of the clusters
png(filename = "Figures/Heatmap_clusters_names_WIP.png", width=2500, height=1500)
x.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10
DoHeatmap(x, features = top10$gene, size = 12, angle = 90) + NoLegend()
dev.off()


#Renaming identitites

cells <- c("Arterial_EC", 
                  "Venous_EC",
                  "Col13+Fibroblast",
                  "CapEC",
                  "Col14+Fibroblast",
                  "AT2",
                  "AT1",
                  "MSC1",
                  "MSC2")

names(cells) <- levels(x)
x3 <- RenameIdents(x, cells)

# Stash cell identity classes
x[["cells"]] <- Idents(object = x3)
rm(x3)

Idents(object = x) <- "cells"
levels(x)

table(x@meta.data$samples, x@meta.data$cells)

#Calculating percentages of each cell type
cell_type_percentages <- x@meta.data %>%
  group_by(samples, cells) %>%
  summarise(cell_type_count = n()) %>%
  group_by(samples) %>%
  mutate(total_cells = sum(cell_type_count),
         cell_type_percentage = (cell_type_count / total_cells) * 100) %>%
  ungroup()

#Column plot
png(filename = "Figures/Cell_percentages_split.png", width=600, height=1000)
ggplot(cell_type_percentages, aes(fill = cells, y = cell_type_percentage, x = samples)) +
  geom_col(position = "stack") +
  labs(x = "Sample", y = "Cell Type Percentage", title = "Cell Type Percentages by Sample") +
  geom_text(aes(label = paste0(round(cell_type_percentage, 1), "%")),
            position = position_stack(vjust = 0.5), color = "white", size = 3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(cells ~ ., scales = "free_y")
dev.off()

#Stacked column plot
png(filename = "Figures/Cell_percentages_split_stacked.png", width=800, height=1000)
ggplot(cell_type_percentages, aes(fill = cells, y = cell_type_percentage, x = samples)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(x = "Sample", y = "Cell Type Percentage", title = "Cell Type Percentages by Sample") +
  geom_text(aes(label = paste0(round(cell_type_percentage, 1), "%")),
            position = position_fill(vjust = 0.5), color = "white", size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

Idents(object = x) <- "cells"
#Umap plot with cells names
png(filename = "Figures/UMP_plot_with_cells.png", width=1200, height=800)
DimPlot(x, reduction = "umap", label = FALSE, pt.size = 1, label.size = 4.5, repel = TRUE, raster=FALSE)
dev.off()


Idents(object = x) <- "samples"


#Umap plot with samples names
png(filename = "Figures/UMP_plot_with_samples_WIP.png", width=1600, height=600)
DimPlot(x, reduction = "umap", group.by = 'samples', label = FALSE, label.size = 4.5, repel = TRUE, raster=FALSE)
dev.off()


#Saving the Seurat object
saveRDS(x, file = "x_names.rds")

####################################################################

#Dot plot [Note: can add: cols = c("blue", "red"), split.by = "orig.ident"] # can add: plot1 + coord_flip()
EC.markers.to.plot <- c("Kdr", 
                     "Vwf", 
                     "Vegfc",
                     "Car4",
                     "Fibin")

png(filename = "Figures/Dot_plot_EC_markers.png", width=800, height=600)
dot <- DotPlot(x, features = rev(EC.markers.to.plot), dot.scale = 8) + RotatedAxis()
dot + coord_flip()
dev.off()
rm(dot)

x <- factor(x$cells,levels=c("Aterial_EC","Venous_EC","CapEC","Col13+Fibroblast","Col14+Fibroblast","MSC1","MSC2","AT1","AT2"))

all.markers.to.plot <- c("Kdr", "Vwf", "Vegfc", "Car4", "Fibin",
                     "Cd34", "Pdgfra", "Col13a1", "Col14a1", "Pdgfrb",
                     "Lgr5", "Lgr6", "Aqp5", "Ager", "Sftpc",
                     "Epcam")
png(filename = "Figures/Dot_plot_combo_markers_WIP.png", width=1800, height=1200)
dot <- DotPlot(x, features = rev(all.markers.to.plot), dot.scale = 25) + RotatedAxis()
dot + coord_flip()
dev.off()
rm(dot)


#Cell cycle genes import
m.s.genes <- c("Mcm5", "Pcna", "Tyms", "Fen1", "Mcm2", "Mcm4", "Rrm1", "Ung", "Gins2", "Mcm6", "Cdca7", "Dtl", "Prim1", "Uhrf1", "Mlf1ip", "Hells", "Rfc2", "Rpa2", "Nasp", "Ras51ap1", "Gmnn", "Wdr76", "Slbp",  "Ccne2", "Ubr7", "Pold3", "Msh2", "Atad2", "Rad51", "Rrm2", "Cdc45", "Cdc6", "Exo1", "Tipin", "Dscc1", "Blm", "CaspP8ap2", "Usp1", "Clspn", "Pola1", "Chaf1b", "Brip1", "E2f8")
m.g2m.genes <- c("Hmgb2", "Cdk1"  ,"Nusap1" ,"Ube2c", "Birc5" ,"Tpx2", "Top2a", "Ndc80" ,"Cks2" , "Nuf2"  , "Cks1b" , "Mki67" , "Tmpo", "Cenpf" , "Tacc3" , "Fam64a",  "Smc4", "Ccnb2", "Ckap2l","Ckap2", "Aurkb" , "Bub1"  , "Kif11" ,"Anp32e" , "Tubb4b", "Gtse1" ,"Kif20b" ,"Hjurp"  , "Cdca3" ,"Hn1", "Cdc20" ,"Ttk"  , "Cdc25c", "Kif2c" , "Rangap1" , "Ncapd2", "Dlgap5", "Cdca2" , "Cdca8" ,"Ect2" , "Kif23" , "Hmmr"  , "Aurka" , "Psrc1" , "Anln", "Lbr",  "Ckap5" , "Cenpe" ,"Ctcf" , "Nek2"  , "G2e3"  , "Gas2l3", "Cbx5"  , "Cenpa")

x <- CellCycleScoring(x, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(x[[]])

#Umap plot with samples name
png(filename = "Figures/UMP_plot_with_cell_cycle_split_condition.png", width=1600, height=800)
DimPlot(x, reduction = "umap", group.by = "Phase", split.by = 'samples', raster=FALSE)
dev.off()



####################################################################


#Saving the Seurat object
saveRDS(x, file = "x_names_cc.rds")




