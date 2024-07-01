#####CAP######

x <- readRDS(file = "x_names.rds")
Idents(object = x) <- "cells"
levels(x)

endo_clusters <- c("cap_EC")
CAP <- subset(x, idents = endo_clusters)
arterialEC <- c("endo1")
ART <- subset(x, idents = arterialEC)
table(CAP$cells)
levels(CAP)

Idents(object = CAP) <- "samples"
levels(CAP)
head(CAP@meta.data)

x2 <- FindNeighbors(CAP, dims = 1:4)
x2 <- FindClusters(CAP, resolution =0.2)

EC_clusters <- c("2")
endo2 <- subset(x2, idents = EC_clusters)

x2 <- RunUMAP(x2, dims = 1:10)

x.markers <- FindAllMarkers(x2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)


# Find differentially expressed features between Vehilcel and naive # might add: min.diff.pct = 0.25, max.cells.per.ident = 200
latent_vars <- c("nCount_RNA", "percent.mt", "S.Score")
CAP_WTI_WTI_Vehicle_vs_Naive <- FindMarkers(CAP, ident.1 = "WTI_Vehicle", ident.2 = "Naive",
                                           min.pct = 0.25, 
                                           logfc.threshold = log(0.5),
                                           test.use = "MAST",
                                           latent.vars = latent_vars)

#Saving results
write.xlsx(CAP_WTI_WTI_Vehicle_vs_Naive, "Cap_WTI_Vehicle_vs_Naive.xlsx", rowNames = TRUE)

##volcano of differential between TPOm and Vehicle 

alpha <- 0.05
fc <- 0.5
df <- rownames_to_column(CAP_WTI_Vehicle_vs_Naive, var = "gene")
significant_genes <- df[df$avg_log2FC > fc & df$p_val < alpha |
                          df$avg_log2FC < -fc & df$p_val < alpha, ]
volcano_plot <- ggplot(df, aes(x = avg_log2FC, y = -log10(p_val))) +
  geom_point(aes(color = abs(avg_log2FC) > fc & p_val < alpha), size = 2) +
  scale_color_manual(name = "Legend",
                     values = c("black", "red"),
                     labels = c("Not Significant", "Significant")) +
  geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-fc, fc), linetype = "dashed", color = "blue") +
  labs(x = "avg. log2 Fold Change", y = "-log10(p val)") +
  geom_text_repel(data = significant_genes,
                  aes(label = gene), color = "red", size = 6, max.overlaps = 15) +
  theme_minimal()+
  theme(panel.background = element_rect(fill = "white"))

ggsave("Figures/Volcano_Plot_CAP_WTI_Vehicle_vs_Naive_WIP.png", plot = volcano_plot, width = 2500, height = 1500, units = "px")
rm(volcano_plot)


# Set the log2 fold change threshold
alpha <- 0.05
fc <- 0.5

# Filter for significant genes with a meaningful log2 fold change and no NAs
filtered_genes <- CAP_WTI_Vehicle_vs_Naive[
  CAP_WTI_Vehicle_vs_Naive$p_val < alpha & abs(CAP_WTI_Vehicle_vs_Naive$avg_log2FC) > fc, 
]

# Separate up-regulated and down-regulated genes
log2FC_t <- 0.5
up <- filtered_genes[filtered_genes$avg_log2FC > log2FC_t, ]
down <- filtered_genes[filtered_genes$avg_log2FC < -log2FC_t, ]

# Display separately up-regulated and down-regulated genes
print(head(up))

print(head(down))

#Creating a list of gene names
up_l <- rownames(up)

BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)

#Select genes 
gene_set <- up_l
pathway_db <- "org.Mm.eg.db"  # Use appropriate organism database

# Perform gene overrepresentation analysis
enrich_result_WTI_Vehicle_vs_Naive_up <- enrichGO(gene = gene_set,
                                 OrgDb = pathway_db,
                                 keyType = "SYMBOL",
                                 ont = "BP",  # Biological Process ontology
                                 pAdjustMethod = "BH")

# Extract the top 10 pathways
top_pathways_WTI_Vehicle_vs_Naive_up <- head(enrich_result_WTI_Vehicle_vs_Naive_up@result, 7)

head(top_pathways_WTI_Vehicle_vs_Naive_up)

# Set larger bottom margin to make more space for X-axis labels
# Adjust the bottom margin (default is c(5,4,4,2) + 0.1)

# Save the plot as PNG with increased size
png("Figures/EnrichGO-CAP_WTI_Vehicle_vs_Naive_up_thesis.png", width = 800, height = 600)  # Adjust width and height as needed
par(mar = c(4, 22, 4, 1) + 0.1, mgp = c(1, 3, 0))
barplot(rev(-log10(top_pathways_WTI_Vehicle_vs_Naive_up$pvalue)),
        names.arg = rev(top_pathways_WTI_Vehicle_vs_Naive_up$Description),
        horiz = T,
        las = 1, cex.names = 1,
        main = "Vehicle vs Naive Upregulated Pathways",
      
        ylab = "-log10(pvalue)",
        col = "lightpink",
        cex.axis = 0.8)

dev.off()  # Close the PNG device

write.xlsx(enrich_result_CAP_up@result, "EnrighGO-CAP_WTI_Vehicle_vs_Naive_up.xlsx", rowNames = FALSE)

#Creating a list of gene names
down_l <- rownames(down)

#Select genes 
gene_set <- down_l
pathway_db <- "org.Mm.eg.db"  # Use appropriate organism database

# Perform gene overrepresentation analysis
enrich_result_WTI_Vehicle_vs_Naive_down <- enrichGO(gene = gene_set,
                                   OrgDb = pathway_db,
                                   keyType = "SYMBOL",
                                   ont = "BP",  # Biological Process ontology
                                   pAdjustMethod = "BH")

# Extract the top 10 pathways
top_pathways_WTI_Vehicle_vs_Naive_down <- head(enrich_result_WTI_Vehicle_vs_Naive_down@result, 7)


# Set larger bottom margin to make more space for X-axis labels
# Adjust the bottom margin (default is c(5,4,4,2) + 0.1)

# Save the plot as PNG with increased size
png("Figures/EnrighGO-CAP_WTI_Vehicle_vs_Naive_down_thesis.png", width = 800, height = 600)  # Adjust width and height as needed
par(mar = c(4, 22, 4, 2) + 0.1, mgp = c(1, 3, 0))
barplot(rev(-log10(top_pathways_WTI_Vehicle_vs_Naive_down$pvalue)),
        names.arg = rev(top_pathways_WTI_Vehicle_vs_Naive_down$Description),
        horiz = T,
        las = 1, cex.names = 1,
        main = "Vehicle vs Naive Downregulated Pathways",
   
        ylab = "-log10(P-value)",
        col = "lightblue2",
        cex.axis = 0.8)

dev.off()  # Close the PNG device

write.xlsx(enrich_result_WTI_Vehicle_vs_Naive_down@result, "EnrighGO-CAP_WTI_Vehicle_vs_Naive_down.xlsx", rowNames = FALSE)




##############################################

# Find differentially expressed features between TPOm and vehicle # might add: min.diff.pct = 0.25, max.cells.per.ident = 200
latent_vars <- c("nCount_RNA", "percent.mt", "S.Score")
CAP_WTI_TPOm_vs_WTI_Vehicle <- FindMarkers(CAP, ident.1 = "WTI_TPOm", ident.2 = "WTI_Vehicle",
                                            min.pct = 0.25, 
                                            logfc.threshold = log(0.5),
                                            test.use = "MAST",
                                            latent.vars = latent_vars)
##volcano of differential between TPOm and Vehicle 

alpha <- 0.05
fc <- 0.5
df <- rownames_to_column(CAP_WTI_TPOm_vs_WTI_Vehicle, var = "gene")
significant_genes <- df[df$avg_log2FC > fc & df$p_val < alpha |
                          df$avg_log2FC < -fc & df$p_val < alpha, ]
volcano_plot <- ggplot(df, aes(x = avg_log2FC, y = -log10(p_val))) +
  geom_point(aes(color = abs(avg_log2FC) > fc & p_val < alpha), size = 2) +
  scale_color_manual(name = "Legend",
                     values = c("black", "red"),
                     labels = c("Not Significant", "Significant")) +
  geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-fc, fc), linetype = "dashed", color = "blue") +
  labs(x = "avg. log2 Fold Change", y = "-log10(p val)") +
  geom_text_repel(data = significant_genes,
                  aes(label = gene), color = "red", size = 6, max.overlaps = 15) +
  theme_minimal()+
  theme(panel.background = element_rect(fill = "white"))

ggsave("Figures/Volcano_Plot_CAP_WTI_TPOm_vs_WTI_Vehicle_WIP.png", plot = volcano_plot, width = 2500, height = 1500, units = "px")
rm(volcano_plot)
# Set the log2 fold change threshold
alpha <- 0.05
fc <- 0.25

# Filter for significant genes with a meaningful log2 fold change and no NAs
filtered_genes <- CAP_WTI_TPOm_vs_WTI_Vehicle[
  CAP_WTI_TPOm_vs_WTI_Vehicle$p_val < alpha & abs(CAP_WTI_TPOm_vs_WTI_Vehicle$avg_log2FC) > fc, 
]

# Separate up-regulated and down-regulated genes
log2FC_t <- 0.25
up <- filtered_genes[filtered_genes$avg_log2FC > log2FC_t, ]
down <- filtered_genes[filtered_genes$avg_log2FC < -log2FC_t, ]

# Display separately up-regulated and down-regulated genes
print(head(up))

print(head(down))

#Creating a list of gene names
up_l <- rownames(up)

BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)

#Select genes 
gene_set <- up_l
pathway_db <- "org.Mm.eg.db"  # Use appropriate organism database

# Perform gene overrepresentation analysis
enrich_result_CAP_up <- enrichGO(gene = gene_set,
                          OrgDb = pathway_db,
                          keyType = "SYMBOL",
                          ont = "BP",  # Biological Process ontology
                          pAdjustMethod = "BH")

# Extract the top pathways of interest
top_pathways_CAP_up <- head(enrich_result_CAP_up@result, 7)
View(enrich_result_CAP_up@result)

# Define the specific values you want to extract
specific_values_up <- c("phosphatidylinositol-mediated signaling", "regulation of autophagy", "regulation of endothelial cell migration", "small GTPase mediated signal transduction", "lipid transport", "regulation of Ras protein signal transduction", "platelet-derived growth factor receptor signaling pathway")


filtered_df <- subset(enrich_result_CAP_up, Description %in% specific_values_up)

# Extract rows corresponding to the specific values
filtered_df_up <- enrich_result_CAP_up[enrich_result_CAP_up$Description %in% specific_values_up, ]

# View the filtered data frame
print(filtered_df)


# Set larger bottom margin to make more space for X-axis labels
# Adjust the bottom margin (default is c(5,4,4,2) + 0.1)

# Save the plot as PNG with increased size
png("Figures/EnrichGO-CAP_WTI_TPOm_vs_WTI_Vehicle_up_.25fc.png", width = 800, height = 600)  # Adjust width and height as needed
par(mar = c(4, 22, 4, 1) + 0.1, mgp = c(1, 3, 0))
barplot(rev(-log10(top_pathways_CAP_up$pvalue)),
        names.arg = rev(top_pathways_CAP_up$Description),
        horiz = T,
        las = 1, cex.names = 1,
        main = "TPOm vs Vehicle Upregulated Pathways",
        xlim = c(0,4),
        ylab = "-log10(pvalue)",
        col = "lightpink",
        cex.axis = 0.8)

dev.off()  # Close the PNG device

write.xlsx(enrich_result_CAP_up@result, "EnrighGO-CAP_WTI_TPOm_vs_WTI_Vehicle_up_.25fc.xlsx", rowNames = FALSE)

#Creating a list of gene names
down_l <- rownames(down)

#Select genes 
gene_set <- down_l
pathway_db <- "org.Mm.eg.db"  # Use appropriate organism database

# Perform gene overrepresentation analysis
enrich_result_CAP_down <- enrichGO(gene = gene_set,
                          OrgDb = pathway_db,
                          keyType = "SYMBOL",
                          ont = "BP",  # Biological Process ontology
                          pAdjustMethod = "BH")

# Extract the top 10 pathways
top_pathways <- head(enrich_result_CAP_down@result, 10)

# Define the specific values you want to extract
specific_values_down <- c("cytoplasmic translation", "protein folding", "regulation of endopeptidase activity", "response to heat", "tumor necrosis factor superfamily cytokine production", "G1/S transition of mitotic cell cycle", "regulation of leukocyte cell-cell adhesion")

# Extract rows corresponding to the specific values
filtered_df_down <- enrich_result_CAP_down[enrich_result_CAP_down$Description %in% specific_values_down, ]

# View the filtered data frame
print(filtered_df_down)


# Set larger bottom margin to make more space for X-axis labels
# Adjust the bottom margin (default is c(5,4,4,2) + 0.1)

# Save the plot as PNG with increased size
png("Figures/EnrighGO-CAP_WTI_TPOm_vs_WTI_Vehicle_down_.25fc.png", width = 800, height = 600)  # Adjust width and height as needed
par(mar = c(4, 22, 4, 2) + 0.1, mgp = c(1, 3, 0))
barplot(rev(-log10(filtered_df_down$pvalue)),
        names.arg = rev(filtered_df_down$Description),
        horiz = T,
        las = 1, cex.names = 1,
        main = "TPOm vs Vehicle Downregulated Pathways",
        xlim = c(0, 10),
        ylab = "-log10(P-value)",
        col = "lightblue2",
        cex.axis = 0.8)

dev.off()  # Close the PNG device

write.xlsx(enrich_result_CAP_down@result, "EnrighGO-CAP_WTI_TPOm_vs_WTI_Vehicle_down_.25fc.xlsx", rowNames = FALSE)

# Set the log2 fold change threshold
alpha <- 0.05
fc <- 0.25

# Filter for significant genes with a meaningful log2 fold change and no NAs
filtered_genes <- CAP_WTI_TPOm_vs_WTI_Vehicle[
  CAP_WTI_TPOm_vs_WTI_Vehicle$p_val < alpha & abs(CAP_WTI_TPOm_vs_WTI_Vehicle$avg_log2FC) > fc, 
]

# Separate up-regulated and down-regulated genes
log2FC_t <- 0.5
up <- filtered_genes[filtered_genes$avg_log2FC > log2FC_t, ]
down <- filtered_genes[filtered_genes$avg_log2FC < -log2FC_t, ]

# Display separately up-regulated and down-regulated genes
print(head(up))

print(head(down))

#Creating a list of gene names
up_l <- rownames(up)

BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)

#Select genes 
gene_set <- up_l
pathway_db <- "org.Mm.eg.db"  # Use appropriate organism database

# Perform gene overrepresentation analysis
enrich_result_CAP_up <- enrichGO(gene = gene_set,
                          OrgDb = pathway_db,
                          keyType = "SYMBOL",
                          ont = "BP",  # Biological Process ontology
                          pAdjustMethod = "BH")

# Extract the top 10 pathways
top_pathways_CAP_up <- head(enrich_result_CAP_up@result, 10)


# Set larger bottom margin to make more space for X-axis labels
# Adjust the bottom margin (default is c(5,4,4,2) + 0.1)

# Save the plot as PNG with increased size
png("Figures/EnrichGO-CAP_WTI_TPOm_vs_WTI_Vehicle_up_WIP.png", width = 1500, height = 1500)  # Adjust width and height as needed
par(mar = c(4, 22, 4, 1) + 0.1, mgp = c(1, 3, 0))
barplot(rev(-log10(top_pathways_CAP_up$pvalue)),
        names.arg = rev(top_pathways_CAP_up$Description),
        horiz = T,
        las = 1, cex.names = 1,
        main = "TPOm vs Vehicle Upregulated Pathways",
        xlim = c(0,4),
        ylab = "-log10(pvalue)",
        col = "lightpink",
        cex.axis = 0.8)

dev.off()  # Close the PNG device

write.xlsx(enrich_result_CAP_up@result, "EnrighGO-CAP_WTI_TPOm_vs_WTI_Vehicle_up.xlsx", rowNames = FALSE)

#Creating a list of gene names
down_l <- rownames(down)

#Select genes 
gene_set <- down_l
pathway_db <- "org.Mm.eg.db"  # Use appropriate organism database

# Perform gene overrepresentation analysis
enrich_result_CAP_down <- enrichGO(gene = gene_set,
                          OrgDb = pathway_db,
                          keyType = "SYMBOL",
                          ont = "BP",  # Biological Process ontology
                          pAdjustMethod = "BH")

# Extract the top 10 pathways
top_pathways <- head(enrich_result_CAP_down@result, 10)


# Set larger bottom margin to make more space for X-axis labels
# Adjust the bottom margin (default is c(5,4,4,2) + 0.1)

# Save the plot as PNG with increased size
png("Figures/EnrighGO-CAP_WTI_TPOm_vs_WTI_Vehicle_down_WIP.png", width = 1500, height = 1500)  # Adjust width and height as needed
par(mar = c(4, 22, 4, 2) + 0.1, mgp = c(1, 3, 0))
barplot(rev(-log10(top_pathways$pvalue)),
        names.arg = rev(top_pathways$Description),
        horiz = T,
        las = 1, cex.names = 1,
        main = "TPOm vs Vehicle Downregulated Pathways",
        xlim = c(0, 10),
        ylab = "-log10(P-value)",
        col = "lightblue2",
        cex.axis = 0.8)

dev.off()  # Close the PNG device

write.xlsx(enrich_result_CAP_down@result, "EnrighGO-CAP_WTI_TPOm_vs_WTI_Vehicle_down.xlsx", rowNames = FALSE)




##Do the clustering between CAP cells

x2 <- FindNeighbors(CAP, dims = 1:4)
x2 <- FindClusters(CAP, resolution =0.2)

png(filename = "Figures/UMP_plot_with_CAP_split_WIP.png", width=1500, height=1500)
DimPlot(x2, pt.size = 6, label = TRUE, label.size = 20, reduction = "umap", raster=FALSE)

dev.off()



# find markers for every cluster compared to all remaining cells, report only the positive ones
x.markers <- FindAllMarkers(x2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)

#Saving markers for each cluster
selected_markers <- print(x.markers %>% group_by(cluster) %>% slice_max(n = 15, order_by = avg_log2FC), n = 300)

write.csv(selected_markers, "selected_VehiclevsNaive_markers.csv")


Idents(object = x2) <- "seurat_clusters"
levels(x2)

#Renaming identitites

CAPcells <- c("CapEC1", 
           "CapEC2",
           "CapEC3",
           "CapEC4")

names(CAPcells) <- levels(x2)
x3 <- RenameIdents(x2, CAPcells)

# Stash cell identity classes
x2[["CAPcells"]] <- Idents(object = x3)
rm(x3)

Idents(object = x2) <- "CAPcells"
levels(x2)

x2 <- ScaleData(x2, features = rownames(x2))
#Heatmap of the clusters
png(filename = "Figures/Heatmap_clusters_names_CAP_WIP.png", width=5, height=6, units = "in", res = 300)
x.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(x2, features = top10$gene, label = TRUE, size = 3, angle = 90)
dev.off()


#####Calculate CAP cell type

cell_type_percentages <- x2@meta.data %>%
  group_by(samples, CAPcells) %>%
  summarise(cell_type_count = n()) %>%
  group_by(samples) %>%
  mutate(total_cells = sum(cell_type_count),
         cell_type_percentage = (cell_type_count / total_cells) * 100) %>%
  ungroup()

#Column plot
png(filename = "Figures/CAPCell_percentages_split_WIP.png", width=1000, height=1000)
ggplot(cell_type_percentages, aes(fill = CAPcells, y = cell_type_percentage, x = samples)) +
  geom_col(position = "stack") +
  labs(x = "Sample", y = "Cell Type Percentage", title = "Cell Type Percentages by Sample") +
  geom_text(aes(label = paste0(round(cell_type_percentage, 1), "%")),
            position = position_stack(vjust = 0.5), color = "white", size = 14) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        axis.text.y = element_text(size = 11))
dev.off()

#Stacked column plot
png(filename = "Figures/CAPCell_percentages_split_stacked.png", width=800, height=1000)
ggplot(cell_type_percentages, aes(fill = CAPcells, y = cell_type_percentage, x = samples)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(x = "Sample", y = "Cell Type Percentage", title = "Cell Type Percentages by Sample") +
  geom_text(aes(label = paste0(round(cell_type_percentage, 1), "%")),
            position = position_fill(vjust = 0.5), color = "white") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

Idents(object = CAPEC1) <- "samples"
Idents(object = CAPEC2) <- "samples"
Idents(object = CAPEC3) <- "samples"
Idents(object = CAPEC4) <- "samples"

features <- c("Hspa1a", "Hspa1b")

endo_clusters <- c("CapEC1")
CAPEC1 <- subset(x2, idents = endo_clusters)
table(CAP$cells)
levels(CAPEC1)

endo_clusters <- c("CapEC2")
CAPEC2 <- subset(x2, idents = endo_clusters)

levels(CAPEC2)

endo_clusters <- c("CapEC3")
CAPEC3 <- subset(x2, idents = endo_clusters)

levels(CAPEC3)

endo_clusters <- c("CapEC4")
CAPEC4 <- subset(x2, idents = endo_clusters)

levels(CAPEC4)

png(filename = "Figures/Violin_plot_CAP_Hspa1a1_CapeC1.png", width=375, height=400)

VlnPlot(CAP, features = c("Inhbb"), cols = c('gray', 'deeppink', 'cyan4'), pt.size = 0, flip = TRUE)

VlnPlot(CAP, features = c("Inhba"), cols = c('gray', 'deeppink', 'cyan4'), pt.size = 0, flip = TRUE)

VlnPlot(CAP, features = c("Acvr1"), cols = c('gray', 'deeppink', 'cyan4'), pt.size = 0, flip = TRUE)

VlnPlot(CAP, features = c("Gdf2"), cols = c('gray', 'deeppink', 'cyan4'), pt.size = 0, flip = TRUE)

VlnPlot(CAP, features = c("Gdf15"), cols = c('gray', 'deeppink', 'cyan4'), pt.size = 0, flip = TRUE)

VlnPlot(CAP, features = c("Fstl3"), cols = c('gray', 'deeppink', 'cyan4'), pt.size = 0, flip = TRUE)

dev.off()

png(filename = "Figures/Violin_plot_CAP_Hspa1a1_CapeC1.png", width=375, height=400)
VlnPlot(CAPEC1, features = c("Hspa1a"), cols = c('gray', 'deeppink', 'cyan4'), pt.size = 0, flip = TRUE)
dev.off()

png(filename = "Figures/Violin_plot_CAP_Hspa1a1_CapEC2.png", width=375, height=400)
VlnPlot(CAPEC2, features = c("Hspa1a"), cols = c('gray', 'deeppink', 'cyan4'), pt.size = 0, flip = TRUE)
dev.off()

png(filename = "Figures/Violin_plot_CAP_Hspa1a1_CapEC3.png", width=375, height=400)
VlnPlot(CAPEC3, features = c("Hspa1a"), cols = c('gray', 'deeppink', 'cyan4'), pt.size = 0, flip = TRUE)
dev.off()

png(filename = "Figures/Violin_plot_CAP_Hspa1a1_CapEC4.png", width=375, height=400)
VlnPlot(CAPEC4, features = c("Hspa1a"), cols = c('gray', 'deeppink', 'cyan4'), pt.size = 0, flip = TRUE)
dev.off()

png(filename = "Figures/Violin_plot_CAP_Hspa1a.png", width=600, height=400)
VlnPlot(CAP, features = c("Hspa1a"), cols = c('gray', 'deeppink', 'cyan4'), pt.size = 0, flip = TRUE)
dev.off()

png(filename = "Figures/Violin_plot_CAP_Hspa1b_WIP.png", width=600, height=400)
VlnPlot(CAP, features = c("Hspa1b"), cols = c('gray', 'deeppink', 'cyan4'), pt.size = 0, flip = TRUE)
dev.off()

png(filename = "Figures/Vln_plot_CAPcells_Hspa1a_WIP.png", width=300, height=250)
FeaturePlot(x2, features = c("Hspa1a"), pt.size = 0, raster=FALSE)
dev.off()

png(filename = "Figures/Vln_plot_CAPcells_Klf2.png", width=300, height=250)
VlnPlot(CAP, features = c("Klf2"), pt.size = 0, raster=FALSE)
dev.off()

png(filename = "Figures/Violin_plot_CAP_Itgb5.png", width=400, height=400)
VlnPlot(CAP, features = c("Itgb5"), pt.size = 0, raster=FALSE)
dev.off()

png(filename = "Figures/Violin_plot_CAP_Flt1.png", width=400, height=400)
VlnPlot(CAP, features = c("Flt1"), pt.size = 0, raster=FALSE)
dev.off()


png(filename = "Figures/Violin_plot_CAP_Asf3.png", width=400, height=400)
VlnPlot(CAP, features = c("Atf3"), pt.size = 0, raster=FALSE)
dev.off()

png(filename = "Figures/Violin_plot_CAP_Hspa1a1_thesis.png", width=375, height=400)
VlnPlot(CAP, features = c("Hspa1a"), cols = c('gray', 'deeppink', 'cyan4'), pt.size = 0, flip = TRUE)
dev.off()

png(filename = "Figures/Violin_plot_CAP_Hspa1b_thesis.png", width=375, height=400)
VlnPlot(CAP, features = c("Hspa1b"), cols = c('gray', 'deeppink', 'cyan4'), pt.size = 0, flip = TRUE)
dev.off()

png(filename = "Figures/Vln_plot_CAPcells_Hspa1a_thesis.png", width=350, height=400)
VlnPlot(x2, features = c("Hspa1a"), pt.size = 0, raster=FALSE)
dev.off()

png(filename = "Figures/Vln_plot_CAPcells_Hspa1b_thesis.png", width=350, height=400)
VlnPlot(x2, features = c("Hspa1b"), pt.size = 0, raster=FALSE)
dev.off()

View(CAP@meta.data)
