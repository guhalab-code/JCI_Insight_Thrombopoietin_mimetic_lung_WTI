####subsetting and defining CAP cells ####

x <- readRDS(file = "x_names.rds")

endo_clusters <- c("cap_EC")
CAP <- subset(x, idents = endo_clusters)
table(CAP$cells)
levels(CAP)

Idents(object = CAP) <- "samples"
levels(CAP)
head(CAP@meta.data)

##Do the clustering between CAP cells

x2 <- FindNeighbors(CAP, dims = 1:4)
x2 <- FindClusters(CAP, resolution =0.2)

features <- c("Hspa1a", "Hspa1b")

x2 <- RunUMAP(x2, dims = 1:10)

png(filename = "Figures/UMP_plot_with_Cap_split_WIP.png", width=3600, height=1200)
DimPlot(x2, pt.size = 10, reduction = "umap", split.by = "samples", raster=FALSE)
dev.off()

png(filename = "Figures/Feature_plot_CAP_Hspa1a.png", width=1000, height=250)
FeaturePlot(x2, features = c("Hspa1a"), pt.size = 1, split.by = "samples", raster=FALSE)
dev.off()

png(filename = "Figures/Feature_plot_CAP_Cxcl1.png", width=1000, height=250)
FeaturePlot(x2, features = c("Mcp1"), pt.size = 2, split.by = "samples", raster=FALSE)
dev.off()


png(filename = "Figures/Heatmap_clusters_names_thesis.png", width=1000, height=1500)
x.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(x2, features = top10$gene, size = 12, angle = 90, ) + NoLegend()
dev.off()

# find markers for every cluster compared to all remaining cells, report only the positive ones
x.markers <- FindAllMarkers(x2, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.5)

View(x.markers)

# rename clusters
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

selected_markers <- print(x.markers %>% group_by(cluster) %>% slice_max(n = 15, order_by = avg_log2FC), n = 300)

write.csv(selected_markers, "selected_CAP_markers.csv")
VlnPlot(x2, features = "Cdkn1a", group.by = "samples")

png(filename = "Figures/Violin_plot_CAP_cdkn1a.png", width=400, height=400)
VlnPlot(x2, features = c("Cdkn1a"), pt.size = 0, group.by = "CAPcells", raster=FALSE)
dev.off()

# Set the log2 fold change threshold
alpha <- 0.05
fc <- 0.25

# Filter for significant genes with a meaningful log2 fold change and no NAs
filtered_genes <- x.markers[
  x.markers$p_val < alpha & abs(x.markers$avg_log2FC) > fc, 
]

# Separate up-regulated and down-regulated genes
log2FC_t <- 0.25
up <- filtered_genes[filtered_genes$avg_log2FC > log2FC_t, ]

#separate out each subcluster from dataframe
CapEC1 <- up[up$cluster == "0",]
CapEC2 <- up[up$cluster == "1",]
CapEC3 <- up[up$cluster == "2",]
CapEC4 <- up[up$cluster == "3",]

#Creating a list of gene names
up_CapEC1 <- rownames(CapEC1)
up_CapEC2 <- rownames(CapEC2)
up_CapEC3 <- rownames(CapEC3)
up_CapEC4 <- rownames(CapEC4)


#Select genes 
gene_set1 <- up_CapEC1
gene_set2 <- up_CapEC2
gene_set3 <- up_CapEC3
gene_set4 <- up_CapEC4

pathway_db <- "org.Mm.eg.db"  # Use appropriate organism database

# Perform gene overrepresentation analysis
enrich_result_CapEC1 <- enrichGO(gene = gene_set1,
                                 OrgDb = pathway_db,
                                 keyType = "SYMBOL",
                                 ont = "BP",  # Biological Process ontology
                                 pAdjustMethod = "BH")

# Extract the pathways of interest
top_pathways_CapEC1 <- head(enrich_result_CapEC1@result, 10)


View(enrich_result_CapEC1@result)

write.xlsx(enrich_result_CapEC1@result, "EnrighGO-CAPEC1.xlsx", rowNames = FALSE)

# Save the plot as PNG with increased size
png("Figures/EnrichGO-CapEC1.png", width = 0, height = 600)  # Adjust width and height as needed
par(mar = c(4, 22, 4, 1) + 0.1, mgp = c(1, 3, 0))
barplot(rev(-log10(top_pathways_CapEC1$pvalue)),
        names.arg = rev(top_pathways_CapEC1$Description),
        horiz = T,
        las = 1, cex.names = 1,
        main = "CapEC1",
       
        ylab = "-log10(pvalue)",
        col = "red1",
        cex.axis = 0.8)

dev.off()  # Close the PNG device

# Perform gene overrepresentation analysis
enrich_result_CapEC2 <- enrichGO(gene = gene_set2,
                                 OrgDb = pathway_db,
                                 keyType = "SYMBOL",
                                 ont = "BP",  # Biological Process ontology
                                 pAdjustMethod = "BH")

# Extract the top 10 pathways
top_pathways_CapEC2 <- head(enrich_result_CapEC2@result, 7)

View(enrich_result_CapEC2@result)

write.xlsx(enrich_result_CapEC2@result, "EnrighGO-CAPEC2.xlsx", rowNames = FALSE)

# Save the plot as PNG with increased size
png("Figures/EnrichGO-CapEC2.png", width = 550, height = 600)  # Adjust width and height as needed
par(mar = c(4, 22, 4, 1) + 0.1, mgp = c(1, 3, 0))
barplot(rev(-log10(top_pathways_CapEC2$pvalue)),
        names.arg = rev(top_pathways_CapEC2$Description),
        horiz = T,
        las = 1, cex.names = 1,
        main = "CapEC2",
        
        ylab = "-log10(pvalue)",
        col = "red1",
        cex.axis = 0.8)

dev.off()  # Close the PNG device

# Perform gene overrepresentation analysis
enrich_result_CapEC3 <- enrichGO(gene = gene_set3,
                                 OrgDb = pathway_db,
                                 keyType = "SYMBOL",
                                 ont = "BP",  # Biological Process ontology
                                 pAdjustMethod = "BH")

# Extract the top 10 pathways
top_pathways_CapEC3 <- head(enrich_result_CapEC3@result, 7)
View(enrich_result_CapEC3@result)
write.xlsx(enrich_result_CapEC3@result, "EnrighGO-CAPEC3.xlsx", rowNames = FALSE)

# Save the plot as PNG with increased size
png("Figures/EnrichGO-CapEC3.png", width = 550, height = 600)  # Adjust width and height as needed
par(mar = c(4, 22, 4, 1) + 0.1, mgp = c(1, 3, 0))
barplot(rev(-log10(top_pathways_CapEC3$pvalue)),
        names.arg = rev(top_pathways_CapEC3$Description),
        horiz = T,
        las = 1, cex.names = 1,
        main = "CapEC3",
        
        ylab = "-log10(pvalue)",
        col = "red1",
        cex.axis = 0.8)

dev.off()  # Close the PNG device

# Perform gene overrepresentation analysis
enrich_result_CapEC4 <- enrichGO(gene = gene_set4,
                                 OrgDb = pathway_db,
                                 keyType = "SYMBOL",
                                 ont = "BP",  # Biological Process ontology
                                 pAdjustMethod = "BH")

# Extract the top 10 pathways
top_pathways_CapEC4 <- head(enrich_result_CapEC4@result, 7)

View(enrich_result_CapEC4@result)

write.xlsx(enrich_result_CapEC4@result, "EnrighGO-CAPEC4.xlsx", rowNames = FALSE)

# Save the plot as PNG with increased size
png("Figures/EnrichGO-CapEC4.png", width = 550, height = 600)  # Adjust width and height as needed
par(mar = c(4, 22, 4, 1) + 0.1, mgp = c(1, 3, 0))
barplot(rev(-log10(top_pathways_CapEC4$pvalue)),
        names.arg = rev(top_pathways_CapEC4$Description),
        horiz = T,
        las = 1, cex.names = 1,
        main = "CapEC4",
        
        ylab = "-log10(pvalue)",
        col = "red1",
        cex.axis = 1)

dev.off()  # Close the PNG device

#Cell cycle genes import
m.s.genes <- c("Mcm5", "Pcna", "Tyms", "Fen1", "Mcm2", "Mcm4", "Rrm1", "Ung", "Gins2", "Mcm6", "Cdca7", "Dtl", "Prim1", "Uhrf1", "Mlf1ip", "Hells", "Rfc2", "Rpa2", "Nasp", "Ras51ap1", "Gmnn", "Wdr76", "Slbp",  "Ccne2", "Ubr7", "Pold3", "Msh2", "Atad2", "Rad51", "Rrm2", "Cdc45", "Cdc6", "Exo1", "Tipin", "Dscc1", "Blm", "CaspP8ap2", "Usp1", "Clspn", "Pola1", "Chaf1b", "Brip1", "E2f8")
m.g2m.genes <- c("Hmgb2", "Cdk1"  ,"Nusap1" ,"Ube2c", "Birc5" ,"Tpx2", "Top2a", "Ndc80" ,"Cks2" , "Nuf2"  , "Cks1b" , "Mki67" , "Tmpo", "Cenpf" , "Tacc3" , "Fam64a",  "Smc4", "Ccnb2", "Ckap2l","Ckap2", "Aurkb" , "Bub1"  , "Kif11" ,"Anp32e" , "Tubb4b", "Gtse1" ,"Kif20b" ,"Hjurp"  , "Cdca3" ,"Hn1", "Cdc20" ,"Ttk"  , "Cdc25c", "Kif2c" , "Rangap1" , "Ncapd2", "Dlgap5", "Cdca2" , "Cdca8" ,"Ect2" , "Kif23" , "Hmmr"  , "Aurka" , "Psrc1" , "Anln", "Lbr",  "Ckap5" , "Cenpe" ,"Ctcf" , "Nek2"  , "G2e3"  , "Gas2l3", "Cbx5"  , "Cenpa")

x2 <- CellCycleScoring(x2, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(x2[[]])

#Umap plot with samples name
png(filename = "Figures/UMP_plot_with_cell_cycle_split_condition_thesis.png", width=1600, height=800)
DimPlot(x2, reduction = "umap", group.by = "Phase", raster=FALSE)
dev.off()

#Column plot
png(filename = "Figures/Cell_cycle_split.png", width=600, height=1000)
ggplot(x2, aes(fill = cells, y = cell_type_percentage, x = samples)) +
  geom_col(position = "stack") +
  labs(x = "Sample", y = "Cell Type Percentage", title = "Cell Type Percentages by Sample") +
  geom_text(aes(label = paste0(round(cell_type_percentage, 1), "%")),
            position = position_stack(vjust = 0.5), color = "white", size = 3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(cells ~ ., scales = "free_y")
dev.off()




count <- sum(x$samples == "Naive")
print(count)

count <- sum(x$samples == "WTI_Vehicle")
print(count)

count <- sum(x$samples == "WTI_TPOm")
print(count)
