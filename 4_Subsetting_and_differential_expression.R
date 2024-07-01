setwd('/Users/jebenglish/Dropbox (EinsteinMed)/Jeb-scRNAseq_code')

#Importing libraries
library(Seurat)
library(SeuratWrappers)
library(patchwork)
library(dplyr)
library(openxlsx)
library(readxl)
library(ggplot2)  
library(dplyr)
library(ggrepel)
library(tibble)
library(MAST)
library(SCPA)
library(msigdbr)
library(cowplot)

set.seed(1234)

#Load .rds file saved after last analysis
x <- readRDS(file = "x_names.rds")

levels(x)
head(x@meta.data)
table(x$cells)


#Changing id to cluster numbers
Idents(object = x) <- "samples"
levels(x)
WTI <- c("WTI_Vehicle")
WTIcells <- subset(x, idents = WTI)
WTIcells

#####################
#Subsetting AT
AT_clusters <- c("AT1", "AT2")
AT <- subset(x, idents = AT_clusters)
table(AT$cells)
levels(AT)

Idents(object = AT) <- "samples"
levels(AT)
head(AT@meta.data)

# Find differentially expressed features between naive and vehicle # might add: min.diff.pct = 0.25, max.cells.per.ident = 200
latent_vars <- c("nCount_RNA", "percent.mt", "S.Score")
AT_WTI_TPOm_vs_WTI_Vehicle <- FindMarkers(AT, ident.1 = "WTI_TPOm", ident.2 = "WTI_Vehicle",
                                           min.pct = 0.25, 
                                           logfc.threshold = log(0.5),
                                           test.use = "MAST",
                                           latent.vars = latent_vars)
head(MSC_WTI_TPOm_vs_WTI_Vehicle)

#Volcano plot
alpha <- 0.05
fc <- 0.5
df <- rownames_to_column(AT_WTI_TPOm_vs_WTI_Vehicle, var = "gene")
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
                  aes(label = gene), color = "red", size = 3, max.overlaps = 10) +
  theme_minimal()+
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"))

ggsave("Figures/Volcano_Plot_AT_WTI_TPOm_vs_WTI_Vehicle.png", plot = volcano_plot, width = 2200, height = 1800, units = "px")
rm(volcano_plot)

######A##########
pathways <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP") %>%
  format_pathways()

# Filter cells based on metadata
AT_Naive <- subset(AT, subset = samples == "Naive")
# Access assay data for the filtered cells
AT_Naive_data <- as.matrix(GetAssayData(object = AT_Naive, assay = "RNA"))

# Filter cells based on metadata
AT_WTI_TPOm <- subset(AT, subset = samples == "WTI_TPOm")
# Access assay data for the filtered cells
AT_WTI_TPOm_data <- as.matrix(GetAssayData(object = AT_WTI_TPOm, assay = "RNA"))

# Filter cells based on metadata
AT_WTI_Vehicle <- subset(AT, subset = samples == "WTI_Vehicle")
# Access assay data for the filtered cells
AT_WTI_Vehicle_data <- as.matrix(GetAssayData(object = AT_WTI_Vehicle, assay = "RNA"))

#Comparing pathways
NaivevVeh_AT_out <- compare_pathways(samples = list(AT_WTI_Vehicle_data, AT_Naive_data),
                                       pathways = pathways)

TPOmvVeh_AT_out <- compare_pathways(samples = list(AT_WTI_TPOm_data, AT_WTI_Vehicle_data),
                                       pathways = pathways)

head(TPOmvVeh_AT_out, 15)

# Filter pathways based on qval # c(1, 5)
filtered_pathways2 <- TPOmvVeh_AT_out[1:15, ]
filtered_pathways2$FC[filtered_pathways2$FC >100] <- 100
filtered_pathways2$FC[filtered_pathways2$FC < -100] <- -100
filtered_pathways2$qval[filtered_pathways2$qval > 7] <- 7

# Create a new variable for pathway order based on FC values
filtered_pathways2$Pathway_order <- with(filtered_pathways2, reorder(Pathway, FC))


#Generating heatmap
png(filename = "Figures/Heatmap_AT_WTI_TPOm vs WTI_Vehicle.png", width=600, height=800)
ggplot(filtered_pathways2, aes(x = 1, y = Pathway_order, fill = FC)) + 
  geom_tile() + 
  ylab("") +
  xlab("") +
  scale_x_continuous(breaks = NULL) +
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red",
                       breaks = c(-5, -2.5, 0, 2.5, 5),
                       limits = c(-5, 5)) +
  geom_point(aes(x = 1, y = Pathway_order, size = qval)) +
  scale_size(limits = c(1, 2)) +
  theme(axis.text.y = element_text(size = 9),
        panel.background = element_rect(fill = "white"))
dev.off()


#####################
#Subsetting MSCs
MSC_clusters <- c("MSCs", "fibroblasts_col14", "fibroblasts_col13", "fibroblasts_3")
MSC <- subset(x, idents = MSC_clusters)
table(MSC$cells)
levels(MSC)

Idents(object = MSC) <- "samples"
levels(MSC)
head(MSC@meta.data)

# Find differentially expressed features between naive and vehicle # might add: min.diff.pct = 0.25, max.cells.per.ident = 200
latent_vars <- c("nCount_RNA", "percent.mt", "S.Score")
MSC_WTI_TPOm_vs_WTI_Vehicle <- FindMarkers(MSC, ident.1 = "WTI_TPOm", ident.2 = "WTI_Vehicle",
                                            min.pct = 0.25, 
                                            logfc.threshold = log(0.5),
                                            test.use = "MAST",
                                            latent.vars = latent_vars)
head(MSC_WTI_TPOm_vs_WTI_Vehicle)

#Volcano plot
alpha <- 0.05
fc <- 0.5
df <- rownames_to_column(MSC_WTI_TPOm_vs_WTI_Vehicle, var = "gene")
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
                  aes(label = gene), color = "red", size = 3, max.overlaps = 10) +
  theme_minimal()+
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"))

ggsave("Figures/Volcano_Plot_MSC_WTI_TPOm_vs_WTI_Vehicle.png", plot = volcano_plot, width = 2200, height = 1800, units = "px")
rm(volcano_plot)

# Filter cells based on metadata
MSC_Naive <- subset(MSC, subset = samples == "Naive")
# Access assay data for the filtered cells
MSC_Naive_data <- as.matrix(GetAssayData(object = MSC_Naive, assay = "RNA"))

# Filter cells based on metadata
MSC_WTI_TPOm <- subset(MSC, subset = samples == "WTI_TPOm")
# Access assay data for the filtered cells
MSC_WTI_TPOm_data <- as.matrix(GetAssayData(object = MSC_WTI_TPOm, assay = "RNA"))

# Filter cells based on metadata
MSC_WTI_Vehicle <- subset(MSC, subset = samples == "WTI_Vehicle")
# Access assay data for the filtered cells
MSC_WTI_Vehicle_data <- as.matrix(GetAssayData(object = MSC_WTI_Vehicle, assay = "RNA"))

#Comparing pathways
MSC_out <- compare_pathways(samples = list(MSC_WTI_TPOm_data, MSC_WTI_Vehicle_data),
                             pathways = pathways)

VehvNaive_MSC_out <- compare_pathways(samples = list(MSC_WTI_Vehicle_data, MSC_Naive_data),
                                   pathways = pathways)

head(MSC_out, 15)

# Filter pathways based on qval # c(1, 5)
filtered_pathways2 <- MSC_out[1:15, ]
filtered_pathways2$FC[filtered_pathways2$FC >100] <- 100
filtered_pathways2$FC[filtered_pathways2$FC < -100] <- -100
filtered_pathways2$qval[filtered_pathways2$qval > 7] <- 7

# Create a new variable for pathway order based on FC values
filtered_pathways2$Pathway_order <- with(filtered_pathways2, reorder(Pathway, FC))


#Generating heatmap
png(filename = "Figures/Heatmap_MSC_WTI_TPOm vs WTI_Vehicle.png", width=600, height=800)
ggplot(filtered_pathways2, aes(x = 1, y = Pathway_order, fill = FC)) + 
  geom_tile() + 
  ylab("") +
  xlab("") +
  scale_x_continuous(breaks = NULL) +
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red",
                       breaks = c(-5, -2.5, 0, 2.5, 5),
                       limits = c(-5, 5)) +
  geom_point(aes(x = 1, y = Pathway_order, size = qval)) +
  scale_size(limits = c(1, 2)) +
  theme(axis.text.y = element_text(size = 9),
        panel.background = element_rect(fill = "white"))
dev.off()

#################################MSC###################################################
#Subsetting ECs
endo_clusters <- c("endo1", "endo2", "cap_EC")
endo <- subset(x, idents = endo_clusters)
table(endo$cells)
levels(endo)

Idents(object = endo) <- "samples"
levels(endo)
head(endo@meta.data)


# Find differentially expressed features between naive and vehicle # might add: min.diff.pct = 0.25, max.cells.per.ident = 200
latent_vars <- c("nCount_RNA", "percent.mt", "S.Score")
endo_WTI_TPOm_vs_WTI_Vehicle <- FindMarkers(endo, ident.1 = "WTI_TPOm", ident.2 = "WTI_Vehicle",
                                   min.pct = 0.25, 
                                   logfc.threshold = log(0.5),
                                   test.use = "MAST",
                                   latent.vars = latent_vars)
head(endo_WTI_TPOm_vs_WTI_Vehicle)

#Volcano plot
alpha <- 0.05
fc <- 1
df <- rownames_to_column(endo_WTI_TPOm_vs_WTI_Vehicle, var = "gene")
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
                  aes(label = gene), color = "red", size = 3, max.overlaps = 10) +
  theme_minimal()+
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"))

ggsave("Figures/Volcano_Plot_endo_WTI_TPOm_vs_WTI_Vehicle.png", plot = volcano_plot, width = 2200, height = 1800, units = "px")
rm(volcano_plot)

#Saving results
write.xlsx(EC_TPOm_vs_Vehicle, "endo_WTI_TPOm_vs_WTI_Vehicle.xlsx", rowNames = TRUE)

######A##########
pathways <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP") %>%
  format_pathways()

# Filter cells based on metadata
endo_Naive <- subset(endo, subset = samples == "Naive")
# Access assay data for the filtered cells
endo_Naive_data <- as.matrix(GetAssayData(object = endo_Naive, assay = "RNA"))

# Filter cells based on metadata
endo_WTI_TPOm <- subset(endo, subset = samples == "WTI_TPOm")
# Access assay data for the filtered cells
endo_WTI_TPOm_data <- as.matrix(GetAssayData(object = endo_WTI_TPOm, assay = "RNA"))

# Filter cells based on metadata
endo_WTI_Vehicle <- subset(endo, subset = samples == "WTI_Vehicle")
# Access assay data for the filtered cells
endo_WTI_Vehicle_data <- as.matrix(GetAssayData(object = endo_WTI_Vehicle, assay = "RNA"))

#Comparing pathways
NaivevVeh_endo_out <- compare_pathways(samples = list(endo_WTI_Vehicle_data, endo_Naive_data),
                        pathways = pathways)

NaivevVeh_endo_out <- compare_pathways(samples = list(endo_WTI_TPOm_data, endo_WTI_Vehicle_data),
                                       pathways = pathways)

#subsetting pathways of interest
NaivevVeh_endo_out <- mutate(NaivevVeh_endo_out, Comparison = 'Vehicle_vs_Naive')
conditions <- c("GOBP_RESPONSE_TO_TRANSFORMING_GROWTH_FACTOR_BETA", 
                "GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_SIGNALING_PATHWAY",
                "GOBP_REGULATION_OF_RESPONSE_TO_CYTOKINE_STIMULUS", 
                "GOBP_RESPONSE_TO_INTERLEUKIN_1", 
                "GOBP_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_BY_P53_CLASS_MEDIATOR", 
                "GOBP_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_RESPONSE_TO_DNA_DAMAGE", 
                "GOBP_NEGATIVE_REGULATION_OF_GROWTH")
selected_rows <- NaivevVeh_endo_out[NaivevVeh_endo_out$Pathway ==  "GOBP_CELLULAR_RESPONSE_TO_CHEMICAL_STRESS", 
                                    NaivevVeh_endo_out$Pathway ==  "GOBP_RESPONSE_TO_OXIDATIVE_STRESS",
                                      NaivevVeh_endo_out$Pathway == "GOBP_RESPONSE_TO_TRANSFORMING_GROWTH_FACTOR_BETA",
                                      NaivevVeh_endo_out$Pathway == "GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_SIGNALING_PATHWAY",
                                      NaivevVeh_endo_out$Pathway == "GOBP_REGULATION_OF_RESPONSE_TO_CYTOKINE_STIMULUS",
                                    ]
selected_rows <- NaivevVeh_endo_out[NaivevVeh_endo_out$Pathway %in% c("GOBP_RESPONSE_TO_TRANSFORMING_GROWTH_FACTOR_BETA", 
                                                                      "GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_SIGNALING_PATHWAY",
                                                                      "GOBP_REGULATION_OF_RESPONSE_TO_CYTOKINE_STIMULUS",
                                                                      "GOBP_CELLULAR_RESPONSE_TO_CHEMICAL_STRESS", 
                                                                      "GOBP_RESPONSE_TO_OXIDATIVE_STRESS"), ]

endo_out <- mutate(endo_out, Comparison = 'TPOm_vs_Vehicle')
selected_rows2 <- endo_out[endo_out$Pathway %in% c("GOBP_RESPONSE_TO_TRANSFORMING_GROWTH_FACTOR_BETA", 
                                                   "GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_SIGNALING_PATHWAY",
                                                   "GOBP_REGULATION_OF_RESPONSE_TO_CYTOKINE_STIMULUS",
                                                   "GOBP_CELLULAR_RESPONSE_TO_CHEMICAL_STRESS", 
                                                   "GOBP_RESPONSE_TO_OXIDATIVE_STRESS"), ]

## combine selected rows of interest
combined_df <- rbind(selected_rows2, selected_rows)
#change names of pathways so
combined_df <- 
png(filename = "Figures/Heatmap_endo_all_comparison.png", width=600, height=400)
ggplot(combined_df, aes(x = Comparison, y = Pathway, fill = FC)) + 
  geom_tile() + 
  ylab("") +
  xlab("") +
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red",
                       breaks = c(-5, -2.5, 0, 2.5, 5),
                       limits = c(-7, 7)) +
  geom_point(aes(x = Comparison, y = Pathway, size = qval)) +
  scale_size(limits = c(.5, 3)) +
  theme(axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 11),
        panel.background = element_rect(fill = "white"))
dev.off()

head(endo_out, 15)

# Filter pathways based on qval # c(1, 5)
filtered_pathways2 <- NaivevVeh_endo_out[1:15, ]
filtered_pathways2$FC[filtered_pathways2$FC >100] <- 100
filtered_pathways2$FC[filtered_pathways2$FC < -100] <- -100
filtered_pathways2$qval[filtered_pathways2$qval > 7] <- 7

# Create a new variable for pathway order based on FC values
filtered_pathways2$Pathway_order <- with(filtered_pathways2, reorder(Pathway, FC))

#Generating heatmap
png(filename = "Figures/Heatmap_endo_WTI_TPOm vs WTI_Vehicle.png", width=600, height=800)
ggplot(combined_df, aes(x = Comparison, y = Pathway, fill = FC)) + 
  geom_tile() + 
  ylab("") +
  xlab("") +
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red",
                       breaks = c(-5, -2.5, 0, 2.5, 5),
                       limits = c(-5, 5)) +
  geom_point(aes(x = 1, y = Pathway, size = qval)) +
  scale_size(limits = c(1, 3)) +
  theme(axis.text.y = element_text(size = 9),
        panel.background = element_rect(fill = "white"))
dev.off()

write.xlsx(scpa_out, "endo_WTI_TPOm vs WTI_Vehicle.xlsx", rowNames = FALSE)

###########################END####################

