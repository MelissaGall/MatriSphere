# Clean, annotate and extract data from MC38 scRNA-seq

# Libraries ----

library(Seurat)
library(SeuratDisk)
library(openxlsx)
library(Matrix.utils)
library(SingleCellExperiment)
library(ggplot2)
library(loupeR)
library(dplyr)

# Load data ----

mc38_obj <- LoadH5Seurat("data/Seurat/WT_KO_MC38.h5seurat", verbose = F) # scRNA-seq file
annotation <- read.table("data/Seurat/seurat_clusters_annotation.csv", sep = ",", header = T) # custom annotation

# Filter ----

mc38_obj <- subset(mc38_obj, subset = mitoPercent < 10 & nFeature_RNA > 200 & nFeature_RNA < 7500)
mc38_obj <- FindNeighbors(mc38_obj, reduction = "pca", verbose = FALSE)
mc38_obj <- FindClusters(mc38_obj, verbose = FALSE, resolution = 0.5)
mc38_obj <- subset(mc38_obj, subset = Class == "WT.MC38") # Keep only WT (remove KO cells)

# Add annotation and remove unknown cluster ----

mc38_obj <- AddMetaData(mc38_obj, annotation)
mc38_obj@meta.data$Subtype <- gsub("-.*", "", mc38_obj@meta.data$seurat_clusters)

mc38_obj@meta.data <- mc38_obj@meta.data %>%
  mutate(Type = case_when(Subtype %in% c("MC38", "Prolif MC38", "Quiescent MC38") ~ "MC38",
                          Subtype %in% c("Endothelial cell", "Fibroblasts", "Pericyte") ~ "Stromal cells",
                          Subtype == "Unk" ~ "Unk",
                          TRUE ~ "Immune cells"))

mc38_obj <- subset(mc38_obj, subset = Type != "Unk")

# Pseudo-bulk ----

## Pseudo bulk for subclusters ----

sce <- as.SingleCellExperiment(mc38_obj)
groups <- data.frame(colData(sce)[, c("seurat_clusters")])
aggr_counts <- aggregate.Matrix(t(counts(sce)), groupings = groups, fun = "sum") 
write.table(t(aggr_counts), "data/pseudobulk_MC38_subclusters.txt", quote = F, sep = "\t", row.names = T, col.names = T)

## Pseudo bulk for subtypes ----

sce <- as.SingleCellExperiment(mc38_obj)
groups <- data.frame(colData(sce)[, c("Subtype")])
aggr_counts <- aggregate.Matrix(t(counts(sce)), groupings = groups, fun = "sum") 
write.table(t(aggr_counts), "data/pseudobulk_MC38_subtypes.txt", quote = F, sep = "\t", row.names = T, col.names = T)

## Pseudo bulk for types ----

sce <- as.SingleCellExperiment(mc38_obj)
groups <- data.frame(colData(sce)[, c("Type")])
aggr_counts <- aggregate.Matrix(t(counts(sce)), groupings = groups, fun = "sum") 
write.table(t(aggr_counts), "data/pseudobulk_MC38_types.txt", quote = F, sep = "\t", row.names = T, col.names = T)


# Save file ----

SaveH5Seurat(mc38_obj, "data/Seurat/WT_MC38_filtered_annotated", overwrite = T)

# Save cLoupe file ----

create_loupe_from_seurat(
  mc38_obj,
  output_dir = "data",
  output_name = "WT_MC38_filtered_annotated_cloupe",
  dedup_clusters = FALSE,
  executable_path = NULL,
  force = TRUE
)