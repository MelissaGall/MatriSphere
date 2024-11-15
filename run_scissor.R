# Run Scissor on scRNA-seq using mouse RNA-seq

# Libraries ----

library(Scissor)
library(Seurat)
library(SeuratDisk)
library(SeuratObject)
library(dplyr)
library(ggplot2)

# Prepare scRNA seq data ----

mc38_obj <- LoadH5Seurat("data/Seurat/WT_MC38_filtered_annotated.h5seurat")
mc38_obj_sub <- subset(mc38_obj, subset = Type == "MC38") # Subset to MC38 cells only


# Prepare bulk data for Scissor ----

mouse_rnaseq <- read.table("data/mouse_rawcounts.csv", header = TRUE, row.names = 1, sep = ",")
mouse_rnaseq[] <- sapply(mouse_rnaseq, as.integer)
mouse_rnaseq <- mouse_rnaseq %>% select(contains("MC38")) # Keep only MC38 cells
mouse_rnaseq <- mouse_rnaseq[, !(colnames(mouse_rnaseq) %in% c("CA_MC38_6", "ECM_MC38_2"))] # Remove outliers
keep <- rowSums(mouse_rnaseq >= 1) >= 3 # Remove low counts rows
mouse_rnaseq <- mouse_rnaseq[keep,]
rownames(mouse_rnaseq) <- make.names(gsub(".*_","", rownames(mouse_rnaseq)), unique = TRUE)
phenotype <- c(rep(0, 5), rep(1, 5))
tag <- c("CA", "ECM")

# Run Scissor ----

scissor_info <- Scissor(as.matrix(mouse_rnaseq), 
                  mc38_obj_sub,
                  phenotype, 
                  alpha = 0.05,
                  tag = tag,
                  family = "binomial",
                  Save_file = 'data/scissor.RData')

# Extract results ----

coefs <- scissor_info$Coefs
names(coefs) <- colnames(mc38_obj)
scissor_select <- rep(0, ncol(mc38_obj))
names(scissor_select) <- colnames(mc38_obj)
scissor_select[scissor_info$Scissor_pos] <- 1
scissor_select[scissor_info$Scissor_neg] <- 2

# Add results to scRNA-seq object ----

mc38_obj <- AddMetaData(mc38_obj, metadata = scissor_select, col.name = "scissor")
mc38_obj <- AddMetaData(mc38_obj, metadata = coefs, col.name = "coef2")

# Plot ----

pdf("plots/Fig_5L_Scissor.pdf")

DimPlot(mc38_obj, 
        group.by = 'scissor', 
        cols = c('grey','indianred1','royalblue'),
        na.value = "grey",
        pt.size = 0.6,
        order = c(2, 1)) + 
  theme(aspect.ratio = 1)

FeaturePlot(mc38_obj,
            features = "coef2",
            order = T, 
            pt.size = 0.6) + 
  scale_color_gradient2(low = "blue",
                        mid = "grey", 
                        high = "red",
                        na.value = "grey") +
  theme(aspect.ratio = 1) + 
  ggtitle("Coefficient")

dev.off()