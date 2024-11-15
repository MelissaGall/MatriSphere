# Generate scRNA-seq plots

# Libraries ----

library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(openxlsx)
library(ComplexHeatmap)
library(circlize)

# Load data ----

mc38_obj <- LoadH5Seurat("data/Seurat/WT_MC38_filtered_annotated.h5seurat") # scRNAseq data
matrisome <- read.xlsx("data/MJB_CRC Matrisome Analysis.xlsx") # Matrisome gene sets
matrisome_2 <- read.xlsx("data/Editable Naba 2012 Matrisome List.xlsx", sheet = 2, startRow = 2)
  
# Dimplots ----

pdf("plots/Fig_1K_MC38_WT_Tumor_annotated.pdf", width = 10)

DimPlot(mc38_obj, group.by = "seurat_clusters", label = T) + theme(aspect.ratio = 1, legend.position = "None")
DimPlot(mc38_obj, group.by = "Subtype", label = T) + theme(aspect.ratio = 1, legend.position = "None")
DimPlot(mc38_obj, group.by = "Type") + 
  theme(aspect.ratio = 1) +
  scale_color_manual(values = c("MC38" = "#619cff", 
                                "Immune cells" = "#f8766d",
                                "Stromal cells" = "#00b935"),
                     labels = c("MC38" = "Tumor cells", 
                                "Immune cells" = "Immune system",
                                "Stromal cells" = "Stromal cells")) 
dev.off()

# FeaturePlots of markers of interest ----

markers <- read.table("data/markers_to_investigate.csv")
markers <- markers$V1

for (marker in markers){
  jpeg(paste0("plots/Fig_S4_S5_WT_MC38_markers_/", marker, ".jpeg"))
  print(FeaturePlot(mc38_obj, features = marker, cols = c( "gray","red")))
  dev.off()
}

# Z-score dotplot ----

## Define function to calculate zscore ----

calculate_zscore  <- function(object, gene_list, column_name){
  
  ### Get expression data
  
  expression_data <- GetAssayData(object, assay = "RNA", slot = "scale.data")
  expression_data <- data.frame(t(as.matrix(expression_data[rownames(expression_data) %in% gene_list,]))) # Keep only data for genes in gene set
  expression_data <- scale(expression_data, center = T)
  
  ### Computre z-score for each cell
  
  mean_exp <- mean(expression_data, na.rm = T)
  sd_exp <- sd(expression_data, na.rm = T) 
  gene_exp_mean <- as.vector(apply(expression_data, 1, mean, na.rm = T)) # Calculate mean for each row (each cell)
  expression_data <- data.frame(cbind(expression_data, gene_exp_mean)) # Add the means in a new column
  expression_data$zscore <- (expression_data$gene_exp_mean - mean_exp) / sd_exp # Compute z-score
  zcore_df <- data.frame(expression_data[, "zscore"])
  rownames(zcore_df) <- rownames(expression_data)
  colnames(zcore_df) <- column_name # Add meaninful column name (= gene set name)
   
  ### Add to metadata
  
  object <- AddMetaData(object, zcore_df)
  
  return(object)

}

## Get genes list  ----

collagens <- matrisome$Gene.symbol[matrisome$Category == "Collagens"]
glycoproteins <- matrisome$Gene.symbol[matrisome$Category == "ECM Glycoproteins"]
proteoglycans <- matrisome$Gene.symbol[matrisome$Category == "Proteoglycans"]
ECMAP <- matrisome$Gene.symbol[matrisome$Category == "ECM-affiliated Proteins"]
ECM_reg <- matrisome$Gene.symbol[matrisome$Category == "ECM Regulators"]
SF <- matrisome$Gene.symbol[matrisome$Category == "Secreted Factors"]

## Calculate zscore for each gene set  ----

mc38_obj <- calculateZscore(mc38_obj, collagens, "Collagens")
mc38_obj <- calculateZscore(mc38_obj, glycoproteins, "ECM.Glycoproteins")
mc38_obj <- calculateZscore(mc38_obj, proteoglycans, "Proteoglycans")
mc38_obj <- calculateZscore(mc38_obj, ECMAP, "ECM.affiliated.Proteins")
mc38_obj <- calculateZscore(mc38_obj, ECM_reg, "ECM.Regulators")
mc38_obj <- calculateZscore(mc38_obj, SF, "Secreted.Factors")

## Generate DotPlots  ----

features <- c("Collagens", "ECM.Glycoproteins", "Proteoglycans", "ECM.affiliated.Proteins", "ECM.Regulators", "Secreted.Factors")

dplot_1 <- DotPlot(mc38_obj, features = features, group.by = "Type", scale = F) +
  scale_colour_gradient2(low = "#0492C2", mid = "#FAFAD2", high = "#D21404", midpoint = 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  coord_flip()

dplot_2 <- DotPlot(mc38_obj, features = features, group.by = "Subtype", scale = F) +
  scale_colour_gradient2(low = "#0492C2", mid = "#FAFAD2", high = "#D21404") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  coord_flip()

## Save 

pdf("plots/F1L_gene_sets_dotplots.pdf", width = 8)
print(dplot_1)
print(dplot_2)
dev.off()

# Core matrisome gene sets heatmap

## Create gene list

core <- matrisome[matrisome$Division == "Core matrisome",]
genes <- core$Gene.symbol

## Extract average expression for each cell type

counts <- AverageExpression(mc38_obj, group.by = "Type", assays = "RNA", slot = "counts")
counts <- counts$RNA
counts <- as.matrix(counts[rownames(counts) %in% genes, ])
counts <- rescale(counts)

## Extract meta data ----

matrisome_sub <- matrisome_2[matrisome_2$Gene.symbol %in% rownames(counts),]
matrisome_sub <- matrisome_sub[match(rownames(counts), matrisome_sub$Gene.symbol), ]
hp_meta <- data.frame("Genes" = rownames(counts),
                      "Category" = matrisome_sub$Category)

## Prepare annotation and color scale ----

split <- factor(hp_meta$Category)

col_fun <- colorRamp2(c(0, 0.5, 1), c("#040380", "white", "#d3323d"))

row_ha <- rowAnnotation(Category = hp_meta$Category,
                       show_annotation_name = F,
                       col = list(Category = c("ECM Glycoproteins" ="forestgreen",
                                               "Collagens" ="orange", 
                                               "Proteoglycans" = "cyan")))
## Plot ----

pdf("plots/Fig_S3_MJB_Core_matrisome_raw_exp_3columns.pdf", height = 12.5)

Heatmap(counts, 
        cluster_columns = FALSE, 
        name = "Average expression", 
        col = col_fun,
        row_split  = split, 
        show_parent_dend_line = FALSE, 
        cluster_row_slices  = TRUE, 
        row_names_side = "left",
        right_annotation = row_ha,
        column_names_centered = F,
        show_column_dend  = F,
        column_names_rot  = 45,
        show_row_dend = F,
        row_title_gp = gpar(fontsize = 0),
        show_row_names = T,
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8),
        width = ncol(counts)*unit(6, "mm"), 
        height =  nrow(counts) * unit(2.5, "mm"),
        row_title  = NULL,
        column_title = "Core matrisome genes - MJB list"
)

dev.off()
