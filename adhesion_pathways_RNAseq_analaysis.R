# Mouse RNA-seq adhesion pathway analysis

# Libraries ----

library(fgsea)
library(msigdbr)
library(dplyr)
library(stringr)
library(DESeq2)
library(ggplot2)
library(ggVennDiagram)
library(ComplexHeatmap)

# Prepare pathways (extract adhesion-related pathways) ----

## Mouse ----

all_gene_sets <- msigdbr(species = "Mus musculus")
all_gene_sets <- data.frame(all_gene_sets)
adhesion_sets <- all_gene_sets[grepl("adhesion", all_gene_sets$gs_name, ignore.case = TRUE),]
msigdbr_list_mouse <- split(x = adhesion_sets$gene_symbol, f = adhesion_sets$gs_name)

## Human ----

all_gene_sets <- msigdbr(species = "Homo sapiens")
all_gene_sets <- data.frame(all_gene_sets)
adhesion_sets <- all_gene_sets[grepl("adhesion", all_gene_sets$gs_name, ignore.case = TRUE),]
msigdbr_list_human <- split(x = adhesion_sets$gene_symbol, f = adhesion_sets$gs_name)

# Load data ----

mouse_data <- read.table("data/ccbr1223_mouse_rawcounts.csv", header = TRUE, row.names = 1, sep = ",")
mouse_data[] <- sapply(mouse_data, as.integer)

human_data <- read.table("data/ccbr1223_human_rawcounts.csv", header = TRUE, row.names = 1, sep = ",")
human_data[] <- sapply(human_data, as.integer)

# Prepare data ----

## For MC38  ----

mc38_data <- mouse_data %>% select(contains("MC38"))
mc38_data <- mc38_data[, !(colnames(mc38_data) %in% c("CA_MC38_6", "ECM_MC38_2"))] # Remove outliers
rownames(mc38_data) <- make.names(gsub(".*_","", rownames(mc38_data)), unique = TRUE)
mc38_data <- mc38_data[order(rownames(mc38_data)),]

metadata_mc38 <- data.frame(row.names = colnames(mc38_data), 
                            "Sample" = colnames(mc38_data),
                            "Group" = c(rep("CA", 5), rep("ECM", 5)))

dds_mc38 <- DESeqDataSetFromMatrix(countData = mc38_data,
                                   colData = metadata_mc38,
                                   design = ~ Group)

keep <- rowSums(counts(dds_mc38) >= 10) >= 3 # Remove low count rows
dds_mc38 <- dds_mc38[keep,]
dds_mc38 <- DESeq(dds_mc38)

results_mc38 <- results(dds_mc38, contrast = c("Group", "ECM", "CA"), alpha =  0.05)
results_mc38 <- results_mc38[order(results_mc38$log2FoldChange, decreasing = TRUE), ]

ranks_mc38 <- results_mc38$log2FoldChange
names(ranks_mc38) <- rownames(results_mc38)
fgsea_results_mc38 <- fgsea(pathways = msigdbr_list_mouse, stats = ranks_mc38, minSize = 15, maxSize = 500)
fgsea_sign_mc38 <- fgsea_results_mc38[fgsea_results_mc38$padj < 0.05,]

## For CT26 ----

ct26_data <- mouse_data %>% select(contains("CT26"))
rownames(ct26_data) <- make.names(gsub(".*_","", rownames(ct26_data)), unique = TRUE)
ct26_data <- ct26_data[order(rownames(ct26_data)),]

metadata_ct26 <- data.frame(row.names = colnames(ct26_data), 
                            "Sample" = colnames(ct26_data),
                            "Group" = c(rep("CA", 6), rep("ECM", 6)))

dds_ct26 <- DESeqDataSetFromMatrix(countData = ct26_data,
                                   colData = metadata_ct26,
                                   design = ~ Group)

keep <- rowSums(counts(dds_ct26) >= 10) >= 3
dds_ct26 <- dds_ct26[keep,]
dds_ct26 <- DESeq(dds_ct26)

results_ct26 <- results(dds_ct26, contrast = c("Group", "ECM", "CA"), alpha =  0.05)
results_ct26 <- results_ct26[order(results_ct26$log2FoldChange, decreasing = TRUE), ]

ranks_ct26 <- results_ct26$log2FoldChange
names(ranks_ct26) <- rownames(results_ct26)
fgsea_results_ct26 <- fgsea(pathways = msigdbr_list_mouse, stats = ranks_ct26, minSize = 15, maxSize = 500)
fgsea_sign_ct26 <- fgsea_results_ct26[fgsea_results_ct26$padj < 0.05,]

## For HT29 ----

ht29_data <- human_data %>% select(contains("HT29"))
rownames(ht29_data)<- make.names(gsub(".*_","", rownames(ht29_data)), unique = TRUE)
ht29_data <- ht29_data[order(rownames(ht29_data)),]

metadata_ht29 <- data.frame(row.names = colnames(ht29_data), 
                            "Sample" = colnames(ht29_data),
                            "Group" = c(rep("ECM", 5), rep("CA", 6), "ECM"))

dds_ht29 <- DESeqDataSetFromMatrix(countData = ht29_data,
                                   colData = metadata_ht29,
                                   design = ~ Group)

keep <- rowSums(counts(dds_HT29) >= 10) >= 3
dds_ht29 <- dds_ht29[keep,]
dds_ht29 <- DESeq(dds_ht29)

results_ht29 <- results(dds_ht29, contrast = c("Group", "ECM", "CA"), alpha =  0.05)
results_ht29 <- results_ht29[order(results_ht29$log2FoldChange, decreasing = TRUE), ]

ranks_ht29 <- results_ht29$log2FoldChange
names(ranks_ht29) <- rownames(results_ht29)
fgsea_results_ht29 <- fgsea(pathways = msigdbr_list_human, stats = ranks_ht29, minSize = 15, maxSize = 500)
fgsea_results_ht29 <- fgsea_results_ht29[fgsea_results_ht29$padj < 0.05,]

# Leading edge Venn ----

## Get intersection ----

MC38_le <- unique(unlist(fgsea_sign_mc38$leadingEdge))
CT26_le <- unique(unlist(fgsea_sign_ct26$leadingEdge))
common_le <- intersect(MC38_le, CT26_le)
MC38_le <- toupper(MC38_le)
CT26_le <- toupper(CT26_le)
HT29_le <- toupper(unique(unlist(fgsea_sign_ht29$leadingEdge)))

## Plot ----

pdf("Plots/Fig_S24a_adhesion_pathways_leading_edge_Venn.pdf")

ggVennDiagram(x = list("MC38" = MC38_le, "CT26" = CT26_le, "HT29" = HT29_le), color = 1, lwd = 0.7) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")

dev.off()

# Heatmaps ----

## Get counts ----

### MC38 ----

norm_counts_MC38 <- counts(dds_mc38, normalized = TRUE)
selected_genes_MC38 <- norm_counts_MC38[common_le, ]
scaled_counts_MC38 <- t(scale(t(selected_genes_MC38)))

### CT26 ----

norm_counts_CT26 <- counts(dds_ct26, normalized = TRUE)
selected_genes_CT26 <- norm_counts_CT26[common_le, ]
scaled_counts_CT26 <- t(scale(t(selected_genes_CT26)))

### HT29 ----

norm_counts_HT29 <- counts(dds_ht29, normalized = TRUE)
HT29_genes <- toupper(common_le)
HT29_genes <- HT29_genes[HT29_genes %in% rownames(norm_counts_HT29)]
selected_genes_HT29 <- norm_counts_HT29[HT29_genes, ]
scaled_counts_HT29 <- t(scale(t(selected_genes_HT29)))

## Plot ----

pdf("Plots/Fig_S24_to_S26_adhesion_pathways_leading_edge_heatmaps.pdf")

Heatmap(scaled_counts_MC38,
        name = "Expression",
        show_row_names = TRUE,
        show_column_names = TRUE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(title = "Expression", 
                                    title_position = "topcenter"),
        column_title = "MC38")

Heatmap(scaled_counts_CT26,
        name = "Expression",
        show_row_names = TRUE,
        show_column_names = TRUE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(title = "Expression", 
                                    title_position = "topcenter"),
        column_title = "CT26")

Heatmap(scaled_counts_HT29,
        name = "Expression",
        show_row_names = TRUE,
        show_column_names = TRUE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(title = "Expression", 
                                    title_position = "topcenter"),
        column_title = "HT29")

dev.off()