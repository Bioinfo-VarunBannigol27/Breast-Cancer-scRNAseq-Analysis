#!/usr/bin/env Rscript
# ============================================================
# Single-Cell RNA-seq Analysis of Breast Cancer Tumor Tissue
# ============================================================
#
# Dataset : Wu et al. 2021, Nature Genetics (GSE176078)
# Platform: 10x Genomics Chromium
# Cells   : 100,064 raw -> 98,593 after QC
# Genes   : 27,719
#
# Pipeline: QC -> Normalization -> Feature Selection -> PCA ->
#           UMAP -> Clustering -> Marker ID -> Cell Type Annotation
#
# Author  : [your name]
# Tools   : R 4.5, Seurat v5, dplyr, ggplot2, presto
# ============================================================

# ---- 0. SETUP ----------------------------------------------

dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)

library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)
library(presto)   # speeds up FindAllMarkers via Wilcoxon rank-sum test

options(timeout = 3600)

# ---- 1. DOWNLOAD & LOAD RAW DATA ----------------------------
# Source: GEO accession GSE176078 (Wu et al. 2021, Nature Genetics)

download.file(
  url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176078/suppl/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz",
  destfile = "data/BRCA_scRNA.tar.gz",
  mode = "wb"
)
untar("data/BRCA_scRNA.tar.gz", exdir = "data/")

matrix_data <- readMM("data/Wu_etal_2021_BRCA_scRNASeq/count_matrix_sparse.mtx")
genes       <- read.delim("data/Wu_etal_2021_BRCA_scRNASeq/count_matrix_genes.tsv",
                           header = FALSE, stringsAsFactors = FALSE)
barcodes    <- read.delim("data/Wu_etal_2021_BRCA_scRNASeq/count_matrix_barcodes.tsv",
                           header = FALSE, stringsAsFactors = FALSE)

rownames(matrix_data) <- genes$V1
colnames(matrix_data) <- barcodes$V1

# ---- 2. CREATE SEURAT OBJECT --------------------------------

brca <- CreateSeuratObject(
  counts       = matrix_data,
  project      = "BreastCancer",
  min.cells    = 3,     # gene must appear in >= 3 cells
  min.features = 200    # cell must express >= 200 genes
)

print(brca)
# An object of class Seurat
# 27719 features across 100064 samples within 1 assay

# Attach published cell-type metadata (for downstream validation only —
# annotation in this project was derived independently from clustering + markers)
metadata_df <- read.csv("data/Wu_etal_2021_BRCA_scRNASeq/metadata.csv", row.names = 1)
brca <- AddMetaData(brca, metadata = metadata_df)

# ---- 3. QUALITY CONTROL -------------------------------------

brca[["percent.mt"]] <- PercentageFeatureSet(brca, pattern = "^MT-")

# QC metric distributions before filtering
plot_qc_before <- VlnPlot(
  brca,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3, pt.size = 0.1
)
ggsave("plots/QC_before_filtering.png", plot_qc_before, width = 12, height = 6, dpi = 300)

# Relationship between metrics
plot1 <- FeatureScatter(brca, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend()
plot2 <- FeatureScatter(brca, feature1 = "nCount_RNA", feature2 = "percent.mt") + NoLegend()
ggsave("plots/QC_scatter_plots.png", plot1 + plot2, width = 12, height = 5, dpi = 300)

summary(brca$nFeature_RNA)
summary(brca$nCount_RNA)
summary(brca$percent.mt)

cat("Cells before QC filtering:", ncol(brca), "\n")

# Filtering thresholds (derived from the distributions above)
brca <- subset(
  brca,
  subset = nFeature_RNA > 200 &
           nFeature_RNA < 6000 &
           percent.mt   < 20
)

cat("Cells after QC filtering:", ncol(brca), "\n")
# 100,064 -> 98,593 cells retained

plot_qc_after <- VlnPlot(
  brca, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3, pt.size = 0.1
)
ggsave("plots/QC_after_filtering.png", plot_qc_after, width = 12, height = 5, dpi = 300)

# NOTE ON DOUBLET REMOVAL:
# DoubletFinder was deliberately skipped at this cell count (~100k).
# It requires constructing artificial doublets and re-running PCA/UMAP,
# which roughly doubles memory and compute load -- a documented
# trade-off for this single-sample, single-machine analysis.

# ---- 4. NORMALIZATION ---------------------------------------

brca <- NormalizeData(
  brca,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# ---- 5. FEATURE SELECTION (Highly Variable Genes) -----------

brca <- FindVariableFeatures(
  brca,
  selection.method = "vst",
  nfeatures = 2000
)

top10_variable <- head(VariableFeatures(brca), 10)
print(top10_variable)

plot_hvg <- VariableFeaturePlot(brca)
plot_hvg_labeled <- LabelPoints(plot = plot_hvg, points = top10_variable, repel = TRUE)
ggsave("plots/Highly_Variable_Genes.png", plot_hvg_labeled, width = 8, height = 6, dpi = 300)

# ---- 6. SCALING -----------------------------------------------

brca <- ScaleData(
  brca,
  features = VariableFeatures(brca),
  vars.to.regress = c("nCount_RNA", "percent.mt")  # remove technical variation
)

# ---- 7. PCA -----------------------------------------------------

brca <- RunPCA(
  brca,
  features = VariableFeatures(brca),
  npcs = 50,
  verbose = TRUE
)

print(brca[["pca"]], dims = 1:5, nfeatures = 5)

pca_plot <- DimPlot(brca, reduction = "pca") + NoLegend()
ggsave("plots/PCA_plot.png", pca_plot, width = 7, height = 6, dpi = 300)

elbow_plot <- ElbowPlot(brca, ndims = 50)
ggsave("plots/PCA_elbow_plot.png", elbow_plot, width = 6, height = 5, dpi = 300)

DimHeatmap(brca, dims = 1:9, cells = 500, balanced = TRUE)

n_pcs <- 30   # chosen from the elbow plot inflection point

# ---- 8. UMAP ----------------------------------------------------

brca <- RunUMAP(brca, dims = 1:n_pcs, verbose = TRUE)

umap_basic <- DimPlot(brca, reduction = "umap")
ggsave("plots/UMAP_basic.png", umap_basic, width = 8, height = 6, dpi = 300)

# ---- 9. CLUSTERING (Louvain, graph-based) -----------------------

brca <- FindNeighbors(brca, dims = 1:n_pcs, k.param = 20)
brca <- FindClusters(brca, resolution = 0.5, algorithm = 1)

n_clusters <- length(unique(brca$seurat_clusters))
cat("Number of clusters found:", n_clusters, "\n")  # 36
print(table(brca$seurat_clusters))

umap_clusters <- DimPlot(
  brca, reduction = "umap", group.by = "seurat_clusters",
  label = TRUE, label.size = 5, repel = TRUE
) + NoLegend()
ggsave("plots/UMAP_clusters_numbered.png", umap_clusters, width = 9, height = 8, dpi = 300)

# ---- 10. MARKER GENE IDENTIFICATION -----------------------------
# Wilcoxon rank-sum test (recommended by Su et al. 2022 review;
# low false-positive rate, robust at large sample sizes)

brca <- JoinLayers(brca)                  # required for Seurat v5 multi-layer objects
brca_small <- subset(brca, downsample = 200)  # downsample per cluster for tractable DE testing

all_markers <- FindAllMarkers(
  brca_small,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.30,
  test.use = "wilcox",
  verbose = TRUE
)

write.csv(all_markers, "results/all_cluster_markers.csv")

top5_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
print(top5_markers)

# ---- 11. MARKER VALIDATION (FeaturePlots + DotPlot) -------------

cancer_markers_plot <- FeaturePlot(
  brca, features = c("EPCAM", "KRT8", "KRT18"),
  ncol = 3, cols = c("lightgrey", "red")
)
ggsave("plots/Cancer_cell_markers.png", cancer_markers_plot, width = 15, height = 5, dpi = 300)

immune_markers_plot <- FeaturePlot(
  brca, features = c("CD3D", "CD79A", "CD68", "GNLY"),
  ncol = 4, cols = c("lightgrey", "red")
)
ggsave("plots/Immune_markers.png", immune_markers_plot, width = 20, height = 5, dpi = 300)

stromal_markers_plot <- FeaturePlot(
  brca, features = c("COL1A1", "PECAM1", "ACTA2"),
  ncol = 3, cols = c("lightgrey", "red")
)
ggsave("plots/Stromal_markers.png", stromal_markers_plot, width = 15, height = 5, dpi = 300)

dot_plot <- DotPlot(
  brca,
  features = c("EPCAM", "KRT8", "CD3D", "CD8A", "CD79A", "MS4A1",
               "CD68", "CD163", "GNLY", "NKG7", "COL1A1", "FAP",
               "PECAM1", "VWF")
) + RotatedAxis()
ggsave("plots/Dotplot_markers.png", dot_plot, width = 14, height = 8, dpi = 300)

# ---- 12. CELL TYPE ANNOTATION -----------------------------------
# Clusters merged into 9 major populations based on canonical
# marker expression (Step 11) cross-checked against literature.

new_cluster_ids <- c(
  "0"="T-cells","1"="T-cells","8"="T-cells","28"="T-cells","34"="T-cells","17"="T-cells","19"="T-cells",
  "3"="Myeloid","14"="Myeloid","30"="Myeloid","32"="Myeloid",
  "9"="B-cells","31"="B-cells","10"="B-cells","25"="B-cells",
  "24"="Plasmablasts",
  "4"="CAFs","5"="CAFs",
  "29"="PVL","2"="PVL",
  "13"="Endothelial","18"="Endothelial",
  "15"="Normal Epithelial","22"="Normal Epithelial",
  "6"="Cancer Epithelial","7"="Cancer Epithelial","11"="Cancer Epithelial",
  "12"="Cancer Epithelial","16"="Cancer Epithelial","20"="Cancer Epithelial",
  "21"="Cancer Epithelial","23"="Cancer Epithelial","26"="Cancer Epithelial",
  "27"="Cancer Epithelial","33"="Cancer Epithelial","35"="Cancer Epithelial"
)

brca <- RenameIdents(brca, new_cluster_ids)
brca$celltype_major_custom <- Idents(brca)

umap_annotated <- DimPlot(
  brca, reduction = "umap", label = TRUE, label.size = 4,
  repel = TRUE, pt.size = 0.5
) + NoLegend() +
  ggtitle("Breast Cancer - Cell Types") +
  theme(plot.title = element_text(hjust = 0.5, size = 16))

ggsave("plots/UMAP_annotated_celltypes.png", umap_annotated, width = 12, height = 8, dpi = 300)

# ---- 13. SAVE FINAL OBJECT --------------------------------------

saveRDS(brca, file = "results/brca_seurat_annotated.rds")

cat("\n=== Analysis complete ===\n")
cat("Final cell count:", ncol(brca), "\n")
cat("Clusters identified:", n_clusters, "\n")
cat("Annotated cell types:", length(unique(brca$celltype_major_custom)), "\n")
