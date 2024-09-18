library(dplyr)
library(Seurat)
library(patchwork)
set.seed = 123
Wmic_trans <- Read10X_h5("filtered_feature_bc_matrix.h5")
Wmic_snRNA <- CreateSeuratObject(counts = Wmic_trans, project = "Wmic_snRNA", min.cells = 3, min.features = 200)
-------------------
#Getting the stats
raw_counts_Wmic <- GetAssayData(Wmic_snRNA, assay = "RNA", slot = "counts")
unique_genes_Wmic <- rowSums(raw_counts_Wmic > 0)
number_of_unique_genes_Wmic <- sum(unique_genes_Wmic > 0)
print(paste("Number of unique genes identified:", number_of_unique_genes_Wmic))
umis_per_cell_Wmic <- Wmic_snRNA@meta.data$nCount_RNA
summary(umis_per_cell_Wmic)
genes_per_cell_Wmic <- Wmic_snRNA@meta.data$nFeature_RNA
summary(genes_per_cell_Wmic)
-------------------
Wmic_snRNA@assays
DefaultAssay(Wmic_snRNA) <- 'RNA'
Wmic_snRNA@active.assay
Wmic_snRNA[["percent.mt"]] <- PercentageFeatureSet(Wmic_snRNA, pattern = "^MT-")
Wmic_snRNA[["percent.cp"]] <- PercentageFeatureSet(Wmic_snRNA, pattern = "^CP-")
head(Wmic_snRNA@meta.data)
VlnPlot(Wmic_snRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.cp"), ncol = 4)
FeatureScatter(Wmic_snRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
Wmic_snRNA <- subset(Wmic_snRNA, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA > 200 & percent.cp < 20 & percent.mt < 1)

Wmic_snRNA <- NormalizeData(Wmic_snRNA, normalization.method = "LogNormalize", scale.factor = 10000)
Wmic_snRNA <- FindVariableFeatures(Wmic_snRNA, selection.method = "vst", nfeatures = 2000)
all.genes.Wmic <- rownames(Wmic_snRNA)
Wmic_snRNA <- ScaleData(Wmic_snRNA, features = all.genes.Wmic)
top10_Wmic <- head(VariableFeatures(Wmic_snRNA), 10)
plot_Wmic <- VariableFeaturePlot(Wmic_snRNA)
LabelPoints(plot = plot_Wmic, points = top10_Wmic, repel = TRUE)

Wmic_snRNA <- RunPCA(Wmic_snRNA, features = VariableFeatures(object = Wmic_snRNA))
ElbowPlot(Wmic_snRNA)
Wmic_snRNA <- FindNeighbors(Wmic_snRNA, dims = 1:16)
Wmic_snRNA <- FindClusters(Wmic_snRNA, resolution = 0.5)
Wmic_snRNA <- RunUMAP(Wmic_snRNA, dims = 1:16)
DimPlot(Wmic_snRNA, reduction = "umap")
--------------------
#Doubletfinder
library(DoubletFinder)
homotypic.prop.Wmic <- modelHomotypic(Wmic_snRNA@meta.data$seurat_clusters)
doublet_ratio_Wmic <- nrow(Wmic_snRNA@meta.data) / 100000
nExp_poi_Wmic <- round(doublet_ratio_Wmic * nrow(Wmic_snRNA@meta.data))
nExp_poi.adj.Wmic <- round(nExp_poi_Wmic * (1 - homotypic.prop.Wmic))
sweep.res.list.Wmic <- paramSweep(Wmic_snRNA, PCs = 1:16, sct = FALSE)
sweep.stats.Wmic <- summarizeSweep(sweep.res.list.Wmic, GT = FALSE)
bcmvn0 <- find.pK(sweep.stats.Wmic)
pk_value_Wmic <- as.numeric(as.character(bcmvn0[which.max(bcmvn0$BCmetric), "pK"]))
Wmic_snRNA <- doubletFinder(Wmic_snRNA, PCs = 1:16, pN = 0.25, pK = pk_value_Wmic, nExp = nExp_poi.adj.Wmic, reuse.pANN = FALSE, sct = FALSE)
head(Wmic_snRNA@meta.data)
Wmic_snRNA <- subset(Wmic_snRNA, subset = DF.classifications_0.25_0.05_634 == "Singlet")
-------------------
Wmic_snRNA <- NormalizeData(Wmic_snRNA, normalization.method = "LogNormalize", scale.factor = 10000)
Wmic_snRNA <- FindVariableFeatures(Wmic_snRNA, selection.method = "vst", nfeatures = 2000)
all.genes.Wmic.df <- rownames(Wmic_snRNA)
Wmic_snRNA <- ScaleData(Wmic_snRNA, features = all.genes.Wmic.df)
Wmic_snRNA <- RunPCA(Wmic_snRNA, features = VariableFeatures(object = Wmic_snRNA))
ElbowPlot(Wmic_snRNA)
Wmic_snRNA <- FindNeighbors(Wmic_snRNA, dims = 1:16)
Wmic_snRNA <- FindClusters(Wmic_snRNA, resolution = 0.5)
Wmic_snRNA <- RunUMAP(Wmic_snRNA, dims = 1:16)
DimPlot(Wmic_snRNA, reduction = "umap")
markers_wolffia_snRNA_final <- FindAllMarkers(Wmic_snRNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
-------------------
#Top markers
markers_wolffia_snRNA_final %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)
-------------------
#No. of nuclei in each cluster
table(Cluster = Wmic_snRNA@active.ident, Batch = Wmic_snRNA$orig.ident)

