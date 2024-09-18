#Processing the Wolffia microscopica dataset
library(dplyr)
library(Seurat)
library(patchwork)
set.seed = 123
Wmic_trans_merged <- Read10X_h5("filtered_feature_bc_matrix.h5")
Wmic_snRNA_merged <- CreateSeuratObject(counts = Wmic_trans_merged, project = "Wmic_snRNA_merged", min.cells = 3, min.features = 200)
Wmic_snRNA_merged@assays
DefaultAssay(Wmic_snRNA_merged) <- 'RNA'
Wmic_snRNA_merged@active.assay
Wmic_snRNA_merged[["percent.mt"]] <- PercentageFeatureSet(Wmic_snRNA_merged, pattern = "^MT-")
Wmic_snRNA_merged[["percent.cp"]] <- PercentageFeatureSet(Wmic_snRNA_merged, pattern = "^CP-")
head(Wmic_snRNA_merged@meta.data)
VlnPlot(Wmic_snRNA_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.cp"), ncol = 4)
FeatureScatter(Wmic_snRNA_merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
Wmic_snRNA_merged <- subset(Wmic_snRNA_merged, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA > 200 & percent.cp < 20 & percent.mt < 1)
Wmic_snRNA_merged <- NormalizeData(Wmic_snRNA_merged, normalization.method = "LogNormalize", scale.factor = 10000)
Wmic_snRNA_merged <- FindVariableFeatures(Wmic_snRNA_merged, selection.method = "vst", nfeatures = 2000)
all.genes.Wmic.merged <- rownames(Wmic_snRNA_merged)
Wmic_snRNA_merged <- ScaleData(Wmic_snRNA_merged, features = all.genes.Wmic.merged)
top10_Wmic_merged <- head(VariableFeatures(Wmic_snRNA_merged), 10)
plot_Wmic_merged <- VariableFeaturePlot(Wmic_snRNA_merged)
LabelPoints(plot = plot_Wmic_merged, points = top10_Wmic_merged, repel = TRUE)
Wmic_snRNA_merged <- RunPCA(Wmic_snRNA_merged, features = VariableFeatures(object = Wmic_snRNA_merged))
ElbowPlot(Wmic_snRNA_merged)
Wmic_snRNA_merged <- FindNeighbors(Wmic_snRNA_merged, dims = 1:16)
Wmic_snRNA_merged <- FindClusters(Wmic_snRNA_merged, resolution = 0.5)
Wmic_snRNA_merged <- RunUMAP(Wmic_snRNA_merged, dims = 1:16)
DimPlot(Wmic_snRNA_merged, reduction = "umap")

library(DoubletFinder)
homotypic.prop.Wmic.merged <- modelHomotypic(Wmic_snRNA_merged@meta.data$seurat_clusters)
doublet_ratio_Wmic_merged <- nrow(Wmic_snRNA_merged@meta.data) / 100000
nExp_poi_Wmic_merged <- round(doublet_ratio_Wmic_merged * nrow(Wmic_snRNA_merged@meta.data))
nExp_poi.adj.Wmic.merged <- round(nExp_poi_Wmic_merged * (1 - homotypic.prop.Wmic.merged))
sweep.res.list.Wmic.merged <- paramSweep(Wmic_snRNA_merged, PCs = 1:16, sct = FALSE)
sweep.stats.Wmic.merged <- summarizeSweep(sweep.res.list.Wmic.merged, GT = FALSE)
bcmvnM <- find.pK(sweep.stats.Wmic.merged)
pk_value_Wmic_merged <- as.numeric(as.character(bcmvnM[which.max(bcmvnM$BCmetric), "pK"]))
Wmic_snRNA_merged <- doubletFinder(Wmic_snRNA_merged, PCs = 1:16, pN = 0.25, pK = pk_value_Wmic_merged, nExp = nExp_poi.adj.Wmic.merged, reuse.pANN = FALSE, sct = FALSE)
Wmic_snRNA_merged <- subset(Wmic_snRNA_merged, subset = DF.classifications_0.25_0.05_634 == "Singlet")


#Replacing the Wolffia microscopica gene IDs with Arabidopsis thaliana orthologs
mapping_microscopica <- read.table("wmic_AT_all_unique_header.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
colnames(mapping_microscopica) <- c("original_gene_id", "arabidopsis_ortholog")
gene_mapping_microscopica <- setNames(mapping_microscopica$arabidopsis_ortholog, mapping_microscopica$original_gene_id)
current_gene_ids_microscopica <- rownames(Wmic_snRNA_merged)
new_gene_ids_microscopica <- gene_mapping_microscopica[current_gene_ids_microscopica]
valid_genes_microscopica <- !is.na(new_gene_ids_microscopica)
filtered_Wmic_snRNA_merged <- Wmic_snRNA_merged[valid_genes_microscopica, ]
new_gene_ids_microscopica <- new_gene_ids_microscopica[valid_genes_microscopica]
rownames(filtered_Wmic_snRNA_merged) <- new_gene_ids_microscopica
-----------
#Processing the Wolffia australiana dataset (Denyer et al., 2024)
library(dplyr)
library(Seurat)
library(patchwork)
set.seed = 123
expression_matrix <- read.table("MERGED_expression-raw.txt")
expression_matrix <- as.matrix(expression_matrix)
Waus_snRNA_merged <- CreateSeuratObject(counts = expression_matrix, project = "Waus_snRNA_merged")
Waus_snRNA_merged@assays
DefaultAssay(Waus_snRNA_merged) <- 'RNA'
Waus_snRNA_merged@active.assay
Waus_snRNA_merged <- NormalizeData(Waus_snRNA_merged, normalization.method = "LogNormalize", scale.factor = 10000)
Waus_snRNA_merged <- FindVariableFeatures(Waus_snRNA_merged, selection.method = "vst", nfeatures = 2000)
all.genes.Waus.merged <- rownames(Waus_snRNA_merged)
Waus_snRNA_merged <- ScaleData(Waus_snRNA_merged, features = all.genes.Waus.merged)
top10_Waus_merged <- head(VariableFeatures(Waus_snRNA_merged), 10)
plot_Waus_merged <- VariableFeaturePlot(Waus_snRNA_merged)
LabelPoints(plot = plot_Waus_merged, points = top10_Waus_merged, repel = TRUE)
Waus_snRNA_merged <- RunPCA(Waus_snRNA_merged, features = VariableFeatures(object = Waus_snRNA_merged))
ElbowPlot(Waus_snRNA_merged)
Waus_snRNA_merged <- FindNeighbors(Waus_snRNA_merged, dims = 1:15)
Waus_snRNA_merged <- FindClusters(Waus_snRNA_merged, resolution = 0.5)
Waus_snRNA_merged <- RunUMAP(Waus_snRNA_merged, dims = 1:15)
DimPlot(Waus_snRNA_merged, reduction = "umap")

library(DoubletFinder)
homotypic.prop.Waus.merged <- modelHomotypic(Waus_snRNA_merged@meta.data$seurat_clusters)
doublet_ratio_Waus_merged <- nrow(Waus_snRNA_merged@meta.data) / 100000
nExp_poi_Waus_merged <- round(doublet_ratio_Waus_merged * nrow(Waus_snRNA_merged@meta.data))
nExp_poi.adj.Waus.merged <- round(nExp_poi_Waus_merged * (1 - homotypic.prop.Waus.merged))
sweep.res.list.Waus.merged <- paramSweep(Waus_snRNA_merged, PCs = 1:15, sct = FALSE)
sweep.stats.Waus.merged <- summarizeSweep(sweep.res.list.Waus.merged, GT = FALSE)
bcmvnMG <- find.pK(sweep.stats.Waus.merged)
pk_value_Waus_merged <- as.numeric(as.character(bcmvnMG[which.max(bcmvnMG$BCmetric), "pK"]))
Waus_snRNA_merged <- doubletFinder(Waus_snRNA_merged, PCs = 1:15, pN = 0.25, pK = pk_value_Waus_merged, nExp = nExp_poi.adj.Waus.merged, reuse.pANN = FALSE, sct = FALSE)
Waus_snRNA_merged <- subset(Waus_snRNA_merged, subset = DF.classifications_0.25_0.09_281 == "Singlet")

#Replacing the Wolffia australiana gene IDs with Arabidopsis thaliana orthologs
mapping_australiana <- read.table("waus_AT_all_unique_header.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
colnames(mapping_australiana) <- c("original_gene_id", "arabidopsis_ortholog")
gene_mapping_australiana <- setNames(mapping_australiana$arabidopsis_ortholog, mapping_australiana$original_gene_id)
current_gene_ids_australiana <- rownames(Waus_snRNA_merged)
new_gene_ids_australiana <- gene_mapping_australiana[current_gene_ids_australiana]
valid_genes_australiana <- !is.na(new_gene_ids_australiana)
filtered_Waus_snRNA_merged <- Waus_snRNA_merged[valid_genes_australiana, ]
new_gene_ids_australiana <- new_gene_ids_australiana[valid_genes_australiana]
rownames(filtered_Waus_snRNA_merged) <- new_gene_ids_australiana

---------------

#Integrated data set of W. microscopica and W. australiana
wmic_waus <- merge(filtered_Wmic_snRNA_merged, 
              y = filtered_Waus_snRNA_merged, 
              add.cell.ids = c("Wmic", "Waus"), 
              project = "Wolffia_merged")
wmic_waus@assays
DefaultAssay(wmic_waus) <- 'RNA'
wmic_waus@active.assay
wmic_waus <- NormalizeData(wmic_waus, normalization.method = "LogNormalize", scale.factor = 10000)
wmic_waus <- FindVariableFeatures(wmic_waus, selection.method = "vst", nfeatures = 2000)
wmic_waus.genes <- rownames(wmic_waus)
wmic_waus <- ScaleData(wmic_waus, features = wmic_waus.genes)
wmic_waus <- RunPCA(wmic_waus, features = VariableFeatures(object = wmic_waus))
ElbowPlot(wmic_waus)
wmic_waus <- FindNeighbors(wmic_waus, dims = 1:15)
wmic_waus <- FindClusters(wmic_waus, resolution = 0.5)
wmic_waus <- RunUMAP(wmic_waus, dims = 1:15)
DimPlot(wmic_waus, reduction = "umap", group.by="orig.ident")

#Batch correction
merged_wolffia_list <- lapply(c(filtered_Wmic_snRNA_merged, filtered_Waus_snRNA_merged1), function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features_merged_wolffia <- SelectIntegrationFeatures(object.list = merged_wolffia_list, nfeatures=3000)
merged_wolffia_anchors <- FindIntegrationAnchors(object.list = merged_wolffia_list, anchor.features = features_merged_wolffia)
merged_wolffia_combined <- IntegrateData(anchorset = merged_wolffia_anchors)
merged_wolffia_combined@assays
DefaultAssay(merged_wolffia_combined) <- "integrated"
merged_wolffia_combined <- ScaleData(merged_wolffia_combined, verbose = FALSE)
merged_wolffia_combined <- RunPCA(merged_wolffia_combined, npcs = 50, verbose = FALSE)
merged_wolffia_combined <- RunUMAP(merged_wolffia_combined, reduction = "pca", dims = 1:50)
merged_wolffia_combined <- FindNeighbors(merged_wolffia_combined, dims = 1:50)
merged_wolffia_combined <- FindClusters(merged_wolffia_combined, resolution = 0.5)
DimPlot(merged_wolffia_combined, reduction = "umap", group.by = "orig.ident")
DimPlot(merged_wolffia_combined, reduction = "umap", label = TRUE)
DimPlot(merged_wolffia_combined, reduction = "umap", split.by = "orig.ident")
markers_combined <- FindAllMarkers(merged_wolffia_combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
