library(Seurat)
library(Signac)
library(Matrix)
library(dplyr)

wolffia_multiome_data <- Read10X_h5("filtered_feature_bc_matrix.h5")
-----------------------
#Create a GRanges object
library(GenomicRanges)
library(data.table)

gtf_file_multi <- "Wmic_reference_genome1.gtf"
gtf_data_multi <- fread(gtf_file_multi, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(gtf_data_multi) <- c("seqnames", "source", "type", "start", "end", "score", "strand", "frame", "attributes")
extract_attributes <- function(attribute_string, attribute_name) {
  sub(paste0(".*", attribute_name, " \"([^\"]+)\".*"), "\\1", attribute_string)
}
gtf_data_multi[, transcript_id := extract_attributes(attributes, "transcript_id")]
gtf_data_multi[, gene_name := extract_attributes(attributes, "gene_name")]
gtf_data_multi[, gene_id := extract_attributes(attributes, "gene_id")]
gtf_data_multi[, gene_biotype := extract_attributes(attributes, "gene_biotype")]
gr_multi <- GRanges(seqnames = gtf_data_multi$seqnames,
               ranges = IRanges(start = gtf_data_multi$start, end = gtf_data_multi$end),
               strand = gtf_data_multi$strand,
               tx_id = gtf_data_multi$transcript_id,
               gene_name = gtf_data_multi$gene_name,
               gene_id = gtf_data_multi$gene_id,
               gene_biotype = gtf_data_multi$gene_biotype,
               type = gtf_data_multi$type)
print(gr_multi)
-----------------------
seurat_wolffia_multi <- CreateSeuratObject(counts = wolffia_multiome_data$`Gene Expression`, assay = "RNA", project = "wolffia_multi")
seurat_wolffia_multi[['ATAC']] <- CreateChromatinAssay(counts = wolffia_multiome_data$`Peaks`, annotation = gr_multi, fragments = "atac_fragments.tsv.gz", sep = c(":", "-"))
seurat_wolffia_multi <- PercentageFeatureSet(seurat_wolffia_multi, pattern = "^MT-", col.name = "percent.mt", assay = "RNA")
seurat_wolffia_multi <- PercentageFeatureSet(seurat_wolffia_multi, pattern = "^CP-", col.name = "percent.cp", assay = "RNA")
seurat_wolffia_multi <- NucleosomeSignal(seurat_wolffia_multi, assay = "ATAC")
seurat_wolffia_multi <- TSSEnrichment(seurat_wolffia_multi, assay = "ATAC")
VlnPlot(seurat_wolffia_multi, features = c("nFeature_RNA", "percent.mt", "percent.cp", "nFeature_ATAC", "TSS.enrichment", "nucleosome_signal"), ncol = 6, pt.size = 0)
seurat_wolffia_multi <- subset(seurat_wolffia_multi, subset = nFeature_RNA > 100 & nFeature_RNA < 3000 & percent.mt < 1 & percent.cp < 20 & nFeature_ATAC > 100 & nFeature_ATAC < 18000 & TSS.enrichment > 1 & nucleosome_signal < 1)

DefaultAssay(seurat_wolffia_multi) <- "RNA"
seurat_wolffia_multi <- NormalizeData(seurat_wolffia_multi) %>%
  FindVariableFeatures(nfeatures = 3000) %>%
  ScaleData() %>%
  RunPCA(npcs = 50) %>%
  RunUMAP(dims = 1:20, reduction.name = "umap_rna", reduction.key = "UMAPRNA_")

seurat_wolffia_multi <- FindNeighbors(seurat_wolffia_multi,
                        reduction = "umap_rna",
                        dims = 1:ncol(Embeddings(seurat_wolffia_multi,"umap_rna"))) %>%
                        FindClusters(resolution = 0.5)
p1 <- DimPlot(seurat_wolffia_multi,
              group.by="RNA_snn_res.0.5",
              reduction="umap_rna", label=T) & NoAxes() & NoLegend()
p1

DefaultAssay(seurat_wolffia_multi) <- "ATAC"
seurat_wolffia_multi <- FindTopFeatures(seurat_wolffia_multi, min.cutoff = 50)
seurat_wolffia_multi <- RunTFIDF(seurat_wolffia_multi, method = 1)
seurat_wolffia_multi <- RunSVD(seurat_wolffia_multi, n = 50)
p1 <- ElbowPlot(seurat_wolffia_multi, ndims = 30, reduction="lsi")
p2 <- DepthCor(seurat_wolffia_multi, n = 30)
p1 | p2
seurat_wolffia_multi <- RunUMAP(seurat_wolffia_multi,
                  reduction = "lsi",
                  dims = 2:30,
                  reduction.name = "umap_atac",
                  reduction.key = "UMAPATAC_")
p1 <- DimPlot(seurat_wolffia_multi,
              group.by = "orig.ident",
              reduction = "umap_atac") & NoAxes()
p1
seurat_wolffia_multi <- FindMultiModalNeighbors(seurat_wolffia_multi,
                                  reduction.list = list("pca", "lsi"),
                                  dims.list = list(1:ncol(Embeddings(seurat_wolffia_multi,"pca")),
                                                   1:ncol(Embeddings(seurat_wolffia_multi,"lsi"))),
                                  modality.weight.name = c("RNA.weight","ATAC.weight"),
                                  verbose = TRUE)
seurat_wolffia_multi <- RunUMAP(seurat_wolffia_multi, nn.name = "weighted.nn", assay = "RNA")
seurat_wolffia_multi <- FindClusters(seurat_wolffia_multi, graph.name = "wsnn", resolution = 0.5)
UMAPPlot(seurat_wolffia_multi, group.by = "orig.ident") & NoAxes()
UMAPPlot(seurat_wolffia_multi, group.by = "wsnn_res.0.5", label=T) 

-------------------
#Marker genes
DefaultAssay(seurat_wolffia_multi) <- "RNA"
markers_wolffia_multiome_final <- FindAllMarkers(seurat_wolffia_multi, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers_wolffia_multiome_final %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)
-------------------
#Plotting expression of top marker genes
DefaultAssay(seurat_wolffia_multi) <- "RNA"
markers.to.plot <- c("Wolffia04G001880", "Wolffia07G002780", "Wolffia08G000040", "Wolffia13G001990", "Wolffia13G005280", "Wolffia14G005420", "Wolffia16G003700")
DotPlot(seurat_wolffia_multi, features = markers.to.plot, cols = c("skyblue", "red", "lightgrey", "black", "darkblue", "green", "darkgreen", "yellow", "darkred", "maroon", "magenta", "brown", "violet", "orange", "purple", "blue"), dot.scale = 8, split.by = "celltype") + RotatedAxis()



