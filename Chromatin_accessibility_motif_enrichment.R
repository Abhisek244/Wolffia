#Linking marker genes to nearby peaks
library(Biostrings)
library(GenomicRanges)
library(data.table)
custom_genome_fasta <- "Wmic_reference_genome.fasta"
custom_genome <- readDNAStringSet(custom_genome_fasta, format = "fasta")
seurat_wolffia_multi <- RegionStats(seurat_wolffia_multi, genome = custom_genome)
seurat_wolffia_multi <- LinkPeaks(
  object = seurat_wolffia_multi,
  peak.assay = "ATAC",
  expression.assay = "RNA",
  genes.use = top_markers_ct$feature
)
-------------------------
#Coverageplots
new_cluster_ids <- c("Cluster 0", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6",  "Cluster 7", "Cluster 8", "Cluster 9", "Cluster 10", "Cluster 11", "Cluster 12", "Cluster 13", "Cluster 14", "Cluster 15")
names(new_cluster_ids) <- levels(seurat_wolffia_multi)
seurat_wolffia_multi <- RenameIdents(seurat_wolffia_multi, new_cluster_ids)
seurat_wolffia_multi$celltype <- Idents(seurat_wolffia_multi)
print(head(seurat_wolffia_multi$celltype))
CoveragePlot(seurat_wolffia_multi,
                   region = "Wolffia04G002500",
                   features = "Wolffia04G002500",
                   group.by = "celltype",
                   extend.upstream = 1000,
                   extend.downstream = 1000)
--------------------------
#TF binding motif enrichment analysis
library(TFBSTools)
library(JASPAR2020)
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'plants', all_versions = FALSE)
)
df_pfm <- data.frame(t(sapply(pfm, function(x)
  c(id=x@ID, name=x@name, symbol=ifelse(!is.null(x@tags$symbol),x@tags$symbol,NA)))))

library(Biostrings)
custom_genome2 <- readDNAStringSet("Wmic_reference_genome.fasta")
DefaultAssay(seurat_wolffia_multi) <- "ATAC"
seurat_wolffia_multi <- AddMotifs(seurat_wolffia_multi, genome = custom_genome2, pfm = pfm)

open_peaks <- AccessiblePeaks(seurat_wolffia_multi)
peaks_matched <- MatchRegionStats(meta.feature = seurat_wolffia_multi[['ATAC']]@meta.features[open_peaks, ],
                                  query.feature = seurat_wolffia_multi[['ATAC']]@meta.features[top_peaks_ct$feature, ], n = 30269)

motif_enrichment_endo <- FindMotifs(seurat_wolffia_multi,
                                    features = top_peaks_ct$feature[top_peaks_ct$group == "0"],
                                    background = peaks_matched) %>%
  mutate(symbol = setNames(ifelse(is.na(df_pfm$symbol), df_pfm$name, df_pfm$symbol), df_pfm$id)[motif]) %>%
  mutate(padj = p.adjust(pvalue, method="BH"))
enriched_motif_endo <- motif_enrichment_endo %>%
  dplyr::filter(padj < 0.01 & fold.enrichment > 0.25) %>%
  top_n(-4, wt = padj)
print(head(enriched_motif_endo))
p1 <- MotifPlot(seurat_wolffia_multi, motifs = enriched_motif_endo$motif[1:2], ncol=2)
p1
