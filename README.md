This repository contains the codes used for single-nuclei RNA-Seq (snRNA-Seq) and single-nuclei ATAC-Seq (snATAC-Seq) data analysis of Wolffia microscopica

**snRNA-Seq analysis**

CellRanger.sh - Commands for mapping the snRNA-Seq reads onto W. microscopica reference genome and constructing the gene-count matrix

snRNA_analysis.R - R script to process the gene-count matrix (constructed using CellRanger.sh) using Seurat to construct the single-nuclei clusters

Monocle3.R - R script to perform developmental trajectory and pseudotime analysis using the single-nuclei clusters obtained using snRNA_analysis.R

Wolffia_Integrated.R - R script to merge the snRNA-Seq dataset of W. microscopica and W. australiana for a comparative analysis

----------------
**#Multiome analysis**
CellRanger_arc.sh - Commands for mapping both snRNA-Seq and snATAC-Seq reads onto W. microscopica reference genome and constructing the gene-count matrix

snRNA_snATAC_multiome_analysis.R - R script to process the gene-count matrix (constructed using CellRanger_arc.sh) using Seurat to construct the single-nuclei clusters

Chromatin_accessibility_motif_enrichment.R - R script for chromatin accessibility profile analysis and TF-binding motif enrichment analysis
