library(Seurat)
library(monocle3)
counts <- GetAssayData(object = Wmic_snRNA, assay = "RNA", layer = "counts")
cell.metadata <- Wmic_snRNA[[]]
feature.metadata <- Wmic_snRNA[["RNA"]][[]]
feature.metadata$gene_short_name <- rownames(x = feature.metadata)
cds <- new_cell_data_set(counts,
                         cell_metadata = cell.metadata,
                         gene_metadata = feature.metadata)
cds <- preprocess_cds(cds, num_dim = 16)
cds <- align_cds(cds, alignment_group = "orig.ident")
cds <- reduce_dimension(cds)
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "seurat_clusters", group_label_size = 4)
cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition", label_cell_groups = FALSE)
cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "partition",
           label_groups_by_cluster=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE)
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)
pr_graph_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=1)
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))
-----------------
#Gene modules
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=1e-2)
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=partitions(cds)[colnames(cds)])
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("Partition ", colnames(agg_mat))
pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=6)
