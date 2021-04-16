# Find gene modules based on a cell's position in UMAP space versus spatial position
# Analysis for Fig S30

# Load startup packages ---------------------------------------------------
suppressPackageStartupMessages({ 
  library(tidyverse)
  library(ggplot2)
  library(monocle3)
  library(ggplot2)
  space_directory = "/Volumes/GoogleDrive/My Drive/sciSpace/"
  setwd(dir=space_directory)
  
  # Pass TRUE if you want to see progress output on some of Monocle 3's operations
  DelayedArray:::set_verbose_block_processing(TRUE)
  # Passing a higher value will make some computations faster but use more memory. 
  # Adjust with caution!
  options(DelayedArray.block.size=1000e6)
  options(future.globals.maxSize= 891289600)
  set.seed(42)
})

spatial_cds = 
  readRDS(file = "Submission_Data/E14_slides/RDS_intermediates/Notebook6_spatial_cds_anatomy.RDS")

slide4e_cells = 
  spatial_cds %>%
  colData() %>%
  as.data.frame() %>%
  filter(max_slide_id == "slide_4E") %>%
  pull(Cell)

slide4e_cds = 
  spatial_cds[,slide4e_cells]

slide4e_cds = 
  slide4e_cds %>%
  estimate_size_factors() %>%
  detect_genes() %>%
  preprocess_cds() %>%
  reduce_dimension(preprocess_method ="PCA",
                   max_components = 2,
                   reduction_method = "UMAP",
                   umap.fast_sgd = F) %>%
  cluster_cells(resolution = 1e-2,
                random_seed = 42)

# Add UMAP coordinates to metadata
colData(slide4e_cds)$umap1 = 
  reducedDim(slide4e_cds,
             type = "UMAP")[,1]

colData(slide4e_cds)$umap2 = 
  reducedDim(slide4e_cds,
             type = "UMAP")[,2]

# Add UMAP clusters to metadata
colData(slide4e_cds)$umap_cluster = 
  slide4e_cds@clusters[["UMAP"]][["clusters"]]

# Perform graph_test on UMAP coordinates this function 
# calculates autocorrelation on the nearest neighbor graph
pr_graph_test_UMAP =
  graph_test(slide4e_cds, 
             neighbor_graph="knn",
             reduction_method = "UMAP",
             cores = 8)

slide4e_cds_umap = slide4e_cds
slide4e_cds_space = slide4e_cds


# Move spatial coordinates to the UMAP slot
reducedDim(x = slide4e_cds_space,
           type = "UMAP") <-
  matrix(cbind(colData(slide4e_cds_space)$coords.x2,
               -colData(slide4e_cds_space)$coords.x1), 
         ncol=2)


slide4e_cds_space = cluster_cells(slide4e_cds_space,
                    resolution = 1e-2,
                    random_seed = 42,
                    reduction_method = "UMAP")

# Perform graph_test on spatial coordinates this function 
# calculates autocorrelation on the nearest neighbor graph
pr_graph_test_spatial =
  graph_test(slide4e_cds_space, 
             neighbor_graph="knn",
             reduction_method = "UMAP",
             cores = 8)


# Get genes that have significant autocorrelation by both measures
genes_of_interest = 
  union(pr_graph_test_UMAP %>% filter(q_value < 0.05) %>% pull(id),
            pr_graph_test_spatial %>% filter(q_value < 0.05) %>% pull(id))

# umap gene modules
gene_module_df_umap <- 
  find_gene_modules(slide4e_cds_umap[genes_of_interest,], 
                    resolution=1e-2)

# spatial gene modules
gene_module_df_space <- 
  find_gene_modules(slide4e_cds_space[genes_of_interest,], 
                    resolution=1e-2)

# Supplemental Figure S30A ------------------------------------------------
module_df_joint = 
  inner_join(gene_module_df_umap %>% dplyr::select(id, UMAP = module),
             gene_module_df_space %>% dplyr::select(id, Spatial = module),
             by = "id")

module_df_joint %>%
  gather(key = "computation",
         value = "module",
         -id) %>%
  ggplot() +
  geom_violin(aes(x = computation,
                   y = as.numeric(module)),
              fill = "grey80") +
  theme_classic() +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 6)) +
  ylab("Genes per Module") +
  xlab("Embedding")+
  ggsave("Figures/Figure_Components/Supplemental_spatial_modules_heatmap/modules_sizes.pdf",
         height = 2,
         width = 1.5)


gene_module_df_space %>%
  left_join(rowData(spatial_cds) %>%
              as.data.frame(),
            by = "id") %>%
  filter(grepl("Hox",gene_short_name))

# Calculate the adjusted rand index for the two groupings
mclust::adjustedRandIndex(module_df_joint$UMAP,module_df_joint$Spatial)

cell_group_df_umap <- 
  tibble::tibble(cell=row.names(colData(slide4e_cds_space)),
                 cell_group=colData(slide4e_cds_space)$cluster)

# Supplemental Figure S30B ------------------------------------------------
agg_mat_umap <- 
  aggregate_gene_expression(slide4e_cds_umap, 
                            gene_module_df_umap, 
                            cell_group_df_umap)

row.names(agg_mat_umap) <- stringr::str_c("Module ", row.names(agg_mat_umap))
colnames(agg_mat_umap) <- stringr::str_c("Cluster ", colnames(agg_mat_umap))



cluster_annotations = 
  colData(slide4e_cds_space) %>%
  as.data.frame() %>%
  mutate(cluster  = paste("Cluster ", cluster,sep = "")) %>%
  dplyr::select(cluster,
                final_cluster_label) %>%
  group_by(cluster,final_cluster_label) %>%
  add_tally() %>%
  ungroup() %>%
  group_by(cluster) %>%
  distinct() %>%
  arrange(-n) %>%
  top_n(n = 1, wt = n) %>%
  distinct() %>%
  ungroup() %>%
  column_to_rownames(var = "cluster") %>%
  dplyr::select(final_cluster_label)



colors =c("Cardiac muscle lineages" = "#FF3333",
          "Chondrocytes" = "#E1B0CF",
          "Choroid Plexus"  = "#00A6A6",
          "Connective Tissue Progenitors" = "#EBCF00",
          "Developing Gut" = "#0050D4",
          "Endothelial Cells" = "#64E8A6",
          "Epithelial Cells" = "#EB7700",
          "Erythroid Lineage" = "#DCE4C9",
          "Fibroblast"  = "#790009",
          "Glial Cells" = "#8172E7",
          "Hepatocytes" = "#8538E8",
          "Lateral Plate Mesoderm" = "#086661",
          "Myocytes" = "#00D9F7",
          "Neuron" = "#B4D900",
          "OPCs" = "#E7E041",
          "Peripheral Neuron" = "#00B562",
          "Radial glia" = "#FF00C3",
          "Schwann Cells" = "#E1AF52",
          "Testis Cells" = "black",
          "White Blood Cells" = "#3288BD")

cluster_colors = list(final_cluster_label = colors)

# Supplemental Figure S30A
pheatmap::pheatmap(agg_mat_umap, 
                   clustering_method="ward.D2",
                   fontsize=6,
                   cluster_rows = T,
                   cluster_cols = T,
                   treeheight_row = 0,
                   treeheight_col = 0,
                   show_rownames = F,
                   show_colnames = F,
                   annotation_col = cluster_annotations,
                   annotation_colors = cluster_colors,
                   annotation_names_col = F,
                   annotation_legend = F,
                   legend = F,
                   height = 2,
                   width = 2,
                   filename = "Figures/Figure_Components/Supplemental_spatial_modules_heatmap/umap_modules_clusters.pdf")

# Supplemental Figure S30C ------------------------------------------------

agg_mat_umap_groups_space <- 
  aggregate_gene_expression(slide4e_cds_umap, 
                            gene_module_df_space, 
                            cell_group_df_umap)


row.names(agg_mat_umap_groups_space) <- stringr::str_c("Module ", row.names(agg_mat_umap_groups_space))
colnames(agg_mat_umap_groups_space) <- stringr::str_c("Cluster ", colnames(agg_mat_umap_groups_space))

pheatmap::pheatmap(agg_mat_umap_groups_space, 
                   clustering_method="ward.D2",
                   fontsize=6,
                   cluster_rows = T,
                   cluster_cols = T,
                   treeheight_row = 0,
                   treeheight_col = 0,
                   show_rownames = F,
                   show_colnames = F,
                   annotation_col = cluster_annotations,
                   annotation_colors = cluster_colors,
                   annotation_names_col = F,
                   annotation_legend = F,
                   legend = F,
                   height = 2,
                   width = 2,
                   filename = "Figures/Figure_Components/Supplemental_spatial_modules_heatmap/spatial_modules_clusters.pdf")

pheatmap::pheatmap(agg_mat_umap_groups_space, 
                   clustering_method="ward.D2",
                   fontsize=6,
                   cluster_rows = T,
                   cluster_cols = T,
                   treeheight_row = 0,
                   treeheight_col = 0,
                   show_rownames = F,
                   show_colnames = F,
                   annotation_col = cluster_annotations,
                   annotation_colors = cluster_colors,
                   annotation_names_col = F,
                   annotation_legend = F,
                   legend = T,
                   height = 1,
                   width = 1,
                   filename = "Figures/Figure_Components/Supplemental_spatial_modules_heatmap/spatial_modules_clusters_legend.pdf")


# Supplemental Figure S30D ------------------------------------------------

cell_group_df_space <- 
  tibble::tibble(cell=row.names(colData(slide4e_cds_space)),
                 cell_group=interaction(colData(slide4e_cds_space)$Row,
                                        colData(slide4e_cds_space)$Col))

agg_mat_space <- 
  aggregate_gene_expression(slide4e_cds_space, 
                            gene_module_df_space, 
                            cell_group_df_space) 



spatial_annotations = 
  data.frame(row.names = interaction(colData(slide4e_cds_space)$Row,
                                     colData(slide4e_cds_space)$Col) %>%
               unique(),
             spatial_coordinate = interaction(colData(slide4e_cds_space)$Row,
                                              colData(slide4e_cds_space)$Col) %>%
               unique())

spatial_annotations = 
  colData(slide4e_cds_space) %>%
  as.data.frame() %>%
  mutate(spatial_coordinate  = interaction(Row,Col)) %>%
  dplyr::select(spatial_coordinate,
                final_cluster_label) %>%
  group_by(spatial_coordinate,
           final_cluster_label) %>%
  add_tally() %>%
  ungroup() %>%
  distinct() %>%
  group_by(spatial_coordinate) %>%
  top_n(n = 1, wt = n) %>%
  sample_n(size = 1) %>%
  column_to_rownames(var = "spatial_coordinate") %>%
  dplyr::select(-n)

  

coordinate_colors = list(final_cluster_label = colors)

pheatmap::pheatmap(agg_mat_space, 
                   clustering_method="ward.D2",
                   fontsize=6,
                   annotation_col = spatial_annotations,
                   cluster_rows = T,
                   cluster_cols = T,
                   treeheight_row = 0,
                   treeheight_col = 0,
                   annotation_colors = coordinate_colors,
                   show_rownames = F,
                   show_colnames = F,
                   annotation_names_col = F,
                   annotation_legend = F,
                   height = 2,
                   width = 5.5,
                   legend = F,
                   filename = "Figures/Figure_Components/Supplemental_spatial_modules_heatmap/spatial_modules_positions.pdf")


pheatmap::pheatmap(agg_mat_space, 
                   clustering_method="ward.D2",
                   fontsize=6,
                   annotation_col = spatial_annotations,
                   cluster_rows = T,
                   cluster_cols = T,
                   treeheight_row = 0,
                   treeheight_col = 0,
                   annotation_colors = coordinate_colors,
                   show_rownames = F,
                   show_colnames = F,
                   annotation_names_col = F,
                   annotation_legend = F,
                   height = 1,
                   width = 1,
                   legend = T,
                   filename = "Figures/Figure_Components/Supplemental_spatial_modules_heatmap/spatial_modules_positions_legend.pdf")


  