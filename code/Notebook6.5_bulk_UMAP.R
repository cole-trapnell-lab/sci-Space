# Single cell data is valuable because it gives you a sense of what the functional
# unit of biology (the cell) is doing. In this notebook we aggregated the nuclei
# recovered in the sci-Space data set per position and asked how this leads to 
# destruction the structure present in single cell data

# Load startup packages ---------------------------------------------------
suppressPackageStartupMessages({ 
  library(tidyverse)
  library(ggplot2)
  library(tictoc)
  library(monocle3)
  library(tidymodels)
  library(furrr)
  library(ggplot2)
  library(automap)
  library(spatstat)
  library(imager)
  library(vec2dtransf)
  library(sp)
  library(sf)
  library(pheatmap)
  
  space_directory = "/Volumes/GoogleDrive/My Drive/sciSpace/"
  setwd(dir=space_directory)
  cc.genes <- readRDS("Submission_Data/bin/cc.genes.mouse.RDS")
  
  # Pass TRUE if you want to see progress output on some of Monocle 3's operations
  DelayedArray:::set_verbose_block_processing(TRUE)
  # Passing a higher value will make some computations faster but use more memory. 
  # Adjust with caution!
  options(DelayedArray.block.size=1000e6)
  options(future.globals.maxSize= 891289600)
  set.seed(42)
})

# Read in sci-Space data
spatial_cds = 
  readRDS(file = "Submission_Data/E14_slides/RDS_intermediates/Notebook6_spatial_cds_anatomy.RDS")

# Calculate values for the text, mean nuclei per position and s.d.
colData(spatial_cds) %>%
  as.data.frame() %>%
  group_by(top_spot,max_slide_id) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  summarise(mean_n = mean(n),
            sd_n = sqrt(var(n)))

# Aggregate positions by summing the gene expression counts from nuclei in 
# each position
aggregated_counts_by_spot =
  colData(spatial_cds) %>%
  as.data.frame() %>% 
  group_by(top_spot,max_slide_id) %>%
  filter(n()>1) %>%
  nest() %>%
  mutate(aggregated_counts = 
           purrr::map(.x = data,
                      .f = function(coldata_subset,
                                    cds){
                        counts(cds)[,coldata_subset$Cell] %>%
                          Matrix::rowSums()
                      },spatial_cds))


aggregated_count_matrix = 
  do.call(cbind, aggregated_counts_by_spot$aggregated_counts)

# Give each column a unique name
colnames(aggregated_count_matrix) = 
  paste(aggregated_counts_by_spot$top_spot,
        "_",
        aggregated_counts_by_spot$max_slide_id,
        sep ="")

# Look at how many spots have one nucleus mapping vs. multiple nuclei
colData(spatial_cds) %>% 
  as.data.frame() %>%   
  group_by(top_spot,max_slide_id) %>%
  mutate(cells_per_spot =
           ifelse(n() > 1,
                  "multiple",
                  "single")) %>%
  group_by(max_slide_id,cells_per_spot) %>%
  summarise(n = n()) %>%
  as.data.frame()

single_cells_in_spot = 
  colData(spatial_cds) %>%
  as.data.frame() %>%
  group_by(top_spot,max_slide_id) %>%
  filter(n() == 1) %>%
  pull(Cell)

colData(spatial_cds) %>% 
  as.data.frame() %>%   
  group_by(top_spot,max_slide_id) %>% 
  summarise(n = n()) %>%
  pull()


# Make an aggregated count matrix
single_cell_count_mat = 
  counts(spatial_cds[,single_cells_in_spot]) %>%
  as.matrix()

colnames(single_cell_count_mat) = 
  paste(colData(spatial_cds[,single_cells_in_spot])$top_spot,
        "_",
        colData(spatial_cds[,single_cells_in_spot])$max_slide_id,
        sep = "")

joint_count_mat =
  cbind(single_cell_count_mat,
        aggregated_count_matrix)


# Make an aggregated column metadata dataframe
coldata = 
  data.frame(row.names = colnames(joint_count_mat) %>% as.character(),
             position = colnames(joint_count_mat) %>% as.character(),
             single_cell = colnames(joint_count_mat) %in% colnames(single_cell_count_mat))

print(paste(dim(spatial_cds)[2],
            " single cells",
            sep = ""))

print(
  paste(dim(spatial_cds)[2] - sum(coldata$single_cell),
        " cells aggregated over ",
        sum(!coldata$single_cell),
        " positions",
        sep = ""))

rowdata = 
  rowData(spatial_cds) %>%
  as.data.frame() 

identical(rownames(rowdata) %>% as.character(),
          rownames(joint_count_mat) %>% as.character())

joint_cds = 
  new_cell_data_set(expression_data = joint_count_mat,
                    cell_metadata = coldata,
                    gene_metadata = rowdata)

# Make sure there is slide_id column
colData(joint_cds)$max_slide_id = 
  stringr::str_split_fixed(colData(joint_cds)$position %>% as.character(),
                           "_",
                           3)[,3]

# UMIs per position
colData(joint_cds)$n.umi = 
  Matrix::colSums(counts(joint_cds))

# genes captured per position
colData(joint_cds)$n.genes = 
  Matrix::colSums(counts(joint_cds) >= 1)

# saveRDS(object = joint_cds,
#         file = "Submission_Data/E14_slides/RDS_intermediates/Notebook6.5_joint_bulk_CDS.RDS")
# joint_cds = 
#   readRDS(file = "Submission_Data/E14_slides/RDS_intermediates/Notebook6.5_joint_bulk_CDS.RDS")

# Calculate per position statistics ---------------------------------------

colData(joint_cds) %>%
  as.data.frame() %>%
  dim()

colData(joint_cds) %>%
  as.data.frame() %>%
  summarise(mean_umis = mean(n.umi),
            sd_umis = sqrt(var(n.umi)),
            mean_genes = mean(n.genes),
            sd_genes = sqrt(var((n.genes))))


colData(joint_cds) %>%
  as.data.frame() %>%
  group_by(max_slide_id) %>%
  summarise(med_genes = median(n.genes))

table_1 = 
  colData(joint_cds) %>%
  as.data.frame() %>%
  group_by(max_slide_id) %>%
  summarise(genes_per_position = mean(n.genes),
            umis_per_position = mean(n.umi),
            number_of_positions = n())

table_2 = 
  colData(spatial_cds) %>%
  as.data.frame() %>% 
  group_by(top_spot, max_slide_id) %>%
  summarise(cells_per_position = n()) %>%
  group_by(max_slide_id) %>%
  summarise(avg_cells_per_position = mean(cells_per_position))

table_3 =
  colData(spatial_cds) %>%
  as.data.frame() %>% 
  group_by(max_slide_id) %>%
  summarise(total_cells = n())

# Data frame used to make Fig S14
left_join(table_1,
           table_2) %>%
  left_join(table_3)
  

# Dimensionality reduction of aggregated data -----------------------------

# Run PCA, UMAP and clustering
bulked_cds = 
  joint_cds[,!colData(joint_cds)$single_cell] %>%
  estimate_size_factors() %>%
  preprocess_cds(num_dim = 100) %>%
  reduce_dimension() %>%
  cluster_cells(resolution = 1e-3)

# Add cluster information to the aggregated CDS
colData(bulked_cds)$bulk_clusters = 
  clusters(bulked_cds)

# Add UMAP embedding information to the aggregated CDS
colData(bulked_cds)$umap1 = reducedDim(x = bulked_cds,type = "UMAP")[,1]
colData(bulked_cds)$umap2 = reducedDim(x = bulked_cds,type = "UMAP")[,2]

colData(bulked_cds) %>%
  as.data.frame() %>%
  group_by()

# How many clusters from the single cell UMAP (Fig 1C) are contained within each spot?
colData(bulked_cds)$clusters_per_position = 
  colData(spatial_cds) %>%
  as.data.frame() %>% 
  mutate(position = 
           paste(top_spot,
                 max_slide_id,
                 sep = "_"),
         cluster = as.character(cluster)) %>%
  filter(position %in% colData(bulked_cds)$position) %>%
  dplyr::select(position, cluster) %>%
  arrange(position,cluster) %>%
  distinct() %>%
  group_by(position) %>%
  summarise(n = n()) %>%
  right_join(colData(bulked_cds) %>%
               as.data.frame(),
             by = "position") %>%
  pull(n)

# How many cell types from the single cell UMAP (Fig 1C) are contained within each spot?
colData(bulked_cds)$celltypes_per_position = 
  colData(spatial_cds) %>%
  as.data.frame() %>% 
  mutate(position = 
           paste(top_spot,
                 max_slide_id,
                 sep = "_")) %>%
  filter(position %in% colData(bulked_cds)$position) %>%
  dplyr::select(position, final_cluster_label) %>%
  arrange(position, final_cluster_label) %>%
  distinct() %>%
  group_by(position) %>%
  summarise(n = n()) %>%
  right_join(colData(bulked_cds) %>%
               as.data.frame(),
             by = "position") %>%
  pull(n)


# What are the labels of the single cell type spots?
colData(bulked_cds)$single_celltypes = 
  colData(spatial_cds) %>%
  as.data.frame() %>% 
  mutate(position = 
           paste(top_spot,
                 max_slide_id,
                 sep = "_")) %>%
  filter(position %in% colData(bulked_cds)$position) %>%
  dplyr::select(position, final_cluster_label) %>%
  distinct() %>%
  group_by(position) %>%
  filter(n() == 1) %>%
  right_join(colData(bulked_cds) %>%
               as.data.frame(),
             by = "position") %>%
  pull(final_cluster_label)


# Supplementary Figure 21A,B ----------------------------------------------

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


# Fig S21A
ggplot() +
  geom_point(data = 
               colData(bulked_cds) %>%
               as.data.frame(),
             aes(x = umap1,
                 y = umap2),
             color = "black",
             stroke = 0,
             size = 0.75) +
  geom_point(data = 
               colData(bulked_cds) %>%
               as.data.frame() %>%
               mutate(one_cluster = ifelse(clusters_per_position == 1,
                                           T,
                                           NA)),
             aes(x = umap1,
                 y = umap2,
                 color = one_cluster),
             stroke = 0,
             size = 0.65) +
  scale_color_manual(values = c("red"),na.value="grey90") +
  theme_void() +
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Supplement_bulk_aggregation/aggregated_cells_by_sing_cluster.png",
         bg = "transparent",
         dpi =450,
         height = 2.5,
         width = 2.5)

# Fig S21B
ggplot() +
  geom_point(data = 
               colData(bulked_cds) %>%
               as.data.frame(),
             aes(x = umap1,
                 y = umap2),
             color = "black",
             stroke = 0,
             size = 0.75) +
  geom_point(data = 
               colData(bulked_cds) %>%
               as.data.frame(),
             aes(x = umap1,
                 y = umap2,
                 color = single_celltypes),
             stroke = 0,
             size = 0.65) +
  scale_color_manual(values = colors,na.value="grey90") +
  theme_void() +
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Supplement_bulk_aggregation/aggregated_cells_by_celltype.png",
         bg = "transparent",
         dpi =450,
         height = 2.5,
         width = 2.5)


# Supplementary Figure S21C -----------------------------------------------

# Sample the same number of single cells to make single cell UMAP
sampled_cells =
  colData(spatial_cds) %>%
  as.data.frame() %>% 
  group_by(top_spot,max_slide_id) %>%
  filter(n()>1) %>%
  sample_n(size = 1) %>%
  pull(Cell) %>%
  as.character()


sampled_cds = 
  spatial_cds[,c(sampled_cells)] %>%
  estimate_size_factors() %>%
  preprocess_cds(num_dim = 100) %>%
  reduce_dimension() %>%
  cluster_cells(resolution = 1e-3)

colData(sampled_cds)$sampled_clusters = 
  clusters(sampled_cds)

colData(sampled_cds)$umap1 = reducedDim(x = sampled_cds,type = "UMAP")[,1]
colData(sampled_cds)$umap2 = reducedDim(x = sampled_cds,type = "UMAP")[,2]

colData(sampled_cds) %>%
  as.data.frame() %>%
  ggplot() +
  geom_point(aes(x = umap1,
                 y = umap2),
             color = "black",
             stroke = 0,
             size = 0.5) +
  geom_point(aes(x = umap1,
                 y = umap2,
                 color = final_cluster_label),
             stroke = 0,
             size = 0.4) +
  scale_color_manual(values = colors,na.value="grey80") +
  theme_void() +
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Supplement_bulk_aggregation/sampled_cells_by_celltype.png",
         bg = "transparent",
         dpi =450,
         height = 2.5,
         width = 2.5)
