# Process single-cell transcriptomes from E14 embryos: Perform dimensionality reduction and
# subclustering on the collected cells.

# Load startup packages ---------------------------------------------------
suppressPackageStartupMessages({ 
  library(tidyverse)
  library(ggplot2)
  library(tidyr)
  library(viridis)
  library(ggridges)
  library(randomcoloR)
  library(ggrepel)
  library(monocle3)
  library(garnett)
  space_directory = "~/Google Drive File Stream/My Drive/sciSpace/"
  setwd(dir=space_directory)
  
  source("Submission_Data/bin/cell_cycle.R")
  cc.genes <- readRDS("Submission_Data/bin/cc.genes.mouse.RDS")
  # Pass TRUE if you want to see progress output on some of Monocle 3's operations
  DelayedArray:::set_verbose_block_processing(TRUE)
  # Passing a higher value will make some computations faster but use more memory. 
  # Adjust with caution!
  options(DelayedArray.block.size=1000e6)
})


# Set a seed to make umap and other non-deterministic steps consistent
set.seed(seed = 42)

spatial_cds = readRDS("Submission_Data/E14_slides/RDS_intermediates/Notebook1_E14_spatial_CDS.RDS")


# Incorporate scrublet scores ---------------------------------------------
scrublet_score = 
  read.table(file = "Submission_Data/E14_slides/RDS_intermediates/scrublet_doublet_scores.txt",
             header = F)[,1]

colData(spatial_cds)$scrublet_score = 
  scrublet_score

ggplot(colData(spatial_cds) %>%
         as.data.frame()) +
  geom_histogram(aes(x = scrublet_score),
                 fill = "grey90",
                 color = "black",
                 size = 0.25) +
  geom_vline(xintercept = 0.15,
             color = "red") +
  monocle3:::monocle_theme_opts() +
  ggsave("Submission_Data/E14_slides/Auxillary_Figures/Notebook2_scrublet_scores.pdf",
         height = 2,
         width = 2)

# Remove cells recieving a doublet score of greater than 0.15
colData(spatial_cds)$scrublet_doublet = 
  colData(spatial_cds)$scrublet_score >= 0.15

sum(colData(spatial_cds)$scrublet_doublet)

colData(spatial_cds) %>%
  as.data.frame() %>%
  group_by(scrublet_doublet) %>%
  summarise(n = n())

# Keep cells that aren't doublets
spatial_cds = spatial_cds[,!colData(spatial_cds)$scrublet_doublet]

# > dim(spatial_cds)
# [1]  52636 122278

# Calculate cell size factor and perform PCA ----------------------------
spatial_cds = 
  spatial_cds %>%
  estimate_size_factors() %>%
  estimate_cell_cycle(g1s_markers = cc.genes$s.genes,
                      g2m_markers = cc.genes$g2m.genes) %>%
  preprocess_cds(method = "PCA",
                 num_dim = 100)

plot_pc_variance_explained(cds = spatial_cds) +
  ggsave(filename = "Submission_Data/E14_slides/Auxillary_Figures/Notebook3_scree_plot.pdf",
         height = 4,
         width = 4)


# Perform UMAP dimensionality reduction -----------------------------------

colData(spatial_cds)$experiment =
  ifelse(colData(spatial_cds)$sample == 1,
         "experiment_1",
         "experiment_2")

colData(spatial_cds)$log.n.umi =
  log10(colData(spatial_cds)$n.umi) 
 
# Perform data integration based on prepared batch -- coded in the sample column
spatial_cds = 
  align_cds(spatial_cds,
            residual_model_formula_str = "~log.n.umi + sample",
            alignment_group = "sample")

# UMAP using a fixed seed and umap.fast_sgd = F
spatial_cds =
  reduce_dimension(spatial_cds,
                   reduction_method = "UMAP",
                   umap.metric = "cosine",
                   umap.fast_sgd = F)

# Add UMAP coordinates to CDS metadata
spatial_cds$umap1 = reducedDim(spatial_cds, type = "UMAP")[,1]
spatial_cds$umap2 = reducedDim(spatial_cds, type = "UMAP")[,2]

# Perform louvain clustering -----------------------------------
spatial_cds = cluster_cells(spatial_cds,
                            random_seed = 42)

# Add cluster and partition information to CDS metadata
colData(spatial_cds)$cluster = clusters(spatial_cds)
colData(spatial_cds)$partition = partitions(spatial_cds)

# Subcluster partitions from clusters  ------------------------------------

# Iteratively subcluster all the partitions
sub_clustered_partitions =
  lapply(X = partitions(spatial_cds) %>% 
         unique(),
       FUN = function(p){
         cds_p = spatial_cds[,colData(spatial_cds)$partition == p]
         
         # Re-estimate size_factors
         cds_p = 
           estimate_size_factors(cds_p)
         cds_p = 
           detect_genes(cds_p, 
                        min_expr = 0.1)
         # Run PCA prior to initalize dimensionality reduction
         cds_p = 
           preprocess_cds(cds = cds_p,
                          method = "PCA",
                          num_dim = 100)
         # Perform UMAP dimensionality reduction
         cds_p = 
           align_cds(cds_p,
                     residual_model_formula_str = "~log.n.umi",
                     alignment_group = "sample")
         
         
         cds_p = 
           reduce_dimension(cds = cds_p,
                            max_components = 2,
                            reduction_method = "UMAP",
                            preprocess_method = 'PCA',
                            umap.fast_sgd = F)
         
         # Recluster these cells
         cds_p = cluster_cells(cds_p,
                               random_seed = 42,
                               resolution = 5e-4)
         
         plot_cells(cds_p) +
           ggsave(filename = 
                    paste0("Submission_Data/E14_slides/Auxillary_Figures/",
                         p,
                         "_partition_subcluster.pdf",
                         sep = ""),
                  height = 3,
                  width = 3)
         
         # Append a sub_cluster relating to the orginal partition
         colData(cds_p)$sub_cluster = 
           paste0(colData(cds_p)$partition,
                  ".",
                  clusters(cds_p,
                           reduction_method = "UMAP"))
         
         colData(cds_p) %>%
           as.data.frame() %>%
           dplyr::select(Cell, sub_cluster,cluster)
         
       })

sub_clustered_partitions =
  do.call(rbind,
          sub_clustered_partitions)


# Check to see if all cells have been assigned a subcluster
sum(sub_clustered_partitions$Cell %in% colData(spatial_cds)$Cell) == dim(spatial_cds)[2]

colData(spatial_cds)$sub_cluster =
  colData(spatial_cds) %>% 
  as.data.frame() %>%
  left_join(sub_clustered_partitions,
            by = c("Cell", "cluster")) %>%
  pull(sub_cluster)

colData(spatial_cds)$sub_cluster %>%
  unique %>%
  length()

is.na(colData(spatial_cds)$sub_cluster) %>% sum()


# Write out resulting RDS and mtx files -----------------------------------------------
saveRDS(object = spatial_cds,
        file = "Submission_Data/E14_slides/RDS_intermediates/Notebook3_E14_spatial_CDS.RDS")


dim(spatial_cds)
