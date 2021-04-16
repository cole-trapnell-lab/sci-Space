# How do transcriptional clusters found in La Manno et. al. Developing Mouse Brain Atlas
# correspond to spatial positions in the sci-Space data?

# This dataset can be downloaded at -- http://mousebrain.org/downloads.html

# G. La Manno, K. Siletti, A. Furlan, D. Gyllborg, E. Vinsland, 
# C. M. Langseth, I. Khven, A. Johnsson, M. Nilsson, P. Lönnerberg, 
# S. Linnarsson, Molecular architecture of the developing mouse brain. 
# Cold Spring Harbor Laboratory (2020), p. 2020.07.02.184051.

# https://www.biorxiv.org/content/10.1101/2020.07.02.184051v1


# Load startup packages ---------------------------------------------------
suppressPackageStartupMessages({ 
  library(tidyverse)
  library(viridis)
  library(ggpubr)
  library(ggrepel)
  library(pheatmap)
  library(monocle3)
  library(FNN)
  
  space_directory = "/Volumes/GoogleDrive/My Drive/sciSpace/"
  setwd(dir=space_directory)
  source("Submission_Data/bin/hotspot_functions.R")
  
  # Pass TRUE if you want to see progress output on some of Monocle 3's operations
  DelayedArray:::set_verbose_block_processing(TRUE)
  # Passing a higher value will make some computations faster but use more memory. 
  # Adjust with caution!
  options(DelayedArray.block.size=1000e6)
  set.seed(42)
})


# Load in cds objects -----------------------------------------------------
all_image_data = 
  readRDS("Submission_Data/E14_slides/RDS_intermediates/Notebook0_1_images_and_transformations.RDS")


spatial_cds = 
  readRDS(file = "Submission_Data/E14_slides/RDS_intermediates/Notebook6.01_spatial_cds_anatomy.RDS")


laManno_joint_cds = 
  readRDS("Submission_Data/E14_slides/RDS_intermediates/sci_space_cortex_neuron_e13p5_e14p5_joint_embedding.RDS")


laManno_coldata = 
  read.table(file = "Submission_Data/E14_slides/RDS_intermediates/LaManno_coldata.tsv",
             header = T,
             row.names = 1,
             sep = "\t")


# Supplemental Figure S22 — Anatomical Region -----------------------------

scaling_df = 
  data.frame(max_slide_id = c("slide_1D",
                              "slide_1E",
                              "slide_1F",
                              "slide_1G",
                              "slide_2G",
                              "slide_2H",
                              "slide_3D",
                              "slide_3F",
                              "slide_3G",
                              "slide_3H",                              
                              "slide_4A",
                              "slide_4D",
                              "slide_4E"),
             x = c(1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1),
             y = c(-1,-1,-1,-1,-1,-1,1,-1,-1,-1,-1,1,1))

jittered_positions = 
  colData(spatial_cds) %>%
  as.data.frame() %>%
  left_join(scaling_df,
            by = "max_slide_id") %>%
  filter(anatomical_annotation == "Cortex",
         partition == "2",
         !is.na(lamanno_Tissue)) %>%
  mutate(lamanno_Tissue = ifelse(grepl(pattern = "Midbrain",
                                       lamanno_Tissue),
                                 "Midbrain",
                                 lamanno_Tissue))

jittered_positions =  
  jittered_positions %>%
  mutate(jittered_y = coords.x1 * y + rnorm(n =dim(jittered_positions)[1],
                                            mean = 0,
                                            sd = 10),
         jittered_x = coords.x2 * x + rnorm(n =dim(jittered_positions)[1],
                                            mean = 0,
                                            sd = 10))


ggplot(jittered_positions) +
  geom_point(aes(y = jittered_y,
                 x = jittered_x),
             color = "black",
             stroke = 0,
             size = 0.75)+
  
  geom_point(aes(y = jittered_y,
                 x = jittered_x,
                 color = lamanno_Tissue),
             stroke = 0,
             size = 0.65)+
  facet_wrap(~slide_id, scales = "free") +
  theme_void() +
  scale_color_brewer(palette = "Spectral") +
  theme(legend.position = "none")+
  ggsave("Figures/Figure_Components/Supplement_brain_trajectory/brain_spatial_transferred_positions.png",
         dpi = 300,
         bg = "transparent")

# How do UMAP clusters map? ---------------------------------------------------------

colData(laManno_joint_cds)$joint_partition = partitions(laManno_joint_cds)
colData(laManno_joint_cds)$joint_cluster = clusters(laManno_joint_cds)

joint_coldata = 
  colData(laManno_joint_cds) %>%
  as.data.frame() %>%
  left_join(laManno_coldata,
            by = c("Cell" = "CellID",
                   "Age")) %>%
  mutate(dataset = ifelse(sample == "LaManno",
                          "LaManno",
                          "sci-Space"))

ggplot(joint_coldata) +
  geom_point(data = joint_coldata %>%
               dplyr::select(-Age),
             aes(x = umap1,
                 y = umap2),
             stroke = 0,
             size = 0.6) +
  geom_point(aes(x = umap1,
                 y = umap2,
                 color = Age),
             stroke = 0,
             size = 0.525) +
  scale_color_manual(values = viridis(n = 5)[3:5]) +
  facet_wrap(~dataset,ncol = 1) +
  theme_void() +
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Supplement_LaManno_Integration/time_point.png",
         height = 6,
         width = 3,
         bg = "transparent")

# Filter out non-neural clusters  
neural_partitions = 
  c("1","2","3","6","9","10","12","16","17","18","19")

laManno_joint_cds_neural_partitions =
  colData(laManno_joint_cds)  %>%
  as.data.frame() %>%
  filter(joint_partition %in% neural_partitions) %>%
  pull(Cell)


laManno_joint_cds = 
  laManno_joint_cds[,laManno_joint_cds_neural_partitions]

# Filter out non-neural cells explicitly
neural_cells_lamanno = 
  laManno_joint_cds %>%
  colData() %>%
  as.data.frame() %>%
  filter(sample == "LaManno") %>%
  left_join(laManno_coldata,
            by = c("Cell" = "CellID")) %>%
  filter(Class %in% c("Neuron",
                      "Glia",
                      "OPCs",
                      "Radial glia")) %>%
  pull(Cell)

sci_space_cells =
  laManno_joint_cds %>%
  colData() %>%
  as.data.frame() %>%
  filter(sample != "LaManno") %>%
  pull(Cell)
  
laManno_joint_cds = 
  laManno_joint_cds[,union(neural_cells_lamanno,sci_space_cells)]

metadata_for_transfer =
  laManno_joint_cds %>%
  colData() %>%
  as.data.frame() %>%
  left_join(laManno_coldata,
            by = c("Cell" = "CellID")) %>%
  left_join(colData(spatial_cds) %>%
              as.data.frame() %>%
              dplyr::select(everything(),
                            -contains("umap"),
                            -sample,
                            -Size_Factor),
            by = "Cell")


lamanno_coldata =
  metadata_for_transfer %>%
  filter(sample == "LaManno")

lamanno_coldata$Clusters = 
  lamanno_coldata$Clusters %>%
  as.character()

spatial_coldata = 
  metadata_for_transfer %>%
  filter(sample != "LaManno")

# Get nearest neighbors
nn = 
  FNN::get.knnx(data =  
                  lamanno_coldata %>%
                  dplyr::select(umap1, umap2) %>%
                  as.matrix(),
                query = spatial_coldata %>%
                  dplyr::select(umap1, umap2) %>%
                  as.matrix(),
                k = 5)


## Use nearest neighbor indices to get top 5 labels
nn_Cluster = 
  sapply(X = seq(1,dim(nn$nn.index)[2]),
         FUN = function(X){
           lamanno_coldata$Clusters[nn$nn.index[,X]]
         }
  ) %>%
  as.data.frame()

rownames(nn_Cluster) =
  spatial_coldata$Cell

nn_Cluster_majority = 
  sapply(X = seq(1,dim(nn_Cluster)[1]),
         FUN = function(X){
           tabulated_labels = 
             nn_Cluster[X,] %>%
             unlist() %>%
             table() %>%
             sort(decreasing = T)
           
           tabulated_labels[1] %>%
             names()
         })

nn_Cluster_fraction =
  sapply(X = seq(1,dim(nn_Cluster)[1]), 
         FUN = function(Y){
           sum(nn_Cluster[Y,] == nn_Cluster_majority[Y])/5
         })

cluster_df = 
  data.frame(Cell = spatial_coldata$Cell,
           nn_Cluster = nn_Cluster_majority %>% as.character(),
           fraction = nn_Cluster_fraction,
           nearest_neighbor = nn_Cluster[,1]) 

cluster_df$nn_Cluster = 
  ifelse(cluster_df$fraction == 0.2,
         cluster_df$nearest_neighbor  %>% as.character(),
         cluster_df$nn_Cluster  %>% as.character())

cluster_df$nn_Cluster = as.numeric(as.character(cluster_df$nn_Cluster))

spatial_coldata = 
  spatial_coldata %>%
  left_join(cluster_df,
            by = "Cell")


# Supplemental Figure S23 -------------------------------------------------

spatial_coldata = 
  spatial_coldata %>%
  left_join(scaling_df,
            by = "max_slide_id") %>%
  mutate(jittered_y = coords.x1 * y + rnorm(n =dim(spatial_coldata)[1],
                                            mean = 0,
                                            sd = 10),
         jittered_x = coords.x2 * x + rnorm(n =dim(spatial_coldata)[1],
                                            mean = 0,
                                            sd = 10))

ggplot() +
  geom_point(data = spatial_coldata %>%
               filter(max_slide_id == "slide_4E") %>%
               group_by(nn_Cluster) %>%
               filter(n() > 5 ) %>%
               ungroup() %>%
               dplyr::select(-nn_Cluster),
             aes(x = jittered_x,
                 y = jittered_y),
             color = "grey80",
             size = 0.35,
             stroke = 0) +
  geom_point(data = spatial_coldata %>%
               filter(max_slide_id == "slide_4E") %>%
               group_by(nn_Cluster) %>%
               filter(n() >= 5 ) %>%
               ungroup(),
             aes(x = jittered_x,
                 y = jittered_y),
             color = "black",
             size = .6,
             stroke = 0) +
  geom_point(data = spatial_coldata %>%
               filter(max_slide_id == "slide_4E") %>%
               group_by(nn_Cluster) %>%
               filter(n() >= 5 ) %>%
               ungroup(),
             aes(x = jittered_x,
                 y = jittered_y),
             color = "red",
             size = 0.45,
             stroke = 0) +
  facet_wrap(~nn_Cluster) +
  theme_void() +
  ggsave("Figures/Figure_Components/Supplement_LaManno_Integration/slide_4e.png",
         height = 6,
         width = 6)


ggplot() +
  geom_point(data = spatial_coldata %>%
               filter(max_slide_id == "slide_3F") %>%
               group_by(nn_Cluster) %>%
               filter(n() > 5 ) %>%
               ungroup() %>%
               dplyr::select(-nn_Cluster),
             aes(x = jittered_x,
                 y = jittered_y),
             color = "grey80",
             size = 0.35,
             stroke = 0) +
  geom_point(data = spatial_coldata %>%
               filter(max_slide_id == "slide_3F") %>%
               group_by(nn_Cluster) %>%
               filter(n() >= 5 ) %>%
               ungroup(),
             aes(x = jittered_x,
                 y = jittered_y),
             color = "black",
             size = .6,
             stroke = 0) +
  geom_point(data = spatial_coldata %>%
               filter(max_slide_id == "slide_3F") %>%
               group_by(nn_Cluster) %>%
               filter(n() >= 5 ) %>%
               ungroup(),
             aes(x = jittered_x,
                 y = jittered_y),
             color = "red",
             size = 0.45,
             stroke = 0) +
  facet_wrap(~nn_Cluster) +
  theme_void() +
  ggsave("Figures/Figure_Components/Supplement_LaManno_Integration/slide_3f.png",
         height = 6,
         width = 6)


# Calculare the spatial statistic to see whether a cluster is significantly enriched -------------------------
cell_subset =
  spatial_coldata %>% 
  group_by(nn_Cluster) %>%
  filter(n() >= 5) %>% 
  pull(Cell)

localG_cds = 
  spatial_cds[,cell_subset]

colData(localG_cds)$nn_Cluster = 
  colData(localG_cds) %>%
  as.data.frame() %>%
  left_join(spatial_coldata %>%
              dplyr::select(Cell,
                            nn_Cluster),
            by = "Cell") %>%
  pull(nn_Cluster)

localG_cds = 
  localG_cds %>%
  estimate_size_factors() %>%
  preprocess_cds(num_dim = 50) %>%
  reduce_dimension() %>%
  cluster_cells()

# replace UMAP coordinates with the spatial coordinates
reducedDim(x = localG_cds,
           type = "UMAP") <-
  matrix(cbind(colData(localG_cds)$coords.x2, 
               colData(localG_cds)$coords.x1), 
         ncol=2)

localG_df =
  lapply(colData(localG_cds)$max_slide_id %>% unique, 
         FUN = function(curr_slide){
           cds_subset = localG_cds[,colData(localG_cds)$max_slide_id == curr_slide]
           
           df = colData(cds_subset) %>% 
             as.data.frame()
           
           lw <- monocle3:::calculateLW(cds = cds_subset, 
                                        k = 10,
                                        verbose = FALSE,
                                        neighbor_graph = "knn",
                                        reduction_method = "UMAP")
           
           
           wc <- spdep::spweights.constants(lw, zero.policy = TRUE, adjust.n = TRUE)
           
           # apply local G calculation for every variable
           label = "nn_Cluster"
           label_df = lapply(unique(df[,label]), function(x) {
             
             lg <- calculateLocalG(df, lw, wc, label,x)
             
             cells_to_select = 
               df %>%
               filter(nn_Cluster == x) %>%
               pull(Cell)
             
             lg = 
               lg %>%
               filter(Cell %in% cells_to_select)
             
             return(inner_join(df,
                               lg,
                               by = "Cell"))
           })
           
           label_df =
             do.call(rbind,
                     label_df)
           
           label_df = 
             getPval(df = label_df,
                     method = "BH",
                     tail = "twosided")
           
           label_df
         })




localG_df =
  do.call(what = rbind,
          args = localG_df)

sig_clusters =
  localG_df %>%
  group_by(nn_Cluster,
           max_slide_id) %>%
  mutate(sig = (median(qval,na.rm = T) < 0.01)) %>%
  group_by(max_slide_id) %>%
  dplyr::select(nn_Cluster,
                max_slide_id,
                sig) %>%
  distinct() %>%
  summarise(percent_spatial = sum(sig,na.rm = T)/n(),
            n = sum(sig,na.rm = T))

sig_clusters

mean(sig_clusters$percent_spatial)
