# Add nearest neighbor and garnett annotations to the CDS along with manual annotations. 
# Plot concordance between the two measures and reconcile cells that are found at a very low 
# percentage in subclusters

# Load startup packages ---------------------------------------------------
suppressPackageStartupMessages({ 
  library(tidyverse)
  library(ggplot2)
  library(tidyr)
  library(viridis)
  library(ggridges)
  library(RColorBrewer)
  library(ggrepel)
  library(monocle3)
  library(garnett)
  library(pheatmap)
  
  space_directory = "~/Google Drive File Stream/My Drive/sciSpace/"
  setwd(dir=space_directory)
  source("Submission_Data/bin/cell_cycle.R")
  cc.genes <- readRDS("Submission_Data/bin/cc.genes.mouse.RDS")
  
  # Pass TRUE if you want to see progress output on some of Monocle 3's operations
  DelayedArray:::set_verbose_block_processing(TRUE)
  # Passing a higher value will make some computations faster but use more memory. 
  # Adjust with caution!
  options(DelayedArray.block.size=1000e6)
  
  # Set a seed to make umap and other non-deterministic steps consistent
  set.seed(seed = 42)
  
})

spatial_cds = readRDS("Submission_Data/E14_slides/RDS_intermediates/Notebook3_E14_spatial_CDS.RDS")

# Read in the nearest neighbor assignments
alignment_coldata = 
  read.table(file = "Submission_Data/E14_slides/RDS_intermediates/alignment_coldata_seurat.tsv",
             sep = "\t",
             header = T) 


# Incorporate nearest neighbor annotations --------------------------------
colData(spatial_cds)$nn_label = 
  colData(spatial_cds) %>%
  as.data.frame() %>%
  left_join(alignment_coldata %>%
              dplyr::select(Cell,
                    nn_majority),
            by = "Cell") %>%
  pull(nn_majority) %>%
  as.character()

colData(spatial_cds)$nn_label_original = colData(spatial_cds)$nn_label
colData(spatial_cds)$lamanno_Class[colData(spatial_cds)$nn_label == "Endothelial cells"] <- "Endothelial Cells"


# Read in the nearest neighbor assignments from La Manno et. al.
alignment_coldata_lamanno = 
  read.table(file = "Submission_Data/E14_slides/RDS_intermediates/alignment_coldata_lamanno.tsv",
             sep = "\t",
             header = T) 

alignment_coldata_lamanno = 
  alignment_coldata_lamanno %>%
  filter(!(nn_class_majority %in% c("Blood","Vascular")))

colData(spatial_cds)$lamanno_Class = 
  colData(spatial_cds) %>%
  as.data.frame() %>%
  left_join(alignment_coldata_lamanno %>%
              dplyr::select(Cell,
                            nn_class_majority),
            by = "Cell") %>%
  pull(nn_class_majority) %>%
  as.character()

# Rename cell types for consistency
colData(spatial_cds)$lamanno_Class[colData(spatial_cds)$lamanno_Class == "Schwann"] <- "Schwann Cells"
colData(spatial_cds)$lamanno_Class[colData(spatial_cds)$lamanno_Class == "Glia"] <- "Glial Cells"
colData(spatial_cds)$lamanno_Class[colData(spatial_cds)$lamanno_Class == "Choroid plexus"] <- "Choroid Plexus"


# Transfer labels from the DMBA -------------------------------------------
colData(spatial_cds)$lamanno_Punchcard = 
  colData(spatial_cds) %>%
  as.data.frame() %>%
  left_join(alignment_coldata_lamanno %>%
              dplyr::select(Cell,
                            nn_punchcard_majority),
            by = "Cell") %>%
  pull(nn_punchcard_majority) %>%
  as.character()

colData(spatial_cds)$lamanno_Tissue = 
  colData(spatial_cds) %>%
  as.data.frame() %>%
  left_join(alignment_coldata_lamanno %>%
              dplyr::select(Cell,
                            nn_tissue_majority),
            by = "Cell") %>%
  pull(nn_tissue_majority) %>%
  as.character()


colData(spatial_cds)$nn_label =
  ifelse(is.na(colData(spatial_cds)$lamanno_Class),
         colData(spatial_cds)$nn_label,
         colData(spatial_cds)$lamanno_Class)


# Calculate if a label is in the minority in a cluster --------------------

minors = 
  colData(spatial_cds) %>% 
  as.data.frame() %>%
  group_by(sub_cluster,
           nn_label) %>%
  summarise(n = n()) %>%
  group_by(sub_cluster) %>%
  mutate(fraction = n/sum(n)) %>%
  dplyr::select(sub_cluster,
                nn_label,
                fraction,
                n)%>% 
  filter(fraction < 0.25) %>%
  dplyr::select(sub_cluster, 
                nn_label) %>%
  mutate(minor = T)


# Get the majority of each cluster
majors = 
  colData(spatial_cds) %>% 
  as.data.frame() %>%
  mutate(sub_cluster = as.character(sub_cluster)) %>%
  group_by(sub_cluster,
           nn_label) %>%
  summarise(n = n()) %>%
  group_by(sub_cluster) %>%
  mutate(fraction = n/sum(n)) %>%
  dplyr::select(sub_cluster,
                cluster_majority = nn_label,
                fraction) %>%
  top_n(n = 1) %>%
  group_by(sub_cluster) %>%
  sample_n(size = 1)


# If a cell is in the minority label in that cluster, assign it to the 
# majority label. This is basically using the cluster to smooth the labels
colData(spatial_cds)$minor = 
  colData(spatial_cds) %>% 
  as.data.frame() %>%
  left_join(minors,
            by = c("sub_cluster","nn_label")) %>%
  pull(minor)


colData(spatial_cds)$major = 
  colData(spatial_cds) %>% 
  as.data.frame() %>%
  left_join(majors,
            by = c("sub_cluster")) %>% 
  pull(cluster_majority) 


colData(spatial_cds)$cluster_extended_nnlabel =
  ifelse(is.na(colData(spatial_cds)$minor),
          colData(spatial_cds)$nn_label,
          colData(spatial_cds)$major)


# Read in manually refined annotations ------------------------------------
refined_annotations = 
  read.table(file = "Submission_Data/E14_slides/RDS_intermediates/manual_cell_annotations_editted.txt",
             sep = "\t",
             header = T,
             colClasses = rep(x = "character",10),
             fill = NA) %>%
  dplyr::select(cluster,
                manual_annotation,
                manual_annotation_2) %>% 
  distinct() %>%
  mutate(cluster = as.character(cluster))

colData(spatial_cds)$manual_annotation = 
  colData(spatial_cds) %>%
  as.data.frame() %>%
  mutate(cluster = as.character(cluster)) %>%
  left_join(refined_annotations,
            by = "cluster") %>%
  pull(manual_annotation) %>%
  as.character()

colData(spatial_cds)$manual_annotation_2 = 
  colData(spatial_cds) %>%
  as.data.frame() %>%
  mutate(cluster = as.character(cluster)) %>%
  left_join(refined_annotations,
            by = "cluster") %>%
  pull(manual_annotation_2) %>%
  as.character()

# The final cluster label is used throughout the study
colData(spatial_cds)$final_cluster_label = 
  ifelse(!is.na(colData(spatial_cds)$manual_annotation),
         colData(spatial_cds)$manual_annotation,
         colData(spatial_cds)$lamanno_Class)

colData(spatial_cds)$final_cluster_label = 
  ifelse(!is.na(colData(spatial_cds)$manual_annotation),
         colData(spatial_cds)$manual_annotation,
         colData(spatial_cds)$lamanno_Class)

colData(spatial_cds)$final_cluster_label = 
  ifelse(colData(spatial_cds)$cluster == "34",
         colData(spatial_cds)$cluster_extended_nnlabel,
         colData(spatial_cds)$final_cluster_label)

majors_2 = 
  colData(spatial_cds) %>% 
  as.data.frame() %>%
  mutate(cluster = as.character(cluster)) %>%
  group_by(cluster,
           final_cluster_label) %>%
  summarise(n = n()) %>%
  group_by(cluster) %>%
  mutate(fraction = n/sum(n)) %>%
  dplyr::select(cluster,
                cluster_majority = final_cluster_label,
                fraction) %>%
  top_n(n = 1) %>%
  group_by(cluster) %>%
  sample_n(size = 1)

colData(spatial_cds)$final_cluster_label = 
  colData(spatial_cds) %>%
  as.data.frame() %>%
  mutate(cluster = as.character(cluster)) %>%
  left_join(majors_2 %>%
              dplyr::rename(cluster_majority2 = cluster_majority),
            by = "cluster") %>%
  mutate(final_cluster_label = ifelse(is.na(final_cluster_label),
                                      cluster_majority2,
                                      final_cluster_label)) %>%
  pull(final_cluster_label)

colData(spatial_cds)$final_cluster_label[colData(spatial_cds)$final_cluster_label == "Endothelial cells"] <- "Endothelial Cells"


# Supplemental Figure 9B --------------------------------------------------
ggplot() +
  geom_point(data = 
               rbind(colData(spatial_cds) %>% as.data.frame() %>% mutate(order = "1"),
                     colData(spatial_cds) %>% as.data.frame() %>% mutate(order = "2")) %>%
               arrange(cluster,order) %>%
               mutate(cluster = ifelse(order ==2,
                                       cluster %>% as.character(),
                                       NA)),
             aes(x = umap1,
                 y = umap2,
                 color = cluster,
                 size = order),
             stroke = 0) +
  scale_size_manual(values = c("1" = 0.35,
                               "2" = 0.25)) +
  theme_void() +
  theme(legend.position = "none") + 
  scale_color_manual(values = cluster_colors,
                     na.value="black") +
  ggsave("Figures/Figure_Components/Supplement_experiment_QC/umap_by_cluster.png",
         dpi = 300,
         height = 2.5,
         width = 2.5)

# Figure 1C --------------------------------------------------


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

ggplot() +
  geom_point(data = 
               rbind(colData(spatial_cds) %>% as.data.frame() %>% mutate(order = "1"),
                     colData(spatial_cds) %>% as.data.frame() %>% mutate(order = "2")) %>%
               arrange(final_cluster_label,order) %>%
               mutate(final_cluster_label = ifelse(order ==2,
                                                   final_cluster_label %>% as.character(),
                                                   NA)),
             aes(x = umap1,
                 y = umap2,
                 color = final_cluster_label,
                 size = order),
             stroke = 0) +
  scale_size_manual(values = c("1" = 0.45,
                               "2" = 0.35)) +
  theme_void() +
  theme(legend.position = "none",
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_color_manual(values = colors,
                     na.value="black")+  
  ggsave("Figures/Figure_Components/Figure1/umap_nnlabel_unknown.png", 
         height = 4.5, 
         width = 4.5,
         dpi = 600)


# Marker Free Garnett Annotations -----------------------------------------

# Read in classifications
label_free_garnett_calls =
  read.table("Submission_Data/E14_slides/RDS_intermediates/garnett_classifications.tsv",
             sep = "\t",
             header = T) %>%
  dplyr::select(Cell,
                cell_type,
                cluster_ext_type)

# Add classifications to metadata
colData(spatial_cds)$garnett_cell_type =
  colData(spatial_cds) %>%
  as.data.frame() %>%
  left_join(label_free_garnett_calls %>%
              dplyr::select(Cell,
                            cell_type),
            by = "Cell") %>%
  pull(cell_type) %>%
  as.character()

colData(spatial_cds)$garnett_cluster_ext_type =
  colData(spatial_cds) %>%
  as.data.frame() %>%
  left_join(label_free_garnett_calls %>%
              dplyr::select(Cell,
                            cluster_ext_type),
            by = "Cell") %>%
  pull(cluster_ext_type) %>%
  as.character()


# Make a matrix to compare the two annotation strategies
matrix_for_heatmap =
  colData(spatial_cds) %>%
  as.data.frame() %>%
  group_by(garnett_cluster_ext_type) %>%
  add_tally(name = "num_garnett_cell_type") %>%
  ungroup() %>%
  group_by(garnett_cluster_ext_type, nn_label_original) %>%
  mutate(percent_nn_in_garnett = n()/num_garnett_cell_type) %>%
  dplyr::select(nn_label_original, garnett_cluster_ext_type, percent_nn_in_garnett) %>%
  drop_na() %>%
  distinct() %>% 
  spread(key =nn_label_original, value =  percent_nn_in_garnett, fill = 0) 

matrix_for_heatmap =
  matrix_for_heatmap %>%
  tibble::column_to_rownames(var = "garnett_cluster_ext_type") %>%
  as.matrix()


# Supplementary Figure 11B --------------------------------------------------

pheatmap(matrix_for_heatmap,
         cellwidth = 4,
         cellheight = 4,
         treeheight_row = 0,
         treeheight_col = 0,
         cluster_rows = F,
         cluster_cols = F,
         fontsize_row = 4,
         fontsize_col = 4,legend = F,
         filename = "Figures/Figure_Components/Supplement_experiment_QC/percent_nn_in_garnett.pdf",
         color = viridis(option = "viridis", n = 40)) 

# Write out resulting CDS -------------------------------------------------
saveRDS(object = spatial_cds,
        file = "Submission_Data/E14_slides/RDS_intermediates/Notebook5_E14_spatial_CDS.RDS")

