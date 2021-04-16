# Align sci-Space data containing E14 embryos to 
# entire mouse organogenesis cell atlas (MOCA) dataset.  

# All computations were run on a computing clustering with sizable amount of
# physical memory owing to the size of the MOCA dataset.

# Downloading MOCA dataset -
# https://oncoscape.v3.sttrcancer.org/atlas.gs.washington.edu.mouse.rna/downloads


# Load startup packages ---------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyverse)
  library(ggrepel)
  library(FNN)
  library(devtools)
  library(monocle3)
  
  DelayedArray:::set_verbose_block_processing(TRUE)
  options(DelayedArray.block.size=1000e7)
  cc.genes = readRDS("/net/trapnell/vol1/home/sanjays/projects/Space/Submission/bin/cc.genes.mouse.RDS")
  
})

setwd(dir="~/projects/Space/published_datasets/Cao_et_al")


append_umap_coordinates = function(cds){
  colData(cds)$umap1 = reducedDim(x = cds,
                                  type = "UMAP")[,1]
  colData(cds)$umap2 = reducedDim(x = cds,
                                  type = "UMAP")[,2]
  return(cds)
}


# Read in Cao -- monocle3 object ------------------------------------------
cds_cao = readRDS("Cao_2019_cds_cleaned_monocle3.RDS")

# Working off of revised annotations 
new_annotations =
  read.table(file = "cell_annotate_20200119.csv",
             sep = ",",
             header = T,
             stringsAsFactors = F)


# Pull in necessary column data -------------------------------------------
colData(cds_cao)$Main_trajectory =
  colData(cds_cao) %>%
  as.data.frame() %>%
  left_join(new_annotations %>% 
              dplyr::select(sample,
                            Main_trajectory),
            by = "sample") %>%
  pull(Main_trajectory.y)


colData(cds_cao)$Sub_trajectory_name =
  colData(cds_cao) %>%
  as.data.frame() %>%
  left_join(new_annotations %>% 
              dplyr::select(sample,
                            Sub_trajectory_name),
            by = "sample") %>%
  pull(Sub_trajectory_name.y)

colData(cds_cao)$Main_cell_type =
  colData(cds_cao) %>%
  as.data.frame() %>%
  left_join(new_annotations %>% 
              dplyr::select(sample,
                            Main_cell_type),
            by = "sample") %>%
  pull(Main_cell_type)

# Read in sci-Space data ------------------------------------------

cds_space = readRDS("/net/trapnell/vol1/home/sanjays/projects/Space/Submission/E14_slides/RDS_intermediates/Notebook3_E14_spatial_CDS.RDS")

# Get the intersecting set of genes between the two datasets
intersecting_genes = 
  intersect(rowData(cds_cao)$gene_id %>%
              as.character(),
            rowData(cds_space)$id)

# Subset MOCA CDS for intersecting genes
cds_cao = cds_cao[intersecting_genes,]

# Select a subset of columns for a combined metadata
coldata_cao = 
  colData(cds_cao) %>%
  as.data.frame() %>%
  dplyr::select(Cell = sample,
                stage = day) %>%
  dplyr::mutate(sample = "Cao")

rownames(coldata_cao) = coldata_cao$Cell

# Subset sci-Space CDS for intersecting genes
cds_space = cds_space[intersecting_genes,]

# Select a subset of columns for a combined metadata
coldata_cds_space = 
  colData(cds_space) %>%
  as.data.frame() %>%
  dplyr::select(Cell, sample) %>%
  mutate(stage = "14")

rownames(coldata_cds_space) = coldata_cds_space$Cell 

# Construct a joint cds
joint_cds =
  new_cell_data_set(expression_data = cbind(counts(cds_cao), 
                                          counts(cds_space)),
                  cell_metadata = rbind(coldata_cao, coldata_cds_space),
                  gene_metadata = rowData(cds_space) %>% 
                    as.data.frame())


# Process joint MOCA/Space CDS --------------------------------------------
colData(joint_cds)$n.umi.log = 
  log10(Matrix::colSums(counts(joint_cds)))

colData(joint_cds)$experiment =
  ifelse(colData(joint_cds)$sample == "Cao",
         "Cao",
         ifelse(colData(joint_cds)$sample == "1",
                "exp1",
                "exp2"))
joint_cds = 
  joint_cds %>%
  estimate_size_factors() %>%
  preprocess_cds(num_dim = 100) %>%
  # Regress out the effect of molecules per cell and perform MNN alignment
  # based on the dataset
  align_cds(residual_model_formula_str = "~n.umi.log",
            alignment_group = "experiment") %>%
  reduce_dimension(verbose = T)

joint_cds = 
  cluster_cells(cds = joint_cds,
                reduction_method = "UMAP",
                verbose = T)

# Add UMAP coordinates to CDS metadata
joint_cds =
  append_umap_coordinates(joint_cds)

# Add clustering done in UMAP space to CDS metadata
colData(joint_cds)$umap_cluster  = 
  clusters(joint_cds, 
           reduction_method="UMAP")

# Add partition from graph abstraction to CDS metadata
colData(joint_cds)$partition = partitions(joint_cds)


# Create a dataframe for plotting purposes
coldata = 
  colData(joint_cds) %>%
  as.data.frame()

# Re-level stages to be in developmental order 
colData(joint_cds)$stage =
  factor(colData(joint_cds)$stage,
            levels = c("9.5",
                       "10.5",
                       "11.5",
                       "12.5",
                       "13.5",
                       "14"))

# pull in annotations from MOCA coldata and Space coldata
coldata = 
  coldata %>%
  left_join(colData(cds_cao) %>%
              as.data.frame() %>% 
              dplyr::select(Cell = sample,
                            Main_trajectory,
                            Sub_trajectory_name,
                            Main_cell_type),
            by = "Cell") 

# Figure 1 — Panel B ------------------------------------------------------

# Choose a set of trajectories to highlight
highlighted_trajectory = 
  c( "Neural tube and notochord trajectory",
     "Mesenchymal trajectory",
     "Epithelial trajectory",
     "Endothelial trajectory",
     "Haematopoiesis trajectory",
     "Hepatocyte trajectory",
     "Neural crest 1",
     "Neural crest 2",
     "Hepatocyte trajectory")

label_position_df = 
  coldata %>%
  filter(Main_trajectory %in% highlighted_trajectory) %>% 
  group_by(Main_trajectory) %>%
  summarise(med_umap1 = median(umap1),
         med_umap2 = median(umap2))

ggplot(coldata) +
  geom_point(data = coldata,
             aes(x = umap1,
                 y = umap2,
                 color = stage),
             size = 0.15,
             stroke = 0) +
  geom_label_repel(data = label_position_df,
                   aes(x = med_umap1,
                       y = med_umap2,
                       label = Main_trajectory),
                  size = 1.5,
                  label.size = .15,
                  alpha = 0.65,
                  seed = 42,
                  box.padding = 0.35,
                  point.padding = 0.25) +
  geom_label_repel(data = label_position_df,
                  aes(x = med_umap1,
                      y = med_umap2,
                      label = Main_trajectory),
                  size = 1.5,
                  label.size = .15,
                  fill = NA,
                  seed = 42,
                  box.padding = 0.35,
                  point.padding = 0.25) +
  theme_void() +
  scale_color_viridis_d() +
  theme(legend.position = "none",
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  ggsave("figures/umap_aligned_all_cells_space68_stage.png", 
         height = 2.5, 
         width = 2.5,
         dpi = 600)

saveRDS(object = joint_cds, 
        file = "bin/joint_space_cao_cds.RDS")