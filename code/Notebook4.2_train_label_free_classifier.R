# Using annotations in the mouse organogenesis cell atlas (MOCA) dataset, we trained
# a garnett classifier.

# https://github.com/cole-trapnell-lab/garnett/tree/marker_free
# marker_free branch
# All computations were run on a computing clustering 

# Load startup packages ---------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyverse)
  library(ggrepel)
  library(FNN)
  library(devtools)
  library(monocle3)
  #load_all("~/bin/monocle3_9becd94f60930c2a9b51770e3818c194dd8201eb/monocle3/")
  load_all("~/bin/garnett/")
  DelayedArray:::set_verbose_block_processing(TRUE)
  options(DelayedArray.block.size=1000e7)
  
})

setwd(dir="~/projects/Space/published_datasets/")


append_umap_coordinates = function(cds){
  colData(cds)$umap1 = reducedDim(x = cds,
                                  type = "UMAP")[,1]
  colData(cds)$umap2 = reducedDim(x = cds,
                                  type = "UMAP")[,2]
  return(cds)
}

# Read in Cao -- monocle3 object ------------------------------------------
# Read in Cao data that is a monocle3 object
cds_cao = readRDS("Cao_et_al/Cao_2019_cds_cleaned_monocle3.RDS")

new_annotations =
  read.table(file = "Cao_et_al/cell_annotate_20200119.csv",
             sep = ",",
             header = T,
             stringsAsFactors = F)

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


cds_cao_e13.5 = 
  cds_cao[,colData(cds_cao)$day == "13.5"]

# Train a classifier on this data -----------------------------------------
cao_e13.5_classifier = 
  train_cell_classifier(cds = cds_cao_e13.5,
                        marker_file = "Cao_et_al/garnet_label_free.marker_file",
                        db="none")




# Read in sci-Space data ------------------------------------------
cds_space = readRDS("/net/trapnell/vol1/home/sanjays/projects/Space/Submission/E14_slides/RDS_intermediates/Notebook3_E14_spatial_CDS.RDS")


# Use trained classifier on sci-space data --------------------------------
cds_space = 
  classify_cells(cds_space, 
                 cao_e13.5_classifier,
                 db = "none",
                 cluster_extend = TRUE)


# Save result -------------------------------------------------------------
write.table(x = colData(cds_space) %>%
              as.data.frame(),
            file = "garnett_classifications.tsv",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)
