# Mapping cells to segemented ROIs. Some of the organs and their boundaries were
# discernible from the DAPI stained section and aided through immunostaining of 
# the adjacent section. In this notebook cardiomyocytes mapping outside the heart
# and hepatocytes outside of the liver are marked for removal

# Load startup packages ---------------------------------------------------
suppressPackageStartupMessages({ 
  library(tidyverse)
  library(ggplot2)
  library(tidyr)
  library(viridis)
  library(purrr)
  library(OpenImageR)
  library(spatstat)
  library(imager)
  library(vec2dtransf)
  library(sp)
  library(sf)
  library(monocle3)
  library(magrittr)
  
  space_directory = "/Volumes/GoogleDrive/My Drive/sciSpace/"
  setwd(dir=space_directory)
  source("Submission_Data/bin/hotspot_functions.R")
  
  # Pass TRUE if you want to see progress output on some of Monocle 3's operations
  DelayedArray:::set_verbose_block_processing(TRUE)
  # Passing a higher value will make some computations faster but use more memory. 
  # Adjust with caution!
  options(DelayedArray.block.size=1000e6)
  
})

load("Submission_Data/E14_slides/RDS_intermediates/Notebook5_initial_mapping.RData")

all_image_data = 
  readRDS("Submission_Data/E14_slides/RDS_intermediates/Notebook0_1_images_and_transformations.RDS")

updated_coldata = 
  readRDS(file = "Submission_Data/E14_slides/RDS_intermediates/Notebook5.4_updated_coldata.RDS")


# Read in the outlines of the organs and creat polygons -------------------
all_celltypes = 
  readRDS("Submission_Data/E14_slides/RDS_intermediates/Notebook5.2_celltype_hulls.RDS")

cells_to_remove =
  updated_coldata %>%
  group_by(final_cluster_label,
           max_slide_id) %>%
  nest() %>%
  inner_join(all_celltypes,
             by = c("max_slide_id" = "slide_id",
                    "final_cluster_label" = "label")) %>%
  filter(final_cluster_label %in% c("Cardiac muscle lineages",
                                       "Hepatocytes")) %>%
  mutate(in_mask = purrr::map2(.x = data,
                               .y = celltype_polygon,
                               .f = function(coldata_subset,
                                             this_mask){
                                 
                                 coldata_subset %>%
                                   dplyr::select(coords.x1,
                                                 coords.x2) %>%
                                   as.matrix() %>%
                                   SpatialPoints() %>%
                                   over(as(object = this_mask,
                                           Class = "Spatial")) %>%
                                   is.na() %>%
                                   not()
                                 
                               })
  ) %>%
  tidyr::unnest(cols = c(data,in_mask)) %>%
  dplyr::select(Cell,
                max_slide_id,
                final_cluster_label,
                in_mask) %>%
  filter(!in_mask) %>%
  pull(Cell)


spatial_cds = 
  readRDS(file = "Submission_Data/E14_slides/RDS_intermediates/Notebook5_E14_spatial_CDS.RDS")

rownames(updated_coldata) = updated_coldata$Cell

updated_coldata = 
  updated_coldata[as.character(colnames(spatial_cds)),]

identical(rownames(updated_coldata),
          colnames(spatial_cds))

for (i in (setdiff(colnames(updated_coldata),colnames(colData(spatial_cds))))){
  colData(spatial_cds)[[i]] = updated_coldata[[i]]
} 

colData(spatial_cds)$cell_to_remove = 
  colData(spatial_cds)$Cell %in% cells_to_remove

# spatial_cds = 
#   spatial_cds[,!(colnames(spatial_cds) %in% cells_to_remove)]

head(colData(spatial_cds))

saveRDS(object = spatial_cds,
        file = "Submission_Data/E14_slides/RDS_intermediates/Notebook5.5_E14_spatial_CDS.RDS")






