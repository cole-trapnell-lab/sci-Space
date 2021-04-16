# Adding in anatomical annotations from manually segmented using the 
# Allen Brain Atlas into 1 of 6 regions: Pallium, Sub Pallium, Thalamus, Hypothalamus, Midbrain and Hindbrain

# This notebook takes those segmentations and records which cells land in each segmentation

# Load startup packages ---------------------------------------------------
suppressPackageStartupMessages({ 
  library(tidyverse)
  library(ggplot2)
  library(tidyr)
  library(viridis)
  library(ggridges)
  library(purrr)
  library(spatstat)
  library(imager)
  library(vec2dtransf)
  library(sp)
  library(sf)
  library(ggrepel)
  library(pheatmap)
  library(monocle3)
  
  space_directory = "~/Google Drive File Stream/My Drive/sciSpace/"
  setwd(dir=space_directory)
  source("Submission_Data/bin/hotspot_functions.R")
  
  # Pass TRUE if you want to see progress output on some of Monocle 3's operations
  DelayedArray:::set_verbose_block_processing(TRUE)
  # Passing a higher value will make some computations faster but use more memory. 
  # Adjust with caution!
  options(DelayedArray.block.size=1000e6)
  
})

all_image_data = 
  readRDS("Submission_Data/E14_slides/RDS_intermediates/Notebook0_1_images_and_transformations.RDS")


spatial_cds = 
  readRDS(file = "Submission_Data/E14_slides/RDS_intermediates/Notebook6_spatial_cds_anatomy.RDS")


# Read in the outlines of the organs and create polygons -------------------

# Read in contours that create sf spatial map objects
# Note points are in the the raw pixel space prior to flipping the y axis and scaling
# Need to join the countours with the image_locations_and_sizes data frame

brain_hulls_path = "Submission_Data/E14_slides/Images/brain_hulls/"
brain_hulls_coordinates = list.files(path = brain_hulls_path)
brain_hulls_coordinates =
  brain_hulls_coordinates[grepl(x = brain_hulls_coordinates,
                              pattern = ".csv")]

all_brain_hulls = list()
for (file_name in brain_hulls_coordinates) {
  curr_slide =
    stringr::str_sub(string = file_name,
                     start = 1,
                     end = 8)
  
  label =
    stringr::str_sub(string = file_name,
                     start = 10) %>%
    stringr::str_replace_all(".csv",
                             "") %>%
    stringr::str_replace_all("-",
                             " ")
  
  all_brain_hulls[[file_name]] = c(file_name, curr_slide, label)
  
}


all_brain_hulls =
  do.call(rbind, all_brain_hulls) %>%
  as.data.frame()

colnames(all_brain_hulls) = c("file_name",
                            "slide_id",
                            "label")

# Create a matrix from coordinates that can be used to make a polygon
all_brain_hulls =
  all_brain_hulls %>%
  group_by(file_name,
           slide_id,
           label) %>%
  nest() %>%
  mutate(celltype_data_frame =
           purrr::map2(.x = data,
                       .y = file_name,
                       .f = function(subset,
                                     this_file){
                         df = read.table(file = paste0(brain_hulls_path,
                                                       this_file,
                                                       sep = ""),
                                         sep = ",")[,c(3,4)] %>%
                           as.data.frame()
                         colnames(df) = c("x","y")
                         
                         df
                         
                       })) %>%
  dplyr::select(-data)


# Read in image locations -------------------------------------------------

landmarks_path = "Submission_Data/E14_slides/Images/Alignment_Landmarks/"

image_locations_and_sizes =
  read.table(paste0(landmarks_path,
                    "filenames_and_sizes.csv",
                    sep = ""),
             sep = ",",
             col.names = c("sybr_path",
                           "dapi_path",
                           "dapi_sybr_path",
                           "slide_id",
                           "starting_width_image",
                           "starting_height_image",
                           "ending_width_image",
                           "ending_height_image",
                           "starting_width_grid",
                           "starting_height_grid",
                           "ending_width_grid",
                           "ending_height_grid"))

# Scale and organs just as the images were scaled 
all_brain_hulls =
  all_brain_hulls %>%
  left_join(image_locations_and_sizes,
            by = "slide_id") %>%
  group_by(slide_id,
           file_name,
           label) %>%
  nest() %>%
  mutate(brain_hull_polygon =
           purrr::map(.x = data,
                      .f = function(subset){
                        curr_df =
                          subset$celltype_data_frame %>%
                          as.data.frame()
                        
                        mat =
                          data.frame(
                            x_scaled = curr_df$x / (subset$starting_width_image/subset$ending_width_image),
                            y_scaled = subset$ending_height_image - (curr_df$y/ (subset$starting_height_image/subset$ending_height_image))
                          ) %>% as.matrix()
                        
                        mat = rbind(mat,
                                    mat[1,])
                        
                        st_polygon(list(mat))
                        
                      })) %>%
  ungroup() %>%
  dplyr::select(slide_id,
                label,
                brain_hull_polygon)


all_brain_hulls$slide_id = 
  all_brain_hulls$slide_id %>%
  stringr::str_replace(pattern = "S",
                       replacement = "s")


organs_to_test = all_brain_hulls$label %>% unique() %>% as.character()

# Get the cells that overlap each annotated organ
annotation_mapped_cells = 
  lapply(X = seq(1,dim(all_brain_hulls)[1]),
         function(X){
           this.slide = all_brain_hulls$slide_id[X]
           
           this.organ = all_brain_hulls$label[X]
           this.polygon = all_brain_hulls$brain_hull_polygon[X][[1]]
           
           cells_in_this_section =
             spatial_cds %>%
             colData() %>%
             as.data.frame() %>%
             filter(max_slide_id == this.slide) %>%
             filter(!is.na(coords.x1),
                    !is.na(coords.x2)) 
           
           
           cell_positions = 
             cells_in_this_section %>%
             dplyr::select(coords.x1,coords.x2) %>%
             as.matrix()
           
           
           cell_positions = 
             SpatialPoints(cell_positions)

           polygon_hull =
             this.polygon %>%
             as(Class = "Spatial")
           
           mask = 
             sp::over(cell_positions,
                      polygon_hull)
           
           cells_in_this_section$in_organ = 
             !(is.na(mask))
           
           cells_in_this_section_and_organ = 
             cells_in_this_section %>%
             filter(in_organ) %>%
             mutate(organ = this.organ) %>%
             dplyr::select(Cell,
                           max_slide_id,
                           final_cluster_label,
                           organ)
           
           cells_in_this_section_and_organ
         })

annotation_mapped_cells = 
  do.call(rbind, annotation_mapped_cells)

# Some cells mapped to two anatomical locations due to overlapping segmentation
# Choose one of the two annotations for these cells
annotation_mapped_cells =
  annotation_mapped_cells %>%
  group_by(Cell) %>%
  sample_n(1)

colData(spatial_cds)$brain_region = 
  colData(spatial_cds) %>%
  as.data.frame() %>%
  left_join(annotation_mapped_cells %>%
              dplyr::select(Cell,
                            organ),
            by = "Cell") %>%
  pull(organ)


# Write out files ---------------------------------------------------------

saveRDS(object = spatial_cds,
        file = "Submission_Data/E14_slides/RDS_intermediates/Notebook6.01_spatial_cds_anatomy.RDS")
