# Position refinement 1. ROI for each cell type based on broad area of expectation
# Cell is moved if found outside the ROI and there is enough sequencing support to 
# move the cell -- The ratio of the top convolved spot in the ROI is within 5 fold 
# of the top convolved spot outside the ROI

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

  # Pass TRUE if you want to see progress output on some of Monocle 3's operations
  DelayedArray:::set_verbose_block_processing(TRUE)
  # Passing a higher value will make some computations faster but use more memory. 
  # Adjust with caution!
  options(DelayedArray.block.size=1000e6)
  
})

load("Submission_Data/E14_slides/RDS_intermediates/Notebook5_initial_mapping.RData")

all_image_data = 
  readRDS("Submission_Data/E14_slides/RDS_intermediates/Notebook0_1_images_and_transformations.RDS")


spatial_cds = 
  readRDS(file = "Submission_Data/E14_slides/RDS_intermediates/Notebook5_E14_spatial_CDS.RDS")


updated_coldata = 
  readRDS(file = "Submission_Data/E14_slides/RDS_intermediates/Notebook5.1_updated_coldata.RDS")


updated_coldata = 
  colData(spatial_cds) %>%
  as.data.frame() %>%
  left_join(updated_coldata,
            by = "Cell")

# Read in the outlines of the organs and creat polygons -------------------

# Read in contours that create sf spatial map objects
# Note points are in the the raw pixel space prior to flipping the y axis and scaling
# Need to join the countours with the image_locations_and_sizes data frame

# Build ROIs
celltypes_path = "Submission_Data/E14_slides/Images/celltype_hulls/"
celltypes_coordinates = list.files(path = celltypes_path)
celltypes_coordinates =
  celltypes_coordinates[grepl(x = celltypes_coordinates,
                          pattern = ".csv")]

all_celltypes = list()
for (file_name in celltypes_coordinates) {
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
  
  all_celltypes[[file_name]] = c(file_name, curr_slide, label)

}


all_celltypes =
  do.call(rbind, all_celltypes) %>%
  as.data.frame()

colnames(all_celltypes) = c("file_name",
                        "slide_id",
                        "label")

# Create a matrix from coordinates that can be used to make a polygon
all_celltypes =
  all_celltypes %>%
  group_by(file_name,
           slide_id,
           label) %>%
  nest() %>%
  mutate(celltype_data_frame =
           purrr::map2(.x = data,
                       .y = file_name,
                      .f = function(subset,
                                    this_file){
                        df = read.table(file = paste0(celltypes_path,
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
all_celltypes =
  all_celltypes %>%
  left_join(image_locations_and_sizes,
            by = "slide_id") %>%
  group_by(slide_id,
           file_name,
           label) %>%
  nest() %>%
  mutate(celltype_polygon =
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
                celltype_polygon)


all_celltypes$slide_id = 
  all_celltypes$slide_id %>%
  stringr::str_replace(pattern = "S",
                       replacement = "s")

  
saveRDS(all_celltypes,
        file = "Submission_Data/E14_slides/RDS_intermediates/Notebook5.2_celltype_hulls.RDS")

cells_to_remap =
  updated_coldata %>%
  group_by(final_cluster_label,
           max_slide_id) %>%
  nest() %>%
  inner_join(all_celltypes,
             by = c("max_slide_id" = "slide_id",
                    "final_cluster_label" = "label")) %>%
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

#16,145 cells 


# Iterate through Cells and find the top combination of spot and sector --------
remap_list = list()
for (i in seq(1,4)){
  
  current_cells = 
    intersect(rownames(sector.hashes.list[[i]]), 
              rownames(spot.hashes.list[[i]])) %>%
    intersect(updated_coldata %>%
                pull(Cell)) %>%
    intersect(cells_to_remap) %>%
    as.character()
  
  print(paste0("Starting iteration ",
               i,
               " now",
               sep =""))
  
  temp = 
    lapply(current_cells, 
           function(x,
                    coldata_to_test,
                    sector_mat,
                    spot_mat,
                    image_data,
                    mask){
             
             this_cell = as.character(x)
             this_slide_id = 
               coldata_to_test[which(coldata_to_test$Cell == this_cell),"max_slide_id"] %>%
               as.character()
             
             # Get the spots which are in the possible universe of spots
             this_slide_layout = 
               image_data[which(image_data$slide_id == this_slide_id),
                          "transformed_oligo_layout"][[1]][[1]]
             
             this_celltype = 
               coldata_to_test[which(coldata_to_test$Cell == this_cell),"final_cluster_label"] %>%
               as.character()
             
             this_mask = 
               mask[intersect(which(mask$slide_id == this_slide_id),
                              which(mask$label == this_celltype)), "celltype_polygon"][[1]][[1]]
             
             this_slide_layout$in_mask = 
               this_slide_layout %>%
               dplyr::select(coords.x1,
                             coords.x2) %>%
               as.matrix() %>%
               SpatialPoints() %>%
               over(as(object = this_mask,
                       Class = "Spatial")) %>%
               is.na() %>%
               not()
             
             
             matching_spots = 
               this_slide_layout %>%
               filter(in_hull,
                      in_mask) %>%
               left_join(data.frame(value_spot =as.numeric(spot_mat[this_cell,]),
                                    hash = colnames(spot_mat)), 
                         by = "hash")

             mat = 
               matching_spots %>%
               dplyr::select(Row,Col,value_spot)
             
             mat = Matrix::sparseMatrix(i = mat$Row, 
                                        j = mat$Col, 
                                        x = mat$value_spot,
                                        dims = c(84,84)) %>%
               Matrix::as.matrix()
             
             res =
               OpenImageR::convolution(mat,
                           kernel = gaussian_blur_3x3,
                           mode = "same")
             
             matching_spots = 
               reshape2::melt(res,varnames = c("Row","Col")) %>%
               dplyr::select(everything(),
                             blurred_spot_value = value) %>%
               right_join(matching_spots,
                          by = c("Row",
                                 "Col"))
               
             # Choose the top spot 
             top_spot = 
               matching_spots %>%
               arrange(desc(blurred_spot_value * value_spot)) %>%
               mutate(Cell = this_cell) %>%
               head(n = 1)  
             
             return(top_spot)
           }, 
           updated_coldata,
           sector.hashes.list[[i]],
           spot.hashes.list[[i]],
           all_image_data,
           all_celltypes)
  
  remap_list[[i]] = 
    do.call(rbind,temp)
}


# Make a data frame from the result
remap_df.list = list()
for( i in seq(1,4)){
  remap_df.list[[i]] =
    data.frame(row.names = remap_list[[i]]$Cell,
               Cell = remap_list[[i]]$Cell,
               top_spot.mask = as.character(remap_list[[i]]$Oligo),
               Row.mask = as.character(remap_list[[i]]$Row),
               Col.mask = as.character(remap_list[[i]]$Col),
               blurred_spot_value.mask = remap_list[[i]]$blurred_spot_value,
               coords.x1.mask = remap_list[[i]]$coords.x1,
               coords.x2.mask = remap_list[[i]]$coords.x2)
  }

remapped_spots =
  do.call(rbind,remap_df.list)


# Only write positions that are within 5 fold of the top spot outside the ROI
new_positions_post_mask =
  inner_join(updated_coldata,remapped_spots) %>% 
  filter((blurred_spot_value/blurred_spot_value.mask <= 5))

# dim(new_positions_post_mask)
# [1] 11614    46

new_positions_post_mask$top_spot = new_positions_post_mask$top_spot.mask %>% as.character()
new_positions_post_mask$Row = new_positions_post_mask$Row.mask %>% as.character()
new_positions_post_mask$Col = new_positions_post_mask$Col.mask%>% as.character()
new_positions_post_mask$coords.x1 = new_positions_post_mask$coords.x1.mask
new_positions_post_mask$coords.x2 = new_positions_post_mask$coords.x2.mask
   
new_positions = new_positions_post_mask[,colnames(updated_coldata)]             

updated_coldata =
  updated_coldata %>% 
  filter(!(Cell %in% new_positions_post_mask$Cell)) %>%
  rbind(new_positions)

saveRDS(object = updated_coldata,
        file = "Submission_Data/E14_slides/RDS_intermediates/Notebook5.2_updated_coldata.RDS")

