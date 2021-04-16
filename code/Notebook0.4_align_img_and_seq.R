# Imaging and Sequencing Registration: Constructing data structures that hold each image 
# and map files that correspond to the outline or hull of the embryo section or contours  
# of the section with cavities denoted as holes in the constructed polygon.

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
  library(monocle3)
  
  space_directory = "/Volumes/GoogleDrive/My Drive/sciSpace/"
  setwd(dir=space_directory)

  # Pass TRUE if you want to see progress output on some of Monocle 3's operations
  DelayedArray:::set_verbose_block_processing(TRUE)
  # Passing a higher value will make some computations faster but use more memory. 
  # Adjust with caution!
  options(DelayedArray.block.size=1000e6)
  
})


# Read in alignment landmarks ---------------------------------------------

# Landmarks: Choosen pairs of corresponding points that map from an image pixel space to
# an idealized SYBR grid point image. {Submission_Data/E14_slides/Images/Alignment_Landmarks/barcode_grid.jpg}
# Landmarks were chosen on original imaged pixel space {Submission_Data/E14_slides/Images/Slide_XX/}. 

landmarks_path = "Submission_Data/E14_slides/Images/Alignment_Landmarks/"
alignment_coordinates = list.files(path = landmarks_path)
alignment_coordinates = 
  alignment_coordinates[grepl(x = alignment_coordinates,
                              pattern = "landmarks")]

# Read in landmarks that link SYBR waypoints on ideal grid to the image
all_landmarks = list()
for (file_name in alignment_coordinates) {
  landmarks = read.table(file = paste0(landmarks_path,
                                       file_name),
                         sep = ",",
                         col.names = 
                           c("point_id", 
                             "dummy",
                             "source_x",
                             "source_y",
                             "target_x",
                             "target_y")) %>%
    dplyr::select(-dummy)
  
  landmarks$slide_id = 
    file_name %>%
    stringr::str_sub(start = 11,
                     end = 18)
  
  all_landmarks[[file_name]] = landmarks
  
}

all_landmarks = 
  do.call(rbind, all_landmarks)

# Read in images and their size in pixels ---------------------------------

# Read in locations and pixel sizes for each of the images
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

# Calculate a scaling factor and affine transformation for the image of each section and the
# idealized SYBR grid image {Submission_Data/E14_slides/Images/Alignment_Landmarks/barcode_grid.jpg}.
all_landmarks = 
  all_landmarks %>%
  left_join(image_locations_and_sizes,
            by = "slide_id") %>%
  mutate(source_x_scaled = source_x / (starting_width_grid/ending_width_grid ),
         source_y_scaled = ending_height_grid - (source_y / (starting_height_grid/ending_height_grid)),
         target_x_scaled = target_x / (starting_width_image/ending_width_image),
         target_y_scaled = ending_height_image - (target_y/ (starting_height_image/ending_height_image)))

affine_trasformations =
  all_landmarks %>%
  group_by(slide_id) %>%
  dplyr::select(slide_id,
                X_Source = source_x_scaled,
                Y_Source = source_y_scaled,
                X_target = target_x_scaled,
                Y_target = target_y_scaled) %>%
  nest() %>% 
  mutate(affine_transformation = 
           purrr::map(.x = data, .f = function(data_subset){
             aff_obj <- AffineTransformation(controlPoints = data_subset)
             calculateParameters(aff_obj)
             aff_obj
           })) %>%
  dplyr::select(-data)

# Read in rescaled images and store in dataframe --------------------------
# Images were rescaled to allow time and space efficeient plotting in figures
rescale_img_path = "Submission_Data/E14_slides/Images/rescaled_images/"

# Read in the DAPI image, SYBR image and the combo
rescaled_images =
  all_landmarks %>%
  dplyr::select(sybr_path,
                dapi_path,
                dapi_sybr_path,
                slide_id) %>%
  distinct() %>%
  group_by(slide_id) %>%
  nest() %>%
  mutate(dapi_image = 
           purrr::map(.x = data,
                      .f = function(data_subset){
                        image_cimg =
                          load.image(paste0(rescale_img_path,
                                            data_subset$dapi_path,
                                            sep = ""))
                        image_cimg = cimg2im(image_cimg)
                        as.data.frame(image_cimg)
                      }),
         dapi_sybr_image = 
           purrr::map(.x = data,
                      .f = function(data_subset){
                        image_cimg =
                          load.image(paste0(rescale_img_path,
                                            data_subset$dapi_sybr_path,
                                            sep = ""))
                        image_cimg = cimg2im(image_cimg)
                        as.data.frame(image_cimg)
                      }),
         sybr_image = 
           purrr::map(.x = data,
                      .f = function(data_subset){
                        image_cimg =
                          load.image(paste0(rescale_img_path,
                                            data_subset$sybr_path,
                                            sep = ""))
                        image_cimg = cimg2im(image_cimg)
                        as.data.frame(image_cimg)
                      })) %>%
  dplyr::select(-data)

# Read in coordinates that create sf spatial map objects for hulls --------
# Note points are in the the raw pixel space prior to flipping the y axis and scaling

# Hulls polygons without holes used to encompass the embryo section. This helps define the
# set of points that are considered for mapping cells to the oligo grid 

hulls_path = "Submission_Data/E14_slides/Images/image_hull/"
hulls_coordinates = list.files(path = hulls_path)
hulls_coordinates = 
  hulls_coordinates[grepl(x = hulls_coordinates,
                          pattern = ".csv")]

all_hulls = list()
for (file_name in hulls_coordinates) {
  curr_slide = 
    stringr::str_sub(string = file_name,
                     start = 1,
                     end = 8)
  
  type =
    stringr::str_sub(string = file_name,
                     start = 10) %>%
    stringr::str_replace_all(".csv",
                             "")
  all_hulls[[file_name]] = c(file_name, curr_slide, type)
  
}


all_hulls = 
  do.call(rbind, all_hulls) %>%
  as.data.frame()

colnames(all_hulls) = c("file_name",
                        "slide_id",
                        "type")

all_hulls = 
  all_hulls %>%
  mutate(file_name_dummy = file_name) %>%
  group_by(file_name,
           slide_id,
           type) %>%
  nest() %>%
  mutate(hull_data_frame = 
           purrr::map(.x = data,
                      .f = function(subset){
                        df = read.table(file = paste0(hulls_path,
                                                      subset$file_name_dummy,
                                                      sep = ""),
                                        sep = ",")[,c(3,4)] %>%
                          as.data.frame()
                        colnames(df) = c("x","y")
                        
                        df
                        
                      })) %>%
  dplyr::select(-data)

# All the hulls are traced on the original pixel space. They must also get scaled to match 
# rescaled images
all_hulls = 
  all_hulls %>%
  left_join(image_locations_and_sizes,
            by = "slide_id") %>%
  group_by(slide_id,
           file_name,
           type) %>%
  nest() %>%
  mutate(hull_polygon = 
           purrr::map(.x = data,
                      .f = function(subset){
                        curr_hull_df = 
                          subset$hull_data_frame %>%
                          as.data.frame()
                        
                        mat = 
                          data.frame(
                            x_scaled = curr_hull_df$x / (subset$starting_width_image/subset$ending_width_image),
                            y_scaled = subset$ending_height_image - (curr_hull_df$y/ (subset$starting_height_image/subset$ending_height_image))
                          ) %>% as.matrix()
                        
                        mat = rbind(mat,
                                    mat[1,])
                        
                        st_polygon(list(mat))
                        
                      })) %>%
  ungroup() %>%
  dplyr::select(slide_id,
                hull_polygon)



# Read in contours that create sf spatial map objects ---------------------
# Note points are in the the raw pixel space prior to flipping the y axis and scaling

# Contours are the traced outlines of both the section's outline and any discernable
# cavities inside the organism
contours_path = "Submission_Data/E14_slides/Images/Image_contours/"
contours_coordinates = list.files(path = contours_path)
contours_coordinates = 
  contours_coordinates[grepl(x = contours_coordinates,
                             pattern = ".csv")]

all_contours = list()
for (file_name in contours_coordinates) {
  curr_slide = 
    stringr::str_sub(string = file_name,
                     start = 1,
                     end = 8)
  
  # Get whether we are looking at the contour or the hole for constructing sf polygons
  type =
    stringr::str_sub(string = file_name,
                     start = 10) %>%
    stringr::str_replace_all("_[0-9]?[0-9]",
                             "") %>%
    stringr::str_replace_all(".csv",
                             "")
  all_contours[[file_name]] = c(file_name, curr_slide, type)
  
}

all_contours = 
  do.call(rbind, all_contours) %>%
  as.data.frame()

colnames(all_contours) = c("file_name",
                           "slide_id",
                           "type")

all_contours = 
  all_contours %>%
  mutate(file_name_dummy = file_name) %>%
  group_by(file_name,
           slide_id,
           type) %>%
  nest() %>%
  mutate(hull_data_frame = 
           purrr::map(.x = data,
                      .f = function(subset){
                        df = read.table(file = paste0(contours_path,
                                                      subset$file_name_dummy,
                                                      sep = ""),
                                        sep = ",")[,c(3,4)] %>%
                          as.data.frame()
                        colnames(df) = c("x","y")
                        
                        df
                        
                      })) %>%
  dplyr::select(-data)

# All the contours are traced on the original pixel space. They must also get scaled to match 
# rescaled images
all_contours = 
  all_contours %>%
  left_join(image_locations_and_sizes,
            by = "slide_id") %>%
  group_by(slide_id,
           file_name,
           type) %>%
  nest() %>%
  mutate(hull_matrix_scaled = 
           purrr::map(.x = data,
                      .f = function(subset){
                        curr_hull_df = 
                          subset$hull_data_frame %>%
                          as.data.frame()
                        
                        mat = 
                          data.frame(
                            x_scaled = curr_hull_df$x / (subset$starting_width_image/subset$ending_width_image),
                            y_scaled = subset$ending_height_image - (curr_hull_df$y/ (subset$starting_height_image/subset$ending_height_image))
                          ) %>% as.matrix()
                        
                        mat = rbind(mat,
                                    mat[1,])
                        
                        mat
                      }))

all_contours = 
  all_contours %>%
  ungroup() %>%
  mutate(type = factor(type, levels = c("contour",
                                        "hole"))) %>% 
  arrange(slide_id,
          type) %>%
  group_by(slide_id) %>%
  nest() %>%
  mutate(slide_polygon = 
           purrr::map(.x = data, 
                      .f = function(subset){
                        st_polygon(subset$hull_matrix_scaled)
                      })) %>%
  dplyr::select(-data)



# Incorporate the Oligo Grid and SYBR coordinates -------------------------
# There were a few different layouts used during different interations of sciSpace grid printing
# This dataframe matches each imaged and sequenced slide with the corresponding oligo layout

oligo_plate_layout_1 = 
  read.table(file = "Submission_Data/Oligo_layout/oligo_plate_layout.tsv", 
             sep = '\t', 
             header = TRUE)

oligo_plate_layout_1_shifted_18  = 
  read.table(file = "Submission_Data/Oligo_layout/shifted18_oligo_plate_layout.tsv", 
             sep = '\t', 
             header = TRUE)


oligo_plate_layout_2 = 
  read.table(file = "Submission_Data/Oligo_layout_2/oligo_plate_layout_2.tsv", 
             sep = '\t', 
             header = TRUE)

# Match each slide to the oligo layout used for that slide
oligo_plate_layout_df = 
  tibble(slide_id = c("Slide_1D",
                      "Slide_1E",
                      "Slide_1F",
                      "Slide_1G",
                      "Slide_2G",
                      "Slide_2H",
                      "Slide_3D",
                      "Slide_3F",
                      "Slide_3G",
                      "Slide_3H",
                      "Slide_4A",
                      "Slide_4C",
                      "Slide_4D",
                      "Slide_4E"),
         oligo_plate_layout =
           list(oligo_plate_layout_1,
                oligo_plate_layout_1,
                oligo_plate_layout_1_shifted_18,
                oligo_plate_layout_1_shifted_18,
                oligo_plate_layout_1_shifted_18,
                oligo_plate_layout_1,
                oligo_plate_layout_2,
                oligo_plate_layout_2,
                oligo_plate_layout_2,
                oligo_plate_layout_2,
                oligo_plate_layout_2,
                oligo_plate_layout_2,
                oligo_plate_layout_2,
                oligo_plate_layout_2),
         procedure = 
           c(1,1,1,1,1,1,
             2,2,2,2,2,2,2,2))


# Join all imaging parts and map files into one data frame ----------------
all_image_data =
  affine_trasformations %>%
  left_join(rescaled_images,
            by = "slide_id") %>%
  left_join(all_contours,
            by = "slide_id") %>%
  left_join(all_hulls,
            by = "slide_id") %>%
  left_join(oligo_plate_layout_df,
            by = "slide_id")


# Perform the affine transformation of the spatial grid for every image
all_image_data =
  all_image_data %>%
  group_by(slide_id) %>%
  nest() %>%
  mutate(transformed_oligo_layout = 
           purrr::map(.x = data,
                      .f = function(df){
                        # Move mapping code into here
                        oligo_plate_layout = df$oligo_plate_layout[[1]]
                        oligo_grid = oligo_plate_layout[, c("Col", "Row")]
                        if (df$procedure == 2){
                          oligo_grid$Row = abs(oligo_grid$Row -84)
                        }
                        oligo_grid_coord_sp <- SpatialPoints(oligo_grid)
                        
                        transformed_grid = 
                          applyTransformation(df$affine_transformation[[1]],
                                              oligo_grid_coord_sp) %>%
                          as.data.frame()
                        
                        transformed_grid$Oligo = oligo_plate_layout[,"Oligo"]
                        
                        transformed_grid = 
                          transformed_grid %>%
                          left_join(oligo_plate_layout,
                                    by = "Oligo")
                        
                        transformed_grid
                      })) %>%
  unnest(cols = c(data))



# Record which oligo points in the sci-space grid fall within the outlined hull
all_image_data$transformed_oligo_layout = 
  lapply(X = seq(1:dim(all_image_data)[1]),
         FUN = function(x, 
                        this_image_data){
           
           curr_df = 
             this_image_data$transformed_oligo_layout[[x]] 
           
           oligo_points = 
             curr_df %>%
             dplyr::select(coords.x1,coords.x2) %>%
             as.matrix()
           
           oligo_points = 
             SpatialPoints(oligo_points)
      
           polygon_hull =
             this_image_data$hull_polygon[[x]] %>%
             as(Class = "Spatial")
           
           mask = 
             sp::over(oligo_points,
                      polygon_hull)
           
           curr_df$in_hull = 
             !(is.na(mask))
           
           curr_df
             
           },all_image_data)

# Add in x1-coordinates for scale bars
sbX1 <- c(1200,1100,1300,1200,750,1050,1100,1400,1050,950,1300,900,1000,1300)
all_image_data$sbX1 = sbX1

# Add in x2-coordinates for scale bars
sbX2 <- sbX1+(500/(4*1.816))
all_image_data$sbX2 = sbX2 

# Set in y-coordinates for scale bars
sbY <- c(-2100,-1900,-1750,-1750,-1350,-1700,-1700,-1700,-1650,-1550,-1550,-1250,-1700,-1600)
all_image_data$sbY = sbY

# Set the rownames for all_image_data 
all_image_data$slide_id = 
  stringr::str_replace_all(all_image_data$slide_id,
                           "S",
                           "s")

rownames(all_image_data) = all_image_data$slide_id


# Save result for future use
saveRDS(object = all_image_data,
        file = "Submission_Data/E14_slides/RDS_intermediates/Notebook0_1_images_and_transformations.RDS")


