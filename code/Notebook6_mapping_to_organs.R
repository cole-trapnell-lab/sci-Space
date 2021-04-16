# Mapping cells to segemented ROIs. Some of the organs and their boundaries were
# discernible from the DAPI stained section and aided through immunostaining of 
# the adjacent section. In this notebook those segmented ROIs are transformed to
# match the appropriate image dimensions and cells mapping to those ROIs are 
# quantified.

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
  
  space_directory = "/Volumes/GoogleDrive/My Drive/sciSpace/"
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
  readRDS(file = "Submission_Data/E14_slides/RDS_intermediates/Notebook5.5_E14_spatial_CDS.RDS")



# Read in the outlines of the organs and creat polygons -------------------

# Read in contours that create sf spatial map objects
# Note points are in the the raw pixel space prior to flipping the y axis and scaling
# Need to join the countours with the image_locations_and_sizes data frame

organs_path = "Submission_Data/E14_slides/Images/organ_hulls/"
organs_coordinates = list.files(path = organs_path)
organs_coordinates =
  organs_coordinates[grepl(x = organs_coordinates,
                          pattern = ".csv")]

all_organs = list()
for (file_name in organs_coordinates) {
  curr_slide =
    stringr::str_sub(string = file_name,
                     start = 1,
                     end = 8)

  tissue =
    stringr::str_sub(string = file_name,
                     start = 10) %>%
    stringr::str_replace_all(".csv",
                             "")
  all_organs[[file_name]] = c(file_name, curr_slide, tissue)

}


all_organs =
  do.call(rbind, all_organs) %>%
  as.data.frame()

colnames(all_organs) = c("file_name",
                        "slide_id",
                        "tissue")

# Create a matrix from coordinates that can be used to make a polygon
all_organs =
  all_organs %>%
  mutate(file_name_dummy = file_name) %>%
  group_by(file_name,
           slide_id,
           tissue) %>%
  nest() %>%
  mutate(organ_data_frame =
           purrr::map(.x = data,
                      .f = function(subset){
                        df = read.table(file = paste0(organs_path,
                                                      subset$file_name_dummy,
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
all_organs =
  all_organs %>%
  left_join(image_locations_and_sizes,
            by = "slide_id") %>%
  group_by(slide_id,
           file_name,
           tissue) %>%
  nest() %>%
  mutate(organ_polygon =
           purrr::map(.x = data,
                      .f = function(subset){
                        curr_organ_df =
                          subset$organ_data_frame %>%
                          as.data.frame()

                        mat =
                          data.frame(
                            x_scaled = curr_organ_df$x / (subset$starting_width_image/subset$ending_width_image),
                            y_scaled = subset$ending_height_image - (curr_organ_df$y/ (subset$starting_height_image/subset$ending_height_image))
                          ) %>% as.matrix()

                        mat = rbind(mat,
                                    mat[1,])

                        st_polygon(list(mat))

                      })) %>%
  ungroup() %>%
  dplyr::select(slide_id,
                tissue,
                organ_polygon)

all_organs =
  all_organs %>%
  spread(key = tissue, -slide_id)



# The neural tube (called spinal cord in the paper) was composed of two different polygons
# Create a multipolygon object to make a single object
all_organs$NeuralTube[[2]] = 
  st_multipolygon(x = list(all_organs$NeuralTube1[[2]],
                           all_organs$NeuralTube2[[2]] ))
all_organs$NeuralTube[[4]] = 
  st_multipolygon(x = list(all_organs$NeuralTube1[[4]],
                           all_organs$NeuralTube2[[4]] ))

all_organs =
  all_organs %>%
  dplyr::select(-NeuralTube1,
                -NeuralTube2)

all_organs$slide_id = 
  all_organs$slide_id %>%
  stringr::str_replace(pattern = "S",
                       replacement = "s")

rownames(all_organs) = all_organs$slide_id

# Set the colors for highlighting organs ----------------------------------

colors_for_organs = 
  data.frame(color_value = RColorBrewer::brewer.pal(n = 5, name = "Set1"),
             organ = c("Cortex",
                       "Heart",
                       "Liver",
                       "Lung",
                       "NeuralTube"),
             row.names = c("Cortex",
                           "Heart",
                           "Liver",
                           "Lung",
                           "NeuralTube"))


# Supplemental Figure 15 — Slide 1 Embryo 1 -------------------------------

slide_1D =
  colData(spatial_cds) %>% 
  as.data.frame() %>% 
  filter(max_slide_id == "slide_1D") 

slide_1D %>%
  group_by(final_cluster_label) %>%
  summarise(n = n()) %>%
  arrange(-n) %>% head(n = 20)


slide_1D_image_hull =
  all_image_data$dapi_image[[1]] %>%
  dplyr::select(x,
                y) %>%
  as.matrix() %>%
  SpatialPoints() %>%
  over(as(object = all_image_data$slide_polygon[[1]],
          Class = "Spatial")) %>%
  is.na() %>%
  not()

slide_1D_image_cropped_by_hull =
  all_image_data$dapi_image[[1]][slide_1D_image_hull %>% as.vector,]

colData(spatial_cds) %>%
  as.data.frame() %>%
  filter(max_slide_id == "slide_1D") %>%
  ggplot() +
  geom_point(aes(y = -coords.x1,
                x = coords.x2)) 





ggplot() +
  geom_raster(data = slide_1D_image_cropped_by_hull,
              aes(y = -x,
                  x = y,
                  fill = value),
              alpha = 0.5) +
  geom_sf(data = all_organs$NeuralTube[[1]] * rbind(c(0, -1),c(1,0)),
          fill = colors_for_organs["NeuralTube","color_value"]) +
  geom_sf(data = all_organs$Lung[[1]] * rbind(c(0, -1),c(1,0)),
          fill = colors_for_organs["Lung","color_value"]) +
  geom_sf(data = all_organs$Heart[[1]] * rbind(c(0, -1),c(1,0)),
          fill = colors_for_organs["Heart","color_value"]) +
  geom_sf(data = all_organs$Liver[[1]] * rbind(c(0, -1),c(1,0)), 
          fill = colors_for_organs["Liver","color_value"]) +
  geom_sf(data = all_organs$Cortex[[1]] * rbind(c(0, -1),c(1,0)),
          fill = colors_for_organs["Cortex","color_value"]) +
  geom_sf(data = all_image_data$slide_polygon[[1]] * rbind(c(0, -1),c(1,0)),
          fill = NA)+
  geom_segment(data = 
                 all_image_data[1,],
               aes(x =sbX1,
                   xend = sbX2,
                   y = sbY,
                   yend = sbY)) +
  scale_fill_gradient(low = "white", high = "black") +
  monocle3:::monocle_theme_opts()+
  theme_void() +
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Figure2/slide_1D_annotations.png",
         height = 2.25,
         width = 2,
         dpi = 300)



# Supplemental Figure 15 — Slide 2 Embryo 2 -------------------------------

slide_1E =
  colData(spatial_cds) %>% 
  as.data.frame() %>% 
  filter(max_slide_id == "slide_1E") 


slide_1E_image_hull =
  all_image_data$dapi_image[[2]] %>%
  dplyr::select(x,
                y) %>%
  as.matrix() %>%
  SpatialPoints() %>%
  over(as(object = all_image_data$slide_polygon[[2]],
          Class = "Spatial")) %>%
  is.na() %>%
  not()

slide_1E_image_cropped_by_hull =
  all_image_data$dapi_image[[2]][slide_1E_image_hull %>% as.vector,]


ggplot() +
  geom_raster(data = slide_1E_image_cropped_by_hull,
              aes(y = -x,
                  x = -y,
                  fill = value),
              alpha = 0.5) +
  geom_sf(data = all_organs$NeuralTube[[2]] * rbind(c(0, -1),c(-1,0)),
          fill = colors_for_organs["NeuralTube","color_value"]) +
  geom_sf(data = all_organs$Lung[[2]] * rbind(c(0, -1),c(-1,0)),
          fill = colors_for_organs["Lung","color_value"]) +
  geom_sf(data = all_organs$Heart[[2]] * rbind(c(0, -1),c(-1,0)),
          fill = colors_for_organs["Heart","color_value"]) +
  geom_sf(data = all_organs$Liver[[2]] * rbind(c(0, -1),c(-1,0)), 
          fill = colors_for_organs["Liver","color_value"]) +
  geom_sf(data = all_organs$Cortex[[2]] * rbind(c(0, -1),c(-1,0)),
          fill = colors_for_organs["Cortex","color_value"]) +
  geom_sf(data = all_image_data$slide_polygon[[2]] * rbind(c(0, -1),c(-1,0)),
          fill = NA)+
  geom_segment(data = 
                 all_image_data[2,],
               aes(x = -sbX1 + 900,
                   xend = -sbX2 + 900,
                   y = sbY,
                   yend = sbY)) +
  scale_fill_gradient(low = "white", high = "black") +
  monocle3:::monocle_theme_opts()+
  theme_void() +
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Supplement_Organ_Composition/slide_1E_annotations.png",
         height = 2.25,
         width = 2,
         dpi = 300)


# Supplemental Figure 15 — Slide 3 Embryo 2 -------------------------------

slide_1F =
  colData(spatial_cds) %>% 
  as.data.frame() %>% 
  filter(max_slide_id == "slide_1F") 


slide_1F_image_hull =
  all_image_data$dapi_image[[3]] %>%
  dplyr::select(x,
                y) %>%
  as.matrix() %>%
  SpatialPoints() %>%
  over(as(object = all_image_data$slide_polygon[[3]],
          Class = "Spatial")) %>%
  is.na() %>%
  not()

slide_1F_image_cropped_by_hull =
  all_image_data$dapi_image[[3]][slide_1F_image_hull %>% as.vector,]


ggplot() +
  geom_raster(data = slide_1F_image_cropped_by_hull,
              aes(y = -x,
                  x = -y,
                  fill = value),
              alpha = 0.5) +
  geom_sf(data = all_organs$NeuralTube[[3]] * rbind(c(0, -1),c(-1,0)),
          fill = colors_for_organs["NeuralTube","color_value"]) +
  geom_sf(data = all_organs$Lung[[3]] * rbind(c(0, -1),c(-1,0)),
          fill = colors_for_organs["Lung","color_value"]) +
  geom_sf(data = all_organs$Heart[[3]] * rbind(c(0, -1),c(-1,0)),
          fill = colors_for_organs["Heart","color_value"]) +
  geom_sf(data = all_organs$Liver[[3]] * rbind(c(0, -1),c(-1,0)), 
          fill = colors_for_organs["Liver","color_value"]) +
  geom_sf(data = all_organs$Cortex[[3]] * rbind(c(0, -1),c(-1,0)),
          fill = colors_for_organs["Cortex","color_value"]) +
  geom_sf(data = all_image_data$slide_polygon[[3]] * rbind(c(0, -1),c(-1,0)),
          fill = NA)+
  geom_segment(data = 
                 all_image_data[3,],
               aes(x = -sbX1 + 900,
                   xend = -sbX2 + 900,
                   y = sbY,
                   yend = sbY)) +
  scale_fill_gradient(low = "white", high = "black") +
  monocle3:::monocle_theme_opts()+
  theme_void() +
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Supplement_Organ_Composition/slide_1F_annotations.png",
         height = 2.25,
         width = 2,
         dpi = 300)



# Supplemental Figure 15 — Slide 4 Embryo 1 -------------------------------

slide_1G_image_hull =
  all_image_data$dapi_image[[4]] %>%
  dplyr::select(x,
                y) %>%
  as.matrix() %>%
  SpatialPoints() %>%
  over(as(object = all_image_data$slide_polygon[[4]],
          Class = "Spatial")) %>%
  is.na() %>%
  not()

slide_1G_image_cropped_by_hull =
  all_image_data$dapi_image[[4]][slide_1G_image_hull %>% as.vector,]


ggplot() +
  geom_raster(data = slide_1G_image_cropped_by_hull,
              aes(y = -x,
                  x = y,
                  fill = value),
              alpha = 0.5) +
  geom_sf(data = all_organs$NeuralTube[[4]] * rbind(c(0, -1),c(1,0)),
          fill = colors_for_organs["NeuralTube","color_value"]) +
  geom_sf(data = all_organs$Lung[[4]] * rbind(c(0, -1),c(1,0)),
          fill = colors_for_organs["Lung","color_value"]) +
  geom_sf(data = all_organs$Heart[[4]] * rbind(c(0, -1),c(1,0)),
          fill = colors_for_organs["Heart","color_value"]) +
  geom_sf(data = all_organs$Liver[[4]] * rbind(c(0, -1),c(1,0)), 
          fill = colors_for_organs["Liver","color_value"]) +
  geom_sf(data = all_organs$Cortex[[4]] * rbind(c(0, -1),c(1,0)),
          fill = colors_for_organs["Cortex","color_value"]) +
  geom_sf(data = all_image_data$slide_polygon[[4]] * rbind(c(0, -1),c(1,0)),
          fill = NA)+
  geom_segment(data = 
                 all_image_data[4,],
               aes(x = sbX1,
                   xend = sbX2,
                   y = sbY,
                   yend = sbY)) +
  scale_fill_gradient(low = "white", high = "black") +
  monocle3:::monocle_theme_opts()+
  theme_void() +
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Supplement_Organ_Composition/slide_1G_annotations.png",
         height = 2.25,
         width = 2,
         dpi = 300)


# Supplemental Figure 15 — Slide 5 Embryo 1 -------------------------------

slide_2G_image_hull =
  all_image_data$dapi_image[[5]] %>%
  dplyr::select(x,
                y) %>%
  as.matrix() %>%
  SpatialPoints() %>%
  over(as(object = all_image_data$slide_polygon[[5]],
          Class = "Spatial")) %>%
  is.na() %>%
  not()

slide_2G_image_cropped_by_hull =
  all_image_data$dapi_image[[5]][slide_2G_image_hull %>% as.vector,]


ggplot() +
  geom_raster(data = slide_2G_image_cropped_by_hull,
              aes(y = -x,
                  x = y,
                  fill = value),
              alpha = 0.5) +
  geom_sf(data = all_organs$NeuralTube[[5]] * rbind(c(0, -1),c(1,0)),
          fill = colors_for_organs["NeuralTube","color_value"]) +
  geom_sf(data = all_organs$Lung[[5]] * rbind(c(0, -1),c(1,0)),
          fill = colors_for_organs["Lung","color_value"]) +
  geom_sf(data = all_organs$Heart[[5]] * rbind(c(0, -1),c(1,0)),
          fill = colors_for_organs["Heart","color_value"]) +
  geom_sf(data = all_organs$Liver[[5]] * rbind(c(0, -1),c(1,0)), 
          fill = colors_for_organs["Liver","color_value"]) +
  geom_sf(data = all_organs$Cortex[[5]] * rbind(c(0, -1),c(1,0)),
          fill = colors_for_organs["Cortex","color_value"]) +
  geom_sf(data = all_image_data$slide_polygon[[5]] * rbind(c(0, -1),c(1,0)),
          fill = NA)+
  geom_segment(data = 
                 all_image_data[5,],
               aes(x = sbX1,
                   xend = sbX2,
                   y = sbY,
                   yend = sbY)) +
  scale_fill_gradient(low = "white", high = "black") +
  monocle3:::monocle_theme_opts()+
  theme_void() +
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Supplement_Organ_Composition/slide_2G_annotations.png",
         height = 2.25,
         width = 2,
         dpi = 300)

# Supplemental Figure 15 — Slide 6 Embryo 1 -------------------------------


slide_2H_image_hull =
  all_image_data$dapi_image[[6]] %>%
  dplyr::select(x,
                y) %>%
  as.matrix() %>%
  SpatialPoints() %>%
  over(as(object = all_image_data$slide_polygon[[6]],
          Class = "Spatial")) %>%
  is.na() %>%
  not()

slide_2H_image_cropped_by_hull =
  all_image_data$dapi_image[[6]][slide_2H_image_hull %>% as.vector,]



ggplot() +
  geom_raster(data = slide_2H_image_cropped_by_hull,
              aes(y = -x,
                  x = y,
                  fill = value),
              alpha = 0.5) +
  geom_sf(data = all_organs$NeuralTube[[6]] * rbind(c(0, -1),c(1,0)),
          fill = colors_for_organs["NeuralTube","color_value"]) +
  geom_sf(data = all_organs$Lung[[6]] * rbind(c(0, -1),c(1,0)),
          fill = colors_for_organs["Lung","color_value"]) +
  geom_sf(data = all_organs$Heart[[6]] * rbind(c(0, -1),c(1,0)),
          fill = colors_for_organs["Heart","color_value"]) +
  geom_sf(data = all_organs$Liver[[6]] * rbind(c(0, -1),c(1,0)), 
          fill = colors_for_organs["Liver","color_value"]) +
  geom_sf(data = all_organs$Cortex[[6]] * rbind(c(0, -1),c(1,0)),
          fill = colors_for_organs["Cortex","color_value"]) +
  geom_sf(data = all_image_data$slide_polygon[[6]] * rbind(c(0, -1),c(1,0)),
          fill = NA)+
  geom_segment(data = 
                 all_image_data[6,],
               aes(x = sbX1,
                   xend = sbX2,
                   y = sbY,
                   yend = sbY)) +
  scale_fill_gradient(low = "white", high = "black") +
  monocle3:::monocle_theme_opts()+
  theme_void() +
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Supplement_Organ_Composition/slide_2H_annotations.png",
         height = 2.25,
         width = 2,
         dpi = 300)



# Supplemental Figure 15 — Slide 7 Embryo 1 -------------------------------


slide_3D_image_hull =
  all_image_data$dapi_image[[7]] %>%
  dplyr::select(x,
                y) %>%
  as.matrix() %>%
  SpatialPoints() %>%
  over(as(object = all_image_data$slide_polygon[[7]],
          Class = "Spatial")) %>%
  is.na() %>%
  not()

slide_3D_image_cropped_by_hull =
  all_image_data$dapi_image[[7]][slide_3D_image_hull %>% as.vector,]



ggplot() +
  geom_raster(data = slide_3D_image_cropped_by_hull,
              aes(y = x,
                  x = y,
                  fill = value),
              alpha = 0.5) +
  geom_sf(data = all_organs$Lung[[7]] * rbind(c(0, 1),c(1,0)),
          fill = colors_for_organs["Lung","color_value"]) +
  geom_sf(data = all_organs$Heart[[7]] * rbind(c(0, 1),c(1,0)),
          fill = colors_for_organs["Heart","color_value"]) +
  geom_sf(data = all_organs$Liver[[7]] * rbind(c(0, 1),c(1,0)), 
          fill = colors_for_organs["Liver","color_value"]) +
  geom_sf(data = all_organs$Cortex[[7]] * rbind(c(0, 1),c(1,0)),
          fill = colors_for_organs["Cortex","color_value"]) +
  geom_sf(data = all_image_data$slide_polygon[[7]] * rbind(c(0, 1),c(1,0)),
          fill = NA)+
  geom_segment(data = 
                 all_image_data[7,],
               aes(x = sbX1,
                   xend = sbX2,
                   y = -sbY - 1500,
                   yend = -sbY - 1500)) +
  scale_fill_gradient(low = "white", high = "black") +
  monocle3:::monocle_theme_opts()+
  theme_void() +
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Supplement_Organ_Composition/slide_3D_annotations.png",
         height = 2.25,
         width = 2,
         dpi = 300)



# Supplemental Figure 15 — Slide 8 Embryo 1 -------------------------------


slide_3F_image_hull =
  all_image_data$dapi_image[[8]] %>%
  dplyr::select(x,
                y) %>%
  as.matrix() %>%
  SpatialPoints() %>%
  over(as(object = all_image_data$slide_polygon[[8]],
          Class = "Spatial")) %>%
  is.na() %>%
  not()

slide_3F_image_cropped_by_hull =
  all_image_data$dapi_image[[8]][slide_3F_image_hull %>% as.vector,]



ggplot() +
  geom_raster(data = slide_3F_image_cropped_by_hull,
              aes(y = -x,
                  x = -y,
                  fill = value),
              alpha = 0.5) +
  geom_sf(data = all_organs$Lung[[8]] * rbind(c(0, -1),c(-1,0)),
          fill = colors_for_organs["Lung","color_value"]) +
  geom_sf(data = all_organs$Heart[[8]] * rbind(c(0, -1),c(-1,0)),
          fill = colors_for_organs["Heart","color_value"]) +
  geom_sf(data = all_organs$Liver[[8]] * rbind(c(0, -1),c(-1,0)), 
          fill = colors_for_organs["Liver","color_value"]) +
  geom_sf(data = all_organs$Cortex[[8]] * rbind(c(0, -1),c(-1,0)),
          fill = colors_for_organs["Cortex","color_value"]) +
  geom_sf(data = all_image_data$slide_polygon[[8]] * rbind(c(0, -1),c(-1,0)),
          fill = NA)+
  geom_segment(data = 
                 all_image_data[8,],
               aes(x = -sbX1 + 850,
                   xend = -sbX2 + 850,
                   y = sbY ,
                   yend = sbY )) +
  scale_fill_gradient(low = "white", high = "black") +
  monocle3:::monocle_theme_opts()+
  theme_void() +
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Supplement_Organ_Composition/slide_3F_annotations.png",
         height = 2.25,
         width = 2,
         dpi = 300)



# Supplemental Figure 15 — Slide 9 Embryo 1 -------------------------------


slide_3G_image_hull =
  all_image_data$dapi_image[[9]] %>%
  dplyr::select(x,
                y) %>%
  as.matrix() %>%
  SpatialPoints() %>%
  over(as(object = all_image_data$slide_polygon[[9]],
          Class = "Spatial")) %>%
  is.na() %>%
  not()

slide_3G_image_cropped_by_hull =
  all_image_data$dapi_image[[9]][slide_3G_image_hull %>% as.vector,]



ggplot() +
  geom_raster(data = slide_3G_image_cropped_by_hull,
              aes(y = -x,
                  x = -y,
                  fill = value),
              alpha = 0.5) +
  geom_sf(data = all_organs$Lung[[9]] * rbind(c(0, -1),c(-1,0)),
          fill = colors_for_organs["Lung","color_value"]) +
  geom_sf(data = all_organs$Heart[[9]] * rbind(c(0, -1),c(-1,0)),
          fill = colors_for_organs["Heart","color_value"]) +
  geom_sf(data = all_organs$Liver[[9]] * rbind(c(0, -1),c(-1,0)),
          fill = colors_for_organs["Liver","color_value"]) +
  geom_sf(data = all_organs$Cortex[[9]] * rbind(c(0, -1),c(-1,0)),
          fill = colors_for_organs["Cortex","color_value"]) +
  geom_sf(data = all_image_data$slide_polygon[[9]] * rbind(c(0, -1),c(-1,0)),
          fill = NA)+
  geom_segment(data = 
                 all_image_data[9,],
               aes(x = -sbX1 + 800,
                   xend = -sbX2 + 800,
                   y = sbY + 50,
                   yend = sbY  + 50)) +
  scale_fill_gradient(low = "white", high = "black") +
  monocle3:::monocle_theme_opts()+
  theme_void() +
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Supplement_Organ_Composition/slide_3G_annotations.png",
         height = 2.25,
         width = 2,
         dpi = 300)



# Supplemental Figure 15 — Slide 10 Embryo 1 -------------------------------


slide_3H_image_hull =
  all_image_data$dapi_image[[10]] %>%
  dplyr::select(x,
                y) %>%
  as.matrix() %>%
  SpatialPoints() %>%
  over(as(object = all_image_data$slide_polygon[[10]],
          Class = "Spatial")) %>%
  is.na() %>%
  not()

slide_3H_image_cropped_by_hull =
  all_image_data$dapi_image[[10]][slide_3H_image_hull %>% as.vector,]



ggplot() +
  geom_raster(data = slide_3H_image_cropped_by_hull,
              aes(y = -x,
                  x = -y,
                  fill = value),
              alpha = 0.5) +
  geom_sf(data = all_organs$Lung[[10]] * rbind(c(0, -1),c(-1,0)),
          fill = colors_for_organs["Lung","color_value"]) +
  geom_sf(data = all_organs$Heart[[10]] * rbind(c(0, -1),c(-1,0)),
          fill = colors_for_organs["Heart","color_value"]) +
  geom_sf(data = all_organs$Liver[[10]] * rbind(c(0, -1),c(-1,0)),
          fill = colors_for_organs["Liver","color_value"]) +
  geom_sf(data = all_organs$Cortex[[10]] * rbind(c(0, -1),c(-1,0)),
          fill = colors_for_organs["Cortex","color_value"]) +
  geom_sf(data = all_image_data$slide_polygon[[10]] * rbind(c(0, -1),c(-1,0)),
          fill = NA)+
  geom_segment(data = 
                 all_image_data[10,],
               aes(x = -sbX1 + 750,
                   xend = -sbX2 + 750,
                   y = sbY,
                   yend = sbY)) +
  scale_fill_gradient(low = "white", high = "black") +
  monocle3:::monocle_theme_opts()+
  theme_void() +
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Supplement_Organ_Composition/slide_3H_annotations.png",
         height = 2.25,
         width = 2,
         dpi = 300)



# Supplemental Figure 15 — Slide 11 Embryo 1 -------------------------------


slide_4A_image_hull =
  all_image_data$dapi_image[[11]] %>%
  dplyr::select(x,
                y) %>%
  as.matrix() %>%
  SpatialPoints() %>%
  over(as(object = all_image_data$slide_polygon[[11]],
          Class = "Spatial")) %>%
  is.na() %>%
  not()

slide_4A_image_cropped_by_hull =
  all_image_data$dapi_image[[11]][slide_4A_image_hull %>% as.vector,]



ggplot() +
  # geom_raster(data = slide_4A_image_cropped_by_hull,
  #             aes(y = -x,
  #                 x = -y,
  #                 fill = value),
  #             alpha = 0.5) +
  geom_sf(data = all_organs$Liver[[11]] * rbind(c(0, -1),c(-1,0)),
          fill = colors_for_organs["Liver","color_value"]) +
  geom_sf(data = all_organs$Cortex[[11]] * rbind(c(0, -1),c(-1,0)),
          fill = colors_for_organs["Cortex","color_value"]) +
  geom_sf(data = all_image_data$slide_polygon[[11]] * rbind(c(0, -1),c(-1,0)),
          fill = NA)+
  geom_segment(data = 
                 all_image_data[11,],
               aes(x = -sbX1 + 750,
                   xend = -sbX2 + 750,
                   y = sbY,
                   yend = sbY)) +
  scale_fill_gradient(low = "white", high = "black") +
  monocle3:::monocle_theme_opts()+
  theme_void() +
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Supplement_Organ_Composition/slide_4A_annotations.png",
         height = 2.25,
         width = 2,
         dpi = 300)




# Supplemental Figure 15 — Slide 13 Embryo 1 -------------------------------


slide_4D_image_hull =
  all_image_data$dapi_image[[13]] %>%
  dplyr::select(x,
                y) %>%
  as.matrix() %>%
  SpatialPoints() %>%
  over(as(object = all_image_data$slide_polygon[[13]],
          Class = "Spatial")) %>%
  is.na() %>%
  not()

slide_4D_image_cropped_by_hull =
  all_image_data$dapi_image[[13]][slide_4D_image_hull %>% as.vector,]



ggplot() +
  geom_raster(data = slide_4D_image_cropped_by_hull,
              aes(y = x,
                  x = y,
                  fill = value),
              alpha = 0.5) +
  geom_sf(data = all_organs$Lung[[12]] * rbind(c(0, 1),c(1,0)),
          fill = colors_for_organs["Lung","color_value"]) +
  geom_sf(data = all_organs$Heart[[12]] * rbind(c(0, 1),c(1,0)),
          fill = colors_for_organs["Heart","color_value"]) +
  geom_sf(data = all_organs$Liver[[12]] * rbind(c(0, 1),c(1,0)),
          fill = colors_for_organs["Liver","color_value"]) +
  geom_sf(data = all_organs$Cortex[[12]] * rbind(c(0, 1),c(1,0)),
          fill = colors_for_organs["Cortex","color_value"]) +
  geom_sf(data = all_image_data$slide_polygon[[13]] * rbind(c(0, 1),c(1,0)),
          fill = NA)+
  geom_segment(data = 
                 all_image_data[13,],
               aes(x = sbX1,
                   xend = sbX2,
                   y = -sbY - 1500,
                   yend = -sbY - 1500)) +
  scale_fill_gradient(low = "white", high = "black") +
  monocle3:::monocle_theme_opts()+
  theme_void() +
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Supplement_Organ_Composition/slide_4D_annotations.png",
         height = 2.25,
         width = 2,
         dpi = 300)



# Supplemental Figure 15 — Slide 14 Embryo 1 -------------------------------


slide_4E_image_hull =
  all_image_data$dapi_image[[14]] %>%
  dplyr::select(x,
                y) %>%
  as.matrix() %>%
  SpatialPoints() %>%
  over(as(object = all_image_data$slide_polygon[[14]],
          Class = "Spatial")) %>%
  is.na() %>%
  not()

slide_4E_image_cropped_by_hull =
  all_image_data$dapi_image[[14]][slide_4E_image_hull %>% as.vector,]



ggplot() +
  geom_raster(data = slide_4E_image_cropped_by_hull,
              aes(y = x,
                  x = y,
                  fill = value),
              alpha = 0.5) +
  geom_sf(data = all_organs$NeuralTube[[13]] * rbind(c(0, 1),c(1,0)),
          fill = colors_for_organs["NeuralTube","color_value"]) +
  geom_sf(data = all_organs$Lung[[13]] * rbind(c(0, 1),c(1,0)),
          fill = colors_for_organs["Lung","color_value"]) +
  geom_sf(data = all_organs$Heart[[13]] * rbind(c(0, 1),c(1,0)),
          fill = colors_for_organs["Heart","color_value"]) +
  geom_sf(data = all_organs$Liver[[13]] * rbind(c(0, 1),c(1,0)),
          fill = colors_for_organs["Liver","color_value"]) +
  geom_sf(data = all_organs$Cortex[[13]] * rbind(c(0, 1),c(1,0)),
          fill = colors_for_organs["Cortex","color_value"]) +
  geom_sf(data = all_image_data$slide_polygon[[14]] * rbind(c(0, 1),c(1,0)),
          fill = NA)+
  geom_segment(data = 
                 all_image_data[14,],
               aes(x = sbX1,
                   xend = sbX2,
                   y = -sbY - 1500,
                   yend = -sbY - 1500)) +
  scale_fill_gradient(low = "white", high = "black") +
  monocle3:::monocle_theme_opts()+
  theme_void() +
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Supplement_Organ_Composition/slide_4E_annotations.png",
         height = 2.25,
         width = 2,
         dpi = 300)


# Mapping cells to each ROI -----------------------------------------------


organs_to_test = 
  c("Cortex","Heart","Liver","Lung","NeuralTube")

all_cell_types = 
  colData(spatial_cds)$final_cluster_label %>%
  unique()

# Get the cells that overlap each annotated organ
annotation_mapped_cells = 
  lapply(X = seq(1,dim(all_organs)[1]),
       function(X){
         this.slide = rownames(all_organs)[X]
         print(this.slide)
         
         organ.res =
           lapply(organs_to_test, 
                function(Y){

                  cells_in_this_section =
                    spatial_cds %>%
                    colData() %>%
                    as.data.frame() %>%
                    filter(max_slide_id == this.slide) %>%
                    filter(!is.na(coords.x1),
                           !is.na(coords.x2)) %>%
                    mutate(final_cluster_label =
                             factor(final_cluster_label,levels = all_cell_types))

                  this.organ = Y
                  print(Y)
                  this.polygon =
                    all_organs[which(all_organs$slide_id == this.slide),this.organ][[1]][[1]]
                  
                  if(!is.null(this.polygon)){
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
                  }
         })
         do.call(rbind,
                 organ.res)
       })

annotation_mapped_cells = 
  do.call(rbind, annotation_mapped_cells)

# All cells that are not given an annotation are labeled unannotated
annotation_mapped_cells=
  rbind(
  annotation_mapped_cells,
  colData(spatial_cds) %>%
    as.data.frame() %>%
    filter(!(Cell %in% annotation_mapped_cells$Cell)) %>%
    mutate(organ = "Unannotated") %>%
    dplyr::select(Cell,
                max_slide_id,
                final_cluster_label,
                organ)
  )

# Some cells mapped to two anatomical locations due to overlapping segmentation
# Choose one of the two annotations for these cells
annotation_mapped_cells =
  annotation_mapped_cells %>%
  group_by(Cell) %>%
  sample_n(1)

colData(spatial_cds)$anatomical_annotation = 
  annotation_mapped_cells  %>%
  dplyr::select(Cell,
                organ) %>%
  right_join(colData(spatial_cds) %>%
               as.data.frame(),
             by = "Cell") %>%
  pull(organ)



# Create heatmaps for celltype x organ ------------------------------------

annotation_mapped_cells = 
  colData(spatial_cds) %>%
  as.data.frame() %>%
  dplyr::select(Cell,
                max_slide_id,
                final_cluster_label,
                organ = anatomical_annotation) %>%
  mutate(final_cluster_label = ifelse(final_cluster_label == "Connective Tissue Progenitors",
                                      "CTPs",
                                      final_cluster_label)) %>%
  group_by(final_cluster_label,
           max_slide_id) %>%
  add_tally(name = "total_label_on_slide") %>%
  ungroup() %>%
  mutate(organ = ifelse(organ == "NeuralTube",
                        "Spinal Cord",
                        organ)) %>%
  mutate(organ = as.factor(organ)) %>%
  group_by(organ,
           max_slide_id,
           final_cluster_label) %>%
  mutate(percent_label = n()/total_label_on_slide) %>%
  dplyr::select(percent_label,
                final_cluster_label,
                organ,
                max_slide_id) %>%
  distinct() %>%
  ungroup() 


annotation_mapped_cells$final_cluster_label =
  ifelse(annotation_mapped_cells$final_cluster_label == "Neuron",
         "Neurons",
         annotation_mapped_cells$final_cluster_label)

annotation_mapped_cells$final_cluster_label =
  ifelse(annotation_mapped_cells$final_cluster_label == "Cardiac muscle lineages",
         "Cardiomyocytes",
         annotation_mapped_cells$final_cluster_label)



annotation_mapped_cells$final_cluster_label =
  ifelse(annotation_mapped_cells$final_cluster_label == "Radial glia",
         "Radial Glia",
         annotation_mapped_cells$final_cluster_label)

cell_types_to_highlight = 
  c("Neurons",
    "OPCs",
    "Radial Glia",
    "Choroid Plexus",
    "Glial Cells",
    "Epithelial Cells",
    "CTPs",
    "Hepatocytes",
    "Cardiomyocytes",
    "Erythroid Lineage" ,
    "White Blood Cells")

row_order =
  c("Unannotated",
    "Cortex",
    "Spinal Cord",
    "Heart",
    "Lung",
    "Liver")

# Supplemental Figure 18 — Panel A — Slide 1 ------------------------------

matrix_for_heatmap_1D = 
  annotation_mapped_cells %>%
  mutate(final_cluster_label = as.factor(final_cluster_label)) %>%
  filter(max_slide_id == "slide_1D") %>% 
  dplyr::select(-max_slide_id) %>%
  spread(key =organ, value =  percent_label, fill = 0, drop = F) %>%
  gather(key = "organ",
         value = "percent_label",
         -final_cluster_label) %>%
  spread(key =final_cluster_label, value =  percent_label, fill = 0, drop = F) %>%
  tibble::column_to_rownames(var = "organ") %>%
  as.matrix()
  
matrix_for_heatmap_1D

matrix_for_heatmap_1D[row_order,cell_types_to_highlight] %>%
  as.data.frame() %>%
  round(digits = 2)

pheatmap(t(matrix_for_heatmap_1D[row_order,cell_types_to_highlight]),
         cellwidth = 7,
         cellheight = 7,
         cluster_rows = F,
         cluster_cols = F,
         treeheight_row = 0,
         treeheight_col = 0,
         show_rownames = F,
         fontsize_row = 8,
         fontsize_col = 8,
         filename = "Figures/Figure_Components/Supplement_Organ_Composition/percent_in_anatomy_slide1.pdf",
         color = viridis(option = "inferno", n = 20),
         legend= F) 

# Supplemental Figure 18 — Panel A — Slide 2 ------------------------------

matrix_for_heatmap_1E = 
  annotation_mapped_cells %>%
  mutate(final_cluster_label = as.factor(final_cluster_label)) %>%
  filter(max_slide_id == "slide_1E") %>% 
  dplyr::select(-max_slide_id) %>%
  spread(key =organ, value =  percent_label, fill = 0, drop = F) %>%
  gather(key = "organ",
         value = "percent_label",
         -final_cluster_label) %>%
  spread(key =final_cluster_label, value =  percent_label, fill = 0, drop = F) %>%
  tibble::column_to_rownames(var = "organ") %>%
  as.matrix()

pheatmap(t(matrix_for_heatmap_1E[row_order,cell_types_to_highlight]),
         cellwidth = 7,
         cellheight = 7,
         cluster_rows = F,
         cluster_cols = F,
         treeheight_row = 0,
         treeheight_col = 0,
         show_rownames = F,
         fontsize_row = 8,
         fontsize_col = 8,
         filename = "Figures/Figure_Components/Supplement_Organ_Composition/percent_in_anatomy_slide_2.pdf",
         color = viridis(option = "inferno", n = 20),
         legend= F) 

# Supplemental Figure 18 — Panel A — Slide 6 ------------------------------

matrix_for_heatmap_2H = 
  annotation_mapped_cells %>%
  mutate(final_cluster_label = as.factor(final_cluster_label)) %>%
  filter(max_slide_id == "slide_2H") %>% 
  dplyr::select(-max_slide_id) %>%
  spread(key =organ, value =  percent_label, fill = 0, drop = F) %>%
  gather(key = "organ",
         value = "percent_label",
         -final_cluster_label) %>%
  spread(key =final_cluster_label, value =  percent_label, fill = 0, drop = F) %>%
  tibble::column_to_rownames(var = "organ") %>%
  as.matrix()



pheatmap(t(matrix_for_heatmap_2H[row_order,cell_types_to_highlight]),
         cellwidth = 7,
         cellheight = 7,
         cluster_rows = F,
         cluster_cols = F,
         treeheight_row = 0,
         treeheight_col = 0,
         show_rownames = T,
         fontsize_row = 8,
         fontsize_col = 8,
         filename = "Figures/Figure_Components/Supplement_Organ_Composition/percent_in_anatomy_slide_6.pdf",
         color = viridis(option = "inferno", n = 20),
         legend= F) 

# Supplemental Figure 18 — Panel A — Slide 5 ------------------------------

matrix_for_heatmap_2G = 
  annotation_mapped_cells %>%
  mutate(final_cluster_label = as.factor(final_cluster_label)) %>%
  filter(max_slide_id == "slide_2G") %>% 
  dplyr::select(-max_slide_id) %>%
  spread(key =organ, value =  percent_label, fill = 0, drop = F) %>%
  gather(key = "organ",
         value = "percent_label",
         -final_cluster_label) %>%
  spread(key =final_cluster_label, value =  percent_label, fill = 0, drop = F) %>%
  tibble::column_to_rownames(var = "organ") %>%
  as.matrix() 

pheatmap(t(matrix_for_heatmap_2G[row_order,cell_types_to_highlight]),
         cellwidth = 7,
         cellheight = 7,
         cluster_rows = F,
         cluster_cols = F,
         treeheight_row = 0,
         treeheight_col = 0,
         show_rownames = F,
         fontsize_row = 8,
         fontsize_col = 8,
         filename = "Figures/Figure_Components/Supplement_Organ_Composition/percent_in_anatomy_slide_5.pdf",
         color = viridis(option = "inferno", n = 20),
         legend= F) 

# Supplemental Figure 18 — Panel A — Slide 4 ------------------------------

matrix_for_heatmap_1G = 
  annotation_mapped_cells %>%
  mutate(final_cluster_label = as.factor(final_cluster_label)) %>%
  filter(max_slide_id == "slide_1G") %>% 
  dplyr::select(-max_slide_id) %>%
  spread(key =organ, value =  percent_label, fill = 0, drop = F) %>%
  gather(key = "organ",
         value = "percent_label",
         -final_cluster_label) %>%
  spread(key =final_cluster_label, value =  percent_label, fill = 0, drop = F) %>%
  tibble::column_to_rownames(var = "organ") %>%
  as.matrix() 


pheatmap(t(matrix_for_heatmap_1G[row_order,cell_types_to_highlight]),
         cellwidth = 7,
         cellheight = 7,
         cluster_rows = F,
         cluster_cols = F,
         treeheight_row = 0,
         treeheight_col = 0,
         show_rownames = F,
         fontsize_row = 8,
         fontsize_col = 8,
         filename = "Figures/Figure_Components/Supplement_Organ_Composition/percent_in_anatomy_slide_4.pdf",
         color = viridis(option = "inferno", n = 20),
         legend= F) 



# Supplemental Figure 18 — Panel A — Slide 3 ------------------------------

matrix_for_heatmap_1F = 
  annotation_mapped_cells %>%
  mutate(final_cluster_label = as.factor(final_cluster_label)) %>%
  filter(max_slide_id == "slide_1F") %>% 
  dplyr::select(-max_slide_id) %>%
  spread(key =organ, value =  percent_label, fill = 0, drop = F) %>%
  gather(key = "organ",
         value = "percent_label",
         -final_cluster_label) %>%
  spread(key =final_cluster_label, value =  percent_label, fill = 0, drop = F) %>%
  tibble::column_to_rownames(var = "organ") %>%
  as.matrix() 


pheatmap(t(matrix_for_heatmap_1F[row_order,cell_types_to_highlight]),
         cellwidth = 7,
         cellheight = 7,
         cluster_rows = F,
         cluster_cols = F,
         treeheight_row = 0,
         treeheight_col = 0,
         show_rownames = F,
         fontsize_row = 8,
         fontsize_col = 8,
         filename = "Figures/Figure_Components/Supplement_Organ_Composition/percent_in_anatomy_slide_3.pdf",
         color = viridis(option = "inferno", n = 20),
         legend= F) 


# Supplemental Figure 18 — Panel A — Slide 8 ------------------------------

matrix_for_heatmap_3D = 
  annotation_mapped_cells %>%
  mutate(final_cluster_label = as.factor(final_cluster_label)) %>%
  filter(max_slide_id == "slide_3D") %>% 
  dplyr::select(-max_slide_id) %>%
  spread(key =organ, value =  percent_label, fill = 0, drop = F) %>%
  gather(key = "organ",
         value = "percent_label",
         -final_cluster_label) %>%
  spread(key =final_cluster_label, value =  percent_label, fill = 0, drop = F) %>%
  tibble::column_to_rownames(var = "organ") %>%
  as.matrix() 


pheatmap(t(matrix_for_heatmap_3D[row_order,cell_types_to_highlight]),
         cellwidth = 7,
         cellheight = 7,
         cluster_rows = F,
         cluster_cols = F,
         treeheight_row = 0,
         treeheight_col = 0,
         show_rownames = F,
         fontsize_row = 8,
         fontsize_col = 8,
         filename = "Figures/Figure_Components/Supplement_Organ_Composition/percent_in_anatomy_slide_7.pdf",
         color = viridis(option = "inferno", n = 20),
         legend= F)


matrix_for_heatmap_3F = 
  annotation_mapped_cells %>%
  mutate(final_cluster_label = as.factor(final_cluster_label)) %>%
  filter(max_slide_id == "slide_3F") %>% 
  dplyr::select(-max_slide_id) %>%
  spread(key =organ, value =  percent_label, fill = 0, drop = F) %>%
  gather(key = "organ",
         value = "percent_label",
         -final_cluster_label) %>%
  spread(key =final_cluster_label, value =  percent_label, fill = 0, drop = F) %>%
  tibble::column_to_rownames(var = "organ") %>%
  as.matrix() 



pheatmap(t(matrix_for_heatmap_3F[row_order,cell_types_to_highlight]),
         cellwidth = 7,
         cellheight = 7,
         cluster_rows = F,
         cluster_cols = F,
         treeheight_row = 0,
         treeheight_col = 0,
         show_rownames = F,
         fontsize_row = 8,
         fontsize_col = 8,
         filename = "Figures/Figure_Components/Supplement_Organ_Composition/percent_in_anatomy_slide_8.pdf",
         color = viridis(option = "inferno", n = 20),
         legend= F)



matrix_for_heatmap_3G = 
  annotation_mapped_cells %>%
  mutate(final_cluster_label = as.factor(final_cluster_label)) %>%
  filter(max_slide_id == "slide_3G") %>% 
  dplyr::select(-max_slide_id) %>%
  spread(key =organ, value =  percent_label, fill = 0, drop = F) %>%
  gather(key = "organ",
         value = "percent_label",
         -final_cluster_label) %>%
  spread(key =final_cluster_label, value =  percent_label, fill = 0, drop = F) %>%
  tibble::column_to_rownames(var = "organ") %>%
  as.matrix() 


pheatmap(t(matrix_for_heatmap_3G[row_order,cell_types_to_highlight]),
         cellwidth = 7,
         cellheight = 7,
         cluster_rows = F,
         cluster_cols = F,
         treeheight_row = 0,
         treeheight_col = 0,
         show_rownames = F,
         fontsize_row = 8,
         fontsize_col = 8,
         filename = "Figures/Figure_Components/Supplement_Organ_Composition/percent_in_anatomy_slide_9.pdf",
         color = viridis(option = "inferno", n = 20),
         legend= F)


matrix_for_heatmap_3H = 
  annotation_mapped_cells %>%
  mutate(final_cluster_label = as.factor(final_cluster_label)) %>%
  filter(max_slide_id == "slide_3H") %>% 
  dplyr::select(-max_slide_id) %>%
  spread(key =organ, value =  percent_label, fill = 0, drop = F) %>%
  gather(key = "organ",
         value = "percent_label",
         -final_cluster_label) %>%
  spread(key =final_cluster_label, value =  percent_label, fill = 0, drop = F) %>%
  tibble::column_to_rownames(var = "organ") %>%
  as.matrix() 


pheatmap(t(matrix_for_heatmap_3H[row_order,cell_types_to_highlight]),
         cellwidth = 7,
         cellheight = 7,
         cluster_rows = F,
         cluster_cols = F,
         treeheight_row = 0,
         treeheight_col = 0,
         show_rownames = F,
         fontsize_row = 8,
         fontsize_col = 8,
         filename = "Figures/Figure_Components/Supplement_Organ_Composition/percent_in_anatomy_slide_10.pdf",
         color = viridis(option = "inferno", n = 20),
         legend= F)



matrix_for_heatmap_4A = 
  annotation_mapped_cells %>%
  mutate(final_cluster_label = as.factor(final_cluster_label)) %>%
  filter(max_slide_id == "slide_4A") %>% 
  dplyr::select(-max_slide_id) %>%
  spread(key =organ, value =  percent_label, fill = 0, drop = F) %>%
  gather(key = "organ",
         value = "percent_label",
         -final_cluster_label) %>%
  spread(key =final_cluster_label, value =  percent_label, fill = 0, drop = F) %>%
  tibble::column_to_rownames(var = "organ") %>%
  as.matrix() 


pheatmap(t(matrix_for_heatmap_4A[row_order,cell_types_to_highlight]),
         cellwidth = 7,
         cellheight = 7,
         cluster_rows = F,
         cluster_cols = F,
         treeheight_row = 0,
         treeheight_col = 0,
         show_rownames = F,
         fontsize_row = 8,
         fontsize_col = 8,
         filename = "Figures/Figure_Components/Supplement_Organ_Composition/percent_in_anatomy_slide_11.pdf",
         color = viridis(option = "inferno", n = 20),
         legend= F)


matrix_for_heatmap_4D = 
  annotation_mapped_cells %>%
  mutate(final_cluster_label = as.factor(final_cluster_label)) %>%
  filter(max_slide_id == "slide_4D") %>% 
  dplyr::select(-max_slide_id) %>%
  spread(key =organ, value =  percent_label, fill = 0, drop = F) %>%
  gather(key = "organ",
         value = "percent_label",
         -final_cluster_label) %>%
  spread(key =final_cluster_label, value =  percent_label, fill = 0, drop = F) %>%
  tibble::column_to_rownames(var = "organ") %>%
  as.matrix() 


pheatmap(t(matrix_for_heatmap_4D[row_order,cell_types_to_highlight]),
         cellwidth = 7,
         cellheight = 7,
         cluster_rows = F,
         cluster_cols = F,
         treeheight_row = 0,
         treeheight_col = 0,
         show_rownames = F,
         fontsize_row = 8,
         fontsize_col = 8,
         filename = "Figures/Figure_Components/Supplement_Organ_Composition/percent_in_anatomy_slide_13.pdf",
         color = viridis(option = "inferno", n = 20),
         legend= F)


matrix_for_heatmap_4E = 
  annotation_mapped_cells %>%
  mutate(final_cluster_label = as.factor(final_cluster_label)) %>%
  filter(max_slide_id == "slide_4E") %>% 
  dplyr::select(-max_slide_id) %>%
  spread(key =organ, value =  percent_label, fill = 0, drop = F) %>%
  gather(key = "organ",
         value = "percent_label",
         -final_cluster_label) %>%
  spread(key =final_cluster_label, value =  percent_label, fill = 0, drop = F) %>%
  tibble::column_to_rownames(var = "organ") %>%
  as.matrix() 


pheatmap(t(matrix_for_heatmap_4E[row_order,cell_types_to_highlight]),
         cellwidth = 7,
         cellheight = 7,
         cluster_rows = F,
         cluster_cols = F,
         treeheight_row = 0,
         treeheight_col = 0,
         show_rownames = T,
         fontsize_row = 8,
         fontsize_col = 8,
         filename = "Figures/Figure_Components/Supplement_Organ_Composition/percent_in_anatomy_slide_14.pdf",
         color = viridis(option = "inferno", n = 20),
         legend= F)

# Figure 2 — Panel D ------------------------------------------------------
ggplot() +
  geom_point(data = 
               colData(spatial_cds) %>%
               as.data.frame() %>%
               dplyr::select(-anatomical_annotation),
             aes(x = umap1,
                 y = umap2),
             color = "grey80",
             stroke = 0,
             size = 0.25) +
  geom_point(data = 
               rbind(colData(spatial_cds) %>% as.data.frame() %>% mutate(order = "1"),
                     colData(spatial_cds) %>% as.data.frame() %>% mutate(order = "2")) %>%
               filter(anatomical_annotation != "Unannotated",
                      anatomical_annotation != "Thymus") %>%
               group_by(anatomical_annotation) %>%
               add_tally() %>%
               ungroup() %>%
               arrange(-n,anatomical_annotation,order) %>%
               mutate(anatomical_annotation = ifelse(order ==2,
                                                     anatomical_annotation %>% as.character(),
                                       NA)),
             aes(x = umap1,
                 y = umap2,
                 color = anatomical_annotation,
                 size = order),
             stroke = 0) +
  scale_size_manual(values = c("1" = 0.3,
                               "2" = 0.25)) +
  theme_void() +
  theme(legend.position = "none",
        strip.text = element_blank()) + 
  scale_color_brewer(palette = "Set1",na.value="black") +
  ggsave("Figures/Figure_Components/Figure2/umap_anatomical_annotation.png",
         dpi = 600,
         height = 3.5,
         width = 3.5)





# write out result files --------------------------------------------------


spatial_cds  = spatial_cds[,!colData(spatial_cds)$cell_to_remove]
saveRDS(object = spatial_cds,
        file = "Submission_Data/E14_slides/RDS_intermediates/Notebook6_spatial_cds_anatomy.RDS")




