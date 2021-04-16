# Visulalize the distribution of cell types and genes expression on the 
# embryo section. Subcluster the excitatory neurons and differentiating 
# mesenchyme -- visualize subclusters of the former on the embryo sections

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
  library(magrittr)
  library(monocle3)
  
  space_directory = "~/Google Drive File Stream/My Drive/sciSpace/"
  setwd(dir=space_directory)

  # Pass TRUE if you want to see progress output on some of Monocle 3's operations
  DelayedArray:::set_verbose_block_processing(TRUE)
  # Passing a higher value will make some computations faster but use more memory. 
  # Adjust with caution!
  options(DelayedArray.block.size=1000e6)
  
})

all_image_data = 
  readRDS(file = "Submission_Data/E14_slides/RDS_intermediates/Notebook0_1_images_and_transformations.RDS")

spatial_cds = 
  readRDS(file = "Submission_Data/E14_slides/RDS_intermediates/Notebook5.5_E14_spatial_CDS.RDS")


# Rename labels for display purposes --------------------------------------

colData(spatial_cds)$final_cluster_label =
  ifelse(colData(spatial_cds)$final_cluster_label == "Fibroblast",
         "Meninges",
         colData(spatial_cds)$final_cluster_label)

colData(spatial_cds)$final_cluster_label =
  ifelse(colData(spatial_cds)$final_cluster_label == "Neuron",
         "Neurons",
         colData(spatial_cds)$final_cluster_label)

colData(spatial_cds)$final_cluster_label =
  ifelse(colData(spatial_cds)$final_cluster_label == "Cardiac muscle lineages",
         "Cardiomyocytes",
         colData(spatial_cds)$final_cluster_label)


colData(spatial_cds)$final_cluster_label =
  ifelse(colData(spatial_cds)$final_cluster_label == "Peripheral Neuron",
         "Peripheral Neurons",
         colData(spatial_cds)$final_cluster_label)

colData(spatial_cds)$final_cluster_label =
  ifelse(colData(spatial_cds)$final_cluster_label == "Radial glia",
         "Radial Glia",
         colData(spatial_cds)$final_cluster_label)

colors =c("Cardiomyocytes" = "#FF3333",
          "Chondrocytes" = "#E1B0CF",
          "Choroid Plexus"  = "#00A6A6",
          "Connective Tissue Progenitors" = "#EBCF00",
          "Developing Gut" = "#0050D4",
          "Endothelial Cells" = "#64E8A6",
          "Epithelial Cells" = "#EB7700",
          "Erythroid Lineage" = "#F77C7C",
          "Meninges"  = "#790009",
          "Glial Cells" = "#8172E7",
          "Hepatocytes" = "#8538E8",
          "Lateral Plate Mesoderm" = "#086661",
          "Myocytes" = "#00D9F7",
          "Neurons" = "#B4D900",
          "OPCs" = "#E7E041",
          "Peripheral Neurons" = "#00B562",
          "Radial Glia" = "#FF00C3",
          "Schwann Cells" = "#E1AF52",
          "Testis Cells" = "black",
          "White Blood Cells" = "#3288BD")


# Supplemental Figure 17 — Slide 1 ----------------------------------------

# Segment each image by cropping based on contoured polygons
slide_1D_image_hull =
  all_image_data$dapi_image[[1]] %>%
  dplyr::select(x,
                y) %>%
  as.matrix() %>%
  SpatialPoints() %>%
  over(as(object = all_image_data$slide_polygon[[1]],
          Class = "Spatial")) %>%
  is.na() %>%
  magrittr::not()

slide_1D_image_cropped_by_hull =
  all_image_data$dapi_image[[1]][slide_1D_image_hull %>% as.vector,]



ggplot() +
  geom_raster(data = slide_1D_image_cropped_by_hull,
              aes(y = -x,
                  x = y,
                  fill = value)) +
  # Jitter the data to show overplotted poisions
  geom_jitter(data = 
                colData(spatial_cds) %>% 
                as.data.frame() %>% 
                filter(max_slide_id == "slide_1D") %>%
                filter(final_cluster_label %in% 
                         (colData(spatial_cds) %>% 
                            as.data.frame() %>% 
                            filter(max_slide_id == "slide_1D") %>%
                            group_by(final_cluster_label) %>%
                            summarise(n = n()) %>%
                            arrange(-n) %>%
                            head(n = 16) %>%
                            pull(final_cluster_label))),
             aes(y = -coords.x1, 
                 x = coords.x2,
                 color = final_cluster_label),
             width = 15,
             height = 15,
             stroke = 0,
             size = 0.35) +
  geom_sf(data = all_image_data$slide_polygon[[1]] * rbind(c(0, -1),c(1,0)),
          fill = NA, 
          color = "black", 
          size = .15) +
  geom_segment(data = 
                 all_image_data[1,] %>%
                 mutate(final_cluster_label= "White Blood Cells"),
               aes(x =sbX1,
                   xend = sbX2,
                   y = sbY,
                   yend = sbY),
               size = 0.25) +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  scale_color_manual(values = colors) +
  scale_fill_gradient(low = "white", high = "grey80") +
  facet_wrap(~final_cluster_label, nrow = 4)+
  theme(legend.position = "none",
        strip.text.x = element_text(size = 3))  +
  ggsave("Figures/Figure_Components/Supplement_celltype_facet_space/slide_1D_facet_space.png",
         dpi = 600,
         height  = 4, 
         width = 3)


# Supplemental Figure 17 — Slide 6 ----------------------------------------

# Segment each image by cropping based on contoured polygons
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
                  fill = value)) +
  # Jitter the data to show overplotted poisions
  geom_jitter(data = 
                colData(spatial_cds) %>% 
                as.data.frame() %>% 
                filter(max_slide_id == "slide_2H") %>%
                filter(final_cluster_label %in% 
                         (colData(spatial_cds) %>% 
                            as.data.frame() %>% 
                            filter(max_slide_id == "slide_2H") %>%
                            group_by(final_cluster_label) %>%
                            summarise(n = n()) %>%
                            arrange(-n) %>%
                            head(n = 16) %>%
                            pull(final_cluster_label))),
              aes(y = -coords.x1,
                  x = coords.x2,
                  color = final_cluster_label),
              width = 15,
              height = 15,
              stroke = 0,
              size = 0.35) +
  geom_sf(data = all_image_data$slide_polygon[[6]] * rbind(c(0, -1),c(1,0)),
          fill = NA,
          color = "black",
          size = .15) +
  geom_segment(data = 
                 all_image_data[6,] %>%
                 mutate(final_cluster_label= "White Blood Cells"),
               aes(x = sbX1,
                   xend = sbX2,
                   y = sbY + 50,
                   yend = sbY + 50),
               size = 0.25) +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  scale_color_manual(values = colors) +
  scale_fill_gradient(low = "white", high = "grey80") +
  facet_wrap(~final_cluster_label, ncol = 4)+
  theme(legend.position = "none",
        strip.text.x = element_text(size = 3)) +
  ggsave("Figures/Figure_Components/Supplement_celltype_facet_space/slide_2H_facet_space.png",
         dpi = 600,
         height  = 4,
         width = 3)

# Supplemental Figure 17 — Slide 3 ----------------------------------------
# Segment each image by cropping based on contoured polygons
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
                  fill = value)) +
  # Jitter the data to show overplotted poisions
  geom_jitter(data = 
                colData(spatial_cds) %>% 
                as.data.frame() %>% 
                filter(max_slide_id == "slide_1F") %>%
                filter(final_cluster_label %in% 
                         (colData(spatial_cds) %>% 
                            as.data.frame() %>% 
                            filter(max_slide_id == "slide_1F") %>%
                            group_by(final_cluster_label) %>%
                            summarise(n = n()) %>%
                            arrange(-n) %>%
                            head(n = 16) %>%
                            pull(final_cluster_label))),
              aes(y = -coords.x1,
                  x = -coords.x2,
                  color = final_cluster_label),
              width = 15,
              height = 15,
              stroke = 0,
              size = 0.35) +
  geom_sf(data = all_image_data$slide_polygon[[3]] * rbind(c(0, -1),c(-1,0)),
          fill = NA,
          color = "black",
          size = .15) +
  geom_segment(data = 
                 all_image_data[3,] %>%
                 mutate(final_cluster_label= "White Blood Cells"),
               aes(x = -sbX1 + 900,
                   xend = -sbX2 + 900,
                   y = sbY + 50,
                   yend = sbY + 50),
               size = 0.25) +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  scale_color_manual(values = colors) +
  scale_fill_gradient(low = "white", high = "grey80") +
  facet_wrap(~final_cluster_label, ncol = 4)+
  theme(legend.position = "none",
        strip.text.x = element_text(size = 3)) +
  ggsave("Figures/Figure_Components/Supplement_celltype_facet_space/slide_1F_facet_space.png",
         dpi = 600,
         height  = 4,
         width = 3)

# Supplemental Figure 17 — Slide 4 ----------------------------------------
# Segment each image by cropping based on contoured polygons
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
                  fill = value)) +
  # Jitter the data to show overplotted poisions
  geom_jitter(data = 
                colData(spatial_cds) %>% 
                as.data.frame() %>% 
                filter(max_slide_id == "slide_1G") %>%
                filter(final_cluster_label %in% 
                         (colData(spatial_cds) %>% 
                            as.data.frame() %>% 
                            filter(max_slide_id == "slide_1G") %>%
                            group_by(final_cluster_label) %>%
                            summarise(n = n()) %>%
                            arrange(-n) %>%
                            head(n = 16) %>%
                            pull(final_cluster_label))),
              aes(y = -coords.x1,
                  x = coords.x2,
                  color = final_cluster_label),
              width = 15,
              height = 15,
              stroke = 0,
              size = 0.35) +
  geom_sf(data = all_image_data$slide_polygon[[4]] * rbind(c(0, -1),c(1,0)),
          fill = NA,
          color = "black",
          size = .15) +
  geom_segment(data = 
                 all_image_data[4,] %>%
                 mutate(final_cluster_label= "White Blood Cells"),
               aes(x = sbX1,
                   xend = sbX2,
                   y = sbY + 25,
                   yend = sbY + 25),
               size = 0.25) +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  scale_color_manual(values = colors) +
  scale_fill_gradient(low = "white", high = "grey80") +
  facet_wrap(~final_cluster_label, ncol = 4)+
  theme(legend.position = "none",
        strip.text.x = element_text(size = 3)) +
  ggsave("Figures/Figure_Components/Supplement_celltype_facet_space/slide_1G_facet_space.png",
         dpi = 600,
         height  = 4,
         width = 3)

# Supplemental Figure 17 — Slide 5 ----------------------------------------
# Segment each image by cropping based on contoured polygons
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
                  fill = value)) +
  # Jitter the data to show overplotted poisions
  geom_jitter(data = 
                colData(spatial_cds) %>% 
                as.data.frame() %>% 
                filter(max_slide_id == "slide_2G") %>%
                filter(final_cluster_label %in% 
                         (colData(spatial_cds) %>% 
                            as.data.frame() %>% 
                            filter(max_slide_id == "slide_2G") %>%
                            group_by(final_cluster_label) %>%
                            summarise(n = n()) %>%
                            arrange(-n) %>%
                            head(n = 16) %>%
                            pull(final_cluster_label))),
              aes(y = -coords.x1,
                  x = coords.x2,
                  color = final_cluster_label),
              width = 15,
              height = 15,
              stroke = 0,
              size = 0.35) +
  geom_sf(data = all_image_data$slide_polygon[[5]] * rbind(c(0, -1),c(1,0)),
          fill = NA,
          color = "black",
          size = .15) +
  geom_segment(data = 
                 all_image_data[5,] %>%
                 mutate(final_cluster_label= "White Blood Cells"),
               aes(x = sbX1,
                   xend = sbX2,
                   y = sbY,
                   yend = sbY),
               size = 0.25) +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  facet_wrap(~final_cluster_label, ncol = 4)+
  scale_color_manual(values = colors) +
  scale_fill_gradient(low = "white", high = "grey80") +
  theme(legend.position = "none",
        strip.text.x = element_text(size = 3)) +
  ggsave("Figures/Figure_Components/Supplement_celltype_facet_space/slide_2G_facet_space.png",
         dpi = 600,
         height  = 4,
         width = 3)
 

# Supplemental Figure 17 — Slide 2 ----------------------------------------

slide_1E =
  colData(spatial_cds) %>%
  as.data.frame() %>%
  filter(max_slide_id == "slide_1E") 

slide_1E$in_contour = 
  slide_1E %>%
  dplyr::select(coords.x1,
                coords.x2) %>%
  as.matrix() %>%
  SpatialPoints() %>%
  over(as(object = all_image_data$slide_polygon[[2]],
          Class = "Spatial")) %>%
  is.na() %>%
  not()

# Segment each image by cropping based on contoured polygons
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
                  fill = value)) +
  # Jitter the data to show overplotted poisions
  geom_jitter(data =
                slide_1E %>%
                filter(final_cluster_label %in% 
                         (slide_1E %>%
                            group_by(final_cluster_label) %>%
                            summarise(n = n()) %>%
                            arrange(-n) %>%
                            head(n = 16) %>%
                            pull(final_cluster_label))),
              aes(y = -coords.x1,
                  x = -coords.x2,
                  color = final_cluster_label),
              width = 15,
              height = 15,
              stroke = 0,
              size = 0.35) +
  geom_sf(data = all_image_data$slide_polygon[[2]] * rbind(c(0, -1),c(-1,0)),
          fill = NA,
          color = "black",
          size = .15) +
  geom_segment(data = 
                 all_image_data[2,] %>%
                 mutate(final_cluster_label= "White Blood Cells"),
               aes(x = -sbX1 + 900,
                   xend = -sbX2 + 900,
                   y = sbY,
                   yend = sbY),
               size = 0.25) +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  scale_color_manual(values = colors) +
  scale_fill_gradient(low = "white", high = "grey80") +
  facet_wrap(~final_cluster_label, ncol = 4)+
  theme(legend.position = "none",
        strip.text.x = element_text(size = 3)) +
  ggsave("Figures/Figure_Components/Supplement_celltype_facet_space/slide_1E_facet_space.png",
         dpi = 600,
         height  = 4,
         width = 3)

# Supplemental Figure 17 — Slide 7 ----------------------------------------
# Segment each image by cropping based on contoured polygons
slide_3D_image_hull =
  all_image_data$dapi_image[[7]] %>%
  dplyr::select(x,
                y) %>%
  as.matrix() %>%
  SpatialPoints() %>%
  over(as(object = all_image_data$slide_polygon[[7]],
          Class = "Spatial")) %>%
  is.na() %>%
  magrittr::not()

slide_3D_image_cropped_by_hull =
  all_image_data$dapi_image[[7]][slide_3D_image_hull %>% as.vector,]

ggplot() +
  geom_raster(data = slide_3D_image_cropped_by_hull,
              aes(y = x,
                  x = y,
                  fill = value)) +
  # Jitter the data to show overplotted poisions
  geom_jitter(data = 
                colData(spatial_cds) %>% 
                as.data.frame() %>% 
                filter(max_slide_id == "slide_3D") %>%
                filter(final_cluster_label %in% 
                         (colData(spatial_cds) %>% 
                            as.data.frame() %>% 
                            filter(max_slide_id == "slide_3D") %>%
                            group_by(final_cluster_label) %>%
                            summarise(n = n()) %>%
                            arrange(-n) %>%
                            head(n = 16) %>%
                            pull(final_cluster_label))),
              aes(y = coords.x1, 
                  x = coords.x2,
                  color = final_cluster_label),
              width = 15,
              height = 15,
              stroke = 0,
              size = 0.35) +
  geom_sf(data = all_image_data$slide_polygon[[7]] * rbind(c(0, 1),c(1,0)),
          fill = NA, 
          color = "black", 
          size = .15) +
  geom_segment(data = 
                 all_image_data[7,] %>%
                 mutate(final_cluster_label= "White Blood Cells"),
               aes(x = sbX1,
                   xend = sbX2,
                   y = -sbY - 1500,
                   yend = -sbY - 1500),
               size = 0.25) +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  scale_color_manual(values = colors) +
  scale_fill_gradient(low = "white", high = "grey80") +
  facet_wrap(~final_cluster_label, nrow = 4)+
  theme(legend.position = "none",
        strip.text.x = element_text(size = 3)) +
  ggsave("Figures/Figure_Components/Supplement_celltype_facet_space/slide_3D_facet_space.png",
         dpi = 600,
         height  = 4, 
         width = 3)


# Supplemental Figure 17 — Slide 8 ----------------------------------------
# Segment each image by cropping based on contoured polygons
slide_3F_image_hull =
  all_image_data$dapi_image[[8]] %>%
  dplyr::select(x,
                y) %>%
  as.matrix() %>%
  SpatialPoints() %>%
  over(as(object = all_image_data$slide_polygon[[8]],
          Class = "Spatial")) %>%
  is.na() %>%
  magrittr::not()

slide_3F_image_cropped_by_hull =
  all_image_data$dapi_image[[8]][slide_3F_image_hull %>% as.vector,]

ggplot() +
  geom_raster(data = slide_3F_image_cropped_by_hull,
              aes(y = -x,
                  x = -y,
                  fill = value)) +
  # Jitter the data to show overplotted poisions
  geom_jitter(data = 
                colData(spatial_cds) %>% 
                as.data.frame() %>% 
                filter(max_slide_id == "slide_3F") %>%
                filter(final_cluster_label %in% 
                         (colData(spatial_cds) %>% 
                            as.data.frame() %>% 
                            filter(max_slide_id == "slide_3F") %>%
                            group_by(final_cluster_label) %>%
                            summarise(n = n()) %>%
                            arrange(-n) %>%
                            head(n = 16) %>%
                            pull(final_cluster_label))),              
              aes(y = -coords.x1, 
                  x = -coords.x2,
                  color = final_cluster_label),
              width = 15,
              height = 15,
              stroke = 0,
              size = 0.35) +
  geom_sf(data = all_image_data$slide_polygon[[8]] * rbind(c(0, -1),c(-1,0)),
          fill = NA, 
          color = "black", 
          size = .15) +
  geom_segment(data = 
                 all_image_data[8,] %>%
                 mutate(final_cluster_label= "White Blood Cells"),
               aes(x = -sbX1 + 850,
                   xend = -sbX2 + 850,
                   y = sbY + 25,
                   yend = sbY + 25),
               size = 0.25) +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  scale_color_manual(values = colors) +
  scale_fill_gradient(low = "white", high = "grey70") +
  facet_wrap(~final_cluster_label, nrow = 4)+
  theme(legend.position = "none",
        strip.text.x = element_text(size = 3)) +
  ggsave("Figures/Figure_Components/Supplement_celltype_facet_space/slide_3F_facet_space.png",
         dpi = 600,
         height  = 4, 
         width = 3)


# Supplemental Figure 17 — Slide 9 ----------------------------------------
# Segment each image by cropping based on contoured polygons
slide_3G_image_hull =
  all_image_data$dapi_image[[9]] %>%
  dplyr::select(x,
                y) %>%
  as.matrix() %>%
  SpatialPoints() %>%
  over(as(object = all_image_data$slide_polygon[[9]],
          Class = "Spatial")) %>%
  is.na() %>%
  magrittr::not()

slide_3G_image_cropped_by_hull =
  all_image_data$dapi_image[[9]][slide_3G_image_hull %>% as.vector,]

ggplot() +
  geom_raster(data = slide_3G_image_cropped_by_hull,
              aes(y = -x,
                  x = -y,
                  fill = value)) +
  # Jitter the data to show overplotted poisions
  geom_jitter(data = 
                colData(spatial_cds) %>% 
                as.data.frame() %>% 
                filter(max_slide_id == "slide_3G") %>%
                filter(final_cluster_label %in% 
                         (colData(spatial_cds) %>% 
                            as.data.frame() %>% 
                            filter(max_slide_id == "slide_3G") %>%
                            group_by(final_cluster_label) %>%
                            summarise(n = n()) %>%
                            arrange(-n) %>%
                            head(n = 16) %>%
                            pull(final_cluster_label))),
              aes(y = -coords.x1, 
                  x = -coords.x2,
                  color = final_cluster_label),
              width = 15,
              height = 15,
              stroke = 0,
              size = 0.35) +
  geom_sf(data = all_image_data$slide_polygon[[9]] * rbind(c(0, -1),c(-1,0)),
          fill = NA, 
          color = "black", 
          size = .15) +
  geom_segment(data = 
                 all_image_data[9,] %>%
                 mutate(final_cluster_label= "White Blood Cells"),
               aes(x = -sbX1 + 800,
                   xend = -sbX2 + 800,
                   y = sbY + 25,
                   yend = sbY  + 25),
               size = 0.25) +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  scale_fill_gradient(low = "white", high = "grey80") +
  scale_color_manual(values = colors) +
  facet_wrap(~final_cluster_label, nrow = 4)+
  theme(legend.position = "none",
        strip.text.x = element_text(size = 3)) +
  ggsave("Figures/Figure_Components/Supplement_celltype_facet_space/slide_3G_facet_space.png",
         dpi = 600,
         height  = 4, 
         width = 3)


# Supplemental Figure 17 — Slide 10 ----------------------------------------
# Segment each image by cropping based on contoured polygons
slide_3H_image_hull =
  all_image_data$dapi_image[[10]] %>%
  dplyr::select(x,
                y) %>%
  as.matrix() %>%
  SpatialPoints() %>%
  over(as(object = all_image_data$slide_polygon[[10]],
          Class = "Spatial")) %>%
  is.na() %>%
  magrittr::not()

slide_3H_image_cropped_by_hull =
  all_image_data$dapi_image[[10]][slide_3H_image_hull %>% as.vector,]

ggplot() +
  geom_raster(data = slide_3H_image_cropped_by_hull,
              aes(y = -x,
                  x = -y,
                  fill = value)) +
  # Jitter the data to show overplotted poisions
  geom_jitter(data = 
                colData(spatial_cds) %>% 
                as.data.frame() %>% 
                filter(max_slide_id == "slide_3H") %>%
                filter(final_cluster_label %in% 
                         (colData(spatial_cds) %>% 
                            as.data.frame() %>% 
                            filter(max_slide_id == "slide_3H") %>%
                            group_by(final_cluster_label) %>%
                            summarise(n = n()) %>%
                            arrange(-n) %>%
                            head(n = 16) %>%
                            pull(final_cluster_label))),
              aes(y = -coords.x1, 
                  x = -coords.x2,
                  color = final_cluster_label),
              width = 15,
              height = 15,
              stroke = 0,
              size = 0.35) +
  geom_sf(data = all_image_data$slide_polygon[[10]] * rbind(c(0, -1),c(-1,0)),
          fill = NA, 
          color = "black", 
          size = .15) +
  geom_segment(data = 
                 all_image_data[10,] %>%
                 mutate(final_cluster_label= "White Blood Cells"),
               aes(x = -sbX1 + 750,
                   xend = -sbX2 + 750,
                   y = sbY,
                   yend = sbY),
               size = 0.25) +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  scale_fill_gradient(low = "white", high = "grey80") +
  scale_color_manual(values = colors) +
  facet_wrap(~final_cluster_label, nrow = 4)+
  theme(legend.position = "none",
        strip.text.x = element_text(size = 3)) +
  ggsave("Figures/Figure_Components/Supplement_celltype_facet_space/slide_3H_facet_space.png",
         dpi = 600,
         height  = 4, 
         width = 3)


# Supplemental Figure 17 — Slide 11 ----------------------------------------
# Segment each image by cropping based on contoured polygons
slide_4A_image_hull =
  all_image_data$dapi_image[[11]] %>%
  dplyr::select(x,
                y) %>%
  as.matrix() %>%
  SpatialPoints() %>%
  over(as(object = all_image_data$slide_polygon[[11]],
          Class = "Spatial")) %>%
  is.na() %>%
  magrittr::not()

slide_4A_image_cropped_by_hull =
  all_image_data$dapi_image[[11]][slide_4A_image_hull %>% as.vector,]

ggplot() +
  geom_raster(data = slide_4A_image_cropped_by_hull,
              aes(y = -x,
                  x = -y,
                  fill = value)) +
  #Jitter the data to show overplotted poisions
  geom_jitter(data = 
                colData(spatial_cds) %>% 
                as.data.frame() %>% 
                filter(max_slide_id == "slide_4A") %>%
                filter(final_cluster_label %in% 
                         (colData(spatial_cds) %>% 
                            as.data.frame() %>% 
                            filter(max_slide_id == "slide_4A") %>%
                            group_by(final_cluster_label) %>%
                            summarise(n = n()) %>%
                            arrange(-n) %>%
                            head(n = 16) %>%
                            pull(final_cluster_label))),
              aes(y = -coords.x1, 
                  x = -coords.x2,
                  color = final_cluster_label),
              width = 15,
              height = 15,
              stroke = 0,
              size = 0.35) +
  geom_sf(data = all_image_data$slide_polygon[[11]] * rbind(c(0, -1),c(-1,0)),
          fill = NA, 
          color = "black", 
          size = .15) +
  geom_segment(data = 
                 all_image_data[11,] %>%
                 mutate(final_cluster_label= "White Blood Cells"),
               aes(x = -sbX1 + 750,
                   xend = -sbX2 + 750,
                   y = sbY,
                   yend = sbY),
               size = 0.25) +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  scale_color_manual(values = colors) +
  scale_fill_gradient(low = "white", high = "grey80") +
  facet_wrap(~final_cluster_label, nrow = 4)+
  theme(legend.position = "none",
        strip.text.x = element_text(size = 3)) +
  ggsave("Figures/Figure_Components/Supplement_celltype_facet_space/slide_4A_facet_space.png",
         dpi = 600,
         height  = 4, 
         width = 3)


# Supplemental Figure 17 — Slide 13 ----------------------------------------
# Segment each image by cropping based on contoured polygons
slide_4D_image_hull =
  all_image_data$dapi_image[[13]] %>%
  dplyr::select(x,
                y) %>%
  as.matrix() %>%
  SpatialPoints() %>%
  over(as(object = all_image_data$slide_polygon[[13]],
          Class = "Spatial")) %>%
  is.na() %>%
  magrittr::not()

slide_4D_image_cropped_by_hull =
  all_image_data$dapi_image[[13]][slide_4D_image_hull %>% as.vector,]

ggplot() +
  geom_raster(data = slide_4D_image_cropped_by_hull,
              aes(y = x,
                  x = y,
                  fill = value)) +
  # Jitter the data to show overplotted poisions
  geom_jitter(data = 
                colData(spatial_cds) %>% 
                as.data.frame() %>% 
                filter(max_slide_id == "slide_4D") %>%
                filter(final_cluster_label %in% 
                         (colData(spatial_cds) %>% 
                            as.data.frame() %>% 
                            filter(max_slide_id == "slide_4D") %>%
                            group_by(final_cluster_label) %>%
                            summarise(n = n()) %>%
                            arrange(-n) %>%
                            head(n = 16) %>%
                            pull(final_cluster_label))),
              aes(y = coords.x1, 
                  x = coords.x2,
                  color = final_cluster_label),
              width = 15,
              height = 15,
              stroke = 0,
              size = 0.35) +
  geom_segment(data = 
                 all_image_data[13,] %>%
                 mutate(final_cluster_label= "White Blood Cells"),
               aes(x = sbX1,
                   xend = sbX2,
                   y = -sbY - 1500,
                   yend = -sbY - 1500),
               size = 0.25) +
  
  geom_sf(data = all_image_data$slide_polygon[[13]] * rbind(c(0, 1),c(1,0)),
          fill = NA, 
          color = "black", 
          size = .15) +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  scale_fill_gradient(low = "white", high = "grey80") +
  scale_color_manual(values = colors) +
  facet_wrap(~final_cluster_label, nrow = 4)+
  theme(legend.position = "none",
        strip.text.x = element_text(size = 3)) +
  ggsave("Figures/Figure_Components/Supplement_celltype_facet_space/slide_4D_facet_space.png",
         dpi = 600,
         height  = 4, 
         width = 3)


# Supplemental Figure 17 — Slide 14 ----------------------------------------
# Segment each image by cropping based on contoured polygons
slide_4E_image_hull =
  all_image_data$dapi_image[[14]] %>%
  dplyr::select(x,
                y) %>%
  as.matrix() %>%
  SpatialPoints() %>%
  over(as(object = all_image_data$slide_polygon[[14]],
          Class = "Spatial")) %>%
  is.na() %>%
  magrittr::not()

slide_4E_image_cropped_by_hull =
  all_image_data$dapi_image[[14]][slide_4E_image_hull %>% as.vector,]

ggplot() +
  geom_raster(data = slide_4E_image_cropped_by_hull,
               aes(y = x,
                  x = y,
                   fill = value)) +
  # Jitter the data to show overplotted poisions
  geom_jitter(data = 
                colData(spatial_cds) %>% 
                as.data.frame() %>% 
                filter(max_slide_id == "slide_4E") %>%
                filter(final_cluster_label %in% 
                         (colData(spatial_cds) %>% 
                            as.data.frame() %>% 
                            filter(max_slide_id == "slide_4E") %>%
                            group_by(final_cluster_label) %>%
                            summarise(n = n()) %>%
                            arrange(-n) %>%
                            head(n = 16) %>%
                            pull(final_cluster_label))),
              aes(y = coords.x1, 
                  x = coords.x2,
                  color = final_cluster_label),
              width = 15,
              height = 15,
              stroke = 0,
              size = 0.35) +
  geom_sf(data = all_image_data$slide_polygon[[14]] * rbind(c(0, 1),c(1,0)),
          fill = NA, 
          color = "black", 
          size = .15) +
  geom_segment(data = 
                 all_image_data[14,] %>%
                 mutate(final_cluster_label= "White Blood Cells"),
               aes(x = sbX1,
                   xend = sbX2,
                   y = -sbY - 1500,
                   yend = -sbY - 1500),
               size = 0.25) +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  scale_fill_gradient(low = "white", high = "grey80") +
  scale_color_manual(values = colors) +
  facet_wrap(~final_cluster_label, nrow = 4)+
  theme(legend.position = "none",
        strip.text.x = element_text(size = 3)) +
  ggsave("Figures/Figure_Components/Supplement_celltype_facet_space/slide_4E_facet_space.png",
         dpi = 600,
         height  = 4, 
         width = 3)

# Figure 2 — Panel C — White blood cells -------------------------------------
ggplot() +
  geom_raster(data = slide_4E_image_cropped_by_hull,
              aes(y = x,
                  x = y,
                  fill = value),
              alpha = 0.9) +
  # Jitter the data to show overplotted poisions
  geom_jitter(data = 
                colData(spatial_cds) %>% 
                as.data.frame() %>% 
                filter(max_slide_id == "slide_4E") %>%
                filter(final_cluster_label == "White Blood Cells"),
              aes(y = coords.x1, 
                  x = coords.x2),
              color = "#3288BD",
              width = 15,
              height = 15,
              stroke = 0,
              size = 0.55) +
  geom_sf(data = all_image_data$slide_polygon[[14]] * rbind(c(0, 1),c(1,0)),
          fill = NA, 
          color = "black", 
          size = .125) +
  geom_segment(data = 
                 all_image_data[14,],
               aes(x = sbX1,
                   xend = sbX2,
                   y = -sbY - 1500,
                   yend = -sbY - 1500),
               size = 0.25) +
  monocle3:::monocle_theme_opts() +
  scale_fill_gradient(low = "white", high = "grey70") +
  theme_void() +
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Figure2/Slide_4E_white_blood_cells.png",
         dpi = 600,
         height  = 1.25, 
         width = 1.25)




# Figure 2 — Panel C — Cardiomyocytes -------------------------------------
# Cardiac muscle 
ggplot() +
  geom_raster(data = slide_4E_image_cropped_by_hull,
              aes(y = x,
                  x = y,
                  fill = value),
              alpha = 0.9) +
  # Jitter the data to show overplotted poisions
  geom_jitter(data = 
                colData(spatial_cds) %>% 
                as.data.frame() %>% 
                filter(max_slide_id == "slide_4E") %>%
                filter(final_cluster_label == "Cardiomyocytes"),
              aes(y = coords.x1, 
                  x = coords.x2),
              color = "red",
              width = 15,
              height = 15,
              stroke = 0,
              size = 0.55) +
  geom_sf(data = all_image_data$slide_polygon[[14]] * rbind(c(0, 1),c(1,0)),
          fill = NA, 
          color = "black", 
          size = .125) +
  geom_segment(data = 
                 all_image_data[14,],
               aes(x = sbX1,
                   xend = sbX2,
                   y = -sbY - 1500,
                   yend = -sbY - 1500),
               size = 0.25) +
  monocle3:::monocle_theme_opts() +
  scale_fill_gradient(low = "white", high = "grey70") +
  theme_void() +
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Figure2/Slide_4E_cardios.png",
         dpi = 600,
         height  = 1.25, 
         width = 1.25)

# Figure 2 — Panel C — Meninges -------------------------------------

ggplot() +
  geom_raster(data = slide_4E_image_cropped_by_hull,
              aes(y = x,
                  x = y,
                  fill = value),
              alpha = 0.9) +
  # Jitter the data to show overplotted poisions
  geom_jitter(data = 
                colData(spatial_cds) %>% 
                as.data.frame() %>% 
                filter(max_slide_id == "slide_4E") %>%
                filter(final_cluster_label == "Meninges"),
              aes(y = coords.x1, 
                  x = coords.x2),
              color = "#790009",
              width = 15,
              height = 15,
              stroke = 0,
              size = 0.55) +
  geom_sf(data = all_image_data$slide_polygon[[14]] * rbind(c(0, 1),c(1,0)),
          fill = NA, 
          color = "black", 
          size = .125) +
  geom_segment(data = 
                 all_image_data[14,],
               aes(x = sbX1,
                   xend = sbX2,
                   y = -sbY - 1500,
                   yend = -sbY - 1500),
               size = 0.25) +
  monocle3:::monocle_theme_opts() +
  scale_fill_gradient(low = "white", high = "grey70") +
  theme_void() +
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Figure2/Slide_4E_VLMC.png",
         dpi = 600,
         height  = 1.25, 
         width = 1.25)

# Figure 2 — Panel C — Neurons -------------------------------------

ggplot() +
  geom_raster(data = slide_4E_image_cropped_by_hull,
              aes(y = x,
                  x = y,
                  fill = value),
              alpha = 0.9) +
  # Jitter the data to show overplotted poisions
  geom_jitter(data = 
                colData(spatial_cds) %>% 
                as.data.frame() %>% 
                filter(max_slide_id == "slide_4E") %>%
                filter(final_cluster_label == "Neurons"),
              aes(y = coords.x1, 
                  x = coords.x2),
              color = "#B4D900",
              width = 15,
              height = 15,
              stroke = 0,
              size = 0.55) +
  geom_sf(data = all_image_data$slide_polygon[[14]] * rbind(c(0, 1),c(1,0)),
          fill = NA, 
          color = "black", 
          size = .125) +
  geom_segment(data = 
                 all_image_data[14,],
               aes(x = sbX1,
                   xend = sbX2,
                   y = -sbY - 1500,
                   yend = -sbY - 1500),
               size = 0.25) +
  monocle3:::monocle_theme_opts() +
  scale_fill_gradient(low = "white", high = "grey70") +
  theme_void() +
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Figure2/Slide_4E_Neuron.png",
         dpi = 600,
         height  = 1.25, 
         width = 1.25)

# Figure 2 —Panel E — Virtual in situ -------------------------------------

# pick out marker and all cells in that slide
markers = c("Slc6a3")
cds_1D = spatial_cds[,colData(spatial_cds)$max_slide_id == "slide_1D"]
cds_subset = cds_1D[rowData(cds_1D)$gene_short_name %in% markers,]

# size factor normalize expression of the marker
marker_expr = 
  (Matrix::t(exprs(cds_subset))/size_factors(cds_subset)) %>% 
  as.matrix() %>%
  as.data.frame()

marker_expr[marker_expr == 0] <- NA

colnames(marker_expr) = 
  rowData(cds_1D)[colnames(marker_expr),"gene_short_name"]

marker_expr =
  marker_expr %>%
  rownames_to_column(var = "Cell") %>%
  gather(key = "gene",
         value = "expression",
         -Cell)

# Append coordinates of expression based on the cell's name
coordinates =
  data.frame(Cell = colData(cds_1D)$Cell,
             spatial_1 = colData(cds_1D)$coords.x1,
             spatial_2 = colData(cds_1D)$coords.x2,
             umap1 = colData(cds_1D)$umap1,
             umap2 = colData(cds_1D)$umap2)

marker_expr = 
  left_join(marker_expr, 
            coordinates,
            by = "Cell")

# Plot expression of the gene through the whole section 
ggplot() +
  geom_raster(data = slide_1D_image_cropped_by_hull,
              aes(y = -x,
                  x = y, 
                  fill = value),
              alpha = 0.25) +
  geom_sf(data = all_image_data$slide_polygon[[1]] * rbind(c(0, -1),c(1,0)),
          fill = NA,
          color = "black",
          size = .125) +
  geom_jitter(data = marker_expr %>%
                filter(!is.na(expression)),
              aes(x =  spatial_2,
                  y = -spatial_1,
                  color = expression),
              stroke = 0,
              size =.65,
              width = 30,
              height = 30) +
  scale_color_viridis_c(option = "inferno") +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  scale_fill_gradient(low = "white", high = "grey60") +
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Figure2/Slide_1D_slc6a3.png",
         dpi = 600,
         height  = 1.25, 
         width = 1.25)

# Zoomed in version to visualize expression in the brain
ggplot() +
  geom_raster(data = slide_1D_image_cropped_by_hull,
              aes(y = -x,
                  x = y, 
                  fill = value),
              alpha = 0.25) +
  geom_sf(data = all_image_data$slide_polygon[[1]] * rbind(c(0, -1),c(1,0)),
          fill = NA,
          color = "black",
          size = .15) +
  geom_jitter(data = marker_expr %>%
                filter(!is.na(expression)),
              aes(x =  spatial_2,
                  y = -spatial_1,
                  color = expression),
              stroke = 0,
              size =.85,
              width = 30,
              height = 30) +
  scale_color_viridis_c(option = "inferno") +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  scale_fill_gradient(low = "white", high = "grey60") +
  theme(legend.position = "none") +
  xlim(500,1100) +
  ylim(-1100,-500) +
  ggsave("Figures/Figure_Components/Figure2/Slide_1D_slc6a3_zoom.png",
         dpi = 600,
         height  = 1, 
         width = 1)
    


# Supplemental Figure 22 ----------------------------------------
# Brain Regions from LaManno et. al. 

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

ggplot(colData(spatial_cds) %>%
         as.data.frame() %>%
         left_join(scaling_df,
                   by = "max_slide_id") %>%
         filter(anatomical_annotation == "Cortex",
                !is.na(lamanno_Tissue))) +
  geom_jitter(aes(y = coords.x1 * y,
                  x = coords.x2 * x,
                  color = lamanno_Tissue),
              size = 0.65,
              height = 15,
              width = 15,
              stroke = 0) +
  facet_wrap(~slide_id,scales = "free", nrow = 3) +
  labs(color = "Tissue \n(LaManno et. al.)") +
  monocle3:::monocle_theme_opts()+
  theme_void() +
  theme(strip.text.x = element_text(hjust = 0)) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  ggsave("Figures/Figure_Components/Supplement_LaManno_Integration/brain_region.pdf",
         height = 4,
         width = 7)

ggplot() +
  geom_sf(data = all_image_data$slide_polygon[[14]] *  rbind(c(0, 1),c(1,0))) +
  theme_void() +
  ggsave("Figures/Figure_Components/Supplement_LaManno_Integration/brain_region_example.pdf",
         height = 1,
         width = 1)



