# Visulalize the alignment of the spatial grid on the imaged embryo section

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
  readRDS("Submission_Data/E14_slides/RDS_intermediates/Notebook0_1_images_and_transformations.RDS")

oligo_plate_layout = 
  read.table(file = "Submission_Data/Oligo_layout/oligo_plate_layout.tsv", 
             sep = '\t', 
             header = TRUE)

spatial_cds = 
  readRDS(file = "Submission_Data/E14_slides/RDS_intermediates/Notebook6_spatial_cds_anatomy.RDS")


# Figure 2 — Panel A ------------------------------------------------------

ggplot() +
  geom_raster(data = all_image_data$dapi_image[[3]],
              aes(y = -x,
                  x = -y, 
                  fill = value)) +
  geom_point(data = all_image_data$transformed_oligo_layout[[3]] ,
             aes(y = -coords.x1,
                 x = -coords.x2,
                 color = SYBR),
             stroke = 0,
             size = 1,
             alpha = 0.5) +
  geom_sf(data = all_image_data$slide_polygon[[3]] * rbind(c(0, -1),c(-1,0)) ,
          fill = NA, 
          color = "white", 
          size = .25) +
  scale_color_manual(values = c("black","green")) +
  theme_void() +
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Figure2/Slide_1F_alignment.png",
         dpi = 600,
         height = 3,
         width = 3)


# Supplemental Figure 13 — Panel B -----------------------------------------

# Make plots for supplement showing grid alignment
ggplot() +
  geom_point(data = 
               oligo_plate_layout %>% 
               filter(!SYBR), 
             aes(x = -Col, 
                 y = -Row),
             color = "black", 
             stroke = 0,
             alpha = 0.5,
             size = .65 ) +
  geom_point(data = 
               oligo_plate_layout %>% 
               filter(SYBR), 
             aes(x = -Col, 
                 y = -Row ), 
             color = "green2", 
             size = 1.25, 
             stroke = 0) +
  theme_void() +
  theme(legend.position = "none") + 
  ggsave("Figures/Figure_Components/Supplement_Grid_Alignment/sybr_grid.png",
         height = 3, 
         width = 3, 
         units = "in", 
         dpi = 550)

# Supplemental Figure 13 — Panel C -----------------------------------------

ggplot() +
  geom_raster(data = all_image_data$sybr_image[[3]],
              aes(x = -y,
                  y = -x, 
                  fill = value)) +
  geom_point(data = all_image_data$transformed_oligo_layout[[3]],
             aes(x = -coords.x2,
                 y = -coords.x1,
                 color = SYBR),
             stroke = 0,
             size = 1,
             alpha = 0.5) +
  geom_sf(data = all_image_data$slide_polygon[[3]] * rbind(c(0, -1),c(-1,0)),
          fill = NA, 
          color = "white", 
          size = .25) +
  scale_color_manual(values = c("black","green")) +
  scale_fill_gradient(low = "black", high = "green") +
  theme_void() +
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Supplement_Grid_Alignment/Slide_1F_alignment.png",
         dpi = 600,
         height = 3,
         width = 3)

# Supplemental Figure 13 — Panel D -----------------------------------------

ggplot() +
  geom_raster(data = all_image_data$sybr_image[[3]],
              aes(x = -y,
                  y = -x, 
                  fill = value),
              alpha = 0.5) +
  geom_point(data = all_image_data$transformed_oligo_layout[[3]],
             aes(x = -coords.x2,
                 y = -coords.x1,
                 color = SYBR),
             stroke = 0,
             size = 6,
             alpha = 0.5,
             fill = NA) +
  geom_sf(data = all_image_data$slide_polygon[[3]] * rbind(c(0, -1),c(-1,0)),
          fill = NA, 
          color = "white", 
          size = 1.5) +
  scale_color_manual(values = c("black","green")) +
  scale_fill_gradient(low = "black", high = "green") +
  ylim(-1500,-1200) +
  xlim(-900,-600) +
  theme_void() +
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Supplement_Grid_Alignment/Slide_1F_alignment_zoom.pdf",
         height = 3,
         width = 3)

