# Visulalize the distribution of genes expression on the embryo section

# Load startup packages ---------------------------------------------------
suppressPackageStartupMessages({ 
  library(tidyverse)
  library(ggplot2)
  library(purrr)
  library(sp)
  library(sf)
  library(monocle3)
  
  space_directory = "/Volumes/GoogleDrive/My Drive/sciSpace/"
  setwd(dir=space_directory)

  set.seed(42)
})

spatial_cds = 
  readRDS(file = "Submission_Data/E14_slides/RDS_intermediates/Notebook6_spatial_cds_anatomy.RDS")

all_image_data = 
  readRDS(file = "Submission_Data/E14_slides/RDS_intermediates/Notebook0_1_images_and_transformations.RDS")


# Function to aggregrate marker expression --------------------------------

aggregate_markers = 
  function(markers,cds){
    cds_subset = cds[rowData(cds)$gene_short_name %in% markers,]
    
    marker_expr = 
      Matrix::t(exprs(cds_subset)) %>% 
      as.matrix() %>%
      as.data.frame()

    marker_expr[marker_expr == 0] <- NA

    colnames(marker_expr) = 
      rowData(cds)[colnames(marker_expr),"gene_short_name"]

    marker_expr =
      marker_expr %>%
      rownames_to_column(var = "Cell") %>%
      gather(key = "gene",
             value = "expression",
             -Cell)

    coordinates =
      data.frame(Cell = colData(cds)$Cell,
                 spatial_1 = colData(cds)$coords.x1,
                 spatial_2 = colData(cds)$coords.x2,
                 umap1 = colData(cds)$umap1,
                 umap2 = colData(cds)$umap2)

    marker_expr = 
      left_join(marker_expr, 
                coordinates,
                by = "Cell")

}

# Function to crop an image based on a polygon shape file -----------------

crop_hull = 
  function(im_df,
           slide_name){
    
    image_hull =
      im_df[which(im_df$slide_id == slide_name),"dapi_image"][[1]][[1]] %>%
      dplyr::select(x,
                    y) %>%
      as.matrix() %>%
      SpatialPoints() %>%
      over(as(object = im_df[which(im_df$slide_id == slide_name),"slide_polygon"][[1]][[1]],
              Class = "Spatial")) %>%
      is.na() %>%
      not()
    
    return(im_df[which(im_df$slide_id == slide_name),"dapi_image"][[1]][[1]][image_hull %>% as.vector,])
  }

# Markers to visualize ----------------------------------------------------

markers =c("Pax1",
           "Pax2",
           "Pitx1",
           "Neurod6")
# Supplemental Figure 15 — Slide 1 ----------------------------------------


slide_1D_marker_exprs =
  aggregate_markers(markers,
                  spatial_cds[,colData(spatial_cds)$max_slide_id == "slide_1D"])

ggplot() +
  geom_sf(data = all_image_data$slide_polygon[[1]] * rbind(c(0, -1),c(1,0)),
          fill = "grey90",
          color = "black",
          size = .1) +
  geom_jitter(data = slide_1D_marker_exprs %>%
                filter(expression <= 10),
              aes(x =  spatial_2,
                  y = -spatial_1,
                  color = expression),
              stroke = 0,
              size =.45,
              width = 12.5,
              height = 12.5) +
  geom_jitter(data = slide_1D_marker_exprs %>%
                filter(expression > 10),
              aes(x =  spatial_2,
                  y = -spatial_1),
              color = "red",
              stroke = 0,
              size =.45,
              width = 12.5,
              height = 12.5) +
  geom_segment(data = 
                 all_image_data[1,],
               aes(x =sbX1,
                   xend = sbX2,
                   y = sbY,
                   yend = sbY),
               size = 0.25) +
  scale_color_viridis_c(option = "inferno",limits = c(0,10)) +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  facet_wrap(~gene,ncol = 1) + 
  scale_fill_gradient(low = "white", high = "grey30") +
  theme(legend.position = "none",
        strip.text = element_blank()) +
  ggsave("Figures/Figure_Components/Supplement_ISH_versus_sciSpace/set1_slide_1D.pdf",
         height = 6,
         width = 1.5)

# Supplemental Figure 15 Slide 14 ------------------------------------------


slide_4E_marker_exprs =
  aggregate_markers(markers,
                    spatial_cds[,colData(spatial_cds)$max_slide_id == "slide_4E"])

ggplot() +
  geom_sf(data = all_image_data$slide_polygon[[14]] * rbind(c(0, 1),c(1,0)),
          fill = "grey90",
          color = "black",
          size = .1) +

  geom_jitter(data = slide_4E_marker_exprs %>%
                filter(expression <= 10),
              aes(x =  spatial_2,
                  y = spatial_1,
                  color = expression),
              stroke = 0,
              size =.45,
              width = 12.5,
              height = 12.5) +
  geom_jitter(data = slide_4E_marker_exprs %>%
                filter(expression > 10),
              aes(x =  spatial_2,
                  y = spatial_1),
              color = "red",
              stroke = 0,
              size =.45,
              width = 12.5,
              height = 12.5) +
  geom_segment(data = 
                 all_image_data[14,],
               aes(x = sbX1,
                   xend = sbX2,
                   y = -sbY - 1500,
                   yend = -sbY - 1500),
               size = 0.25) +
  scale_color_viridis_c(option = "inferno",limits = c(0,10)) +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  facet_wrap(~gene,ncol = 1) + 
  scale_fill_gradient(low = "white", high = "grey30") +
  theme(legend.position = "none",
        strip.text = element_blank()) +
  ggsave("Figures/Figure_Components/Supplement_ISH_versus_sciSpace/set1_slide_4E.pdf",
         height = 6,
         width = 1.5)

# Supplemental Figure 15 Slide 11 ------------------------------------------

slide_4A_marker_exprs =
  aggregate_markers(markers,
                    spatial_cds[,colData(spatial_cds)$max_slide_id == "slide_4A"])


ggplot() +
  geom_sf(data = all_image_data$slide_polygon[[11]] * rbind(c(0, -1),c(-1,0)),
          fill = "grey90",
          color = "black",
          size = .1) +
  geom_jitter(data = slide_4A_marker_exprs %>%
                filter(expression <= 10),
              aes(x =  -spatial_2,
                  y = -spatial_1,
                  color = expression),
              stroke = 0,
              size =.45,
              width = 12.5,
              height = 12.5) +
  geom_jitter(data = slide_4A_marker_exprs %>%
                filter(expression > 10),
              aes(x =  -spatial_2,
                  y = -spatial_1),
              color = "red",
              stroke = 0,
              size =.45,
              width = 12.5,
              height = 12.5) +
  geom_segment(data = 
                 all_image_data[11,],
               aes(x = -sbX1 + 750,
                   xend = -sbX2 + 750,
                   y = sbY,
                   yend = sbY),
               size = 0.25) +
  scale_color_viridis_c(option = "inferno",limits = c(0,10)) +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  facet_wrap(~gene,ncol = 1) + 
  theme(legend.position = "none",
        strip.text = element_blank()) +
  ggsave("Figures/Figure_Components/Supplement_ISH_versus_sciSpace/set1_slide_4A.pdf",
         height = 6,
         width = 1.5)

# Supplemental Figure 15 Slide 8 ------------------------------------------


slide_3F_marker_exprs =
  aggregate_markers(markers,
                    spatial_cds[,colData(spatial_cds)$max_slide_id == "slide_3F"])


ggplot() +
  geom_sf(data = all_image_data$slide_polygon[[8]] * rbind(c(0, -1),c(-1,0)),
          fill = "grey90",
          color = "black",
          size = .1) +
  geom_jitter(data = slide_3F_marker_exprs %>%
                filter(expression <= 10),
              aes(x =  -spatial_2,
                  y = -spatial_1,
                  color = expression),
              stroke = 0,
              size =.45,
              width = 12.5,
              height = 12.5) +
  geom_jitter(data = slide_3F_marker_exprs %>%
                filter(expression > 10),
              aes(x =  -spatial_2,
                  y = -spatial_1),
              color = "red",
              stroke = 0,
              size =.45,
              width = 12.5,
              height = 12.5) +
  geom_segment(data = 
                 all_image_data[8,],
               aes(x = -sbX1 + 850,
                   xend = -sbX2 + 850,
                   y = sbY + 25,
                   yend = sbY + 25),
               size = 0.25) +
  scale_color_viridis_c(option = "inferno",limits = c(0,10)) +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  facet_wrap(~gene,ncol = 1) + 
  scale_fill_gradient(low = "white", high = "grey30") +
  theme(legend.position = "none",
        strip.text = element_blank()) +
  ggsave("Figures/Figure_Components/Supplement_ISH_versus_sciSpace/set1_slide_3F.pdf",
         height = 6,
         width = 1.5)

  
