# Position refinement 3. Ask if single cells apart from other cells of a cluster have
# a comparable mapping near other cells of the same type

# Load startup packages ---------------------------------------------------
suppressPackageStartupMessages({ 
  library(tidyverse)
  library(ggplot2)
  library(tidyr)
  library(purrr)
  library(spatstat)
  library(imager)
  library(vec2dtransf)
  library(sp)
  library(sf)
  library(monocle3)
  library(magrittr)
  library(OpenImageR)
  
  space_directory = "~/Google Drive File Stream/My Drive/sciSpace/"
  setwd(dir=space_directory)
  source("Submission_Data/bin/hotspot_functions.R")
  
  # Pass TRUE if you want to see progress output on some of Monocle 3's operations
  DelayedArray:::set_verbose_block_processing(TRUE)
  # Passing a higher value will make some computations faster but use more memory. 
  # Adjust with caution!
  options(DelayedArray.block.size=1000e6)
  
})

load("Submission_Data/E14_slides/RDS_intermediates/Notebook5_initial_mapping.RData")


updated_coldata =
  readRDS(file = "Submission_Data/E14_slides/RDS_intermediates/Notebook5.3_updated_coldata.RDS")


flat_blur_3x3 = 
  matrix(c(1,1,1,
           1,1,1,
           1,1,1),
         nrow = 3)

density_df = 
  updated_coldata %>%
  mutate(Row = as.character(Row),
         Col = as.character(Col)) %>%
  group_by(max_slide_id,
           cluster) %>%
  nest() %>%
  mutate(test = 
           purrr::map(.x = data,
                      .f = function(cd){
                        long = 
                          cd %>% 
                          mutate(Row = as.numeric(Row),
                                 Col = as.numeric(Col)) %>%
                          group_by(Row,
                                   Col) %>%
                          summarise(n = n()) 
                        
                        mat = Matrix::sparseMatrix(i = long$Row,
                                                   j = long$Col,
                                                   x = long$n,
                                                   dims = c(84,84)) %>%
                          Matrix::as.matrix()
                        
                        res =
                          convolution(mat,
                                      kernel = flat_blur_3x3,
                                      mode = "same")
                        
                        res
                      }))


density_df$max = 
  sapply(X = seq(1,dim(density_df)[1]),
         function(X){
           density_df[X,"test"][[1]][[1]] %>% max()
         })

density_df =
  updated_coldata %>%
  mutate(Row = as.character(Row),
         Col = as.character(Col)) %>%
  group_by(cluster, max_slide_id) %>%
  nest() %>%
  drop_na() %>%
  left_join(density_df %>%
              dplyr::select(cluster,
                            max_slide_id,
                            test,
                            max)) %>% 
  filter(max > 3) %>%
  mutate(in_density = 
           purrr:::map2(.x = data,
                        .y = test,
                        .f = function(coldata,
                                      celltype_density){
                          mask = 
                            reshape2::melt(celltype_density,
                                 varnames = c("Row","Col")) %>%
                            filter(value > 3)
                          
                          coldata %>%
                            mutate(Row = as.numeric(Row),
                                   Col = as.numeric(Col)) %>%
                            dplyr::select(Cell,
                                          Row,
                                          Col) %>%
                            left_join(mask,
                                      by = c("Row","Col")) %>%
                            mutate(in_mask = !is.na(value))
                          
                        })) %>%
  unnest(cols = in_density)

density_df$in_mask %>% sum()

cells_to_exclude_from_density_mapping =
  updated_coldata %>%
  filter(final_cluster_label %in% c("Erythroid lineage",
                                       "White blood cells",
                                       "Endothelial cells",
                                       "Fibroblast")) %>%
  pull(Cell)
  
retest_cells =
  density_df %>%
  filter(!in_mask) %>%
  filter(!(Cell %in% cells_to_exclude_from_density_mapping)) %>%
  pull(Cell)

length(retest_cells)

gaussian_blur_3x3 = 
  matrix(c(2,3,2,
           3,5,3,
           2,3,2),
         nrow = 3) * (1/25)

# Iterate through Cells and find the top combination of spot and sector --------
density_list = list()
for (i in seq(1,4)){
  
  current_cells = 
    intersect(rownames(sector.hashes.list[[i]]), 
              rownames(spot.hashes.list[[i]])) %>%
    intersect(retest_cells) %>%
    as.character()
  
  print(length(current_cells))
  
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
                    mask_df){
             
             this_cell = as.character(x)
             this_slide_id = 
               coldata_to_test[which(coldata_to_test$Cell == this_cell),"max_slide_id"] %>%
               as.character()
             
             # Get the spots which are in the possible universe of spots
             this_slide_layout = 
               image_data[which(image_data$slide_id == this_slide_id),
                          "transformed_oligo_layout"][[1]][[1]]
             
             this_cluster = 
               coldata_to_test[which(coldata_to_test$Cell == this_cell) ,"cluster"]
             
             mask = 
               mask_df[which(mask_df$Cell == this_cell),"test"][[1]][[1]] %>%
               reshape2::melt(varnames = c("Row","Col")) %>%
               filter(value > 3)
             

             matching_spots = 
               this_slide_layout %>%
               filter(in_hull) %>%
               left_join(data.frame(value_spot =as.numeric(spot_mat[this_cell,]),
                                    hash = colnames(spot_mat)), 
                         by = "hash")
             
             mat = 
               matching_spots %>%
               dplyr::select(Row,Col,value_spot)
             
             mat = Matrix::sparseMatrix(i = mat$Row, 
                                        j = mat$Col, 
                                        x = mat$value_spot) %>%
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
                                 "Col"))%>%
               inner_join(mask,
                          by = c("Row","Col"))
             
             
             # Choose the top spot 
             top_spot = 
               matching_spots %>%
               arrange(desc(blurred_spot_value)) %>%
               mutate(Cell = this_cell) %>%
               head(n = 1)  
             
             return(top_spot)
           }, 
           updated_coldata,
           sector.hashes.list[[i]],
           spot.hashes.list[[i]],
           all_image_data,
           density_df)
  
  density_list[[i]] = 
    do.call(rbind,temp)
}


# Make a data frame from the result
density.refined.df.list = list()
for( i in seq(1,1)){
  density.refined.df.list[[i]] =
    data.frame(row.names = density_list[[i]]$Cell,
               Cell = density_list[[i]]$Cell,
               top_spot.mask = as.character(density_list[[i]]$Oligo),
               Row.mask = as.character(density_list[[i]]$Row),
               Col.mask = as.character(density_list[[i]]$Col),
               blurred_spot_value.mask = density_list[[i]]$blurred_spot_value,
               coords.x1.mask = density_list[[i]]$coords.x1,
               coords.x2.mask = density_list[[i]]$coords.x2)
}


density.refined.spots =
  do.call(rbind,density.refined.df.list)


new_density_positions = 
  inner_join(density.refined.spots,
          updated_coldata,
            by = "Cell")


new_density_positions$Col.mask = as.character(new_density_positions$Col.mask)
new_density_positions$Row.mask = as.character(new_density_positions$Row.mask)

new_density_positions = 
  new_density_positions %>%
  mutate(dist_spot = sqrt((as.numeric(as.character(Row.mask)) - as.numeric(as.character(Row)))^2 + 
                            (as.numeric(as.character(Col.mask)) - as.numeric(as.character(Col)))^2))



new_density_positions  = 
  new_density_positions %>% 
  filter(blurred_spot_value/blurred_spot_value.mask <= 2.5) %>%
  filter(dist_spot > 4)

# > dim(new_density_positions)
# [1] 1004   47


new_density_positions$top_spot = new_density_positions$top_spot.mask %>% as.character()
new_density_positions$Row = new_density_positions$Row.mask %>% as.character()
new_density_positions$Col = new_density_positions$Col.mask%>% as.character()
new_density_positions$coords.x1 = new_density_positions$coords.x1.mask
new_density_positions$coords.x2 = new_density_positions$coords.x2.mask

new_positions = new_density_positions[,colnames(updated_coldata)]             

updated_coldata =
  updated_coldata %>% 
  filter(!(Cell %in% new_density_positions$Cell)) %>%
  rbind(new_positions)

rownames(updated_coldata) = 
  updated_coldata$Cell



saveRDS(object = updated_coldata,
        file = "Submission_Data/E14_slides/RDS_intermediates/Notebook5.4_updated_coldata.RDS")

