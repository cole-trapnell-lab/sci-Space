# Use hash data to find the top mapped position for each nucleus.
# Briefly this is done by getting the spatial pattern of hash oligos
# Performing a gaussian convolution to take into account spatial background signal
# Map to the top combination of convolved spot / sector combination

# Load startup packages ---------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(Matrix)
  library(OpenImageR)
  library(reshape2)
  library(sp)
  library(vec2dtransf)
  library(monocle3)
  DelayedArray:::set_verbose_block_processing(TRUE)
  options(DelayedArray.block.size=1000e7)
})

space_directory = "/Volumes/GoogleDrive/My Drive/sciSpace/"
setwd(dir=space_directory)

# Read in CDS
spatial_cds = readRDS("Submission_Data/E14_slides/RDS_intermediates/Notebook5_E14_spatial_CDS.RDS")

# Read in hashTable -- contains spatial information
hash_table.list = 
  readRDS("Submission_Data/E14_slides/RDS_intermediates/Notebook0_hash_table.RDS")

# Read in the layout of the grid with respect to spatial positions derived from imaging
all_image_data = 
  readRDS(file = "Submission_Data/E14_slides/RDS_intermediates/Notebook0_1_images_and_transformations.RDS")


# Separate Sector and Spot Oligos ----------------------------------------
sector.hashes.list = 
  lapply(hash_table.list,
         function(x){
           # Make a matrix of sector oligos
           sector_oligo =
             x %>%
             ungroup() %>%
             dplyr::rename(Oligo = oligo) %>%
             filter(axis == 2) %>%
             dplyr::select(Cell, Oligo, count)%>%
             spread(key = Oligo,
                    value = count,
                    fill = 0) 
           
           sector_oligo_cell = sector_oligo$Cell
           sector_oligo = 
             sector_oligo[,2:ncol(sector_oligo)] %>%
             as.matrix()
           rownames(sector_oligo) = sector_oligo_cell
           sector_oligo
         })


# Make a matrix of spot oligos
spot.hashes.list = 
  lapply(hash_table.list,
         function(x){
           spot_oligo =
             x %>%
             ungroup() %>%
             dplyr::rename(Oligo = oligo) %>%  
             filter(axis == 3) %>%
             dplyr::select(Cell, Oligo, count)%>%
             spread(key = Oligo, 
                    value = count,
                    fill = 0) 
           
           spot_oligo_cell= spot_oligo$Cell
           
           spot_oligo = 
             spot_oligo[,2:ncol(spot_oligo)] %>%
             as.matrix()
           rownames(spot_oligo) = spot_oligo_cell
           spot_oligo
         })

# Choose cells have both spot and sector oligos along with inclusion in the CDS
cells_to_test =
  intersect(sapply(sector.hashes.list,
                   FUN = rownames) %>% 
              do.call(what = base::c), 
            sapply(spot.hashes.list,
                   FUN = rownames) %>% 
              do.call(what = base::c)) %>% 
  intersect(colData(spatial_cds)$Cell) 

# No cells dropped
coldata_spatial_cds = 
  colData(spatial_cds) %>%
  as.data.frame() %>%
  filter(Cell %in% cells_to_test)

rownames(coldata_spatial_cds) = 
  coldata_spatial_cds$Cell

# Gaussian blur matrix for the sector oligo
gaussian_blur_7x7 = 
  matrix(c(1,6,15,20,15,6,1,
           6,36,90,120,90,36,6,
           15,90,225,300,225,90,15,
           20,120,300,400,400,120,20,
           15,90,225,300,225,90,15,
           6,36,90,120,90,36,6,
           1,6,15,20,15,6,1),
         nrow = 7) * (1/4096)

# Modified gaussian blur matrix for the spot oligo
gaussian_blur_3x3 = 
  matrix(c(2,3,2,
           3,5,3,
           2,3,2),
         nrow = 3) * (1/25)


# Iterate through Cells and find the top combination of spot and sector --------
top_spots_list = list()
for (i in seq(1,4)){
  
  current_cells = 
    intersect(rownames(sector.hashes.list[[i]]), 
              rownames(spot.hashes.list[[i]])) %>%
    intersect(coldata_spatial_cds$Cell) %>%
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
                    image_data){
             this_cell = as.character(x)
             this_slide_id = 
               coldata_to_test[this_cell,"max_slide_id"] %>%
               as.character()
             
             # Get the spots which are in the possible universe of spots
             this_slide_layout = 
               image_data[which(rownames(image_data) == this_slide_id),
                          "transformed_oligo_layout"][[1]][[1]]
             
             matching_spots = 
               this_slide_layout %>%
               filter(in_hull) %>%
               left_join(data.frame(value_spot =as.numeric(spot_mat[this_cell,]),
                                    hash = colnames(spot_mat)), 
                         by = "hash") %>%
               left_join(data.frame(value_sector =as.numeric(sector_mat[this_cell,]),
                                    sector = colnames(sector_mat)),
                         by = "sector") %>%
               mutate(count_product = value_spot * value_sector)
             
             mat = 
               matching_spots %>%
               dplyr::select(Row,Col,value_spot)
             
             mat = Matrix::sparseMatrix(i = mat$Row, 
                                        j = mat$Col, 
                                        x = mat$value_spot) %>%
               Matrix::as.matrix()
             
             res =
               convolution(mat,
                           kernel = gaussian_blur_3x3,
                           mode = "same")
             
             matching_spots = 
               melt(res,varnames = c("Row","Col")) %>%
               dplyr::select(everything(),
                             blurred_spot_value = value) %>%
               right_join(matching_spots,
                          by = c("Row",
                                 "Col"))
             
             mat2 =
               matching_spots %>%
               dplyr::select(Row,Col,value_sector)
             
             mat2 = Matrix::sparseMatrix(i = mat2$Row, 
                                         j = mat2$Col, 
                                         x = mat2$value_sector) %>%
               Matrix::as.matrix()
             
             res2 =     
               convolution(mat2,
                           kernel = gaussian_blur_7x7,
                           mode = "same")
             
             matching_spots = 
               melt(res2,varnames = c("Row","Col")) %>%
               dplyr::select(everything(),
                             blurred_sector_value = value) %>%
               right_join(matching_spots,
                          by = c("Row",
                                 "Col"))
             
             # Choose the top spot 
             top_spot = 
               matching_spots %>%
               mutate(total_count_product = 
                        blurred_sector_value * blurred_spot_value * value_spot) %>%
               arrange(desc(total_count_product)) %>%
               mutate(Cell = this_cell) %>%
               head(n = 2)
             
             top_spot = 
               top_spot %>%
               mutate(ratio = 
                        max(top_spot$total_count_product)/min(top_spot$total_count_product)) %>%
               arrange(desc(total_count_product)) %>%
               head(n = 1)
             
             
             return(top_spot)
           }, 
           coldata_spatial_cds,
           sector.hashes.list[[i]],
           spot.hashes.list[[i]],
           all_image_data)
  
  top_spots_list[[i]] = 
    do.call(rbind,temp)
}


# Make a data frame from the result
df.list = list()
for( i in seq(1,4)){
  df.list[[i]] =
    data.frame(row.names = top_spots_list[[i]]$Cell,
               Cell = top_spots_list[[i]]$Cell,
               top_spot = as.character(top_spots_list[[i]]$Oligo),
               Row = as.character(top_spots_list[[i]]$Row),
               Col = as.character(top_spots_list[[i]]$Col),
               value_spot = top_spots_list[[i]]$value_spot,
               blurred_spot_value = top_spots_list[[i]]$blurred_spot_value,
               value_sector = top_spots_list[[i]]$value_sector,
               blurred_sector_value = top_spots_list[[i]]$blurred_sector_value,
               count_product = top_spots_list[[i]]$count_product)
}

top_spots =
  do.call(rbind,df.list)

top_spots =
  coldata_spatial_cds %>%
  left_join(top_spots,
            by = "Cell") 


#saveRDS(top_spots, file = "Submission_Data/E14_slides/RDS_intermediates/Notebook1_topspots.RDS")
#top_spots = readRDS("Submission_Data/E14_slides/RDS_intermediates/Notebook1_topspots.RDS")

# Add spot position and pixel coordinates to the cells -------------------------------
# Add oligo layout parameters into the data frame
transformed_oligo_layouts =
  lapply(X = all_image_data$slide_id,
         FUN = function(slide){
           temp = all_image_data[which(all_image_data$slide_id == slide),
                                 "transformed_oligo_layout"][[1]][[1]]
           temp$max_slide_id = slide
           temp %>%
             dplyr::select(max_slide_id,
                           Oligo,
                           Row,
                           Col,
                           coords.x1,
                           coords.x2,
                           SYBR)
         })


updated_coldata = 
  do.call(rbind,transformed_oligo_layouts) %>% 
  mutate(Row = as.character(Row),
         Col = as.character(Col)) %>%
  right_join(top_spots %>%
               mutate(Row = as.character(Row),
                      Col = as.character(Col)),
             by = c("Oligo" = "top_spot",
                    "max_slide_id",
                    "Row" = "Row",
                    "Col" = "Col"))



rownames(updated_coldata) =
  updated_coldata$Cell

saveRDS(updated_coldata,
        file = "Submission_Data/E14_slides/RDS_intermediates/Notebook5.1_updated_coldata.RDS")

save.image("Submission_Data/E14_slides/RDS_intermediates/Notebook5_initial_mapping.RData")
