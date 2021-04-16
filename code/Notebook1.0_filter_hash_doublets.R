# Notebook that removes cells that have too many RNA UMIs -- likely doublets --  
# or contain too few sector (5) or spot (10) hash oligo UMIs

# Load startup packages ---------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(Matrix)
  library(OpenImageR)
  library(reshape2)
  library(sp)
  library(vec2dtransf)
  library(monocle3)
  space_directory = "/Volumes/GoogleDrive/My Drive/sciSpace/"
  setwd(dir=space_directory)
  
  source("data_from_experiments/bin/chiSq_test_functions.R")
  DelayedArray:::set_verbose_block_processing(TRUE)
  options(DelayedArray.block.size=1000e7)
})


# Read in CDS
spatial_cds = readRDS("Submission_Data/E14_slides/RDS_intermediates/Notebook0_E14_spatial_CDS.RDS")

# Read in hashTable -- contains spatial information
hash_table.list = 
  readRDS("Submission_Data/E14_slides/RDS_intermediates/Notebook0_hash_table.RDS")

# Read in the layout of the grid with respect to spatial positions derived from imaging
all_image_data = 
  readRDS(file = "Submission_Data/E14_slides/RDS_intermediates/Notebook0_1_images_and_transformations.RDS")

# Rename slide_id in all_image_data to match column data in CDS
substr(all_image_data$slide_id, 
       start = 1, 
       stop = 1) = 
  stringr::str_to_lower(substr(x = all_image_data$slide_id, 
                               start = 1, 
                               stop = 1))

rownames(all_image_data) = all_image_data$slide_id

# Remove cells that are at the upper tail of the UMI distribution --------------
cells_to_keep = 
  colData(spatial_cds) %>%
  as.data.frame() %>%
  group_by(max_slide_id) %>%
  filter(n.umi  < 2*sqrt(var(n.umi)) + mean(n.umi)) %>%
  pull(Cell) %>%
  as.character() 

spatial_cds = 
  spatial_cds[,cells_to_keep]

# Filter on cells that have at least Hashes per cell ---------------------------------------------------------
passing_cells = list()

for(i in seq(1,length(hash_table.list))){
  passing_axis_3 =
    hash_table.list[[i]] %>%
    filter(Cell %in% colData(spatial_cds)$Cell) %>%
    dplyr::select(axis,Cell, oligo, count) %>%
    filter(axis == 3) %>%
    group_by(Cell,axis) %>%
    top_n(n = 1) %>%  
    filter(count >= 10) %>%
    pull(Cell)
  
  passing_axis_2 = 
    hash_table.list[[i]] %>%
    filter(Cell %in% colData(spatial_cds)$Cell) %>%
    dplyr::select(axis,Cell, oligo, count) %>%
    filter(axis == 2) %>%
    group_by(Cell,axis) %>%
    top_n(n = 1) %>%  
    filter(count >= 5) %>%
    pull(Cell)
  
  passing_cells[[i]] = intersect(passing_axis_2, passing_axis_3)
}

colData(spatial_cds)$sample %>% length()
# [1] 132382

passing_cells = 
  do.call(c,passing_cells)

# > length(passing_cells)
# [1] 132382

spatial_cds = spatial_cds[,passing_cells]
colData(spatial_cds)$sample %>% length()


# Write out RDS that has these cells removed -----------------------

saveRDS(object = spatial_cds,
        file = "Submission_Data/E14_slides/RDS_intermediates/Notebook1_E14_spatial_CDS.RDS")
# Write out mtx file that contains the new slide information

writeMM(obj = counts(spatial_cds),
        file = "Submission_Data/E14_slides/RDS_intermediates/Notebook1_E14_count_matrix.mtx")

