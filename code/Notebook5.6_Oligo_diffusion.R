# Sample mapped cells to look at the signal from adjacent areas
# Scripts used to generate fig S8C-E and associated numbers

# Load startup packages ---------------------------------------------------
suppressPackageStartupMessages({ 
  library(tidyverse)
  library(ggplot2)
  library(tidyr)
  library(viridis)
  library(ggridges)
  library(RColorBrewer)
  library(ggrepel)
  library(monocle3)
  library(garnett)
  library(pheatmap)
  
  space_directory = "/Volumes/GoogleDrive/My Drive/sciSpace/"
  setwd(dir=space_directory)

  # Set a seed to make umap and other non-deterministic steps consistent
  set.seed(seed = 42)
  
})

# Load previously saved files
load("Submission_Data/E14_slides/RDS_intermediates/Notebook5_initial_mapping.RData")

spatial_cds = 
  readRDS(file = "Submission_Data/E14_slides/RDS_intermediates/Notebook5_E14_spatial_CDS.RDS")

updated_coldata = 
  readRDS(file = "Submission_Data/E14_slides/RDS_intermediates/Notebook5.1_updated_coldata.RDS")

updated_coldata = 
  colData(spatial_cds) %>%
  as.data.frame() %>%
  left_join(updated_coldata,
            by = "Cell")

# Function for fitting a gaussian distribution
fitG =
  function(x,y,mu,sig,scale){
    
    f = function(p){
      d = p[3]*dnorm(x,mean=p[1],sd=p[2])
      sum((d-y)^2)
    }
    
    optim(c(mu,sig,scale),f)
  }

# Slide 1 ------------------------------------------------------------

# Sample 1000 cells
sampled_cells_1D = 
  updated_coldata %>% 
  filter(max_slide_id == "slide_1D") %>%
  sample_n(size = 1000) %>%
  pull(Cell)


# Get the oligo mapping for these 1000 cells
temp_1D = 
    lapply(sampled_cells_1D, 
           function(x,
                    coldata_to_test,
                    sector_mat,
                    spot_mat,
                    image_data){
             
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
                                        x = mat$value_spot,
                                        dims = c(84,84)) %>%
               Matrix::as.matrix()
             
             top_spot_col = coldata_to_test[which(coldata_to_test$Cell == x),"Col"] %>% 
               as.character() %>% 
               as.numeric()
             top_spot_row = coldata_to_test[which(coldata_to_test$Cell == x),"Row"] %>% 
               as.character() %>% 
               as.numeric()
             
             sub_mat = mat[(top_spot_row-3):(top_spot_row+3),
                           (top_spot_col-3):(top_spot_col+3)]
             return( sub_mat)
           }, 
           updated_coldata,
           sector.hashes.list[[1]],
           spot.hashes.list[[1]],
           all_image_data)


# Make a plot for a reviewer figure showing signal from adjacent positions
data.frame(`0.95` = sapply(temp_1D,FUN = function(x){ sum(x/max(x) > .95) > 1}),
           `0.90` = sapply(temp_1D,FUN = function(x){ sum(x/max(x) > .90) > 1}),
           `0.75` = sapply(temp_1D,FUN = function(x){ sum(x/max(x) > .75) > 1}),
           `0.50` = sapply(temp_1D,FUN = function(x){ sum(x/max(x) > .50) > 1})) %>%
  gather(key = "Percent Maximum",value = "Present") %>%
  ggplot(aes(x = `Percent Maximum`)) +
  geom_bar(aes(y = 4*(..count..)/sum(..count..),
               fill = Present),
           color = "black",
           size = .25) +
  scale_y_continuous(labels=scales::percent) +
  scale_x_discrete(labels = c("X0.50" = "50% Max",
                              "X0.75" = "75% Max",
                              "X0.90" = "90% Max",
                              "X0.95" = "95% Max")) +
  ylab("Percentage of Cells") +
  xlab("Percentage of Maximum Count") +
  scale_fill_brewer(palette = "Set1") +
  coord_flip() +
  monocle3:::monocle_theme_opts() +
  ggsave("Figures/Figure_Components/Supplement_Grid_Alignment/barplot_percent_max_spot.pdf",
         height = 3,
         width = 4.5)

# Aggregate and average all the positions surrouding the top position
# Make this into a matrix
avg_mat_1D = Reduce("+", temp_1D) / length(temp_1D) 

avg_df_1D = 
  avg_mat_1D %>%
  as.data.frame() 

colnames(avg_df_1D) = seq(1,7) %>% as.character()
row.names(avg_df_1D) = seq(1,7) %>% as.character()

# Supplemental Figure 8C
avg_df_1D %>%
  rownames_to_column(var = "Row") %>%
  gather(key = "Col",
         value = "counts",
         -Row) %>%
  ggplot(aes(x = Row,
         y = Col,
         fill = counts)) +
  geom_tile(color = "black",
            size = 0.25) +
  scale_fill_viridis_c(option = "D",name = "Avg. Count") +
  theme_void() +
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Supplement_Grid_Alignment/heatmap_counts_1D.pdf",
         height = 2,
         width = 2)

# Supplemental Figure 8C legend
avg_df_1D %>%
  rownames_to_column(var = "Row") %>%
  gather(key = "Col",
         value = "counts",
         -Row) %>%
  ggplot(aes(x = Row,
             y = Col,
             fill = counts)) +
  geom_tile(color = "black",
            size = 0.25) +
  scale_fill_viridis_c(option = "D",name = "Avg. Count") +
  theme_void() +
  ggsave("Figures/Figure_Components/Supplement_Grid_Alignment/heatmap_counts_1D_legend.pdf",
         height = 2,
         width = 2)

# Supplemental Figure 8C row marginal
avg_df_1D %>%
  rownames_to_column(var = "Row") %>%
  gather(key = "Col",
         value = "counts",
         -Row) %>%
  group_by(Row) %>%
  summarise(marginal_count = sum(counts)) %>%
  ggplot(aes(x = Row,
             y = marginal_count)) +
  geom_bar(color = "black",
            size = 0.25,
           stat = "identity",
           fill = "grey80") +
  coord_flip() +
  theme_void() +
  ggsave("Figures/Figure_Components/Supplement_Grid_Alignment/heatmap_row_marginal_1D.pdf",
         height =2,
         width = 0.5)

# Supplemental Figure 8C Col marginal
avg_df_1D %>%
  rownames_to_column(var = "Row") %>%
  gather(key = "Col",
         value = "counts",
         -Row) %>%
  group_by(Col) %>%
  summarise(marginal_count = sum(counts)) %>%
  ggplot(aes(x = Col,
             y = marginal_count)) +
  geom_bar(color = "black",
           size = 0.25,
           stat = "identity",
           fill = "grey80") +
  theme_void() +
  ggsave("Figures/Figure_Components/Supplement_Grid_Alignment/heatmap_Col_marginal_1D.pdf",
         height =0.5,
         width = 2)

# Calculate the gaussian fit for the marginal distributions ---------------

fit_slide_1D_rows = 
  fitG(seq(1,7),rowSums(avg_mat_1D),mu = 4, sig = 0.5, scale = 1500)


ggplot() +
  geom_bar(data = avg_df_1D %>%
             rownames_to_column(var = "Row") %>%
             gather(key = "Col",
                    value = "counts",
                    -Row) %>%
             group_by(Row) %>%
             summarise(marginal_count = sum(counts)),
           aes(x = Row,
               y = marginal_count),
           color = "black",
           size = 0.25,
           stat = "identity",
           fill = "grey80") +
  geom_line(data =data.frame(x = seq(1,7,by = 0.1),
                       y = fit_slide_1D_rows$par[3]*dnorm(x = seq(1,7,by = 0.1),fit_slide_1D_rows$par[1],fit_slide_1D_rows$par[2])),
            aes(x = x,
                y = y)) +
  monocle3:::monocle_theme_opts() +
  ylab("Marginal Counts") +
  ggsave("Figures/Figure_Components/Supplement_Grid_Alignment/row_cumulative_gaussian_fit.pdf")
  

fit_slide_1D_cols = fitG(seq(1,7),colSums(avg_mat_1D),mu = 4, sig = 0.5, scale = 1500)


ggplot() +
  geom_bar(data = avg_df_1D %>%
             rownames_to_column(var = "Row") %>%
             gather(key = "Col",
                    value = "counts",
                    -Row) %>%
             group_by(Col) %>%
             summarise(marginal_count = sum(counts)),
           aes(x = Col,
               y = marginal_count),
           color = "black",
           size = 0.25,
           stat = "identity",
           fill = "grey80") +
  geom_line(data =data.frame(x = seq(1,7,by = 0.1),
                             y = fit_slide_1D_cols$par[3]*dnorm(x = seq(1,7,by = 0.1),fit_slide_1D_cols$par[1],fit_slide_1D_cols$par[2])),
            aes(x = x,
                y = y)) +
  monocle3:::monocle_theme_opts() +
  ylab("Marginal Counts") +
  ggsave("Figures/Figure_Components/Supplement_Grid_Alignment/col_cumulative_gaussian_fit.pdf")


# Slide 4  ------------------------------------------------------------

# Sample 1000 cells
sampled_cells_1G = 
  updated_coldata %>% 
  filter(max_slide_id == "slide_1G") %>%
  sample_n(size = 1000) %>%
  pull(Cell)

# Get the oligo mapping for these 1000 cells
temp_1G = 
  lapply(sampled_cells_1G, 
         function(x,
                  coldata_to_test,
                  sector_mat,
                  spot_mat,
                  image_data){
           
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
                                      x = mat$value_spot,
                                      dims = c(84,84)) %>%
             Matrix::as.matrix()
           
           top_spot_col = coldata_to_test[which(coldata_to_test$Cell == x),"Col"] %>% 
             as.character() %>% 
             as.numeric()
           top_spot_row = coldata_to_test[which(coldata_to_test$Cell == x),"Row"] %>% 
             as.character() %>% 
             as.numeric()
           
           sub_mat = mat[(top_spot_row-3):(top_spot_row+3),
                         (top_spot_col-3):(top_spot_col+3)]
           return( sub_mat)
         }, 
         updated_coldata,
         sector.hashes.list[[1]],
         spot.hashes.list[[1]],
         all_image_data)

# Take the average counts surrounding the top position
avg_mat_1G = Reduce("+", temp_1G) / length(temp_1G) 

avg_df_1G = 
  avg_mat_1G %>%
  as.data.frame() 

colnames(avg_df_1G) = seq(1,7) %>% as.character()
row.names(avg_df_1G) = seq(1,7) %>% as.character()

# Supplemental Figure 8D 
avg_df_1G %>%
  rownames_to_column(var = "Row") %>%
  gather(key = "Col",
         value = "counts",
         -Row) %>%
  ggplot(aes(x = Row,
             y = Col,
             fill = counts)) +
  geom_tile(color = "black",
            size = 0.25) +
  scale_fill_viridis_c(option = "D",name = "Total Count") +
  theme_void() +
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Supplement_Grid_Alignment/heatmap_counts_1G.pdf",
         height = 2,
         width = 2)

# Supplemental Figure 8D legend
avg_df_1G %>%
  rownames_to_column(var = "Row") %>%
  gather(key = "Col",
         value = "counts",
         -Row) %>%
  ggplot(aes(x = Row,
             y = Col,
             fill = counts)) +
  geom_tile(color = "black",
            size = 0.25) +
  scale_fill_viridis_c(option = "D",name = "Avg. Count") +
  theme_void() +
  ggsave("Figures/Figure_Components/Supplement_Grid_Alignment/heatmap_counts_1G_legend.pdf",
         height = 2,
         width = 2)

# Supplemental Figure 8D Row Marginal
avg_df_1G %>%
  rownames_to_column(var = "Row") %>%
  gather(key = "Col",
         value = "counts",
         -Row) %>%
  group_by(Row) %>%
  summarise(marginal_count = sum(counts)) %>%
  ggplot(aes(x = Row,
             y = marginal_count)) +
  geom_bar(color = "black",
           size = 0.25,
           stat = "identity",
           fill = "grey80") +
  coord_flip() +
  theme_void() +
  ggsave("Figures/Figure_Components/Supplement_Grid_Alignment/heatmap_row_marginal_1G.pdf",
         height =2,
         width = 0.5)

# Supplemental Figure 8D Column Marginal
avg_df_1G %>%
  rownames_to_column(var = "Row") %>%
  gather(key = "Col",
         value = "counts",
         -Row) %>%
  group_by(Col) %>%
  summarise(marginal_count = sum(counts)) %>%
  ggplot(aes(x = Col,
             y = marginal_count)) +
  geom_bar(color = "black",
           size = 0.25,
           stat = "identity",
           fill = "grey80") +
  theme_void() +
  ggsave("Figures/Figure_Components/Supplement_Grid_Alignment/heatmap_Col_marginal_1G.pdf",
         height =0.5,
         width = 2)

# Gaussian fits for marginal distribution
fit_slide_1G_rows = fitG(seq(1,7),rowSums(avg_mat_1G),mu = 4, sig = 0.5, scale = 600)
fit_slide_1G_cols = fitG(seq(1,7),colSums(avg_mat_1G),mu = 4, sig = 0.5, scale = 600)

# Slide 13  -------------------------------------------------------------

# Sample cells
sampled_cells_4D = 
  updated_coldata %>% 
  filter(max_slide_id == "slide_4D") %>%
  filter(sample == "2") %>%
  sample_n(size = 1000) %>%
  pull(Cell)

# Map cells
temp_4D = 
  lapply(sampled_cells_4D, 
         function(x,
                  coldata_to_test,
                  sector_mat,
                  spot_mat,
                  image_data){
           
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
                                      x = mat$value_spot,
                                      dims = c(84,84)) %>%
             Matrix::as.matrix()
           
           top_spot_col = coldata_to_test[which(coldata_to_test$Cell == x),"Col"] %>% 
             as.character() %>% 
             as.numeric()
           top_spot_row = coldata_to_test[which(coldata_to_test$Cell == x),"Row"] %>% 
             as.character() %>% 
             as.numeric()
           
           sub_mat = mat[(top_spot_row-3):(top_spot_row+3),
                         (top_spot_col-3):(top_spot_col+3)]
           return( sub_mat)
         }, 
         updated_coldata,
         sector.hashes.list[[2]],
         spot.hashes.list[[2]],
         all_image_data)

# Average positions surrounding top position
avg_mat_4D = Reduce("+", temp_4D) / length(temp_4D) 

avg_df_4D = 
  avg_mat_4D %>%
  as.data.frame() 

colnames(avg_df_4D) = seq(1,7) %>% as.character()
row.names(avg_df_4D) = seq(1,7) %>% as.character()

# Supplemental Figure 8E
avg_df_4D %>%
  rownames_to_column(var = "Row") %>%
  gather(key = "Col",
         value = "counts",
         -Row) %>%
  ggplot(aes(x = Row,
             y = Col,
             fill = counts)) +
  geom_tile(color = "black",
            size = 0.25) +
  scale_fill_viridis_c(option = "D",name = "Total Count") +
  theme_void() +
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Supplement_Grid_Alignment/heatmap_counts_4D.pdf",
         height = 2,
         width = 2)


# Supplemental Figure 8E legend
avg_df_4D %>%
  rownames_to_column(var = "Row") %>%
  gather(key = "Col",
         value = "counts",
         -Row) %>%
  ggplot(aes(x = Row,
             y = Col,
             fill = counts)) +
  geom_tile(color = "black",
            size = 0.25) +
  scale_fill_viridis_c(option = "D",name = "Avg. Count") +
  theme_void() +
  ggsave("Figures/Figure_Components/Supplement_Grid_Alignment/heatmap_counts_4D_legend.pdf",
         height = 2,
         width = 2)

# Supplemental Figure 8E Row marginal
avg_df_4D %>%
  rownames_to_column(var = "Row") %>%
  gather(key = "Col",
         value = "counts",
         -Row) %>%
  group_by(Row) %>%
  summarise(marginal_count = sum(counts)) %>%
  ggplot(aes(x = Row,
             y = marginal_count)) +
  geom_bar(color = "black",
           size = 0.25,
           stat = "identity",
           fill = "grey80") +
  coord_flip() +
  theme_void() +
  ggsave("Figures/Figure_Components/Supplement_Grid_Alignment/heatmap_row_marginal_4D.pdf",
         height =2,
         width = 0.5)

# Supplemental Figure 8E Col marginal
avg_df_4D %>%
  rownames_to_column(var = "Row") %>%
  gather(key = "Col",
         value = "counts",
         -Row) %>%
  group_by(Col) %>%
  summarise(marginal_count = sum(counts)) %>%
  ggplot(aes(x = Col,
             y = marginal_count)) +
  geom_bar(color = "black",
           size = 0.25,
           stat = "identity",
           fill = "grey80") +
  theme_void() +
  ggsave("Figures/Figure_Components/Supplement_Grid_Alignment/heatmap_Col_marginal_4D.pdf",
         height =0.5,
         width = 2)

# Gaussian fits for marginal distribution
fit_slide_4D_rows = fitG(seq(1,7),rowSums(avg_mat_4D),mu = 4, sig = 0.5, scale = 400)
fit_slide_4D_cols = fitG(seq(1,7),colSums(avg_mat_4D),mu = 4, sig = 0.5, scale = 400)


# Slide 14 4E -------------------------------------------------------------

# Sample 1000 cells
sampled_cells_4E = 
  updated_coldata %>% 
  filter(max_slide_id == "slide_4E") %>%
  filter(sample == "2") %>%
  sample_n(size = 1000) %>%
  pull(Cell)

# Get the oligo mapping for these 1000 cells
temp_4E = 
  lapply(sampled_cells_4E, 
         function(x,
                  coldata_to_test,
                  sector_mat,
                  spot_mat,
                  image_data){
           
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
                                      x = mat$value_spot,
                                      dims = c(84,84)) %>%
             Matrix::as.matrix()
           
           top_spot_col = coldata_to_test[which(coldata_to_test$Cell == x),"Col"] %>% 
             as.character() %>% 
             as.numeric()
           top_spot_row = coldata_to_test[which(coldata_to_test$Cell == x),"Row"] %>% 
             as.character() %>% 
             as.numeric()
           
           sub_mat = mat[(top_spot_row-3):(top_spot_row+3),
                         (top_spot_col-3):(top_spot_col+3)]
           return( sub_mat)
         }, 
         updated_coldata,
         sector.hashes.list[[2]],
         spot.hashes.list[[2]],
         all_image_data)

# Average positions surrounding top position
avg_mat_4E = Reduce("+", temp_4E) / length(temp_4E) 

avg_df_4E = 
  avg_mat_4E %>%
  as.data.frame() 

colnames(avg_df_4E) = seq(1,7) %>% as.character()
row.names(avg_df_4E) = seq(1,7) %>% as.character()

# Supplemental Figure 8F heatmap 
avg_df_4E %>%
  rownames_to_column(var = "Row") %>%
  gather(key = "Col",
         value = "counts",
         -Row) %>%
  ggplot(aes(x = Row,
             y = Col,
             fill = counts)) +
  geom_tile(color = "black",
            size = 0.25) +
  scale_fill_viridis_c(option = "D",name = "Total Count") +
  theme_void() +
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Supplement_Grid_Alignment/heatmap_counts_4E.pdf",
         height = 2,
         width = 2)


# Supplemental Figure 8F heatmap legend
avg_df_4E %>%
  rownames_to_column(var = "Row") %>%
  gather(key = "Col",
         value = "counts",
         -Row) %>%
  ggplot(aes(x = Row,
             y = Col,
             fill = counts)) +
  geom_tile(color = "black",
            size = 0.25) +
  scale_fill_viridis_c(option = "D",name = "Avg. Count") +
  theme_void() +
  ggsave("Figures/Figure_Components/Supplement_Grid_Alignment/heatmap_counts_4E_legend.pdf",
         height = 2,
         width = 2)

# Supplemental Figure 8F Row marginal
avg_df_4E %>%
  rownames_to_column(var = "Row") %>%
  gather(key = "Col",
         value = "counts",
         -Row) %>%
  group_by(Row) %>%
  summarise(marginal_count = sum(counts)) %>%
  ggplot(aes(x = Row,
             y = marginal_count)) +
  geom_bar(color = "black",
           size = 0.25,
           stat = "identity",
           fill = "grey80") +
  coord_flip() +
  theme_void() +
  ggsave("Figures/Figure_Components/Supplement_Grid_Alignment/heatmap_row_marginal_4E.pdf",
         height =2,
         width = 0.5)

# Supplemental Figure 8F Col marginal
avg_df_4E %>%
  rownames_to_column(var = "Row") %>%
  gather(key = "Col",
         value = "counts",
         -Row) %>%
  group_by(Col) %>%
  summarise(marginal_count = sum(counts)) %>%
  ggplot(aes(x = Col,
             y = marginal_count)) +
  geom_bar(color = "black",
           size = 0.25,
           stat = "identity",
           fill = "grey80") +
  theme_void() +
  ggsave("Figures/Figure_Components/Supplement_Grid_Alignment/heatmap_Col_marginal_4E.pdf",
         height =0.5,
         width = 2)

# Gaussian fits for marginal distribution
fit_slide_4E_rows = fitG(seq(1,7),rowSums(avg_mat_4E),mu = 4, sig = 0.5, scale = 400)
fit_slide_4E_cols = fitG(seq(1,7),colSums(avg_mat_4E),mu = 4, sig = 0.5, scale = 400)

mean(c(fit_slide_1D_cols$par[2],fit_slide_1D_rows$par[2],
     fit_slide_1G_cols$par[2],fit_slide_1G_rows$par[2],
     fit_slide_4D_cols$par[2],fit_slide_4D_rows$par[2],
     fit_slide_4E_cols$par[2],fit_slide_4E_rows$par[2]))



