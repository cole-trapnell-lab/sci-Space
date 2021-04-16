# Calculate how much variance is explained based on the different labels
# that were either measured (spatial position) or assigned (celltype)

# Run as parrallel jobs on a computing cluster. Each model was run 50 times.
# Results collated and plotted at the bottom of this script

library(tidyverse)
library(dplyr)
library(purrr)
library(monocle3)

spatial_cds_to_test = readRDS("/net/trapnell/vol1/home/sanjays/projects/Space/Submission/E14_slides/RDS_intermediates/Notebook8.2_cds_to_test.RDS")

run.id = paste(round(runif(n = 4)* 100 ),collapse = "") 

# Calculate the norm of a vector
norm_vec  =
  function(x) {
    sqrt(sum(x^2))
  }

# Compute the angular distance between two cells
angular_distance =
  function(cell_1, cell_2){
    if(identical(cell_1,cell_2)){
      0
    }
    else{
      (2/pi * acos(cell_1 %*% cell_2 /(norm_vec(cell_1) * norm_vec(cell_2))))^2
    }
  }

# Rescale so that the magnitude of the vector is unit
unit_vector =
  function(vector){
    vector / norm_vec(vector)
  }

# Take a normalized expression count matrix
# Return the unit mean vector for that matrix
unit_mean =
  function(count_matrix){
    Matrix::rowSums(count_matrix) %>%
      unit_vector()
  }



# Takes a cds as input
# Returns a sparse matrix of the same dimensions as output with reads
# sampled from the average cds
sampled_cds = function(cds){
  umis_per_cell =
    Matrix::colSums(counts(cds))
  
  probability_dist =
    Matrix::rowSums(counts(cds))/sum(Matrix::rowSums(counts(cds)))
  
  sampled_cells =
    lapply(X = umis_per_cell,
           FUN = function(n_size){
             rmultinom(n = 1,
                       size = n_size,
                       prob = probability_dist) %>%
               as(Class = "sparseMatrix")
           })
  
  sampled_cells =
    do.call(cbind,
            sampled_cells)
  
  colnames(sampled_cells) = colnames(cds)
  counts(cds) = sampled_cells
  cds
}


compute_conditional_variance =
  function(pdata_df,cds_to_test){
    
    normalalized_cells =
      spatial_cds_to_test[,pdata_df$Cell] %>%
      monocle3:::normalize_expr_data(norm_method = "log",
                                     pseudo_count = 1)
    
    angular_mean = normalalized_cells %>% unit_mean()
    
    ang_dist =
      sapply(seq(1,ncol(normalalized_cells)),function(X,mat,ang_mean){
        
        this.cell = mat[,X] %>% unit_vector()
        
        angular_distance(this.cell, ang_mean)},
        normalalized_cells, angular_mean)
    
    sampled_cells =
      sampled_cds(spatial_cds_to_test[,pdata_df$Cell]) %>%
      monocle3:::normalize_expr_data(norm_method = "log",
                                     pseudo_count = 1)
    
    ang_dist_sparsity =
      sapply(seq(1,ncol(sampled_cells)),
             function(X,mat,ang_mean){
               
               this.cell = mat[,X] %>% unit_vector()
               
               angular_distance(this.cell,ang_mean)},
             sampled_cells,
             angular_mean)
    
    tibble(angular_distance = ang_dist,
           angular_distance_sparsity = ang_dist_sparsity,
           number_of_cells = dim(normalalized_cells)[2])
  }



# Calculate the within group variance for a model that only knows celltype --------
# Get the cells associated with each group and store the cds associated with that group
celltype_model =
  spatial_cds_to_test %>%
  colData() %>%
  as.data.frame() %>%
  group_by(celltype_bin) %>%
  nest() %>% 
  mutate( res = purrr::map(.x = data,.f = compute_conditional_variance,spatial_cds_to_test)) %>%
  dplyr::select(-data) %>%
  unnest(cols = c("res")) %>%
  ungroup() %>%
  mutate(mean_ang_dist = mean(angular_distance),
         mean_ang_dist_sp = mean(angular_distance_sparsity),
         calculation = "celltype",
         permuted = F,
         run.id = run.id)

write.table(celltype_model,
            file = paste("/net/trapnell/vol1/home/sanjays/projects/Space/Submission/E14_slides/RDS_intermediates/variance_decomposition/",
                         "celltype.model_",
                         run.id,
                         ".tsv",
                         sep = ""),
            col.names = T,
            row.names = F,
            sep = "\t",
            quote = F)

# Celltype model with permuted labels
celltype_model_permuted = 
  spatial_cds_to_test %>%
  colData() %>%
  as.data.frame() %>%
  mutate(shuffled.celltype = sample(celltype_bin,replace = F)) %>%
  group_by(shuffled.celltype) %>%
  nest() %>%
  mutate(res = purrr::map( .x = data, .f = compute_conditional_variance ,spatial_cds_to_test)) %>%
  dplyr::select(-data) %>%
  unnest(cols = c("res")) %>%
  ungroup() %>%
  mutate(mean_ang_dist = mean(angular_distance),
         mean_ang_dist_sp = mean(angular_distance_sparsity),
         calculation = "celltype",
         permuted = T,
         run.id = run.id)

write.table(celltype_model_permuted,
            file = paste("/net/trapnell/vol1/home/sanjays/projects/Space/Submission/E14_slides/RDS_intermediates/variance_decomposition/",
                         "celltype.model.permuted_",
                         run.id,
                         ".tsv",
                         sep = ""),
            col.names = T,
            row.names = F,
            sep = "\t",
            quote = F)

# Calculate the within group variance for a model that only knows cluster --------
# Get the cells associated with each group and store the cds associated with that group
cluster_model =
  spatial_cds_to_test %>%
  colData() %>%
  as.data.frame() %>%
  group_by(cluster) %>%
  nest() %>% 
  mutate( res = purrr::map(.x = data,.f = compute_conditional_variance,spatial_cds_to_test)) %>%
  dplyr::select(-data) %>%
  unnest(cols = c("res")) %>%
  ungroup() %>%
  mutate(mean_ang_dist = mean(angular_distance),
         mean_ang_dist_sp = mean(angular_distance_sparsity),
         calculation = "cluster",
         permuted = F,
         run.id = run.id)

write.table(cluster_model,
            file = paste("/net/trapnell/vol1/home/sanjays/projects/Space/Submission/E14_slides/RDS_intermediates/variance_decomposition/",
                         "cluster.model_",
                         run.id,
                         ".tsv",
                         sep = ""),
            col.names = T,
            row.names = F,
            sep = "\t",
            quote = F)

# Celltype model with permuted labels
cluster_model_permuted = 
  spatial_cds_to_test %>%
  colData() %>%
  as.data.frame() %>%
  mutate(shuffled.cluster = sample(cluster,replace = F)) %>%
  group_by(shuffled.cluster) %>%
  nest() %>%
  mutate(res = purrr::map( .x = data, .f = compute_conditional_variance ,spatial_cds_to_test)) %>%
  dplyr::select(-data) %>%
  unnest(cols = c("res")) %>%
  ungroup() %>%
  mutate(mean_ang_dist = mean(angular_distance),
         mean_ang_dist_sp = mean(angular_distance_sparsity),
         calculation = "cluster",
         permuted = T,
         run.id = run.id)

write.table(cluster_model_permuted,
            file = paste("/net/trapnell/vol1/home/sanjays/projects/Space/Submission/E14_slides/RDS_intermediates/variance_decomposition/",
                         "cluster.model.permuted_",
                         run.id,
                         ".tsv",
                         sep = ""),
            col.names = T,
            row.names = F,
            sep = "\t",
            quote = F)

# Calculate the within group variance for a model that only knows spatial position --------
# Get the cells associated with each spatial group and calculate angular_distance to that mean
spatial_model =
  spatial_cds_to_test %>%
  colData() %>%
  as.data.frame() %>%
  group_by(spatial_bin) %>%
  nest() %>%
  mutate(res = purrr::map(.x = data, .f = compute_conditional_variance, spatial_cds_to_test)) %>%
  dplyr::select(-data) %>%
  unnest(cols = c("res")) %>%
  ungroup() %>%
  mutate(mean_ang_dist = mean(angular_distance),
         mean_ang_dist_sp = mean(angular_distance_sparsity),
         calculation = "spatial",
         permuted = F,
         run.id = run.id)


write.table(spatial_model,
            file = paste("/net/trapnell/vol1/home/sanjays/projects/Space/Submission/E14_slides/RDS_intermediates/variance_decomposition/",
                         "spatial.model_",
                         run.id,
                         ".tsv",
                         sep = ""),
            col.names = T,
            row.names = F,
            sep = "\t",
            quote = F)

# Spatial model with permuted labels
spatial_model_permuted =
  spatial_cds_to_test %>%
  colData() %>%
  as.data.frame() %>%
  mutate(shuffled.spatial.bin = sample(spatial_bin,replace = F)) %>%
  group_by(shuffled.spatial.bin) %>%
  nest() %>%
  mutate(res = purrr::map(.x = data, .f = compute_conditional_variance, spatial_cds_to_test)) %>%
  dplyr::select(-data) %>%
  unnest(cols = c("res")) %>%
  ungroup() %>%
  mutate(mean_ang_dist = mean(angular_distance),
         mean_ang_dist_sp = mean(angular_distance_sparsity),
         calculation = "spatial",
         permuted = T,
         run.id = run.id)

write.table(spatial_model_permuted,
            file = paste("/net/trapnell/vol1/home/sanjays/projects/Space/Submission/E14_slides/RDS_intermediates/variance_decomposition/",
                         "spatial.model.permuted_",
                         run.id,
                         ".tsv",
                         sep = ""),
            col.names = T,
            row.names = F,
            sep = "\t",
            quote = F)


ct_spatial_model =
  # Get the cells associated with each spatial group and calculate angular_distance to that mean
  spatial_cds_to_test %>%
  colData() %>%
  as.data.frame() %>%
  group_by(celltype_spatial_bin) %>%
  nest() %>%
  mutate(res = purrr::map(.x = data, .f = compute_conditional_variance, spatial_cds_to_test)) %>%
  dplyr::select(-data) %>%
  unnest(cols = c("res")) %>%
  ungroup() %>%
  mutate(mean_ang_dist = mean(angular_distance),
         mean_ang_dist_sp = mean(angular_distance_sparsity),
         calculation = "spatial_celltype",
         permuted = F,
         run.id = run.id)


write.table(ct_spatial_model,
            file = paste("/net/trapnell/vol1/home/sanjays/projects/Space/Submission/E14_slides/RDS_intermediates/variance_decomposition/",
                         "celltype.spatial.model_",
                         run.id,
                         ".tsv",
                         sep = ""),
            col.names = T,
            row.names = F,
            sep = "\t",
            quote = F)

# Celltype/spatial model with permuted labels
ct_spatial_model_permuted =
  # Get the cells associated with each spatial group and calculate angular_distance to that mean
  spatial_cds_to_test %>%
  colData() %>%
  as.data.frame() %>%
  mutate(shuffled.ct.spatial.bin = sample(celltype_spatial_bin,replace = F)) %>%
  group_by(shuffled.ct.spatial.bin) %>%
  nest() %>%
  mutate(res = purrr::map(.x = data, .f = compute_conditional_variance, spatial_cds_to_test)) %>%
  dplyr::select(-data) %>%
  unnest(cols = c("res")) %>%
  ungroup() %>%
  mutate(mean_ang_dist = mean(angular_distance),
         mean_ang_dist_sp = mean(angular_distance_sparsity),
         calculation = "spatial_celltype",
         permuted = T,
         run.id = run.id)


write.table(ct_spatial_model_permuted,
            file = paste("/net/trapnell/vol1/home/sanjays/projects/Space/Submission/E14_slides/RDS_intermediates/variance_decomposition/",
                         "celltype.spatial.model.permuted_",
                         run.id,
                         ".tsv",
                         sep = ""),
            col.names = T,
            row.names = F,
            sep = "\t",
            quote = F)


# Cluster Spatial Model ---------------------------------------------------
cluster_spatial_model =
  # Get the cells associated with each spatial group and calculate angular_distance to that mean
  spatial_cds_to_test %>%
  colData() %>%
  as.data.frame() %>%
  group_by(cluster_spatial_bin) %>%
  nest() %>%
  mutate(res = purrr::map(.x = data, .f = compute_conditional_variance, spatial_cds_to_test)) %>%
  dplyr::select(-data) %>%
  unnest(cols = c("res")) %>%
  ungroup() %>%
  mutate(mean_ang_dist = mean(angular_distance),
         mean_ang_dist_sp = mean(angular_distance_sparsity),
         calculation = "spatial_cluster",
         permuted = F,
         run.id = run.id)


write.table(cluster_spatial_model,
            file = paste("/net/trapnell/vol1/home/sanjays/projects/Space/Submission/E14_slides/RDS_intermediates/variance_decomposition/",
                         "cluster.spatial.model_",
                         run.id,
                         ".tsv",
                         sep = ""),
            col.names = T,
            row.names = F,
            sep = "\t",
            quote = F)

# Celltype/spatial model with permuted labels
cluster_spatial_model_permuted =
  # Get the cells associated with each spatial group and calculate angular_distance to that mean
  spatial_cds_to_test %>%
  colData() %>%
  as.data.frame() %>%
  mutate(shuffled.cluster.spatial.bin = sample(cluster_spatial_bin,replace = F)) %>%
  group_by(shuffled.cluster.spatial.bin) %>%
  nest() %>%
  mutate(res = purrr::map(.x = data, .f = compute_conditional_variance, spatial_cds_to_test)) %>%
  dplyr::select(-data) %>%
  unnest(cols = c("res")) %>%
  ungroup() %>%
  mutate(mean_ang_dist = mean(angular_distance),
         mean_ang_dist_sp = mean(angular_distance_sparsity),
         calculation = "spatial_cluster",
         permuted = T,
         run.id = run.id)


write.table(cluster_spatial_model_permuted,
            file = paste("/net/trapnell/vol1/home/sanjays/projects/Space/Submission/E14_slides/RDS_intermediates/variance_decomposition/",
                         "cluster.spatial.model.permuted_",
                         run.id,
                         ".tsv",
                         sep = ""),
            col.names = T,
            row.names = F,
            sep = "\t",
            quote = F)



# Calculate a null model where only spatial label is permuted -------------
ct_spatial_model_permuted_celltype_marginal_slide = 
  # To calculat the variance explained for each celltype upon accounting for space,
  # calculate the average variance for a model with that many degrees of freedom where
  # the spatial label is permuted
  spatial_cds_to_test %>%
  colData() %>%
  as.data.frame() %>%
  # Make sure that spatial bins are only shuffled within a celltype, to make sure the
  # degrees of freedom are the same
  group_by(celltype_bin, max_slide_id) %>%
  mutate(shuffled.spatial.bin = sample(spatial_bin,replace = F)) %>%
  mutate(ct.shuffled.spatial.bin = paste(celltype_bin,shuffled.spatial.bin, sep = ".")) %>%
  ungroup() %>%
  group_by(ct.shuffled.spatial.bin) %>%
  nest() %>%
  mutate(res = purrr::map(.x = data, .f = compute_conditional_variance, spatial_cds_to_test)) %>%
  dplyr::select(-data) %>%
  unnest(cols = c("res")) %>%
  ungroup() %>%
  mutate(mean_ang_dist = mean(angular_distance),
         mean_ang_dist_sp = mean(angular_distance_sparsity),
         calculation = "marginal_celltype",
         permuted = T,
         run.id = run.id)



write.table(ct_spatial_model_permuted_celltype_marginal_slide,
            file = paste("/net/trapnell/vol1/home/sanjays/projects/Space/Submission/E14_slides/RDS_intermediates/variance_decomposition/celltype.permuted.spatial.model.marginal_",
                         run.id,
                         ".tsv",
                         sep = ""),
            col.names = T,
            row.names = F,
            sep = "\t",
            quote = F)



# Calculate a null model where only spatial label is permuted -------------
cluster_spatial_model_permuted_celltype_marginal_slide = 
  # To calculat the variance explained for each celltype upon accounting for space,
  # calculate the average variance for a model with that many degrees of freedom where
  # the spatial label is permuted
  spatial_cds_to_test %>%
  colData() %>%
  as.data.frame() %>%
  # Make sure that spatial bins are only shuffled within a celltype, to make sure the
  # degrees of freedom are the same
  group_by(celltype_bin, max_slide_id) %>%
  mutate(shuffled.spatial.bin = sample(spatial_bin,replace = F)) %>%
  mutate(cluster.shuffled.spatial.bin = paste(paste("cluster",cluster,sep = "_"),shuffled.spatial.bin, sep = ".")) %>%
  ungroup() %>%
  group_by(cluster.shuffled.spatial.bin) %>%
  nest() %>%
  mutate(res = purrr::map(.x = data, .f = compute_conditional_variance, spatial_cds_to_test)) %>%
  dplyr::select(-data) %>%
  unnest(cols = c("res")) %>%
  ungroup() %>%
  mutate(mean_ang_dist = mean(angular_distance),
         mean_ang_dist_sp = mean(angular_distance_sparsity),
         calculation = "marginal_cluster",
         permuted = T,
         run.id = run.id)



write.table(cluster_spatial_model_permuted_celltype_marginal_slide,
            file = paste("/net/trapnell/vol1/home/sanjays/projects/Space/Submission/E14_slides/RDS_intermediates/variance_decomposition/cluster.permuted.spatial.model.marginal_",
                         run.id,
                         ".tsv",
                         sep = ""),
            col.names = T,
            row.names = F,
            sep = "\t",
            quote = F)


# Read in computed values and collate results -------------------------------------------------

# Load startup packages ---------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(ggridges)
  library(purrr)
  library(monocle3)
  library(sp)
  library(sf)
  
  space_directory = "/Volumes/GoogleDrive/My Drive/sciSpace/"
  setwd(dir=space_directory)
})

file_paths = paste("Submission_Data/E14_slides/RDS_intermediates/variance_decomposition/",
                   list.files(path = "Submission_Data/E14_slides/RDS_intermediates/variance_decomposition/",pattern = "tsv$"),
                   sep = "")

results = 
  data.frame(file_path = file_paths,
             file_name = list.files(path = "Submission_Data/E14_slides/RDS_intermediates/variance_decomposition/",pattern = "tsv$"))

results = 
  results %>%
  mutate(file_path= as.character(file_path)) %>%
  group_by(file_path) %>%
  nest() %>%
  mutate(var.decomp = purrr:: map(.x =file_path,
                                  .f = function(x){
                                    read.table(x, header = TRUE, sep = '\t')}))

results =
  results %>%
  unnest(cols = "data")

results$model_type = 
  stringr::str_split_fixed(string = results$file_name,pattern = "_",n = 2)[,1]


# Estimate Technical Variance ---------------------------------------------
# A lot of the variance can be attributed to technical factors like sampling
# estimate the technical variation by sampling cells from the mean pseudobulked
# CDS and ask how much variance is present

global_mean = unit_mean(counts(spatial_cds_to_test))

resampled_cds = sampled_cds(spatial_cds_to_test[,test_cells])

technical_variance_estimate =
  lapply(seq(1,ncol(resampled_cds)),
         function(X,
                  mat,
                  shuffled_mat,
                  ang_mean){
           
           this.cell =
             mat[,X] %>%
             unit_vector()
           
           this.shuffled.cell =
             shuffled_mat[,X] %>%
             unit_vector()
           
           data.frame(ang.dist = angular_distance(this.cell, ang_mean),
                      shuffled.ang.dist = angular_distance(this.shuffled.cell, ang_mean))
         },spatial_cds_to_test[,test_cells] %>%
           monocle3:::normalize_expr_data(norm_method = "log",
                                          pseudo_count = 1),
         resampled_cds %>%
           monocle3:::normalize_expr_data(norm_method = "log",
                                          pseudo_count = 1),
         global_mean)

technical_variance_estimate =
  do.call(rbind,
          technical_variance_estimate)

mean(technical_variance_estimate$shuffled.ang.dist/
       technical_variance_estimate$ang.dist) *100



# Supplemental Figure 29  -------------------------------------------------

cols_to_select =
  c("mean_ang_dist",
    "mean_ang_dist_sp",
    "calculation",
    "permuted",
    "run.id")

tmp.list = list()
for(i in seq(1,length(results$var.decomp))) {
  tmp.list[[i]] = 
    results$var.decomp[[i]][,cols_to_select]  %>%
    distinct()}

do.call(rbind,
        tmp.list) %>%
  group_by(calculation,permuted) %>%
  summarise(mean_ang_dist =mean(mean_ang_dist),
            mean_ang_dist_sp = mean(mean_ang_dist_sp)) %>%
  mutate(delta = mean_ang_dist - mean_ang_dist_sp) %>%
  dplyr::select(calculation,permuted,delta) %>%
  spread(permuted,delta) %>%
  mutate(var_explained = 1- `FALSE`/`TRUE`) %>%
  drop_na() %>% 
  filter(!grepl("cluster",calculation)) %>%
  ggplot() +
  geom_bar(aes(x = reorder(calculation,var_explained),
               y = var_explained*100,
               fill = calculation),
           stat = "identity",
           size = 0.5,
           color = "black") +
  monocle3:::monocle_theme_opts() +
  scale_fill_brewer(palette = "Set1") +
  ylab("Variance Explained") +
  xlab("Model") +
  scale_y_continuous(breaks = seq(0,100,10),
                     labels = c(0,paste(seq(10,100,10),
                                        "%",
                                        sep =""))) +
  scale_x_discrete(labels = c("Cell Type",
                              "Spatial",
                              "Cell Type &\nSpatial")) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) +
  ggsave("Figures/Figure_Components/Figure4/Variance_Explained.pdf",
         height = 3,
         width = 3)


# Calculate celltype marginals after accounting for space ---------------------------

full_model =
  results %>% 
  filter(model_type == "celltype.spatial.model") %>%
  pull(var.decomp) %>%
  do.call(what = rbind) %>%
  mutate(celltype = stringr::str_split_fixed(celltype_spatial_bin,
                                             pattern = "\\.",
                                             n = 4)[,1])%>%
  group_by(celltype) %>%
  summarise(mean_ct_angular_distance = mean(angular_distance),
            mean_ct_angular_distance_sparsity = mean(angular_distance_sparsity))

reduced_model =
  results %>% 
  filter(model_type == "celltype.permuted.spatial.model.marginal") %>%
  pull(var.decomp) %>%
  do.call(what = rbind) %>%
  mutate(celltype = stringr::str_split_fixed(ct.shuffled.spatial.bin,
                                             pattern = "\\.",
                                             n = 4)[,1])%>%
  group_by(celltype) %>%
  summarise(permutated_mean_ct_angular_distance = mean(angular_distance),
            permuted_mean_ct_angular_distance_sparsity = mean(angular_distance_sparsity))



full_join(full_model,
          reduced_model,
          by = "celltype") %>%
  mutate(numerator = (mean_ct_angular_distance - mean_ct_angular_distance_sparsity),
         denominator = (permutated_mean_ct_angular_distance - permuted_mean_ct_angular_distance_sparsity),
         var_explained =  1 - numerator / denominator) %>%
  filter(denominator > 0,
         numerator > 0,
         denominator > numerator) %>%
  dplyr::select(var_explained,
                celltype) %>%
  ggplot() +
  geom_bar(aes(x = reorder(celltype,var_explained),
               y = var_explained),
           stat = "identity",
           fill = "grey80",
           color = "black",
           size = 0.25) +
  coord_flip() +
  monocle:::monocle_theme_opts() +
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8)) +
  xlab("") +
  ylab("Variance Explained") +
  ggsave("Figures/Figure_Components/Figure4/variance_explained_celltype.pdf",
         dpi = 600,
         height = 2,
         width = 3)


# Break values out by slide

full_model_slide =
  results %>% 
  filter(model_type == "celltype.spatial.model") %>%
  pull(var.decomp) %>%
  do.call(what = rbind) %>%
  mutate(celltype = stringr::str_split_fixed(celltype_spatial_bin,
                                             pattern = "\\.",
                                             n = 4)[,1],
         slide = stringr::str_split_fixed(celltype_spatial_bin,
                                          pattern = "\\.",
                                          n = 4)[,4])%>%
  group_by(celltype,slide) %>%
  summarise(mean_ct_angular_distance = mean(angular_distance),
            mean_ct_angular_distance_sparsity = mean(angular_distance_sparsity))


reduced_model_slide =
  results %>% 
  filter(model_type == "celltype.permuted.spatial.model.marginal") %>%
  pull(var.decomp) %>%
  do.call(what = rbind) %>%
  mutate(celltype = stringr::str_split_fixed(ct.shuffled.spatial.bin,
                                             pattern = "\\.",
                                             n = 4)[,1],
         slide = stringr::str_split_fixed(ct.shuffled.spatial.bin,
                                          pattern = "\\.",
                                          n = 4)[,4])%>%
  group_by(celltype,slide) %>%
  summarise(permutated_mean_ct_angular_distance = mean(angular_distance),
            permuted_mean_ct_angular_distance_sparsity = mean(angular_distance_sparsity))


# Figure 4 — Panel B ------------------------------------------------------


slide_colors = 
  c("slide_1D" = "#8DC9DA",
    "slide_1E" = "#4FB961",
    "slide_1F" = "#C01E67",
    "slide_1G" = "grey",
    "slide_2G" = "#C0A630",
    "slide_2H" = "#CD3029",
    "slide_3D" = "#E2AFCE",
    "slide_3F" = "#876F50",
    "slide_3G" = "#8538E8",
    "slide_3H" = "#60E7B2",
    "slide_4A" = "#E7E041",
    "slide_4C" = "#EE7733",
    "slide_4D" =  "#0077BB", 
    "slide_4E" = "#E744DC"
  )

spatial_cds =
  readRDS(file = "Submission_Data/E14_slides/RDS_intermediates/Notebook6_spatial_cds_anatomy.RDS")



full_join(full_model_slide,
          reduced_model_slide,
          by = c("celltype","slide")) %>%
  left_join(colData(spatial_cds) %>%
              as.data.frame() %>%
              group_by(max_slide_id,final_cluster_label) %>%
              summarise(n = n()),
            by = c("slide" = "max_slide_id",
                   "celltype" = "final_cluster_label")) %>%
  mutate(numerator = (mean_ct_angular_distance - mean_ct_angular_distance_sparsity),
         denominator = (permutated_mean_ct_angular_distance - permuted_mean_ct_angular_distance_sparsity),
         var_explained =  1 - numerator / denominator) %>% 
  dplyr::select(var_explained,
                celltype,
                slide,
                n) %>%
  filter(n > 200) %>% 
  mutate(var_explained = 
           ifelse(var_explained < 0,
                0,
                var_explained)) %>% 
  filter(var_explained < 1) %>% 
  group_by(celltype,
           slide) %>%
  mutate(var_explained = 100 * var_explained,
         median_var = mean(var_explained)) %>% 
  ungroup() %>% 
  ggplot() +
  geom_boxplot(aes(x = reorder(celltype,median_var),
                   y = var_explained),
               fill = "grey90",
               color = "black",
               size = 0.25,
               outlier.stroke = 0,
               outlier.size = 0) +
  geom_jitter(aes(x = reorder(celltype,median_var),
                 y = var_explained,
                 fill = slide,
                 size = log10(n)/10),
              shape = 21,
              stroke = 0.2,
              color = "black",
              width = 0.15,
              height = 0)+
  coord_flip() +
  monocle:::monocle_theme_opts() +
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        legend.position = "none")  +
  scale_fill_manual(values = slide_colors) +
  scale_size(range = c(0.25, 3)) +
  xlab("") +
  scale_y_continuous(breaks = seq(0,100,10),
                     labels = c(0,paste(seq(10,100,10),
                                        "%",
                                        sep ="")))  +
  ylab("Var. Explained by Spatial Position") +
  ggsave("Figures/Figure_Components/Figure4/variance_explained_celltype_slide.pdf",
         dpi = 600,
         height = 2,
         width = 3)

