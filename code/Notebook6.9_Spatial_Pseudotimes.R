# Code used to generate figure 5 -- We isolate three distinct trajectories 
# in the developing brain and map their spatial locations

# Load startup packages ---------------------------------------------------
suppressPackageStartupMessages({ 
  library(tidyverse)
  library(viridis)
  library(ggpubr)
  library(ggrepel)
  library(pheatmap)
  library(monocle3)
  library(FNN)
  
  space_directory = "/Volumes/GoogleDrive/My Drive/sciSpace/"
  setwd(dir=space_directory)

  # Pass TRUE if you want to see progress output on some of Monocle 3's operations
  DelayedArray:::set_verbose_block_processing(TRUE)
  # Passing a higher value will make some computations faster but use more memory. 
  # Adjust with caution!
  options(DelayedArray.block.size=1000e6)
  set.seed(42)
})

all_image_data = 
  readRDS("Submission_Data/E14_slides/RDS_intermediates/Notebook0_1_images_and_transformations.RDS")


spatial_cds = 
  readRDS(file = "Submission_Data/E14_slides/RDS_intermediates/Notebook6.01_spatial_cds_anatomy.RDS")

# Make a new CNS CDS ------------------------------------------------------

cns_cds = 
  new_cell_data_set(expression_data = counts(spatial_cds[,partitions(spatial_cds) == 2]),
                    cell_metadata = colData(spatial_cds[,partitions(spatial_cds) == 2]),
                    gene_metadata = rowData(spatial_cds[,partitions(spatial_cds) == 2])) %>%
  estimate_size_factors() %>%
  preprocess_cds(method = "PCA") %>%
  align_cds(residual_model_formula_str = "~log.n.umi + sample",
            alignment_group = "sample") %>%
  reduce_dimension(umap.metric = "cosine",
                   umap.fast_sgd = F,
                   max_components = 3) %>%
  cluster_cells(random_seed = 42)

# These cells form the three parallel trajectories
cells_for_traj = names(clusters(cns_cds))[clusters(cns_cds) %in% c(2,3,4,5)]

plot_cells_3d(cns_cds,color_cells_by = "cluster")


# Neural Trajectory CDS analysis ------------------------------------------
trajectory_cds = 
  new_cell_data_set(expression_data = counts(spatial_cds[,cells_for_traj]),
                    cell_metadata = colData(spatial_cds[,cells_for_traj]),
                    gene_metadata = rowData(spatial_cds[,cells_for_traj])) %>%
  estimate_size_factors() %>%
  preprocess_cds(method = "PCA") %>%
  align_cds(residual_model_formula_str = "~log.n.umi + sample",
            alignment_group = "sample") %>%
  reduce_dimension(umap.metric = "cosine",
                   umap.fast_sgd = F,
                   max_components = 2) %>%
  cluster_cells(random_seed = 42)

# Learn a principle graph
trajectory_cds = learn_graph(trajectory_cds)

colData(trajectory_cds)$umap1 = reducedDim(trajectory_cds,type = "UMAP")[,1]
colData(trajectory_cds)$umap2 = reducedDim(trajectory_cds,type = "UMAP")[,2]


# Figure 5 A --------------------------------------------------------------

colors =c("Choroid Plexus"  = "#00A6A6",
          "Fibroblast"  = "#790009",
          "Glial Cells" = "#8172E7",
          "Neuron" = "#B4D900",
          "OPCs" = "#E7E041",
          "Radial glia" = "#FF00C3")

ggplot() +
  geom_point(data =
               colData(trajectory_cds) %>%
               as.data.frame() %>%
               filter(!is.na(brain_region)),
             aes(x = umap1,
                 y = umap2),
             color = "black",
             stroke = 0,
             size = .75) +
  geom_point(data =
               colData(trajectory_cds) %>%
               as.data.frame() %>%
               filter(!is.na(brain_region)),
             aes(x = umap1,
                 y = umap2,
                 color = final_cluster_label),
             stroke = 0,
             size = 0.5) +
  scale_color_manual(values = colors)  +
  theme_void() +
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Figure4/brain_celltype.png",
         height = 1,
         width = 1,
         dpi = 600,
         bg = "transparent")

# Figure 5C ---------------------------------------------------------------

ggplot() +
  geom_point(data =
               colData(trajectory_cds) %>%
               as.data.frame() %>%
               filter(!is.na(brain_region)),
             aes(x = umap1,
                 y = umap2),
             color = "black",
             stroke = 0,
             size = .75) +
  geom_point(data =
               colData(trajectory_cds) %>%
               as.data.frame() %>%
               filter(!is.na(brain_region)),
             aes(x = umap1,
                 y = umap2,
                 color = brain_region),
             stroke = 0,
             size = 0.5) +
  scale_color_manual(values = c("#0758BE",
                                "#4993BF",
                                "#2FAB7C",
                                "#F6E34A",
                                "#F99713",
                                "#DF424C"))  +
  theme_void() +
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Figure3/brain_region.png",
         height = 1,
         width = 1,
         dpi = 600,
         bg = "transparent")


# Learn Pseudotime on UMAP embedding --------------------------------------

# Choose the roots for each trajectory
# Chosen around the radial glia at the base of each path
trajectory_cds = order_cells(trajectory_cds)
plot_cells(trajectory_cds,color_cells_by = "pseudotime")

colData(trajectory_cds)$Pseudotime = pseudotime(trajectory_cds)

# Separate the three trajectories so that pseudotime can be scaled within each
# trajectory
pallium = choose_cells(trajectory_cds)
midbrain = choose_cells(trajectory_cds)

colData(trajectory_cds)$trajectory = 
  ifelse(colData(trajectory_cds)$Cell %in% colnames(pallium),
         "1",
         ifelse(colData(trajectory_cds)$Cell %in% colnames(midbrain),
                "3",
                "2"))
  
colData(trajectory_cds)$scaled_pseudotime =
  colData(trajectory_cds) %>%
  as.data.frame() %>%
  group_by(trajectory) %>%
  mutate(scaled_trajectory = Pseudotime/max((Pseudotime))) %>%
  pull(scaled_trajectory)
  
plot_cells(trajectory_cds,color_cells_by = "scaled_pseudotime")


# Figure 5D ---------------------------------------------------------------

ggplot() +
  geom_point(data =
               colData(trajectory_cds) %>%
               as.data.frame() %>%
               filter(!is.na(brain_region)),
             aes(x = umap1,
                 y = umap2),
             color = "black",
             stroke = 0,
             size = .75) +
  geom_point(data =
               colData(trajectory_cds) %>%
               as.data.frame() %>%
               filter(!is.na(brain_region)),
             aes(x = umap1,
                 y = umap2,
                 color = scaled_pseudotime),
             stroke = 0,
             size = 0.5) +
  scale_color_viridis_c(option = "C")  +
  theme_void() +
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Figure4/brain_pseudotime.png",
         height = 1,
         width = 1,
         dpi = 600,
         bg = "transparent")


ggplot() +
  geom_point(data =
               colData(trajectory_cds) %>%
               as.data.frame() %>%
               filter(!is.na(brain_region)),
             aes(x = umap1,
                 y = umap2,
                 color = scaled_pseudotime),
             stroke = 0,
             size = 0.5) +
  scale_color_viridis_c(option = "C")  +
  theme_void() +
  ggsave("Figures/Figure_Components/Figure4/brain_pseudotime_legend.pdf",
         height = 5,
         width = 5)



# Figure 5E ---------------------------------------------------------------

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

jittered_points = 
  colData(spatial_cds) %>%
  as.data.frame() %>%
  left_join(scaling_df,
            by = "max_slide_id") %>%
  filter(!(Cell %in% colnames(trajectory_cds)),
         anatomical_annotation == "Cortex",
         slide_id %in% c("Slide 8",
                         "Slide 9",
                         "Slide 13",
                         "Slide 11",
                         "Slide 14")) %>%
  dplyr::select(-brain_region)

jittered_points =
  jittered_points %>%
mutate(jittered_y = coords.x1 * y + rnorm(n =dim(jittered_points)[1],
                                          mean = 0,
                                          sd = 10),
       jittered_x = coords.x2 * x + rnorm(n =dim(jittered_points)[1],
                                          mean = 0,
                                          sd = 10))

jittered_points2 = 
  colData(trajectory_cds) %>%
  as.data.frame() %>%
  left_join(scaling_df,
            by = "max_slide_id") %>%
  filter(is.finite(Pseudotime),
         anatomical_annotation == "Cortex",
         final_cluster_label %in% c("Neuron","Radial glia"),
         slide_id %in% c("Slide 8",
                         "Slide 9",
                         "Slide 13",
                         "Slide 11",
                         "Slide 14"))
  
jittered_points2 =
  jittered_points2 %>%
  mutate(jittered_y = coords.x1 * y + rnorm(n =dim(jittered_points2)[1],
                                            mean = 0,
                                            sd = 10),
         jittered_x = coords.x2 * x + rnorm(n =dim(jittered_points2)[1],
                                            mean = 0,
                                            sd = 10))


ggplot() +
  geom_point(data =jittered_points,
              aes(y = jittered_y,
                  x = jittered_x),
              color = "grey90",
              size = .5,
              stroke = 0) +
  geom_point(data = jittered_points2,
             aes(y = jittered_y,
                 x = jittered_x),
             color = "black",
             size = 0.65,
             stroke = 0) +
  geom_point(data = jittered_points2,
              aes(y = jittered_y,
                  x = jittered_x,
                  color = scaled_pseudotime),
              size = 0.5,
              stroke = 0) +
  facet_wrap(~slide_id,scales = "free",ncol =1) +
  monocle3:::monocle_theme_opts()+
  theme_void() +
  theme(legend.position = "none",
        strip.text.x = element_text(size = 6)) +
  scale_color_viridis_c(option = "C") +
  ggsave("Figures/Figure_Components/Figure4/brain_spatial_pseudotime.pdf",
         height = 5.5,
         width = 1.35)


# Supplemental Figure 33  ----------------------------------------------------

ggplot() +
  geom_jitter(data =
                colData(spatial_cds) %>%
                as.data.frame() %>%
                left_join(scaling_df,
                          by = "max_slide_id") %>%
                filter(!(Cell %in% colnames(trajectory_cds)),
                       anatomical_annotation == "Cortex") %>%
                dplyr::select(-brain_region),
              aes(y = coords.x1 * y,
                  x = coords.x2 * x),
              color = "grey70",
              size = 1,
              height = 14,
              width = 14,
              stroke = 0) +
  geom_jitter(data = 
                colData(trajectory_cds) %>%
                as.data.frame() %>%
                left_join(scaling_df,
                          by = "max_slide_id") %>%
                filter(is.finite(Pseudotime),
                       anatomical_annotation == "Cortex",
                       final_cluster_label %in% c("Neuron","Radial glia")),
              aes(y = coords.x1 * y,
                  x = coords.x2 * x,
                  color = scaled_pseudotime),
              size = .75,
              height = 14,
              width = 14,
              stroke = 0) +
  facet_wrap(~slide_id,scales = "free") +
  monocle3:::monocle_theme_opts()+
  theme_void() +
  theme(legend.position = "none",
        strip.text.x = element_text(hjust = 0)) +
  scale_color_viridis_c(option = "C") +
  ggsave("Figures/Figure_Components/Supplement_brain_trajectory/brain_spatial_pseudotime_all.png",
         dpi = 300,
         bg = "transparent")



# Plot gene expression for Figure 5B ----------------------------------------------------

get_umis = 
  function(cds,gene){
    this.gene = rowData(cds)$id[rowData(cds)$gene_short_name == gene]
    expr = counts(cds)[this.gene,] %>% as.vector()
    
    e = 
      colData(cds) %>%
      as.data.frame() %>%
      mutate(umap1 = reducedDim(x = cds, type = "UMAP")[,1],
             umap2 = reducedDim(x = cds, type = "UMAP")[,2],
             expr = expr)
    
    p = 
      ggplot()+ 
      geom_point(data = 
                   e %>%
                   filter(expr == 0),
                 aes(x= umap1,
                     y = umap2),
                 color = "grey80",
                 stroke = 0,
                 size = 0.25,
                 alpha = 0.25) +
      geom_point(data = 
                   e %>%
                   filter(expr != 0) %>%
                   arrange(expr),
                 aes(x= umap1,
                     y = umap2,
                     color = log2(expr)),
                 stroke = 0,
                 size = 0.25) +
      theme_void() +
      scale_color_viridis(option = "D",
                          limits = c(0,6)) +
      theme(legend.position = "none")
      
    return(p)
  }

# Midbrain Markers
get_umis(trajectory_cds, gene = "Pax3") +
  ggsave("Figures/Figure_Components/Figure4/Pax3.png",
         height = 0.5,
         width = 0.5,
         dpi = 600)

get_umis(trajectory_cds, gene = "Pou4f1") +
  ggsave("Figures/Figure_Components/Figure4/Pou4f1.png",
         height = 0.5,
         width = 0.5,
         dpi = 600)

get_umis(trajectory_cds, gene = "Tfap2b") +
  ggsave("Figures/Figure_Components/Figure4/Tfap2b.png",
         height = 0.5,
         width = 0.5,
         dpi = 600)

# GABAergic Markers
get_umis(trajectory_cds, gene = "Gad2") +
  ggsave("Figures/Figure_Components/Figure4/Gad2.png",
         height = 0.5,
         width = 0.5,
         dpi = 600)
get_umis(trajectory_cds, gene = "Isl1") +
  ggsave("Figures/Figure_Components/Figure4/Isl1.png",
         height = 0.5,
         width = 0.5,
         dpi = 600)
get_umis(trajectory_cds, gene = "Lhx6") +
  ggsave("Figures/Figure_Components/Figure4/Lhx6.png",
         height = 0.5,
         width = 0.5,
         dpi = 600)

# Radial Glia
get_umis(trajectory_cds, gene = "Sox9")+
  ggsave("Figures/Figure_Components/Figure4/Sox9.png",
         height = 0.5,
         width = 0.5,
         dpi = 600)

get_umis(trajectory_cds, gene = "Pax6")+
  ggsave("Figures/Figure_Components/Figure4/Pax6.png",
         height = 0.5,
         width = 0.5,
         dpi = 600)

get_umis(trajectory_cds, gene = "Eomes") +
  ggsave("Figures/Figure_Components/Figure4/Eomes.png",
         height = 0.5,
         width = 0.5,
         dpi = 600)

get_umis(trajectory_cds, gene = "Eomes") +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.key.height = ) +
  ggsave("Figures/Figure_Components/Figure4/gene_expression_legend.png",
         height = 2,
         width = 2,
         dpi = 300)


get_umis(trajectory_cds, gene = "Neurod6")+
  ggsave("Figures/Figure_Components/Figure4/Neurod6.png",
         height = 0.5,
         width = 0.5,
         dpi = 600)

get_umis(trajectory_cds, gene = "Emx1")+
  ggsave("Figures/Figure_Components/Figure4/Emx1.png",
         height = 0.5,
         width = 0.5,
         dpi = 600)

get_umis(trajectory_cds, gene = "Tfap2d")+
  ggsave("Figures/Figure_Components/Figure4/Tfap2d.png",
         height = 0.5,
         width = 0.5,
         dpi = 600)


get_umis(trajectory_cds, gene = "Reln")+
  ggsave("Figures/Figure_Components/Supplement_brain_trajectory/Reln.png",
         height = 0.5,
         width = 0.5,
         dpi = 600)


get_umis(trajectory_cds, gene = "Dab1")+
  ggsave("Figures/Figure_Components/Supplement_brain_trajectory/Dab1.png",
         height = 0.5,
         width = 0.5,
         dpi = 600)



get_umis(trajectory_cds, gene = "Crmp1")+
  ggsave("Figures/Figure_Components/Supplement_brain_trajectory/Crmp1.png",
         height = 0.5,
         width = 0.5,
         dpi = 600)

get_umis(trajectory_cds, gene = "Dcx")+
  ggsave("Figures/Figure_Components/Supplement_brain_trajectory/Dcx.png",
         height = 0.5,
         width = 0.5,
         dpi = 600)


# Pseudotime dependent genes ----------------------------------------------

# Pallium pseudotime models -- Isolate pallium cells
pallium_cells_trajectory_1 = 
  trajectory_cds %>%
  colData() %>%
  as.data.frame() %>%
  filter(brain_region == "pallium",
         final_cluster_label %in% c("Neuron",
                                    "Radial glia")) %>%
  pull(Cell)

pallium_trajectory = choose_cells(trajectory_cds[,pallium_cells_trajectory_1])
pallium_trajectory = order_cells(pallium_trajectory)
plot_cells(pallium_trajectory,color_cells_by = "pseudotime")

# Fit a linear model using pseudotime as the predictor
pallium_trajectory_degs = 
  fit_models(pallium_trajectory,
             model_formula_str = "~splines::ns(pseudotime, df=3)")

new_data_pallium = 
  data.frame(
    Size_Factor = 1 ,
    pseudotime = seq(min(pseudotime(pallium_trajectory)),
                     max(pseudotime(pallium_trajectory)), 
                     length.out= 100))

# Use the model to predict the gene expression along this trajectory
predictions_pallium =
  model_predictions(model_tbl = pallium_trajectory_degs, 
                    new_data=new_data_pallium)


coeff_pallium =  
  pallium_trajectory_degs %>%
  filter(status == "OK") %>%
  coefficient_table() 

sig_genes_pallium = 
  coeff_pallium %>% 
  filter(term != "(Intercept)",
         q_value < 0.01) %>% 
  pull(id) %>%
  as.character()


# Sub Pallium Trajectory --------------------------------------------------

sub_pallium_cells_trajectory_2 = 
  trajectory_cds %>%
  colData() %>%
  as.data.frame() %>%
  filter(brain_region == "sub pallium",
         final_cluster_label %in% c("Neuron",
                                    "Radial glia")) %>%
  pull(Cell)

sub_pallium_trajectory = choose_cells(trajectory_cds[,sub_pallium_cells_trajectory_2])
sub_pallium_trajectory = order_cells(sub_pallium_trajectory)


# Fit a linear model using pseudotime as the predictor
sub_pallium_trajectory_degs = 
  fit_models(sub_pallium_trajectory,
             model_formula_str = "~splines::ns(pseudotime, df=3)")

new_data_sub_pallium = 
  data.frame(
    Size_Factor = 1 ,
    pseudotime = seq(min(pseudotime(sub_pallium_trajectory)),
                     max(pseudotime(sub_pallium_trajectory)), 
                     length.out= 100))

# Use the model to predict the gene expression along this trajectory
predictions_sub_pallium =
  model_predictions(model_tbl = sub_pallium_trajectory_degs, 
                    new_data=new_data_sub_pallium)


coeff_sub_pallium =  
  sub_pallium_trajectory_degs %>%
  filter(status == "OK") %>%
  coefficient_table() 

sig_genes_sub_pallium = 
  coeff_sub_pallium %>% 
  filter(term != "(Intercept)",
         q_value < 0.01) %>% 
  pull(id) %>%
  as.character()


# Midbrain Pseudonyme DEGs ------------------------------------------------

midbrain_cells_trajectory_3 = 
  trajectory_cds %>%
  colData() %>%
  as.data.frame() %>%
  filter(brain_region == "midbrain",
         final_cluster_label %in% c("Neuron",
                                    "Radial glia")) %>%
  pull(Cell)

midbrain_trajectory = choose_cells(trajectory_cds[,midbrain_cells_trajectory_3])
midbrain_trajectory = order_cells(midbrain_trajectory)

# Fit a linear model using pseudotime as the predictor
midbrain_trajectory_degs = 
  fit_models(midbrain_trajectory,
             model_formula_str = "~splines::ns(pseudotime, df=3)")

new_data_midbrain = 
  data.frame(
    Size_Factor = 1 ,
    pseudotime = seq(min(pseudotime(midbrain_trajectory)),
                     max(pseudotime(midbrain_trajectory)), 
                     length.out= 100))

# Use the model to predict the gene expression along this trajectory
predictions_midbrain =
  model_predictions(model_tbl = midbrain_trajectory_degs, 
                    new_data=new_data_midbrain)


coeff_midbrain =  
  midbrain_trajectory_degs %>%
  filter(status == "OK") %>%
  coefficient_table() 

sig_genes_midbrain = 
  coeff_midbrain %>% 
  filter(term != "(Intercept)",
         q_value < 0.01) %>% 
  pull(id) %>%
  as.character()


# Write out pseudotime dependent genes for File S3 -----------------------

columns.to.select = 
  c("id","gene_short_name","term","estimate","std_err","test_val","p_value","q_value")

rbind(coeff_pallium %>%
        filter(status == "OK") %>%
        dplyr::select(columns.to.select) %>%
        mutate(trajectory = "Pallium"),
      coeff_sub_pallium %>%
        filter(status == "OK") %>%
        dplyr::select(columns.to.select) %>%
        mutate(trajectory = "Sub Pallium"),
      coeff_midbrain %>%
        filter(status == "OK") %>%
        dplyr::select(columns.to.select) %>%
        mutate(trajectory = "Midbrain")) %>%
  write.table("Submission_Data/E14_slides/RDS_intermediates/Supplemental_File_S3.tsv",
              quote = F,
              sep = "\t",
              col.names = T,
              row.names = F)
      
# Figure 5F - Pseudotime Heatmaps -----------------------

row_center = function(m){
  # Row-center the data.
  m=m[!apply(m,1,sd)==0,]
  m=Matrix::t(scale(Matrix::t(m),center=TRUE))
  m=m[is.na(row.names(m)) == FALSE,]
  m[is.nan(m)] = 0
  return(m)
}

union_all_genes = 
  union(sig_genes_midbrain,
        sig_genes_pallium) %>%
  union(sig_genes_sub_pallium)

print(paste("Significant genes in the pallium trajectory",
            length(sig_genes_pallium %>% unique())))

print(paste("Significant genes in the sub pallium trajectory",
            length(sig_genes_sub_pallium %>% unique())))

print(paste("Significant genes in the midbrain trajectory",
            length(sig_genes_midbrain %>% unique())))

print(paste("Significant genes across the three trajectories",
            length(union_all_genes)))

intersect_all_genes = 
  intersect(sig_genes_midbrain,
        sig_genes_pallium) %>%
  intersect(sig_genes_sub_pallium)

print(paste("Significant genes common to the three trajectories",
            length(intersect_all_genes)))



predictions_intersection = 
  cbind(row_center(predictions_midbrain[intersect_all_genes,]),
        row_center(predictions_sub_pallium[intersect_all_genes,]),
        row_center(predictions_pallium[intersect_all_genes,]))



#predictions_intersection = row_center(predictions_intersection)
predictions_intersection[predictions_intersection > 3] <- 3
predictions_intersection[predictions_intersection < -1] <- -1

bks <- seq(-1,3, by = 0.1)
hmcols <- viridis::viridis(length(bks) - 1,option = "D")

col_gaps_ind = c(dim(new_data_midbrain)[1],dim(new_data_midbrain)[1] +dim(new_data_sub_pallium)[1])

rownames(predictions_intersection) = rowData(trajectory_cds)[rownames(predictions_intersection),"gene_short_name"]
rowSums(is.na(predictions_intersection)) %>% sort()

ph = 
  pheatmap(predictions_intersection,
           cluster_rows = T,
           cluster_cols = F,
           gaps_col = col_gaps_ind,
           treeheight_row = 0,
           treeheight_col = 0,
           show_rownames = F,
           show_colnames = F,
           breaks=bks,
           legend = F,
           color=hmcols,
           cellwidth = 0.15,
           cellheight = 0.15,
           cutree_rows = 4)

# Used as input to geneontology.org for characterization
cutree(ph$tree_row, k = 4)[cutree(ph$tree_row, k = 4) == 3] %>%
  names()%>%
  write.table("~/Desktop/intersected_genes.tsv",sep = "\t",row.names = F,col.names = F, quote = F)

col.annotations = 
  data.frame(row.names = cutree(ph$tree_row,k = 4) %>% names(),
             cut = cutree(ph$tree_row,k = 4) %>% as.character())
annotation_colors = 
  list(cut = c("1" = RColorBrewer::brewer.pal(n = 4,name = "Set1")[1],
         "2" = RColorBrewer::brewer.pal(n = 4,name = "Set1")[2],
         "4" = RColorBrewer::brewer.pal(n = 4,name = "Set1")[4],
         "3" = RColorBrewer::brewer.pal(n = 4,name = "Set1")[3]))                                

pheatmap(predictions_intersection,
           filename = "Figures/Figure_Components/Figure4/pseudotime_heatmap.png",
           cluster_rows = T,
           cluster_cols = F,
           gaps_col = col_gaps_ind,
           treeheight_row = 0,
           treeheight_col = 0,
           show_rownames = F,
           show_colnames = F,
           breaks=bks,
           legend = F,
           color=hmcols,
           cellwidth = 0.15,
           cellheight = 0.15)

pheatmap(predictions_intersection,
         filename = "Figures/Figure_Components/Figure4/pseudotime_heatmap_legend.png",
         cluster_rows = T,
         cluster_cols = F,
         gaps_col = col_gaps_ind,
         treeheight_row = 0,
         treeheight_col = 0,
         show_rownames = F,
         show_colnames = F,
         breaks=bks,
         legend = T,
         color=hmcols,
         cellwidth = 0.05,
         cellheight = 0.45)

# Trajectory Characterization Fig S32 ---------------------------------------------

laManno_coldata = 
  read.table(file = "Submission_Data/E14_slides/RDS_intermediates/LaManno_coldata.tsv",
             header = T,
             row.names = 1,
             sep = "\t")

cluster_age = 
  lamanno_coldata %>%
  group_by(Clusters) %>%
  add_tally() %>%
  group_by(Clusters, Age.x,n) %>%
  summarise(nn = n()) %>%
  mutate(fraction = nn/n,
         Age = stringr::str_sub(string = Age.x,
                                start = 2) %>%
           as.numeric()) %>%
  group_by(Clusters) %>%
  summarise(mean_age = sum(Age*fraction))

# Supplemental Figure S32D
ggplot() +
  geom_point(data =
               colData(trajectory_cds) %>%
               as.data.frame() %>%
               filter(!is.na(brain_region)),
             aes(x = umap1,
                 y = umap2),
             color = "black",
             stroke = 0,
             size = .75) +
  geom_point(data =
               colData(trajectory_cds) %>%
               as.data.frame() %>%
               filter(!is.na(brain_region)) %>%
               left_join(cluster_age,
                         by = c("lamanno_nn_Cluster" = "Clusters")),
             aes(x = umap1,
                 y = umap2,
                 color = mean_age),
             stroke = 0,
             size = 0.5) +
  scale_color_viridis_c(option = "B")  +
  theme_void() +
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Supplement_brain_trajectory/brain_pseudotime_mean_age.png",
         height = 2,
         width = 2,
         dpi = 600,
         bg = "transparent") 

# Supplemental Figure S32D legend
ggplot() +
  geom_point(data =
               colData(trajectory_cds) %>%
               as.data.frame() %>%
               filter(!is.na(brain_region)),
             aes(x = umap1,
                 y = umap2),
             color = "black",
             stroke = 0,
             size = .75) +
  geom_point(data =
               colData(trajectory_cds) %>%
               as.data.frame() %>%
               filter(!is.na(brain_region)) %>%
               left_join(cluster_age,
                         by = c("lamanno_nn_Cluster" = "Clusters")),
             aes(x = umap1,
                 y = umap2,
                 color = mean_age),
             stroke = 0,
             size = 0.5) +
  scale_color_viridis_c(option = "B")  +
  theme_void() +
  theme(legend.position = "right") +
  ggsave("Figures/Figure_Components/Supplement_brain_trajectory/brain_pseudotime_mean_age_legend.png",
         height = 2,
         width = 2,
         dpi = 600,
         bg = "transparent") 

# Supplemental Figure S32E
colData(trajectory_cds) %>%
  as.data.frame() %>%
  filter(!is.na(brain_region)) %>%
  filter(Cell %in% c(colnames(pallium_trajectory),
                     colnames(sub_pallium_trajectory),
                     colnames(midbrain_trajectory))) %>%
  mutate(trajectory2 = ifelse(Cell %in% colnames(pallium_trajectory),
                              "Pallium",
                              ifelse(Cell %in% colnames(midbrain_trajectory),
                                     "Midbrain",
                                     "Sub Pallium"))) %>%
  left_join(cluster_age,
            by = c("lamanno_nn_Cluster" = "Clusters")) %>%
  ggplot(aes(x = mean_age,
             y = scaled_pseudotime)) +
  geom_point(size = 1) +
  geom_smooth(method = lm,
              se = T) +
  stat_cor(size=2, color="grey31",label.x.npc = "left", geom = "label") +
  facet_wrap(~trajectory2) +
  theme_classic() +
  theme(strip.background = element_rect(fill= NA,
                                        color = NA)) +
  ylim(0,1) + 
  xlab("Mean Age (Emb. Day)") +
  ylab("Scaled Pseudotime") +
  ggsave("Figures/Figure_Components/Supplement_brain_trajectory/lm.pdf",
         height = 2.5,
         width = 8)

# Supplemental Figure S32A 
ggplot() +
  geom_point(data =
               colData(trajectory_cds) %>%
               as.data.frame() %>%
               filter(!is.na(brain_region)),
             aes(x = umap1,
                 y = umap2),
             color = "black",
             stroke = 0,
             size = .75) +
  geom_point(data =
               colData(trajectory_cds) %>%
               as.data.frame() %>%
               filter(!is.na(brain_region)),
             aes(x = umap1,
                 y = umap2,
                 color = g1s_score),
             stroke = 0,
             size = 0.5) +
  scale_color_viridis_c(option = "D")  +
  theme_void() +
  theme(legend.position = "none") +
  guides(colour = guide_colorbar(title = "G1S Score")) +
  ggsave("Figures/Figure_Components/Supplement_brain_trajectory/brain_g1s_score.png",
         height = 2,
         width = 2,
         dpi = 600,
         bg = "transparent") 

# Supplemental Figure S32A legend
ggplot() +
  geom_point(data =
               colData(trajectory_cds) %>%
               as.data.frame() %>%
               filter(!is.na(brain_region)),
             aes(x = umap1,
                 y = umap2),
             color = "black",
             stroke = 0,
             size = .75) +
  geom_point(data =
               colData(trajectory_cds) %>%
               as.data.frame() %>%
               filter(!is.na(brain_region)),
             aes(x = umap1,
                 y = umap2,
                 color = g1s_score),
             stroke = 0,
             size = 0.5) +
  scale_color_viridis_c(option = "D")  +
  theme_void() +
  theme(legend.position = "right") +
  guides(colour = guide_colorbar(title = "G1S Score")) +
  ggsave("Figures/Figure_Components/Supplement_brain_trajectory/brain_g1s_score_legend.png",
         height = 2,
         width = 2,
         dpi = 600,
         bg = "transparent") 

# Supplemental Figure S32B 
ggplot() +
  geom_point(data =
               colData(trajectory_cds) %>%
               as.data.frame() %>%
               filter(!is.na(brain_region)),
             aes(x = umap1,
                 y = umap2),
             color = "black",
             stroke = 0,
             size = .75) +
  geom_point(data =
               colData(trajectory_cds) %>%
               as.data.frame() %>%
               filter(!is.na(brain_region)),
             aes(x = umap1,
                 y = umap2,
                 color = g2m_score),
             stroke = 0,
             size = 0.5) +
  scale_color_viridis_c(option = "D")  +
  theme_void() +
  theme(legend.position = "none") +
  guides(colour = guide_colorbar(title = "G2M\nScore")) +
  ggsave("Figures/Figure_Components/Supplement_brain_trajectory/brain_g2m_score.png",
         height = 2,
         width = 2,
         dpi = 600,
         bg = "transparent") 


# Supplemental Figure S32B legend
ggplot() +
  geom_point(data =
               colData(trajectory_cds) %>%
               as.data.frame() %>%
               filter(!is.na(brain_region)),
             aes(x = umap1,
                 y = umap2),
             color = "black",
             stroke = 0,
             size = .75) +
  geom_point(data =
               colData(trajectory_cds) %>%
               as.data.frame() %>%
               filter(!is.na(brain_region)),
             aes(x = umap1,
                 y = umap2,
                 color = g2m_score),
             stroke = 0,
             size = 0.5) +
  scale_color_viridis_c(option = "D")  +
  theme_void() +
  theme(legend.position = "right") +
  guides(colour = guide_colorbar(title = "G2M\nScore")) +
  ggsave("Figures/Figure_Components/Supplement_brain_trajectory/brain_g2m_score_legend.png",
         height = 2,
         width = 2,
         dpi = 600,
         bg = "transparent") 

