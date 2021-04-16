# Look at the spatial auto-correlation of genes in the spatial coordinates that are subset by
# cell type. For the same set of cells measure the auto-correlation of in UMAP space.

# UMAP-autocorrelation captures genes that are intrinsic to a cell state,
# Spatial-autocorrelation should capture genes that are clustered in space
# Genes with high values for both indicates that the cell state is spatially localized

# Load startup packages ---------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(purrr)
  library(magrittr)
  library(sp)
  library(sf)
  library(monocle3)
  library(ggrepel)
  library(ggridges)
  library(purrr)
  library(modelr)
  library(tidyr)
  library(dplyr)

  space_directory = "/Volumes/GoogleDrive/My Drive/sciSpace/"
  setwd(dir=space_directory)
  
  # Pass TRUE if you want to see progress output on some of Monocle 3's operations
  DelayedArray:::set_verbose_block_processing(TRUE)
  # Passing a higher value will make some computations faster but use more memory.
  options(DelayedArray.block.size=1000e6)
  set.seed(42)
})

spatial_cds = readRDS(file = "Submission_Data/E14_slides/RDS_intermediates/Notebook6_spatial_cds_anatomy.RDS")

all_image_data =
  readRDS(file = "Submission_Data/E14_slides/RDS_intermediates/Notebook0_1_images_and_transformations.RDS")

# Loop through cell types and calculate autocorrelation -------------------

collate_spatial_autocorrelation =
  function(coldata_subset,
           original_cds){

    # Subset the CDS for cells of this subset
    cells_to_test = coldata_subset$Cell
    cds_subset = original_cds[,cells_to_test]


    # Process the subsetted CDS and perform dimensionality reduction
    # and clustering
    cds_subset =
      cds_subset %>%
      estimate_size_factors() %>%
      detect_genes() %>%
      preprocess_cds() %>%
      reduce_dimension(max_components = 2,
                       reduction_method = "UMAP",
                       umap.fast_sgd = F,
                       preprocess_method = "PCA") %>%
      cluster_cells(random_seed = 42)

    # Add UMAP coordinates to the metadata
    colData(cds_subset)$umap1 =
      reducedDim(cds_subset,
                 type = "UMAP")[,1]

    colData(cds_subset)$umap2 =
      reducedDim(cds_subset,
                 type = "UMAP")[,2]

    # Add umap clusters to the metadata
    colData(cds_subset)$umap_cluster = cds_subset@clusters[["UMAP"]][["clusters"]]

    # Perform graph_test on UMAP coordinates this function
    # calculates autocorrelation on the nearest neighbor graph
    pr_graph_test_UMAP =
      graph_test(cds_subset,
                 neighbor_graph="knn",
                 reduction_method = "UMAP")

    pr_deg_ids_UMAP =
      row.names(subset(pr_graph_test_UMAP,
                       q_value < 0.05))

    # Move UMAP coordinates to the "tSNE" slot in reduced dimensions
    reducedDim(x = cds_subset,
               type = "tSNE") =
      reducedDim(cds_subset,
                 type = "UMAP")

    # Move clusters to the "tSNE" slot in reduced dimensions
    cds_subset@clusters[["tSNE"]] =
      cds_subset@clusters[["UMAP"]]

    # Check again to make sure things are correct
    identical(rownames(reducedDim(x = cds_subset,
                                  type = "tSNE")),
              colnames(cds_subset) %>% as.character())

    # Move spatial coordinates to the "UMAP" slot in reduced dimensions
    reducedDim(x = cds_subset,
               type = "UMAP") <-
      matrix(cbind(colData(cds_subset)$coords.x1,
                   colData(cds_subset)$coords.x2),
             ncol=2)

    # Define spatial clusters of cells
    cds_subset = cluster_cells(cds_subset,
                        resolution = 1e-2,
                        random_seed = 42,
                        reduction_method = "UMAP")

    # Perform graph_test on spatial coordinates this function
    # calculates autocorrelation on the nearest neighbor graph
    pr_graph_test_spatial =
      graph_test(cds_subset,
                 neighbor_graph="knn",
                 reduction_method = "UMAP")

    # Join the results together
    spatial_and_umap_stats =
      full_join( pr_graph_test_UMAP,
                 pr_graph_test_spatial,
                 by = c("id",
                        "gene_short_name")) %>%
      filter(status.x == "OK",
             status.y == "OK")

    colnames(spatial_and_umap_stats) =
      stringr::str_replace(colnames(spatial_and_umap_stats),
                           pattern = ".x",
                           replacement = "_umap") %>%
      stringr::str_replace(pattern = ".y",
                           replacement = "_spatial")

    spatial_and_umap_stats =
      spatial_and_umap_stats %>%
      left_join(rowData(cds_subset) %>%
                  as.data.frame(),
                by = c("id","gene_short_name"))

  }

# Tally autocorrelation results for every celltype/slide combination that
# has more than 100 cells
autocor_res =
  spatial_cds %>%
  colData() %>%
  as.data.frame() %>%
  filter(!is.na(coords.x1),
         !is.na(coords.x2))%>%
  group_by(max_slide_id,
           final_cluster_label) %>%
  add_tally() %>%
  filter(n > 100) %>%
  nest() %>%
  mutate(morans_I = purrr::map(.x = data,
                               .f = collate_spatial_autocorrelation,
                               spatial_cds)) %>%
  dplyr::select(-data) %>%
  unnest(cols = c("morans_I"))


# Write out results -------------------------------------------------------

saveRDS(object = autocor_res,
        file = "Submission_Data/E14_slides/RDS_intermediates/20201003_spatial_autocorrelation_results.RDS")

# Repeat Run without splitting out by slide -------------------------------

autocor_res_slide =
  spatial_cds %>%
  colData() %>%
  as.data.frame() %>%
  filter(!is.na(coords.x1),
         !is.na(coords.x2))%>%
  group_by(max_slide_id) %>%
  add_tally() %>%
  filter(n > 100) %>%
  nest() %>%
  mutate(morans_I = purrr::map(.x = data,
                               .f = collate_spatial_autocorrelation,
                               spatial_cds)) %>%
  dplyr::select(-data) %>%
  unnest(cols = c("morans_I"))


# Write out results -------------------------------------------------------

saveRDS(object = autocor_res_slide,
        file = "Submission_Data/E14_slides/RDS_intermediates/20201003_spatial_autocorrelation_results_per_slide.RDS")


# File S2 -----------------------------------------------------------------
autocor_res = readRDS("Submission_Data/E14_slides/RDS_intermediates/20201003_spatial_autocorrelation_results.RDS")
autocor_res_slide = readRDS("Submission_Data/E14_slides/RDS_intermediates/20201003_spatial_autocorrelation_results_per_slide.RDS")

slide_key =
  spatial_cds %>%
  colData() %>%
  as.data.frame() %>%
  dplyr::select(slide_id,
                max_slide_id) %>%
  distinct()

autocor_res_subset = 
  autocor_res %>%
  left_join(slide_key) %>%
  ungroup() %>%
  dplyr::select(-max_slide_id,
                cell_type = final_cluster_label,
                slide_id,
                -contains("umap"),
                id,
                gene_short_name,
                morans_I_spatial,
                q_value_spatial) %>%
  mutate(test = "cell_type")

autocor_res_slide_subset = 
  autocor_res_slide %>%
  left_join(slide_key) %>%
  ungroup() %>%
  mutate(cell_type = NA) %>%
  dplyr::select(-max_slide_id,
                cell_type,
                slide_id,
                -contains("umap"),
                id,
                gene_short_name,
                morans_I_spatial,
                q_value_spatial) %>%
  mutate(test = "slide")

rbind(autocor_res_subset,autocor_res_slide_subset) %>%
  write.table(file = "Supplemental_Tables/File_S2_autocorr_res_ct_slide.tsv",
              sep = "\t",
              col.names = T,
              quote = F,
              row.names = F)

autocor_res_slide %>%
  group_by(gene_short_name) %>%
  dplyr::select()



# Figure 3A ---------------------------------------------------------------

autocor_res =
  readRDS("Submission_Data/E14_slides/RDS_intermediates/20201003_spatial_autocorrelation_results.RDS")

autocor_res %>%
  filter(q_value_spatial < 0.001,
         morans_I_spatial > 0.1) %>%
  pull(id) %>%
  unique() %>%
  length


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


# Number of spatially autocorrelated genes per slide and cell type
autocor_res %>%
  filter(q_value_spatial < 0.001,
         num_cells_expressed > 10) %>%
  group_by(final_cluster_label, max_slide_id) %>%
  add_tally() %>%
  dplyr::select(final_cluster_label,
                max_slide_id,
                n) %>%
  distinct() %>%
  group_by(final_cluster_label) %>%
  ggplot() +
  geom_boxplot(aes(x = reorder(final_cluster_label,n),
                   y = log10(n)),
               size = 0.25,
               outlier.stroke = 0,
               outlier.size = 0) +
  geom_jitter(aes(x = reorder(final_cluster_label,n),
                  y = log10(n),
                  color = max_slide_id),
              size = 0.5) +
  monocle3:::monocle_theme_opts() +
  scale_color_manual(values = slide_colors,
                     name = "Slide ID",
                     labels = c("Slide 1",
                                "Slide 2",
                                "Slide 3",
                                "Slide 4",
                                "Slide 5",
                                "Slide 6",
                                "Slide 7",
                                "Slide 8",
                                "Slide 9",
                                "Slide 10",
                                "Slide 11",
                                "Slide 12",
                                "Slide 13",
                                "Slide 14")) +
  theme(axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        #axis.title.x  = element_text(size = 8),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "none") +
  coord_flip() +
  ylab("Spatially Significant Genes") +
  ggsave("Figures/Figure_Components/Figure3/spatial_genes.pdf",
         height = 2.1,
         width = 2.25)

ctp_genes = 
  autocor_res %>%
  filter(q_value_spatial < 0.001,
         num_cells_expressed > 10) %>%
  group_by(final_cluster_label, max_slide_id) %>%
  add_tally() %>%
  dplyr::select(final_cluster_label,
                max_slide_id,
                n) %>%
  distinct() %>%
  filter(final_cluster_label == "Connective Tissue Progenitors")
mean(ctp_genes$n)
sd(ctp_genes$n)


neuron_genes = autocor_res %>%
  filter(q_value_spatial < 0.001,
         num_cells_expressed > 10) %>%
  group_by(final_cluster_label, max_slide_id) %>%
  add_tally() %>%
  dplyr::select(final_cluster_label,
                max_slide_id,
                n) %>%
  distinct() %>%
  filter(final_cluster_label == "Neuron")
mean(neuron_genes$n)
sd(neuron_genes$n)



# 
# # Variance of the Ratio ---------------------------------------------------
# 
# autocor_res_slide = readRDS("Submission_Data/E14_slides/RDS_intermediates/20201003_spatial_autocorrelation_results_per_slide.RDS")
# 
# d = 
#   autocor_res_slide %>%
#   mutate(new_stat = abs(morans_I_spatial/morans_I_umap)) %>%
#   
#   dplyr::select(max_slide_id,
#                 gene_short_name,
#                 new_stat,
#                 num_cells_expressed) %>%
#   group_by(gene_short_name) %>%
#   summarise(mean_new_stat = mean(new_stat),
#             mean_num_cells_expressed = mean(num_cells_expressed),
#             var_new_stat = var(new_stat,na.rm  = T))%>%
#   filter(mean_num_cells_expressed > 10) %>%
#   arrange(var_new_stat)
# 
# ggplot() +
#   geom_point(data = d,
#              aes(x = (mean_num_cells_expressed),
#                  y = (var_new_stat)),
#              size = 1.05,
#              stroke = 0,
#              color = "black") +
#   geom_point(data = d,
#              aes(x = (mean_num_cells_expressed),
#                  y = (var_new_stat)),
#              size = 0.95,
#              stroke = 0,
#              color = "grey80") +
#   geom_density_2d(data = d,
#                   aes(x = mean_num_cells_expressed,
#                       y = var_new_stat),
#                   alpha = 0.75) +
#   scale_y_log10(breaks = c(0.01,0.1,1,10,1e2,1e3,1e4,1e5,1e6),
#                 labels = c("0.01","0.1","1","10","100",
#                            "1,000","10,000","100,000","1,000,000")) +
#   scale_x_log10(breaks = c(1,10,100,1e3,1e4),
#                 labels = c("1","10","100","1,000","10,000")) +
#   xlab("Avg. Cells Expressing Gene") +
#   ylab("Variance Moran's I [Space]/[UMAP]") +
#   monocle3:::monocle_theme_opts() +
#   ggsave("Figures/Figure_Components/Reviewer_Figures/variance_moransI.png",
#          dpi= 450,
#          bg = "transparent",
#          height = 4, 
#          width = 4)
#   
