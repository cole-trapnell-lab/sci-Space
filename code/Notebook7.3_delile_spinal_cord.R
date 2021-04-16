# Downloaded Delile et. al. dataset can be found here -- 
# https://github.com/juliendelile/MouseSpinalCordAtlas/tree/master/output

# Looking to see whether having more cells in the spinal cord, better defines  
# subsets of cells in which the Hox genes are restricted.


# Load startup packages ---------------------------------------------------
suppressPackageStartupMessages({ 
  library(tidyverse)
  library(ggplot2)
  library(tidyr)
  library(viridis)
  library(ggridges)
  library(purrr)
  library(monocle3)
  library(RColorBrewer)
  library(ggrepel)
  
  space_directory = "/Volumes/GoogleDrive/My Drive/sciSpace/"
  setwd(dir=space_directory)

  # Pass TRUE if you want to see progress output on some of Monocle 3's operations
  DelayedArray:::set_verbose_block_processing(TRUE)
  # Passing a higher value will make some computations faster but use more memory. 
  # Adjust with caution!
  options(DelayedArray.block.size=1000e6)
  # Set a seed to make umap and other non-deterministic steps consistent
  set.seed(seed = 42)
  
})

spatial_cds = readRDS("Submission_Data/E14_slides/RDS_intermediates/Notebook6_spatial_cds_anatomy.RDS")


##### Colors
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
colors = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


# Explore Delile et. al. dataset ------------------------------------------
# Make sure that the orginial dataset "makes sense" in time and annotation 

count_matrix_spinal = 
  Matrix::readMM("Published_Data/Delile_2019/MouseSpinalCordAtlas/analysis/spinal_count_matrix.mtx") %>%
  as(Class = "dgCMatrix")

coldata_spinal = read.table("Published_Data/Delile_2019/MouseSpinalCordAtlas/analysis/phenoData_annotated.csv",
                     header = T,
                     sep = "\t",
                     row.names = 1)

# Modify cell names to match expression data
rownames(coldata_spinal) = 
  stringr::str_replace_all(string = rownames(coldata_spinal),pattern = "-",replacement = "\\.")

coldata_spinal$Cell = 
  rownames(coldata_spinal)

rowdata = read.table("Published_Data/Delile_2019/MouseSpinalCordAtlas/analysis/spinal_rowdata.tsv",
                     header = T,
                     sep = "\t",
                     row.names = 1)

colnames(count_matrix_spinal) = 
  rownames(coldata_spinal)
rownames(count_matrix_spinal) =
  rownames(rowdata)

spinal_cord_cds = 
  new_cell_data_set(expression_data = 
                      count_matrix_spinal,
                    cell_metadata = coldata_spinal,
                    gene_metadata = rowdata)


# Clear memory occupied by the count matrix
rm(count_matrix_spinal)

spinal_cord_cds = 
  spinal_cord_cds %>%
  estimate_size_factors() %>%
  detect_genes() %>%
  preprocess_cds(num_dim = 100,
                 method = "PCA") %>%
  reduce_dimension(max_components = 2,
                   umap.fast_sgd = F) %>%
  cluster_cells(reduction_method = "UMAP")

plot_pc_variance_explained(spinal_cord_cds) +
  ggsave("Published_Data/Delile_2019/MouseSpinalCordAtlas/analysis/spinal_cord_scree_plot.pdf")


plot_cells(spinal_cord_cds,
           color_cells_by = "Type_step1",
           cell_stroke = 0,
           cell_size = 0.5) + 
  theme_void()+
  theme(legend.position = "none") +
  ggsave("Published_Data/Delile_2019/MouseSpinalCordAtlas/analysis/umap_Type_step1.pdf",
         height = 3, 
         width = 3)

plot_cells(spinal_cord_cds,
           color_cells_by = "Type_step2",
           cell_stroke = 0,
           cell_size = 0.5) + 
  theme_void()+
  theme(legend.position = "none") +
  ggsave("Published_Data/Delile_2019/MouseSpinalCordAtlas/analysis/umap_Type_step2.pdf",
         height = 3, 
         width = 3)

plot_cells(spinal_cord_cds,
           color_cells_by = "Type_step2",
           cell_stroke = 0,
           cell_size = 0.5) + 
  theme_void()+
  theme(legend.position = "none") +
  ggsave("Published_Data/Delile_2019/MouseSpinalCordAtlas/analysis/umap_Type_step2.pdf",
         height = 3, 
         width = 3)


plot_cells(spinal_cord_cds,
           color_cells_by = "timepoint",
           cell_stroke = 0,
           cell_size = 0.5) + 
  scale_color_viridis_c() +
  theme_void()+
  theme(legend.position = "none") +
  ggsave("Published_Data/Delile_2019/MouseSpinalCordAtlas/analysis/umap_Type_step2.pdf",
         height = 3, 
         width = 3)




saveRDS(object = spinal_cord_cds,
        file = "Published_Data/Delile_2019/MouseSpinalCordAtlas/analysis/spinal_cord_cds.RDS")



# Subcluster and process the neuron annotation ----------------------------
spinal_cord_cds = 
  readRDS(file = "Published_Data/Delile_2019/MouseSpinalCordAtlas/analysis/spinal_cord_cds.RDS")


plot_cells(spinal_cord_cds,
           color_cells_by = "Type_step1")

neurons = 
  colData(spinal_cord_cds) %>%
  as.data.frame() %>%
  filter(Type_step1 == "Neuron",
         timepoint == "13.5") %>%
  pull(Cell)
         
spinal_cord_cds_neuron = 
  spinal_cord_cds[,neurons]

spinal_cord_cds_neuron =
  spinal_cord_cds_neuron %>%
  estimate_size_factors() %>%
  detect_genes() %>%
  preprocess_cds(num_dim = 100,
                 method = "PCA") %>%
  reduce_dimension(max_components = 2,
                   umap.fast_sgd = F) %>%
  cluster_cells(reduction_method = "UMAP")




# Supplementary figure S25- Panel E - Delile, et al. E13.5 spinal cord neurons UMAP embedded  ----------------------------------------

colData(spinal_cord_cds_neuron)$umap1 = 
  reducedDim(spinal_cord_cds_neuron,
             type = "UMAP")[,1]

colData(spinal_cord_cds_neuron)$umap2 = 
  reducedDim(spinal_cord_cds_neuron,
             type = "UMAP")[,2]


spinal_cord_cds_neuron %>%
  colData() %>%
  as.data.frame() %>%
  ggplot() +
  geom_point(aes(x = umap1,
                 y = umap2,
                 color = Neuron_subtypes),
             stroke = 0,
             size = 0.75,
             alpha = 0.75) +
  geom_text_repel(data =
                    spinal_cord_cds_neuron %>%
                    colData() %>%
                    as.data.frame() %>%
                    filter(!is.na(Neuron_subtypes)) %>%
                    group_by(Neuron_subtypes) %>%
                    filter(n() > 50) %>%
                    mutate(med_umap1 = median(umap1),
                           med_umap2 = median(umap2)) %>%
                    filter(Neuron_subtypes != "Null_Neuron") %>%
                    dplyr::select(Neuron_subtypes,
                                   med_umap1,
                                   med_umap2) %>%
                    distinct(),
                  aes(x = med_umap1,
                      y = med_umap2,
                       label = Neuron_subtypes),
                  size = 3,
                  segment.size = 0.1) +
  scale_color_manual(values = colors)+
  theme_void() +
  theme(legend.position = "none")+
  ggsave("Figures/Figure_Components/Supplement_hox_gene_expression/delile_umap_type_step2.png",
         height = 3.5,
         width = 3)
    


# Choose HoxA cluster genes and visualize ------
marker = c("Hoxa3",
           "Hoxa4",
           "Hoxa5",
           "Hoxa7",
           "Hoxa9",
           "Hoxa10")

# Subset CDS
cds_subset = spinal_cord_cds_neuron[rowData(spinal_cord_cds_neuron)$gene_short_name %in% marker,]

# Normalize per gene counts based on size factor
marker_expr =
  (Matrix::t(exprs(cds_subset))/size_factors(cds_subset)) %>%
  as.matrix() %>%
  as.data.frame()

marker_expr[marker_expr == 0] <- NA

colnames(marker_expr) = rowData(spinal_cord_cds_neuron)[colnames(marker_expr),"gene_short_name"]

# Make wide data (matrix) into a long data format
marker_expr =
  marker_expr %>%
  rownames_to_column(var = "Cell") %>%
  gather(key = "gene",
         value = "expression",
         -Cell)


# Put HoxC cluster in order
marker_expr$gene =
  factor(marker_expr$gene,
         levels = c("Hoxa3",
                    "Hoxa4",
                    "Hoxa5",
                    "Hoxa7",
                    "Hoxa9",
                    "Hoxa10"))

coordinates =
  data.frame(Cell = colData(spinal_cord_cds_neuron)$Cell,
             umap1 = colData(spinal_cord_cds_neuron)$umap1,
             umap2 = colData(spinal_cord_cds_neuron)$umap2)

marker_expr =
  left_join(marker_expr,
            coordinates,
            by = "Cell")


# Supplementary figure S25- Panel F - Delile, et al. E13.5 spinal cord neurons HoxA cluster expression ----------------------------------------

ggplot() +
  geom_point(data = marker_expr %>%
               filter(is.na(expression)),
             aes(x = umap1,
                 y = umap2),
             color = "grey80",
             stroke = 0,
             size =.65) +
  geom_point(data = marker_expr %>%
               filter(!is.na(expression)),
             aes(x = umap1,
                 y = umap2,
                 color = expression),
             stroke = 0,
             size =.75) +
  facet_wrap(~gene,
             nrow = 2) +
  scale_color_viridis_c(option = "plasma",
                        "") +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "right",
        legend.key.height = unit(0.15,"in"),
        legend.key.width = unit(0.075,"in"),
        strip.text.x = element_text(size = 10)) +
  ggsave("Figures/Figure_Components/Supplement_hox_gene_expression/Delile_hox_umap.png",
         width = 4.5,
         height = 3)
