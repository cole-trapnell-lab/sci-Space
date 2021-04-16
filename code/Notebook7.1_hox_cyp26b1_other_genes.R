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
  library(RColorBrewer)
  library(sp)
  library(sf)
  library(monocle3)
  library(ggrepel)
  library(ggridges)
  library(viridis)
  library(purrr)
  library(modelr)
  
  space_directory = "/Users/maryregier/Desktop/Full_sciSpace"
  setwd(dir=space_directory)
  source("Submission_Data/bin/cell_cycle.R")
  
  # Pass TRUE if you want to see progress output on some of Monocle 3's operations
  DelayedArray:::set_verbose_block_processing(TRUE)
  # Passing a higher value will make some computations faster but use more memory. 
  options(DelayedArray.block.size=1000e6)
  set.seed(42)
})

# Load datasets ----

spatial_cds = 
  readRDS(file = "Submission_Data/E14_slides/RDS_intermediates/Notebook6_spatial_cds_anatomy.RDS")

all_image_data = 
  readRDS(file = "Submission_Data/E14_slides/RDS_intermediates/Notebook0_1_images_and_transformations.RDS")




### Generate Connective Tissue Progenitor umap subclusters ----

#Subset cells annotated as connective tissue progenitors
cells_to_pull =
  colData(spatial_cds) %>%
  as.data.frame() %>%
  filter(final_cluster_label %in% c("Connective Tissue Progenitors")) %>%
  pull(Cell)

ctp_cds = spatial_cds[,cells_to_pull]

# Process Connective Tissue Progenitors CDS
ctp_cds = 
  ctp_cds %>%
  estimate_size_factors() %>%
  preprocess_cds() %>%
  reduce_dimension(max_components = 2,
                   reduction_method = "UMAP",
                   umap.fast_sgd = F) %>%
  cluster_cells(resolution = 1e-4,
                random_seed = 42)

ctp_cds = detect_genes(ctp_cds)

plot_cells(ctp_cds,
           cell_size = 0.75)

# Add UMAP coordinates to metadata
colData(ctp_cds)$umap1 = 
  reducedDim(ctp_cds,
             type = "UMAP")[,1]

colData(ctp_cds)$umap2 = 
  reducedDim(ctp_cds,
             type = "UMAP")[,2]

# Add UMAP clusters to metadata
colData(ctp_cds)$umap_cluster = ctp_cds@clusters[["UMAP"]][["clusters"]]




slide_names = c("slide_1D","slide_1E","slide_1F","slide_1G","slide_2G","slide_2H","slide_3D","slide_3F","slide_3G","slide_3H","slide_4A","slide_4C","slide_4D","slide_4E")

# set color palette
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
colors = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
colors[7]= "slateblue4"

# Supp fig 25 - Panel A - subclusters in UMAP space --------
ggplot(colData(ctp_cds) %>%
         as.data.frame()) +
  geom_point(aes(x = umap1,
                 y = umap2,
                 color = umap_cluster),
             stroke = 0,
             size =1) +
  monocle3:::monocle_theme_opts() +
  scale_color_manual(values = colors) + 
  theme_void() + 
  theme(legend.position = "none") + 
  ggsave("Figures/Figure_Components/Supplement_connective_tissue_progenitors/ctp_umap.png",
         height = 7,
         width = 7) 

# Map select subclusters to each slide 
# Supp fig 25 - Panel B components - subclusters in real space ----

# Select subclusters and corresponding colors
subclusters = c(1,7,8)
subset_colors = colors[subclusters]

for (ind in 1:14) {
  plot_slide = colData(ctp_cds) %>%
    as.data.frame() %>%
    filter(max_slide_id == slide_names[ind]) %>%
    filter(umap_cluster %in% subclusters)
  plotstring = paste("Figures/Figure_Components/Supplement_connective_tissue_progenitors/ctp_",slide_names[ind],"clusters.png",sep = "")
  
  ggplot() +
    geom_sf(data = all_image_data$slide_polygon[[ind]],
            fill = "grey90",
            color = "grey70",
            size = .15) +
    geom_jitter(data = plot_slide,
                aes(x = coords.x1 ,
                    y = coords.x2,
                    color = umap_cluster),
                stroke = 0,
                size =1,
                width = 25,
                height = 25) +
    monocle3:::monocle_theme_opts() +
    scale_color_manual(values = subset_colors) + 
    theme_void() +
    facet_wrap(~umap_cluster,
               ncol = 1) + 
    theme_void() +
    theme(strip.text.x = element_blank()) +
    theme(legend.position = "none") +
    ggsave(plotstring,         
           height = 4,
           width = 1.5)
}

# Also create reversed order plots for opposing orientation to simplify figure composition
rev_subclusters = c(8,7,1)

rev_subset_colors = colors[rev_subclusters]

for (ind in 1:14) {
  
  plot_slide = colData(ctp_cds) %>%
    as.data.frame() %>%
    filter(max_slide_id == slide_names[ind]) %>%
    filter(umap_cluster %in% rev_subclusters)
  plot_slide$umap_cluster =
    factor(plot_slide$umap_cluster,
           levels = rev_subclusters)
  
  plotstring = paste("Figures/Figure_Components/Supplement_connective_tissue_progenitors/ctp_",slide_names[ind],"clusters_rev.png",sep = "")
  
  ggplot() +
    geom_sf(data = all_image_data$slide_polygon[[ind]],
            fill = "grey90",
            color = "grey70",
            size = .15) +
    geom_jitter(data = plot_slide,
                aes(x = coords.x1 ,
                    y = coords.x2,
                    color = umap_cluster),
                stroke = 0,
                size =1,
                width = 25,
                height = 25) +
    monocle3:::monocle_theme_opts() +
    scale_color_manual(values = rev_subset_colors) + 
    theme_void() +
    facet_wrap(~umap_cluster,
               ncol = 1) + 
    theme_void() +
    theme(strip.text.x = element_blank())+
    theme(legend.position = "none") +
    ggsave(plotstring,         
           height = 4,
           width = 1.5)
}









### Annotate neurons using Delile and La Manno datasets -------

# Add LaManno types to spatial_cds ----
lamanno_type = colData(spatial_cds) %>%
  as.data.frame() %>%
  dplyr::select(Cell,lamanno_Punchcard) %>%
  separate(col = lamanno_Punchcard,into = c("Region_Stage","Kept","Rgl","Type"),sep = "_")



lamanno_type$Type[lamanno_type$Type == "GABA" | lamanno_type$Type == "GABAGly" | lamanno_type$Type == "Inhibitory"]  <- "Inhibitory - Brain" 
lamanno_type$Type[lamanno_type$Type == "NonGABA" | lamanno_type$Type == "NonGABAGly" | lamanno_type$Type == "Excitatory"]  <- "Excitatory - Brain" 
colData(spatial_cds)$lamanno_type = lamanno_type[,5]

# Subcluster and process the Delile neuron annotation -----
spinal_cord_cds = 
  readRDS(file = "Published_Data/Delile_2019/MouseSpinalCordAtlas/analysis/spinal_cord_cds.RDS")



e13.5_cells = 
  colData(spinal_cord_cds) %>%
  as.data.frame() %>% 
  filter(timepoint %in% c("13.5")) %>%
  pull(Cell)

spinal_cord_cds_neuron = 
  spinal_cord_cds[,e13.5_cells]


# Join the Delile and sciSpace spinal cord datasets and transfer labels ---------------------------------------------------

rowData(spatial_cds)$id_2 = 
  stringr::str_split_fixed(rowData(spatial_cds)$id,
                           "\\.",
                           2)[,1]

intersecting_genes = 
  intersect(rowData(spinal_cord_cds)$id %>%
              as.character(),
            rowData(spatial_cds)$id_2)

rownames(spatial_cds) =rowData(spatial_cds)$id_2


cds_space = spatial_cds[intersecting_genes,colData(spatial_cds)$final_cluster_label == "Neuron"]

spinal_cds = spinal_cord_cds_neuron[intersecting_genes,]

identical(rownames(cds_space) %>%
            as.character(),
          rownames(spinal_cds) %>%
            as.character())


joint_coldata = 
  rbind(colData(spinal_cds) %>%
          as.data.frame() %>%
          dplyr::select(Cell,
                        annotation = Type_step2) %>%
          mutate(sample = "Spinal Atlas"),
        colData(cds_space) %>%
          as.data.frame() %>%
          dplyr::select(Cell,
                        annotation = nn_label) %>%
          mutate(sample = "Spatial Atlas")
  )

rownames(joint_coldata) = 
  joint_coldata$Cell


joint_cds =
  new_cell_data_set(expression_data = cbind(counts(spinal_cds), 
                                            counts(cds_space)),
                    cell_metadata = joint_coldata,
                    gene_metadata = rowData(spinal_cds) %>% 
                      as.data.frame())

# Run Seurat CCA + Data Integration + transfer labels 
library(scater)
library(Seurat)

count_mat = assay(joint_cds)
coldata_df = colData(joint_cds) %>% as.data.frame()
coldata_df$Cell = rownames(coldata_df)

cds_seurat = 
  CreateSeuratObject(counts = count_mat,
                     project = "spine",
                     assay = "RNA",
                     meta.data = coldata_df)

# split dataset by experiment
cds.list <- SplitObject(cds_seurat, split.by = "sample")

# normalize data and find variable genes
for (i in 1:length(cds.list)) {
  cds.list[[i]] <- NormalizeData(cds.list[[i]], verbose = FALSE)
  cds.list[[i]] <- FindVariableFeatures(cds.list[[i]], selection.method = "vst", 
                                        nfeatures = 2000, verbose = FALSE)
}

# transfer annotations from Spinal Altas to Spatial Atlas
cds.query <- cds.list[["Spatial Atlas"]]
cds.anchors <- FindTransferAnchors(reference = cds.list[["Spinal Atlas"]], query = cds.query, 
                                        dims = 1:30)
predictions <- TransferData(anchorset = cds.anchors, refdata = cds.list[["Spinal Atlas"]]$annotation, 
                            dims = 1:30)
Cell = rownames(predictions)
delile_annot = cbind(Cell,predictions[1])

# Reclassify Delile subtypes more broadly as excitatory/inhibitory
delile_annot$delile_lamanno_type <- 
  plyr::revalue((delile_annot$predicted.id),
                c("dp4" = 'Progenitors', 
                  "V0" = 'Inhibitory',
                  "V2a" = 'Inhibitory',
                  "V2b" = 'Inhibitory',
                  "Null_Progenitor" = 'Progenitors', 
                  "dl4" = 'Inhibitory',
                  "Null_Neuron" = 'Unannotated',
                  "Outliers" = 'Unannotated',
                  "dl3" = 'Excitatory',
                  "dl1" = 'Excitatory',
                  "dl2" = 'Excitatory',
                  "dl5" = 'Excitatory',
                  "dl6" = 'Excitatory',
                  "V1" = 'Inhibitory',
                  "V3" = 'Inhibitory',
                  "p1" = 'Progenitors',
                  "p3" = 'Progenitors',
                  "p0" = 'Progenitors',
                  "p2" = 'Progenitors',
                  "pMN" = 'Progenitors',
                  "dp6" = 'Progenitors',
                  "dp3" = 'Progenitors',
                  "dp5" = 'Progenitors',
                  "dp2" = 'Progenitors',
                  "dp1" = 'Progenitors'))

delile_annot = delile_annot %>%
  select(-predicted.id)

# Merge Delile labeled neurons with non-neurons
non_neuron = as.data.frame(colData(spatial_cds))%>%
  filter(final_cluster_label != "Neuron")
non_neuron$delile_lamanno_type = NA
delile_labels = as.data.frame(colData(spatial_cds)) %>%
  filter(final_cluster_label == "Neuron")
delile_labels = left_join(delile_labels,delile_annot)
full_annot = rbind(delile_labels,non_neuron)
newcoldata = left_join(as.data.frame(colData(spatial_cds)),full_annot)

# Use LaManno types for the brain anatomical annotation
for (i in 1:length(newcoldata[,1])){
  if (newcoldata[i,41]=="Cortex"){
    newcoldata[i,43] = newcoldata[i,42]
  } 
}

# Reclassify LaManno subtypes more broadly as excitatory/inhibitory
newcoldata$delile_lamanno_type = newcoldata$delile_lamanno_type %>% replace_na("Unannotated")

# Append Delile/LaManno types to spatial_cds
colData(spatial_cds)$delile_lamanno_type = newcoldata$delile_lamanno_type


### Analyze Slide 1 Neurons -------

#Subset cells of interest
cells_to_pull =
  colData(spatial_cds) %>%
  as.data.frame() %>%
  filter(max_slide_id == "slide_1D") %>%
  filter(final_cluster_label %in% c("Neuron")) %>%
  pull(Cell)

cds = spatial_cds[,cells_to_pull]


# Process neuron CDS 
cds = 
  cds %>%
  estimate_size_factors() %>%
  preprocess_cds() %>%
  reduce_dimension(max_components = 2,
                   reduction_method = "UMAP",
                   umap.fast_sgd = F) %>%
  cluster_cells(resolution = 1e-3,
                random_seed = 42)

cds = detect_genes(cds)

plot_cells(cds,
           cell_size = 0.75)

# Add UMAP coordinates to metadata
colData(cds)$umap1 = 
  reducedDim(cds,
             type = "UMAP")[,1]

colData(cds)$umap2 = 
  reducedDim(cds,
             type = "UMAP")[,2]

# Add UMAP clusters to metadata
colData(cds)$umap_cluster = cds@clusters[["UMAP"]][["clusters"]]

# Perform graph_test on UMAP coordinates this function 
# calculates autocorrelation on the nearest neighbor graph
pr_graph_test_UMAP =
  graph_test(cds, 
             neighbor_graph="knn",
             reduction_method = "UMAP",
             cores = 8)


pr_deg_ids_UMAP = 
  row.names(subset(pr_graph_test_UMAP, 
                   q_value < 0.05))

# Move UMAP coordinates to the "tSNE" slot in reduced dimensions
reducedDim(x = cds,
           type = "tSNE") =
  reducedDim(cds, 
             type = "UMAP")

# Move UMAP clusters to the "tSNE" slot in "clusters"
cds@clusters[["tSNE"]] = 
  cds@clusters[["UMAP"]]

# Check again to make sure things are correct
identical(rownames(reducedDim(x = cds,
                              type = "tSNE")),
          colnames(cds) %>% as.character())

# Move spatial coordinates to the UMAP slot
reducedDim(x = cds,
           type = "UMAP") <-
  matrix(cbind(colData(cds)$coords.x2,
               -colData(cds)$coords.x1), 
         ncol=2)


cds = cluster_cells(cds,
                    resolution = 1e-3,
                    random_seed = 42,
                    reduction_method = "UMAP")

# Perform graph_test on spatial coordinates this function 
# calculates autocorrelation on the nearest neighbor graph
pr_graph_test_spatial =
  graph_test(cds, 
             neighbor_graph="knn",
             reduction_method = "UMAP",
             cores = 8)

pr_deg_ids_spatial = 
  row.names(subset(pr_graph_test_spatial, 
                   q_value < 0.05))

# Get the genes in the union of the two sets
genes_of_interest = 
  intersect(pr_deg_ids_UMAP,
            pr_deg_ids_spatial)


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
  left_join(rowData(cds) %>%
              as.data.frame(),
            by = c("id","gene_short_name"))

# Supp fig 25 — Panel A - neurons in real space ------------------------------------------------------
ggplot(colData(cds) %>% 
         as.data.frame()) +
  geom_sf(data = all_image_data$slide_polygon[[1]],
          fill = "grey80",
          alpha = 0.25) +
  geom_jitter(data = colData(cds) %>% 
                as.data.frame() %>%
                filter(delile_lamanno_type == "Unannotated"),
              aes(x = coords.x1,
                  y = coords.x2),
              color = "gray",
              stroke = 0,
              size =1,
              width = 25,
              height = 25) +
  geom_jitter(data = colData(cds) %>% 
                as.data.frame() %>%
                filter(delile_lamanno_type != "Unannotated"),
              aes(x = coords.x1,
                  y = coords.x2,
                  color = delile_lamanno_type),
              stroke = 0,
              size =1,
              width = 25,
              height = 25) +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Supplement_hox_gene_expression/annot_neuron_spatial.pdf",
         height = 2,
         width = 2.25)


# Supp fig 25 — Panel B - neurons in umap space ------------------------------------------------------
ggplot(colData(cds) %>%
         as.data.frame()) +
  geom_point(data = colData(cds) %>% 
                  as.data.frame() %>%
                  filter(delile_lamanno_type == "Unannotated"),
             aes(x = umap1,
                 y = umap2),
             color = "gray",
             stroke = 0,
             size =1) +
    geom_point(data = colData(cds) %>% 
                 as.data.frame() %>%
                 filter(delile_lamanno_type != "Unannotated"),
               aes(x = umap1,
                   y = umap2,
                   color = delile_lamanno_type),
               stroke = 0,
               size =1) +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Supplement_hox_gene_expression/annot_neuron_umap.pdf",
         height = 2,
         width = 2)



# Figure 3 — Panel B — hox gene autocorrelation space vs umap ---------------------------------


hox_genes =
  rowData(cds)$gene_short_name[grepl(x = rowData(cds)$gene_short_name,
                                     pattern = "Hox[a-e][0-9]+$",
                                     ignore.case = F)]
hox_genes
ggplot(spatial_and_umap_stats %>%
         filter(log(morans_I_spatial) > -20)) +
  geom_abline(color = "blue",
              size = .25,
              alpha = 0.3) +
  geom_point(aes(x = log10(morans_I_umap),
                 y = log10(morans_I_spatial)),
             stroke = 0,
             size = 0.5,
             color = "grey80") +
  geom_point(data =spatial_and_umap_stats %>%
               filter(gene_short_name %in% hox_genes),
             aes(x = log10(morans_I_umap),
                 y = log10(morans_I_spatial)),
             stroke = 0,
             size = 1,
             color = "#117733") +
  geom_text_repel(data =spatial_and_umap_stats %>%
                    filter(gene_short_name %in% hox_genes),
                  aes(x = log10(morans_I_umap),
                      y = log10(morans_I_spatial),
                      label = gene_short_name),
                  color = "black",
                  size = 2,box.padding = 0.25,
                  segment.size = 0.25) +
  
  monocle3:::monocle_theme_opts() +
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8)) +
  xlim(-2.5,NA) +
  ylim(-2,NA) +
  xlab("log(Morans I) [UMAP]") +
  ylab("log(Morans I) [Spatial]") +
  ggsave("Figures/Figure_Components/Figure3/neuron_moransI_Hox_comparison.pdf",
         height = 2.35,
         width = 2.35)


# Figure 3 — Panel C - hox spatial morans I boxplot comparison------------------------------------------------------

rowData(cds)$mean_expr = rowSums(counts(cds))/dim(cds)[2]

mean_expr_hox = 
  rowData(cds) %>%
  as.data.frame() %>%
  filter(gene_short_name %in% hox_genes) %>%
  pull(mean_expr) %>%
  mean()


expression_matched_genes =
  rowData(cds) %>%
  as.data.frame() %>%
  filter(!(gene_short_name %in% hox_genes)) %>%
  filter(mean_expr < mean_expr_hox + 0.025,
         mean_expr > mean_expr_hox - 0.025) %>%
  pull(gene_short_name)


spatial_and_umap_stats %>%
  filter(gene_short_name %in% c(expression_matched_genes,hox_genes)) %>%
  mutate(hox = ifelse(gene_short_name %in% hox_genes,"Hox Gene", "Other")) %>%
  filter(log(morans_I_spatial) > -20) %>%
  ggplot() +
  geom_boxplot(aes(y = log10(morans_I_spatial),
                   x = hox,
                   fill = hox),
               size = 0.25,
               color = "black",
               outlier.size = 0.75,
               outlier.stroke = 0) +
  scale_fill_manual(values = c("#117733",
                               "grey90")) +
  monocle3:::monocle_theme_opts() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 6,
                                   angle = 45,
                                   hjust = 1),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8)) +
  ylab("log10(Morans I) [Spatial]") +
  ggsave("Figures/Figure_Components/Figure3/moransI_boxplot_comparison.pdf",
         height = 2.5,
         width = 1)


hox_binary =
  spatial_and_umap_stats %>%
  filter(gene_short_name %in% c(expression_matched_genes,hox_genes)) %>%
  mutate(hox = gene_short_name %in% hox_genes) %>%
  pull(hox)

# Statistical test between Hox genes and non-Hox genes
t.test(x = 
         spatial_and_umap_stats %>%
         filter(gene_short_name %in% c(hox_genes)) %>%
         pull(morans_I_spatial),
       y = spatial_and_umap_stats %>%
         filter(gene_short_name %in% c(expression_matched_genes)) %>%
         pull(morans_I_spatial))


# Supp fig 25 - Panel C - hox UMAP morans I boxplot comparison ----
spatial_and_umap_stats %>%
  filter(gene_short_name %in% c(expression_matched_genes,hox_genes)) %>%
  mutate(hox = ifelse(gene_short_name %in% hox_genes,"Hox Gene", "Other")) %>%
  #filter(log(morans_I_spatial) > -20) %>%
  ggplot() +
  geom_boxplot(aes(y = log10(morans_I_umap),
                   x = hox,
                   fill = hox),
               size = 0.25,
               color = "black",
               outlier.size = 0.75,
               outlier.stroke = 0) +
  scale_fill_manual(values = c("#117733",
                               "grey90")) +
  monocle3:::monocle_theme_opts() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 6,
                                   angle = 45,
                                   hjust = 1),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8)) +
  ylab("log10(Morans I) [UMAP]") +
  ggsave("Figures/Figure_Components/Supplement_hox_gene_expression/moransI_umap_boxplot_comparison.pdf",
         height = 2.5,
         width = 1)


# Statistical test between Hox genes and non-Hox genes
t.test(x = 
         spatial_and_umap_stats %>%
         filter(gene_short_name %in% c(hox_genes)) %>%
         pull(morans_I_umap),
       y = spatial_and_umap_stats %>%
         filter(gene_short_name %in% c(expression_matched_genes)) %>%
         pull(morans_I_umap))


# Figure 3 — Panel D — Hox gene spatial plot ------------------------------

# Choose HoxA cluster genes to visualize

marker = c("Hoxa3",
           "Hoxa4",
           "Hoxa5",
           "Hoxa7",
           "Hoxa9",
           "Hoxa10")

# Subset CDS
cds_subset = cds[rowData(cds)$gene_short_name %in% marker,]

# Normalize per gene counts based on size factor
marker_expr =
  (Matrix::t(exprs(cds_subset))/size_factors(cds_subset)) %>%
  as.matrix() %>%
  as.data.frame()

marker_expr[marker_expr == 0] <- NA

colnames(marker_expr) = rowData(cds)[colnames(marker_expr),"gene_short_name"]

# Make wide data (matrix) into a long data format
marker_expr =
  marker_expr %>%
  rownames_to_column(var = "Cell") %>%
  gather(key = "gene",
         value = "expression",
         -Cell)

# Add spatial coordinates to normalized count values for cells
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

# Put HoxA cluster in order
marker_expr$gene =
  factor(marker_expr$gene,
         levels = c("Hoxa3",
                    "Hoxa4",
                    "Hoxa5",
                    "Hoxa7",
                    "Hoxa9",
                    "Hoxa10"))

# Plot result
ggplot() +
  geom_sf(data = all_image_data$slide_polygon[[1]] * rbind(c(0, -1),c(1,0)),
          fill = "grey90",
          color = "grey70",
          size = .15) +
  geom_jitter(data = marker_expr %>%
                filter(!is.na(expression)),
              aes(x = spatial_2,
                  y = -spatial_1,
                  color = expression),
              stroke = 0,
              size =.55,
              width = 25,
              height = 25) +
  facet_wrap(~gene,
             ncol = 3) +
  scale_color_viridis_c(option = "plasma") +
  scale_fill_gradient(low = "white", high = "grey30") +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Figure3/neuron_hoxA_spatial.pdf",
         height = 2,
         width = 3)




# Supp fig 25 - Panel D - Hox gene UMAP plot----
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
             ncol = 3) +
  scale_color_viridis_c(option = "plasma",
                        "") +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "right",
        legend.key.height = unit(0.15,"in"),
        legend.key.width = unit(0.075,"in"),
        strip.text.x = element_text(size = 10)) +
  ggsave("Figures/Figure_Components/Supplement_hox_gene_expression/neuron_hoxA_umap.pdf",
         width = 3,
         height = 2)



# Pick out other genes from each category to highlight ----

# Some genes are randomly chosen and some are selected
other_genes =
  c(spatial_and_umap_stats %>% 
      filter(q_value_umap < 0.01) %>% 
      mutate(fold_diff = morans_I_spatial/morans_I_umap) %>%
      arrange(desc(fold_diff)) %>%
      dplyr::select(gene_short_name,fold_diff) %>%
      as.data.frame() %>%
      tail(n = 100) %>%
      sample_n(size = 6) %>% 
      pull(gene_short_name))


equal_genes =
  c(spatial_and_umap_stats %>% 
      filter(q_value_umap < 0.01,
             q_value_spatial < 0.01) %>% 
      mutate(fold_diff = morans_I_spatial/morans_I_umap) %>%
      filter(fold_diff > 0.66,
             fold_diff < 1.66,
             morans_I_umap > 0.15) %>%
      dplyr::select(gene_short_name,fold_diff) %>%
      as.data.frame() %>%
      sample_n(size = 6) %>% 
      pull(gene_short_name))

spatial_genes_to_highlight =
  c("Calb1",
    "Esrrb",
    "L3mbtl1",
    "Fmo6",
    "Jam2",
    "Col9a1",
    "Cyp26b1")

# Figure 3 — Panel E - other genes morans I comparison ------------------------------------------------------

ggplot(spatial_and_umap_stats %>%
         filter(log10(morans_I_spatial) > -4,
                q_value_spatial < 0.01)) +
  geom_abline(color = "blue",
              size = .25,
              alpha = 0.25) +
  geom_point(aes(x = log10(morans_I_umap),
                 y = log10(morans_I_spatial)),
             stroke = 0,
             size = 0.5,
             color = "grey80") +
  
  geom_point(data =spatial_and_umap_stats %>%
               filter(gene_short_name %in% other_genes),
             aes(x = log10(morans_I_umap),
                 y = log10(morans_I_spatial)),
             stroke = 0,
             size = 1,
             color = "blue") +
  geom_point(data =spatial_and_umap_stats %>%
               filter(gene_short_name %in% spatial_genes_to_highlight),
             aes(x = log10(morans_I_umap),
                 y = log10(morans_I_spatial)),
             stroke = 0,
             size = 1,
             color = "red") +
  geom_point(data =spatial_and_umap_stats %>%
               filter(gene_short_name %in% equal_genes),
             aes(x = log10(morans_I_umap),
                 y = log10(morans_I_spatial)),
             stroke = 0,
             size = 1,
             color = "black") +
  geom_text_repel(data =spatial_and_umap_stats %>%
                    filter(gene_short_name %in% other_genes|
                             gene_short_name %in% spatial_genes_to_highlight|
                             gene_short_name %in% equal_genes),
                  aes(x = log10(morans_I_umap),
                      y = log10(morans_I_spatial),
                      label = gene_short_name),
                  color = "black",
                  size = 2.5) +
  monocle3:::monocle_theme_opts() + 
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8)) +
  xlab("log(Morans I) [UMAP]") +
  ylab("log(Morans I) [Spatial]") +
  xlim(-2.5,NA) +
  ylim(-2,NA) +
  ggsave("Figures/Figure_Components/Figure3/neuron_other_gene_morans_I.pdf",
         height = 2.35,
         width = 2.35) 




# Figure 3 — Panel F - other genes spatial ------------------------------------------------------

spatial_genes_to_plot =
  c("Calb1",
    "Esrrb",
    "L3mbtl1",
    "Fmo6",
    "Jam2",
    "Col9a1")

marker = spatial_genes_to_plot

# Subset CDS
cds_subset = cds[rowData(cds)$gene_short_name %in% marker,]

# Normalize per gene counts based on size factor
marker_expr =
  (Matrix::t(exprs(cds_subset))/size_factors(cds_subset)) %>%
  as.matrix() %>%
  as.data.frame()

marker_expr[marker_expr == 0] <- NA

colnames(marker_expr) = rowData(cds)[colnames(marker_expr),"gene_short_name"]

# Make wide data (matrix) into a long data format
marker_expr =
  marker_expr %>%
  rownames_to_column(var = "Cell") %>%
  gather(key = "gene",
         value = "expression",
         -Cell)

# Add spatial coordinates to normalized count values for cells
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

# Plot result
ggplot() +
  geom_sf(data = all_image_data$slide_polygon[[1]] * rbind(c(0, -1),c(1,0)),
          fill = "grey90",
          color = "grey70",
          size = .15) +
  geom_jitter(data = marker_expr %>%
                filter(expression > 6),
              aes(x = spatial_2,
                  y = -spatial_1),
                  color = "yellow",
              stroke = 0,
              size =.55,
              width = 25,
              height = 25) +
  geom_jitter(data = marker_expr %>%
                filter(!is.na(expression)) %>%
                filter(expression <= 6),
              aes(x = spatial_2,
                  y = -spatial_1,
                  color = expression),
              stroke = 0,
              size =.55,
              width = 25,
              height = 25) +
  facet_wrap(~gene,
             ncol = 3) +
  scale_color_viridis_c(option = "plasma") +
  scale_fill_gradient(low = "white", high = "grey30") +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Figure3/neuron_other_spatial.pdf",
         height = 2,
         width = 3)


# Supp fig 25 - Panel D - other genes UMAP ----

ggplot() +
  geom_point(data = marker_expr %>%
               filter(is.na(expression)),
             aes(x = umap1,
                 y = umap2),
             color = "grey80",
             stroke = 0,
             size =.65) +
  geom_point(data = marker_expr %>%
               filter(expression > 6),
             aes(x = umap1,
                 y = umap2),
                 color = "yellow",
             stroke = 0,
             size =.75) +
  geom_point(data = marker_expr %>%
               filter(expression <= 6),
             aes(x = umap1,
                 y = umap2,
                color = expression),
             stroke = 0,
             size =.75) +
  facet_wrap(~gene,
             ncol = 6) +
  scale_color_viridis_c(option = "plasma",
                        "") +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "right",
        legend.key.height = unit(0.15,"in"),
        legend.key.width = unit(0.075,"in"),
        strip.text.x = element_text(size = 10)) +
  ggsave("Figures/Figure_Components/Supplement_other_spatial_&RNAscope/neuron_other_umap.pdf",
         width = 6,
         height = 1)



# Supp fig 25 — Panel A&B - Cyp26b1 and Hoxa10 spatial and UMAP ------------------------------------------------------


spatial_genes_to_plot =
  c("Cyp26b1", "Hoxa10")


marker = spatial_genes_to_plot

# Subset CDS
cds_subset = cds[rowData(cds)$gene_short_name %in% marker,]

# Normalize per gene counts based on size factor
marker_expr =
  (Matrix::t(exprs(cds_subset))/size_factors(cds_subset)) %>%
  as.matrix() %>%
  as.data.frame()

marker_expr[marker_expr == 0] <- NA

colnames(marker_expr) = rowData(cds)[colnames(marker_expr),"gene_short_name"]

# Make wide data (matrix) into a long data format
marker_expr =
  marker_expr %>%
  rownames_to_column(var = "Cell") %>%
  gather(key = "gene",
         value = "expression",
         -Cell)

# Add spatial coordinates to normalized count values for cells
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


# Plot result
ggplot() +
  geom_sf(data = all_image_data$slide_polygon[[1]] * rbind(c(0, -1),c(1,0)),
          fill = "grey80",
          alpha = 0.25) +
  geom_jitter(data = marker_expr %>%
                filter(!is.na(expression)),
              aes(x = spatial_2,
                  y = -spatial_1,
                  color = expression),
              stroke = 0,
              size =1,
              width = 25,
              height = 25) +
  facet_wrap(~gene,
             ncol = 2) +
  scale_color_viridis_c(option = "plasma") +
  scale_fill_gradient(low = "white", high = "grey30") +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Supplement_other_spatial_&RNAscope/neuron_cyp26b1_hoxa10_spatial.pdf",
         height = 2,
         width = 3)


ggplot() +
  geom_point(data = marker_expr %>%
               filter(is.na(expression)),
             aes(x = umap1,
                 y = umap2),
             color = "grey80",
             stroke = 0,
             size =.8) +
  geom_point(data = marker_expr %>%
               filter(!is.na(expression)),
             aes(x = umap1,
                 y = umap2,
                 color = expression),
             stroke = 0,
             size =1) +
  facet_wrap(~gene,
             nrow = 1) +
  scale_color_viridis_c(option = "plasma",
                        "") +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "right",
        legend.key.height = unit(0.15,"in"),
        legend.key.width = unit(0.075,"in"),
        strip.text.x = element_text(size = 10)) +
  ggsave("Figures/Figure_Components/Supplement_other_spatial_&RNAscope/neuron_cyp26b1_hoxa10_umap.pdf",
         width = 3,
         height = 2)







### Compare Slide 14 Neurons ------

#Subset cells of interest
cells_to_pull =
  colData(spatial_cds) %>%
  as.data.frame() %>%
  filter(max_slide_id == "slide_4E") %>%
  filter(final_cluster_label %in% c("Neuron")) %>%
  pull(Cell)

cds_4E = spatial_cds[,cells_to_pull]


# Process neuron CDS 
cds_4E = 
  cds_4E %>%
  estimate_size_factors() %>%
  preprocess_cds() %>%
  reduce_dimension(max_components = 2,
                   reduction_method = "UMAP",
                   umap.fast_sgd = F) %>%
  cluster_cells(resolution = 1e-3,
                random_seed = 42)

cds_4E = detect_genes(cds_4E)

plot_cells(cds,
           cell_size = 0.75)

# Add UMAP coordinates to metadata
colData(cds_4E)$umap1 = 
  reducedDim(cds_4E,
             type = "UMAP")[,1]

colData(cds_4E)$umap2 = 
  reducedDim(cds_4E,
             type = "UMAP")[,2]

# Add UMAP clusters to metadata
colData(cds_4E)$umap_cluster = cds_4E@clusters[["UMAP"]][["clusters"]]

# Perform graph_test on UMAP coordinates this function 
# calculates autocorrelation on the nearest neighbor graph
pr_graph_test_UMAP =
  graph_test(cds_4E, 
             neighbor_graph="knn",
             reduction_method = "UMAP",
             cores = 8)


pr_deg_ids_UMAP = 
  row.names(subset(pr_graph_test_UMAP, 
                   q_value < 0.05))

# Move UMAP coordinates to the "tSNE" slot in reduced dimensions
reducedDim(x = cds_4E,
           type = "tSNE") =
  reducedDim(cds_4E, 
             type = "UMAP")

# Move UMAP clusters to the "tSNE" slot in "clusters"
cds_4E@clusters[["tSNE"]] = 
  cds_4E@clusters[["UMAP"]]

# Check again to make sure things are correct
identical(rownames(reducedDim(x = cds_4E,
                              type = "tSNE")),
          colnames(cds_4E) %>% as.character())

# Move spatial coordinates to the UMAP slot
reducedDim(x = cds_4E,
           type = "UMAP") <-
  matrix(cbind(colData(cds_4E)$coords.x2,
               -colData(cds_4E)$coords.x1), 
         ncol=2)


cds_4E = cluster_cells(cds_4E,
                    resolution = 1e-3,
                    random_seed = 42,
                    reduction_method = "UMAP")

# Perform graph_test on spatial coordinates this function 
# calculates autocorrelation on the nearest neighbor graph
pr_graph_test_spatial =
  graph_test(cds_4E, 
             neighbor_graph="knn",
             reduction_method = "UMAP",
             cores = 8)

pr_deg_ids_spatial = 
  row.names(subset(pr_graph_test_spatial, 
                   q_value < 0.05))

# Get the genes in the union of the two sets
genes_of_interest = 
  intersect(pr_deg_ids_UMAP,
            pr_deg_ids_spatial)


spatial_and_umap_stats_4E =
  full_join( pr_graph_test_UMAP,
             pr_graph_test_spatial,
             by = c("id",
                    "gene_short_name")) %>%
  filter(status.x == "OK",
         status.y == "OK")
colnames(spatial_and_umap_stats_4E) = 
  stringr::str_replace(colnames(spatial_and_umap_stats_4E),
                       pattern = ".x",
                       replacement = "_umap") %>%
  stringr::str_replace(pattern = ".y",
                       replacement = "_spatial")


spatial_and_umap_stats_4E = 
  spatial_and_umap_stats_4E %>%
  left_join(rowData(cds_4E) %>%
              as.data.frame(),
            by = c("id","gene_short_name"))

# Supp fig 26 - Panel E - Slide 14 other genes morans I comparison -----


ggplot(spatial_and_umap_stats_4E %>%
         filter(log10(morans_I_spatial) > -4,
                q_value_spatial < 0.01)) +
  geom_abline(color = "blue",
              size = .25,
              alpha = 0.25) +
  geom_point(aes(x = log10(morans_I_umap),
                 y = log10(morans_I_spatial)),
             stroke = 0,
             size = 0.5,
             color = "grey80") +
  
  geom_point(data = spatial_and_umap_stats_4E %>%
               filter(gene_short_name %in% other_genes),
             aes(x = log10(morans_I_umap),
                 y = log10(morans_I_spatial)),
             stroke = 0,
             size = 1,
             color = "blue") +
  geom_point(data = spatial_and_umap_stats_4E %>%
               filter(gene_short_name %in% spatial_genes_to_highlight),
             aes(x = log10(morans_I_umap),
                 y = log10(morans_I_spatial)),
             stroke = 0,
             size = 1,
             color = "red") +
  geom_point(data = spatial_and_umap_stats_4E %>%
               filter(gene_short_name %in% equal_genes),
             aes(x = log10(morans_I_umap),
                 y = log10(morans_I_spatial)),
             stroke = 0,
             size = 1,
             color = "black") +
  geom_text_repel(data = spatial_and_umap_stats_4E %>%
                    filter(gene_short_name %in% other_genes|
                             gene_short_name %in% spatial_genes_to_highlight|
                             gene_short_name %in% equal_genes),
                  aes(x = log10(morans_I_umap),
                      y = log10(morans_I_spatial),
                      label = gene_short_name),
                  color = "black",
                  size = 2.5) +
  monocle3:::monocle_theme_opts() + 
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8)) +
  xlab("log(Morans I) [UMAP]") +
  ylab("log(Morans I) [Spatial]") +
  xlim(-2.5,NA) +
  ylim(-2,NA) +
  ggsave("Figures/Figure_Components/Supplement_hox_gene_expression/4E_neuron_other_gene_morans_I_all.pdf",
         height = 2.35,
         width = 2.35) 


