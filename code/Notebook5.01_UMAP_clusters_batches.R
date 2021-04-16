# Plot supplemental QC figures mainly in Supplemental Figure 9

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
  
  space_directory = "~/Google Drive File Stream/My Drive/sciSpace/"
  setwd(dir=space_directory)
  source("Submission_Data/bin/cell_cycle.R")
  cc.genes <- readRDS("Submission_Data/bin/cc.genes.mouse.RDS")
  
  # Pass TRUE if you want to see progress output on some of Monocle 3's operations
  DelayedArray:::set_verbose_block_processing(TRUE)
  # Passing a higher value will make some computations faster but use more memory. 
  # Adjust with caution!
  options(DelayedArray.block.size=1000e6)
  
  # Set a seed to make umap and other non-deterministic steps consistent
  set.seed(seed = 42)
  
})

spatial_cds = readRDS("Submission_Data/E14_slides/RDS_intermediates/Notebook5_E14_spatial_CDS.RDS")

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


# Supplemental Figure 9 - Panel A ----------------------------------------

colData(spatial_cds)$n.genes = 
  Matrix::colSums(counts(spatial_cds) >= 1)

colData(spatial_cds) %>%
  as.data.frame() %>%
  dplyr::select(n.umi,
                n.genes,
                Cell,
                slide_id) %>%
  gather(key = "type",
         value = "n",
         c(n.umi,n.genes)) %>% 
  ggplot() +
  geom_boxplot(aes(x = slide_id,
                   y = log10(n),
                   fill = type),
               outlier.stroke =  0,
               outlier.size =  1,
               size = 0.25) +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,
                                   size = 6,
                                   hjust = 1),
        axis.text.y = element_text(size = 6)) +
  monocle3:::monocle_theme_opts() +
  xlab("") +
  ylab("Log10(Observations)") +
  ggsave("Figures/Figure_Components/Supplement_experiment_QC/UMIs_per_slide.pdf",
         height = 2,
         width = 2.5) 


# Supplemental Figure 9 Panel C -- cluster proportions --------------------

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
colors = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
cluster_colors = randomcoloR::randomColor(count = 53)

colData(spatial_cds) %>%
  as.data.frame() %>% 
  group_by(max_slide_id) %>%
  add_tally(name = "num_in_slide") %>%
  group_by(slide_id, cluster, num_in_slide) %>%
  summarise(n = n()) %>%
  mutate(percent_slide = n/num_in_slide,
         cluster2 = 
           cluster %>%
           as.character() %>%
           as.numeric()) %>% 
  ggplot() +
  geom_bar(aes(x = slide_id,
               y = percent_slide,
               fill = cluster),
           color = "black",
           size = 0.15,
           stat = "identity") +
  monocle3:::monocle_theme_opts() +
  theme(legend.position = "none") + 
  xlab("")  +
  ylab("Percent Cluster in Slide") +
  theme(axis.text.x = element_text(size = 6,angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 8)) +
  scale_fill_manual(values = cluster_colors) +
  ggsave("Figures/Figure_Components/Supplement_experiment_QC/cluster_barplot.pdf",
         height = 3.0,
         width = 2.5)

# Supplemental Figure 9 Panel B -- cluster proportions --------------------

ggplot() +
  geom_point(data = 
               rbind(colData(spatial_cds) %>% as.data.frame() %>% mutate(order = "1"),
                     colData(spatial_cds) %>% as.data.frame() %>% mutate(order = "2")) %>%
               arrange(cluster,order) %>%
               mutate(cluster = ifelse(order ==2,
                                       cluster %>% as.character(),
                                       NA)),
             aes(x = umap1,
                 y = umap2,
                 color = cluster,
                 size = order),
             stroke = 0) +
  scale_size_manual(values = c("1" = 0.35,
                               "2" = 0.25)) +
  theme_void() +
  theme(legend.position = "none") + 
  scale_color_manual(values = cluster_colors,
                     na.value="black") +
  ggsave("Figures/Figure_Components/Supplement_experiment_QC/umap_by_cluster.png",
         dpi = 300,
         height = 2.5,
         width = 2.5)


# Supplemental Figure 9 Panel D -- cluster proportions --------------------

ggplot() +
  geom_point(data = 
               colData(spatial_cds) %>%
               as.data.frame(),
             aes(x = umap1,
                 y = umap2),
             color = "black",
             size = 0.35,
             stroke = 0) +
  geom_point(data = 
               colData(spatial_cds) %>%
               as.data.frame() %>% 
               mutate(mouse = ifelse(max_slide_id %in% c("slide_1F","slide_1E"),
                                     "emb2",
                                     "emb1")) %>%
               sample_n(size = dim(spatial_cds)[2]),
             aes(x = umap1,
                 y = umap2,
                 color = mouse),
             size = 0.25,
             stroke = 0) +
  theme_void() +
  theme(legend.position = "none") + 
  scale_color_brewer(palette = "Set2") +
  ggsave("Figures/Figure_Components/Supplement_experiment_QC/umap_by_emb.png",
         dpi = 300,
         height = 1.5,
         width = 1.5)


