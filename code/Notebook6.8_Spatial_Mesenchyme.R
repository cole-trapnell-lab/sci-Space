# Map the sub-clusters comprising the connective tissue cell type label
# Used to make Fig 4C

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


# Make a chondrocyte CDS --------------------------------------------------

chondrocyte_cds = 
  new_cell_data_set(expression_data = counts(spatial_cds[,colData(spatial_cds)$final_cluster_label == "Chondrocytes"]),
                    cell_metadata = colData(spatial_cds[,colData(spatial_cds)$final_cluster_label == "Chondrocytes"]),
                    gene_metadata = rowData(spatial_cds[,colData(spatial_cds)$final_cluster_label == "Chondrocytes"])) %>%
  estimate_size_factors() %>%
  preprocess_cds(method = "PCA") %>%
  align_cds(residual_model_formula_str = "~log.n.umi + sample",
            alignment_group = "sample") %>%
  reduce_dimension(umap.metric = "cosine",
                   umap.fast_sgd = F,
                   max_components = 3)
chondrocyte_cds = 
  chondrocyte_cds  %>%
  cluster_cells(resolution = 1e-3)
  
plot_cells_3d(chondrocyte_cds)

colData(chondrocyte_cds)$chondrocyte_cluster = clusters(chondrocyte_cds)

# This scaling data frame just reorients all the coordinates so that all the
# slides are facing the same way

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
  inner_join(scaling_df,
             by = "max_slide_id") 

jittered_points =
  jittered_points %>%
  mutate(jittered_y = coords.x1 * y + rnorm(n =dim(jittered_points)[1],
                                            mean = 0,
                                            sd = 10),
         jittered_x = coords.x2 * x + rnorm(n =dim(jittered_points)[1],
                                            mean = 0,
                                            sd = 10))

jittered_points2 = 
  colData(chondrocyte_cds) %>%
  as.data.frame() %>%
  inner_join(scaling_df,
             by = "max_slide_id")

jittered_points2 =
  jittered_points2 %>%
  mutate(jittered_y = coords.x1 * y + rnorm(n =dim(jittered_points2)[1],
                                            mean = 0,
                                            sd = 15),
         jittered_x = coords.x2 * x + rnorm(n =dim(jittered_points2)[1],
                                            mean = 0,
                                            sd = 15))
# Supplemental Figure S31B ------------------------------------------------
colors =  
  c("gold2","deepskyblue2","chartreuse3","navy","darkgreen","orangered2","plum4","darkgoldenrod4")

ggplot() +
  geom_point(data = jittered_points,
              aes(x = jittered_x,
                  y = jittered_y),
              color = "grey90",
              stroke = 0,
              size = 0.5,
             alpha = 0.5) +
    geom_point(data = jittered_points2 %>% 
                 arrange(as.numeric(as.character(chondrocyte_cluster))),
              aes(x = jittered_x,
                  y = jittered_y,
                  color = chondrocyte_cluster),
              stroke = 0,
              size = 0.85) +
  facet_wrap(~slide_id,scales = "free",nrow =3)  +
  theme_void() + 
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(size=2.5),
                               title = "Subcluster")) +
  ggsave("Figures/Figure_Components/Supplement_Connective_Tissue/chondrocyte_clusters.png",
         dpi = 600,
         height = 4,
         width = 6,
         bg = "transparent")


# Figure 4C ---------------------------------------------------------------

# Figure 4C (left)
ggplot() +
  geom_sf(data = all_image_data$slide_polygon[[6]] * rbind(c(0, -1),c(1,0)),
          fill = NA, 
          color = "black", 
          size = .15) +
  geom_point(data = jittered_points2 %>% 
               filter(slide_id %in% c("Slide 6")) %>%
               arrange(as.character(chondrocyte_cluster)),
             aes(x = jittered_x,
                 y = jittered_y,
                 color = chondrocyte_cluster),
             stroke = 0,
             size = 0.640) +
  theme_void() + 
  scale_color_manual(values = colors) +
  theme(legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size=2.5),
                               title = "Subcluster")) +
  ggsave("Figures/Figure_Components/Figure4/slide6_chondrocyte.png",
         height = 1,
         width = 1,
         dpi = 600)

# Figure 4C (middle)
ggplot() +
  geom_sf(data = all_image_data$slide_polygon[[11]] * rbind(c(0, -1),c(-1,0)),
          fill = NA, 
          color = "black", 
          size = .15) +
  geom_point(data = jittered_points2 %>% 
               filter(slide_id %in% c("Slide 11")) %>%
               arrange(as.character(chondrocyte_cluster)),
             aes(x = jittered_x,
                 y = jittered_y,
                 color = chondrocyte_cluster),
             stroke = 0,
             size = 0.650) +
  theme_void() + 
  scale_color_manual(values = colors) +
  theme(legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size=2.5),
                               title = "Subcluster")) +
  ggsave("Figures/Figure_Components/Figure4/slide11_chondrocyte.png",
         height = 1,
         width = 1,
         dpi = 600)

# Figure 4C (right)
ggplot() +
  geom_sf(data = all_image_data$slide_polygon[[14]] * rbind(c(0, 1),c(1,0)),
          fill = NA, 
          color = "black", 
          size = .15) +
  geom_point(data = jittered_points2 %>% 
               filter(slide_id %in% c("Slide 14")) %>%
               arrange(as.character(chondrocyte_cluster)),
             aes(x = jittered_x,
                 y = jittered_y,
                 color = chondrocyte_cluster),
             stroke = 0,
             size = 0.640) +
  theme_void() + 
  scale_color_manual(values = colors) +
  theme(legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size=2.5),
                               title = "Subcluster")) +
  ggsave("Figures/Figure_Components/Figure4/slide14_chondrocyte.png",
         height = 1,
         width = 1,
         dpi = 600)

# Perform differential gene expression testing ----------------------------
two.set.differential.gene.test = 
  function(cds, set.1.filter, set.2.filter, formal = F, cores = 1,
           covariates = NULL, sample.n.max = NULL) {
    s1.cds = cds[, set.1.filter]
    s2.cds = cds[, set.2.filter]
    
    if (!is.null(sample.n.max)) {
      set.seed(42)
      if (ncol(s1.cds) > sample.n.max)
        s1.cds = s1.cds[, sample(1:ncol(s1.cds), sample.n.max)]
      if (ncol(s2.cds) > sample.n.max)
        s2.cds = s2.cds[, sample(1:ncol(s2.cds), sample.n.max)]
    }
    
    message(paste("# of cells sampled in set 1:", ncol(s1.cds)))
    message(paste("# of cells sampled in set 2:", ncol(s2.cds)))
    
    s1.norm.expr = monocle3:::normalize_expr_data(s1.cds,
                                                  norm_method = "size_only",
                                                  pseudo_count = 0)
    s2.norm.expr = monocle3:::normalize_expr_data(s2.cds,
                                                  norm_method = "size_only",
                                                  pseudo_count = 0)
    
    s1.tpm = Matrix::rowSums(s1.norm.expr)
    s1.tpm = s1.tpm / sum(s1.tpm) * 1000000
    s2.tpm = Matrix::rowSums(s2.norm.expr)
    s2.tpm = s2.tpm / sum(s2.tpm) * 1000000
    
    s1.n.umi = Matrix::rowSums(counts(s1.cds))
    s2.n.umi = Matrix::rowSums(counts(s2.cds))
    
    higher.expr = ifelse(s1.tpm > s2.tpm, "Set 1", "Set 2")
    higher.tpm = ifelse(higher.expr == "Set 1", s1.tpm, s2.tpm)
    
    s1.ratio = s1.tpm / (s2.tpm + 1)
    s2.ratio = s2.tpm / (s1.tpm + 1)
    log2.ratio = ifelse(
      s1.tpm == 0 & s2.tpm == 0, 0, ifelse(
        higher.expr == "Set 1", log2(s1.ratio), log2(s2.ratio)))
    
    s1.n.expr = Matrix::rowSums(counts(s1.cds))
    s2.n.expr = Matrix::rowSums(counts(s2.cds))
    
    res = data.frame(
      gene = fData(cds)$gene_short_name,
      set.1.n.umi = s1.n.umi,
      set.2.n.umi = s2.n.umi,
      set.1.tpm = s1.tpm,
      set.2.tpm = s2.tpm,
      higher.expr = higher.expr,
      log2.ratio = log2.ratio,
      heuristic.score = log2.ratio * log2(higher.tpm+1)
    ) %>% arrange(-heuristic.score)
    
    if (formal) {
      pData(cds)$tmp = ifelse(set.1.filter, 1, ifelse(set.2.filter, 2, NA))
      
      cds.subset = cds[, set.1.filter | set.2.filter]
      cds.subset = estimate_size_factors(cds.subset)
      #cds.subset = estimateDispersions(cds.subset)
      cds.subset = detect_genes(cds.subset,min_expr =  0.1)
      
      expressed.genes = subset(fData(cds.subset), num_cells_expressed >= 5)[, 1]
      message(paste(length(expressed.genes), "genes expressed in at least 5 cells across both sets"))
      message("Computing differential expression p-values")
      
      DEG = fit_models(cds = cds.subset[expressed.genes,], model_formula_str = "~ tmp")
      
      res = inner_join(res,
                       DEG %>% 
                         coefficient_table() %>%
                         select(everything(),
                                gene = gene_short_name),
                       by = "gene")
      
      pData(cds)$tmp = NULL
    }
    
    return(res)
  }

differential_test_results = 
  lapply(X = colData(chondrocyte_cds)$chondrocyte_cluster %>% unique() %>% as.character(),
       FUN = function(cluster){
         two.set.differential.gene.test(chondrocyte_cds,
                                        pData(chondrocyte_cds)$chondrocyte_cluster == cluster,
                                        pData(chondrocyte_cds)$chondrocyte_cluster != cluster,
                                        sample.n.max = 2500,
                                        formal = F) %>%
           mutate(chondrocyte_cluster = cluster)
       })

differential_test_results  = 
  do.call(rbind,differential_test_results)


differential_test_results = 
  differential_test_results %>%
  mutate(log2.ratio2 = log2(set.1.tpm /(set.2.tpm)))

  
# Supplemental Figure S31B ------------------------------------------------
differential_test_results %>% 
  filter(higher.expr == "Set 1",
         chondrocyte_cluster %in% c("6")) %>%
  group_by(chondrocyte_cluster) %>%
  arrange(chondrocyte_cluster,-heuristic.score) %>%
  top_n(n = 25,wt = heuristic.score) %>%
  as.data.frame() %>%
  ggplot() +
  geom_point(aes(x =log2.ratio,
                 y =reorder(gene,log2.ratio))) +
  xlab("Log2 Fold Change") +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 6)) +
  ggsave("Figures/Figure_Components/Supplement_Connective_Tissue/cluster6_genes.pdf",
         dpi = 600,
         height = 2.5,
         width = 2,
         bg = "transparent")


differential_test_results %>% 
  filter(higher.expr == "Set 1",
         chondrocyte_cluster %in% c("3")) %>%
  group_by(chondrocyte_cluster) %>%
  arrange(chondrocyte_cluster,-heuristic.score) %>%
  top_n(n = 25,wt = heuristic.score) %>%
  as.data.frame() %>%
  ggplot() +
  geom_point(aes(x =log2.ratio,
                 y =reorder(gene,log2.ratio))) +
  xlab("Log2 Fold Change") +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 6)) +
  ggsave("Figures/Figure_Components/Supplement_Connective_Tissue/cluster7_genes.pdf",
         dpi = 600,
         height = 2.5,
         width = 2,
         bg = "transparent")

