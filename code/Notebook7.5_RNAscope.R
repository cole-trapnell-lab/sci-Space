# Estimate coexpression of genes by RNAScope HiPlex Assay slide scans

# Load startup packages ---------------------------------------------------
suppressPackageStartupMessages({ 
  library(tidyverse)
  library(ggplot2)
  library(purrr)
  library(sp)
  library(sf)
  library(vec2dtransf)
  library(monocle3)
  library(garnett)
  library(spatstat)
  
  space_directory = "/Users/maryregier/Desktop/Full_sciSpace"
  setwd(dir=space_directory)
  
  # Pass TRUE if you want to see progress output on some of Monocle 3's operations
  DelayedArray:::set_verbose_block_processing(TRUE)
  # Passing a higher value will make some computations faster but use more memory. 
  # Adjust with caution!
  options(DelayedArray.block.size=1000e6)
  set.seed(42)
})

# Set parameters -----
# 'ks' and 'ks2' -> set to 2 <- number of other points within sub-cellular distance -> set to 10 <- um of the test point required for designation as a positive point (ks) and a coexpressing point (ks2)
ks = 2
ks2 = 2
cell = 10



#set transcript gene names and indices
genelist = c("Cyp26b1", "Gad2", "Slc17a6", "Hoxa10" , "Pax6", "Slc17a7", "Lum", "Cldn5")

# set 'assay_ind' for Cyp26b1 -> set to 1 <- and for Hoxa1 -> set to 4 <-
assay_ind = 4
marker_ind = c(2, 3, 5, 6, 7, 8)

# set indices for probe sets
probe_ind1 = c(1:4)  # first set -> set to 1:4 <-
probe_ind2 = c(5:8)  # second set -> set to 5:8 <-

# set slide names and index
# run set to 1 for slide A and 2 for slide B
slide_list = c(168, 170)
slide_ind = 1
###

### coordinates for aligning the second four probes scan to the first four probes scan
aff_coord = matrix(c(0.9973, 0.0004, 0.0010, 0.9995, -1689.8, -1279.5609, .9998, 0.0001, 0.0002, 0.9989, 166.2391, -1530.4), nrow = 6, ncol = 2)
# obtained from QuPath 0.2.3 using: Analyze/Interactive image alignment/Auto-align/Estimate transform
# settings: Registration type -> Affine transform ; Alignment type -> Image intensity ; Pixel size 20
# source: Bankhead, P. et al. (2017). QuPath: Open source software for digital pathology image analysis. Scientific Reports. https://doi.org/10.1038/s41598-017-17204-5
###


coordstring = "Submission_Data/E14_slides/RNAScopeCoord/"

slide = slide_list[slide_ind]

colorlist = c("gray", "blue", "magenta", "gray" , "green", "cyan", "orange", "green")

all_pos = data.frame(x = NA, y = NA, gene = NA)
all_coex = data.frame(x = NA, y = NA, gene = NA)

# Load data -----

# import first four probes 
for (i in probe_ind1){
  gene = genelist[i]
  gene_string = paste(coordstring, slide, "-", gene, "-points.tsv", sep = "")
  gene_coord = read.table(gene_string, header = T)
  
  gene_pp = as.ppp(gene_coord, c(min(gene_coord[,1]),max(gene_coord[,1]),min(gene_coord[,2]),max(gene_coord[,2]))) 
  
  if (i == 4){
    gene_coord$gene = nndist(gene_pp, k=1)<cell
  }
  else {
    gene_coord$gene = nndist(gene_pp, k=ks)<cell
  }
  gene_coord = gene_coord %>% filter(gene_coord$gene == TRUE)
  gene_coord$gene = gene
  
  all_pos = rbind(gene_coord, all_pos)
}

# import probes 5-8 and align to probes 1-4 with affine transformation
for (i in probe_ind2){
  gene = genelist[i]
  gene_string = paste(coordstring, slide, "-", gene, "-points.tsv", sep = "")
  gene_coord = read.table(gene_string, header = T)
  
  gene_pp = as.ppp(gene_coord, c(min(gene_coord[,1]),max(gene_coord[,1]),min(gene_coord[,2]),max(gene_coord[,2]))) 
  
  aff = affine(gene_pp, matrix(c(aff_coord[1:4,slide_ind]),ncol=2))
  aff = as.data.frame(aff)
  
  aff[, 1] = aff[, 1] - aff_coord[5,slide_ind]
  aff[, 2] = aff[, 2] - aff_coord[6,slide_ind]

  gene_coord = aff

  aff_pp = as.ppp(aff, c(min(aff[,1]),max(aff[,1]),min(aff[,2]),max(aff[,2])))
 
  if (assay_ind == 4){
    gene_coord$gene = nndist(aff_pp, k=1)<cell}
  else {
    gene_coord$gene = nndist(aff_pp, k=ks)<cell
  }
  gene_coord = gene_coord %>% filter(gene_coord$gene == TRUE)
  gene_coord$gene = gene
  
  all_pos = rbind(gene_coord, all_pos)
}

# load sciSpace data
spatial_cds =
  readRDS(file = "Submission_Data/E14_slides/RDS_intermediates/Notebook6_spatial_cds_anatomy.RDS")

all_image_data = 
  readRDS(file = "Submission_Data/E14_slides/RDS_intermediates/Notebook0_1_images_and_transformations.RDS")


# ID RNAscope coexpression positions -----

assay_pos = all_pos %>% filter(all_pos$gene == genelist[assay_ind])
assay_pp <- as.ppp(assay_pos, c(min(assay_pos[,1]),max(assay_pos[,1]),min(assay_pos[,2]),max(assay_pos[,2])))

for (i in marker_ind){
  test_pos = all_pos %>% filter(all_pos$gene == genelist[i])
  test_coex = test_pos
  test_pp <- as.ppp(test_pos, c(min(test_pos[,1]),max(test_pos[,1]),min(test_pos[,2]),max(test_pos[,2])))
  
  if (assay_ind == 4){
  test_coex$gene <- nncross(test_pp,assay_pp, k = 1) < cell}
else {
  test_coex$gene <- nncross(test_pp,assay_pp, k = ks2) < cell}
  test_coex = test_coex %>% filter(test_coex$gene == TRUE)
  test_coex$gene = genelist[i]

  all_coex = rbind(test_coex, all_coex)
}

# Make individual marker gene/Cyp26b1 coexpression dataframes 
Lumdp = all_coex %>% filter(all_coex$gene == "Lum")
Gad2dp = all_coex %>% filter(all_coex$gene == "Gad2")
Cldn5dp = all_coex %>% filter(all_coex$gene == "Cldn5")
Pax6dp = all_coex %>% filter(all_coex$gene == "Pax6")
Slc17a6dp = all_coex %>% filter(all_coex$gene == "Slc17a6")
Slc17a7dp = all_coex %>% filter(all_coex$gene == "Slc17a7")


# Fig 3 - Panel G - & Supplementary figure 26 - Panel C - RNAscope all markers coexpression with Cyp26b1 or Hoxa10 -----

plotstring = paste("Figures/Figure_Components/Supplement_other_spatial_&RNAscope/",slide,"_allmarkers_",genelist[assay_ind],"coex.png", sep = "")
point_size = 0.4

ggplot () +
  geom_jitter(data = all_pos,
              aes(x = all_pos[,1],
                  y = all_pos[,2]),
              color = "grey95",
              stroke = 0,
              size = point_size) +
  geom_jitter(data = assay_pos,
             aes(x = assay_pos[,1],
                 y = assay_pos[,2]),
             color = "grey70",
             stroke = 0,
             size = point_size) +
  geom_jitter(data = Lumdp,
              aes(x = Lumdp[,1],
                  y = Lumdp[,2]),
              color = "orange",
              stroke = 0,
              size = point_size) +
  geom_jitter(data = Gad2dp,
             aes(x = Gad2dp[,1],
                 y = Gad2dp[,2]),
             color = "blue",
             stroke = 0,
             size = point_size) +

  geom_jitter(data = Cldn5dp,
             aes(x = Cldn5dp[,1],
                 y = Cldn5dp[,2]),
             color = "red",
             stroke = 0,
             size = point_size) +
  geom_jitter(data = Pax6dp,
             aes(x = Pax6dp[,1],
                 y = Pax6dp[,2]),
             color = "green",
             stroke = 0,
             size = point_size) +
  geom_jitter(data = Slc17a6dp,
             aes(x = Slc17a6dp[,1],
                 y = Slc17a6dp[,2]),
             color = "magenta",
             stroke = 0,
             size = point_size) +
  geom_jitter(data = Slc17a7dp,
             aes(x = Slc17a7dp[,1],
                 y = Slc17a7dp[,2]),
             color = "cyan",
             stroke = 0,
             size = point_size) +

  theme_void() +
  ggsave(plotstring, 
         height = 3,
         width = 2)

# Pull expression and Cyp26b1 coexpression data for Neurons and Endos for slide A ------

if (slide_ind == 1 && assay_ind == 1) {
  
# combine neural markers expression (no Cyp26b1 coexpression )
Slc17a6_pos = all_pos %>% filter(all_pos$gene == "Slc17a6")
Slc17a7_pos = all_pos %>% filter(all_pos$gene == "Slc17a7")
Gad2_pos = all_pos %>% filter(all_pos$gene == "Gad2")

all_excitatory_pos = rbind(Slc17a6_pos, Slc17a7_pos)
all_neuron_pos = rbind(all_excitatory_pos, Gad2_pos)

# pull endothelial mark expression (no Cyp26b1 coexpression)
Cldn5_pos = all_pos %>% filter(all_pos$gene == "Cldn5")

# combine neural markers Cyp26b1 coexpression
all_excitatory_dp = rbind(Slc17a6dp,Slc17a7dp)
all_neuron_dp = rbind(all_excitatory_dp,Gad2dp)

# Supplementary figure 26 - Panel C- All neurons coexpression with Cyp26b1 for slide A -----

plotstring = paste("Figures/Figure_Components/Supplement_slide14_neuron_v_endo/",slide,"_all_neuron_",genelist[assay_ind],"coex.png", sep = "")
ggplot () +

  geom_jitter(data = all_pos,
             aes(x = all_pos[,1],
                 y =  all_pos[,2]),
             color = "grey95",
             stroke = 0,
             size = .4) +
  geom_jitter(data = all_neuron_pos,
             aes(x = all_neuron_pos[,1],
                 y =  all_neuron_pos[,2]),
             color = "light blue",
             stroke = 0,
             size = .6) +
  geom_jitter(data = all_neuron_dp,
             aes(x = all_neuron_dp[,1],
                 y = all_neuron_dp[,2]),
             color = "navy",
             stroke = 0,
             size = .6) +


  theme_void()+
  ggsave(plotstring, 
         height = 3,
         width = 2)

# # Supp fig - Panel F - Endos coexpression with Cyp26b1 for slide A -----
# 
# plotstring = paste("Figures/Figure_Components/Supplement_RNAScope/",slide,"_Cldn5pos_",genelist[assay_ind],"coex.png", sep = "")
# ggplot () +
#   geom_jitter(data = all_pos,
#               aes(x = all_pos[,1],
#                   y =  all_pos[,2]),
#               color = "grey95",
#               stroke = 0,
#               size = .4) +
# 
#   geom_jitter(data = Cldn5_pos,
#               aes(x = Cldn5_pos[,1],
#                   y =  Cldn5_pos[,2]),
#               color = "pink",
#               stroke = 0,
#               size = .6) +
#   geom_jitter(data = Cldn5dp,
#               aes(x = Cldn5dp[,1],
#                   y = Cldn5dp[,2]),
#               color = "red",
#               stroke = 0,
#               size = .6) +
#   
#   
#   theme_void()+
#   ggsave(plotstring, 
#          height = 3,
#          width = 2)
# 
# }
# 


# Functions for aggregrating and plotting sciSpace data ----

# Function to aggregrate marker expression 

aggregate_markers = 
  function(markers,cds){
    cds_subset = cds[rowData(cds)$gene_short_name %in% markers,]
    
    marker_expr = 
      (Matrix::t(exprs(cds_subset))/size_factors(cds_subset)) %>% 
      as.matrix() %>%
      as.data.frame()
    
    marker_expr[marker_expr == 0] <- NA
    
    colnames(marker_expr) = 
      rowData(cds)[colnames(marker_expr),"gene_short_name"]
    
    marker_expr =
      marker_expr %>%
      rownames_to_column(var = "Cell") %>%
      gather(key = "gene",
             value = "expression",
             -Cell)
    
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
    
  }

# Function to crop an image based on a polygon shape file 

crop_hull = 
  function(im_df,
           slide_name){
    
    image_hull =
      im_df[which(im_df$slide_id == slide_name),"dapi_image"][[1]][[1]] %>%
      dplyr::select(x,
                    y) %>%
      as.matrix() %>%
      SpatialPoints() %>%
      over(as(object = im_df[which(im_df$slide_id == slide_name),"slide_polygon"][[1]][[1]],
              Class = "Spatial")) %>%
      is.na() %>%
      not()
    
    return(im_df[which(im_df$slide_id == slide_name),"dapi_image"][[1]][[1]][image_hull %>% as.vector,])
  }



### Slide 1 -------
slide_image =
  crop_hull(all_image_data,
            "slide_1D")

# Find coexpressing cells ----------------------------------------------------

  assay_marker = genelist[assay_ind]
  markers = c(assay_marker)
  
  sciSpace_marker_exprs =
    aggregate_markers(markers,
                      spatial_cds[,colData(spatial_cds)$max_slide_id == "slide_1D"])
  sciSpace_marker_exprs <- na.omit(sciSpace_marker_exprs)
  
for (i in marker_ind){
  
test_marker = genelist[i]
assay_marker = genelist[assay_ind]
markers = c(test_marker, assay_marker)

slide_marker_exprs =
  aggregate_markers(markers,
                    spatial_cds[,colData(spatial_cds)$max_slide_id == "slide_1D"])

slide_marker_exprs <- na.omit(slide_marker_exprs)
slide_marker_exprs_dup <- slide_marker_exprs[duplicated(slide_marker_exprs$Cell),]
coex_string = paste(genelist[i])
slide_marker_exprs_dup[,2] <- coex_string
sciSpace_marker_exprs <- rbind(sciSpace_marker_exprs,slide_marker_exprs_dup)
}

assay_sciSpace = sciSpace_marker_exprs %>% filter(sciSpace_marker_exprs$gene == assay_marker)    
Cyp26b1_sciSpace = sciSpace_marker_exprs %>% filter(sciSpace_marker_exprs$gene == "Cyp26b1")  
Slc17a6_sciSpace = sciSpace_marker_exprs %>% filter(sciSpace_marker_exprs$gene == "Slc17a6")
Slc17a7_sciSpace = sciSpace_marker_exprs %>% filter(sciSpace_marker_exprs$gene == "Slc17a7")
Gad2_sciSpace = sciSpace_marker_exprs %>% filter(sciSpace_marker_exprs$gene == "Gad2")
Cldn5_sciSpace = sciSpace_marker_exprs %>% filter(sciSpace_marker_exprs$gene == "Cldn5")
Pax6_sciSpace = sciSpace_marker_exprs %>% filter(sciSpace_marker_exprs$gene == "Pax6")
Lum_sciSpace = sciSpace_marker_exprs %>% filter(sciSpace_marker_exprs$gene == "Lum")

point_size = .6

# Fig 3 - Panel G - & Supplementary figure 26 - Panel C  — Slide 1 ----------------------------------------

plotstring = paste("Figures/Figure_Components/Supplement_other_spatial_&RNAscope/","sciSpace_1D_",genelist[assay_ind],"coex.png", sep = "")

ggplot () +
  
  geom_sf(data = all_image_data$slide_polygon[[1]] * rbind(c(0, 1),c(1,0)),
          fill = "gray95",
          color = "grey90",
          size = .15) +
  geom_jitter(data = assay_sciSpace,
             aes(x = assay_sciSpace[,5],
                 y = assay_sciSpace[,4]),
             color = "grey70",
             stroke = 0,
             size = point_size,
             width = 15,
             height = 15) +
  geom_jitter(data = Slc17a6_sciSpace,
              aes(x = Slc17a6_sciSpace[,5],
                  y = Slc17a6_sciSpace[,4]),
              color = "magenta",
              stroke = 0,
              size = point_size,
              width = 15,
              height = 15) +
  geom_jitter(data = Gad2_sciSpace,
             aes(x = Gad2_sciSpace[,5],
                 y = Gad2_sciSpace[,4]),
             color = "blue",
             stroke = 0,
             size = point_size,
             width = 15,
             height = 15) +
  geom_jitter(data = Lum_sciSpace,
             aes(x = Lum_sciSpace[,5],
                 y = Lum_sciSpace[,4]),
             color = "orange",
             stroke = 0,
             size = point_size,
             width = 15,
             height = 15) +
  geom_jitter(data = Cldn5_sciSpace,
             aes(x = Cldn5_sciSpace[,5],
                 y = Cldn5_sciSpace[,4]),
             color = "red",
             stroke = 0,
             size = point_size,
             width = 15,
             height = 15) +
  geom_jitter(data = Pax6_sciSpace,
             aes(x = Pax6_sciSpace[,5],
                 y = Pax6_sciSpace[,4]),
             color = "green",
             stroke = 0,
             size = point_size,
             width = 15,
             height = 15) +
  geom_jitter(data = Slc17a7_sciSpace,
             aes(x = Slc17a7_sciSpace[,5],
                 y = Slc17a7_sciSpace[,4]),
             color = "cyan",
             stroke = 0,
             size = point_size,
             width = 15,
             height = 15) +
  
  theme_void() +
  ggsave(plotstring, 
         height = 3,
         width = 2)


### Slide 14 ------


slide_image =
  crop_hull(all_image_data,
            "slide_4E")

# Find coexpressing cells ----------------------------------------------------

assay_marker = genelist[assay_ind]
markers = c(assay_marker)

sciSpace_marker_exprs =
  aggregate_markers(markers,
                    spatial_cds[,colData(spatial_cds)$max_slide_id == "slide_4E"])
sciSpace_marker_exprs <- na.omit(sciSpace_marker_exprs)

for (i in marker_ind){
  
  test_marker = genelist[i]
  assay_marker = genelist[assay_ind]
  markers = c(test_marker, assay_marker)
  
  slide_marker_exprs =
    aggregate_markers(markers,
                      spatial_cds[,colData(spatial_cds)$max_slide_id == "slide_4E"])
  
  slide_marker_exprs <- na.omit(slide_marker_exprs)
  slide_marker_exprs_dup <- slide_marker_exprs[duplicated(slide_marker_exprs$Cell),]
  coex_string = paste(genelist[i])
  slide_marker_exprs_dup[,2] <- coex_string
  sciSpace_marker_exprs <- rbind(sciSpace_marker_exprs,slide_marker_exprs_dup)
}

assay_sciSpace = sciSpace_marker_exprs %>% filter(sciSpace_marker_exprs$gene == assay_marker)    
Cyp26b1_sciSpace = sciSpace_marker_exprs %>% filter(sciSpace_marker_exprs$gene == "Cyp26b1")  
Slc17a6_sciSpace = sciSpace_marker_exprs %>% filter(sciSpace_marker_exprs$gene == "Slc17a6")
Slc17a7_sciSpace = sciSpace_marker_exprs %>% filter(sciSpace_marker_exprs$gene == "Slc17a7")
Gad2_sciSpace = sciSpace_marker_exprs %>% filter(sciSpace_marker_exprs$gene == "Gad2")
Cldn5_sciSpace = sciSpace_marker_exprs %>% filter(sciSpace_marker_exprs$gene == "Cldn5")
Pax6_sciSpace = sciSpace_marker_exprs %>% filter(sciSpace_marker_exprs$gene == "Pax6")
Lum_sciSpace = sciSpace_marker_exprs %>% filter(sciSpace_marker_exprs$gene == "Lum")

point_size = .6

# Fig 3 - Panel G - & Supplementary figure 26 - Panel C  — Slide 14 ----------------------------------------

plotstring = paste("Figures/Figure_Components/Supplement_other_spatial_&RNAscope/","sciSpace_4E_",genelist[assay_ind],"coex.png", sep = "")

ggplot () +
  
  geom_sf(data = all_image_data$slide_polygon[[14]] * rbind(c(0, 1),c(1,0)),
          fill = "gray95",
          color = "grey90",
          size = .15) +
  geom_jitter(data = assay_sciSpace,
              aes(x = assay_sciSpace[,5],
                  y = assay_sciSpace[,4]),
              color = "grey70",
              stroke = 0,
              size = point_size,
              width = 15,
              height = 15) +
  geom_jitter(data = Slc17a6_sciSpace,
              aes(x = Slc17a6_sciSpace[,5],
                  y = Slc17a6_sciSpace[,4]),
              color = "magenta",
              stroke = 0,
              size = point_size,
              width = 15,
              height = 15) +
  geom_jitter(data = Gad2_sciSpace,
              aes(x = Gad2_sciSpace[,5],
                  y = Gad2_sciSpace[,4]),
              color = "blue",
              stroke = 0,
              size = point_size,
              width = 15,
              height = 15) +
  geom_jitter(data = Lum_sciSpace,
              aes(x = Lum_sciSpace[,5],
                  y = Lum_sciSpace[,4]),
              color = "orange",
              stroke = 0,
              size = point_size,
              width = 15,
              height = 15) +
  geom_jitter(data = Cldn5_sciSpace,
              aes(x = Cldn5_sciSpace[,5],
                  y = Cldn5_sciSpace[,4]),
              color = "red",
              stroke = 0,
              size = point_size,
              width = 15,
              height = 15) +
  geom_jitter(data = Pax6_sciSpace,
              aes(x = Pax6_sciSpace[,5],
                  y = Pax6_sciSpace[,4]),
              color = "green",
              stroke = 0,
              size = point_size,
              width = 15,
              height = 15) +
  geom_jitter(data = Slc17a7_sciSpace,
              aes(x = Slc17a7_sciSpace[,5],
                  y = Slc17a7_sciSpace[,4]),
              color = "cyan",
              stroke = 0,
              size = point_size,
              width = 15,
              height = 15) +
  
  theme_void() +
  ggsave(plotstring, 
         height = 3,
         width = 2)

