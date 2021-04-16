# In vogue right now is the mapping of single cell data to data from spatially aggreagted
# arrays. This analysis aims to argue that single cell resolution is preferable for
# mapping. By using Seurat's integration workflow we ask how well data integration
# can recover the position of nucleus that has been mapped with sci-Space

# Load startup packages ---------------------------------------------------
suppressPackageStartupMessages({ 
  library(tidyverse)
  library(ggplot2)
  library(monocle3)
  library(ggplot2)
  library(sp)
  library(sf)
  library(Seurat)
  library(patchwork)

  space_directory = "/Volumes/GoogleDrive/My Drive/sciSpace/"
  setwd(dir=space_directory)

  # Pass TRUE if you want to see progress output on some of Monocle 3's operations
  DelayedArray:::set_verbose_block_processing(TRUE)
  # Passing a higher value will make some computations faster but use more memory. 
  # Adjust with caution!
  options(DelayedArray.block.size=1000e6)
  options(future.globals.maxSize= 891289600)
  set.seed(42)
})

# Read in single cell data
spatial_cds = 
  readRDS(file = "Submission_Data/E14_slides/RDS_intermediates/Notebook6_spatial_cds_anatomy.RDS")

# Read in previously aggregated sci-Space data
joint_cds = 
  readRDS(file = "Submission_Data/E14_slides/RDS_intermediates/Notebook6.5_joint_bulk_CDS.RDS")

colData(joint_cds)$position = as.character(colData(joint_cds)$position)

colData(joint_cds)$top_spot =
  paste(stringr::str_split_fixed(string = colData(joint_cds)$position,
                                 pattern = "_",
                                 n = 3)[,1],
        stringr::str_split_fixed(string = colData(joint_cds)$position,
                                 pattern = "_",
                                 n = 3)[,2],
        sep = "_")


# Making a Spatial Seurat Object  ---------------------------------------

# This analysis was done on slide 14 (Slide 4E) 
colData(joint_cds) %>%
  as.data.frame() %>%
  filter(max_slide_id == "slide_4E") %>%
  dim()

colData(spatial_cds) %>%
  as.data.frame() %>%
  filter(max_slide_id == "slide_4E") %>%
  dplyr::select(top_spot,max_slide_id,coords.x1,coords.x2) %>%
  arrange(top_spot,coords.x1,coords.x2) %>%
  distinct()



# An applicable note used to load sci-Space data into Seurat

# Load Brain data from 10x genomics

# https://github.com/satijalab/seurat/issues/2790 
# Hi,
# 
# This is an old issue but for anyone else who runs into this and doesn't want to go through Spatial, I also had a hard time finding documentation for spatial seurat construction, although I probably missed it somewhere. (@ Seurat maintainers, I think it could be helpful to add a section in the spatial vignette pointing towards the canonical way of constructing it. Both for 10X/Slideseq (which do have the nice constructors) and also for generic spatial data for other modalities like FISH)
# 
# After reading through some of the constructor files, I think the easiest route to go through for generic XY spatial data is the ReadSlideSeq constructor which takes in a .csv file. But if you already have the coordinates loaded into a variable (say a list of X, Y coordinates and BARCODES (same order as in your CreateSeuratObject constructor just to be safe) then this preliminarily worked for me:
# 
# slide.seq = CreateSeuratObject(counts = COUNTS_MTX, assay="Spatial")
# 
# coord.df = data.frame(x=X, y=Y, stringsAsFactors=FALSE) # (stringsAsFactors only if also have a separate barcodes column)
# rownames(coord.df) = BARCODES
# 
# slide.seq@images$image =  new(
# Class = 'SlideSeq',
# assay = "Spatial",
# key = "image_",
# coordinates = coord.df
# )

slide_4e_bulked_cds = 
  joint_cds[,colData(joint_cds)$max_slide_id == "slide_4E"]

# X position of cells
colData(slide_4e_bulked_cds)$X_pos = 
  colData(slide_4e_bulked_cds) %>%
  as.data.frame() %>% 
  dplyr::select(top_spot,
                max_slide_id) %>%
  left_join(colData(spatial_cds) %>%
              as.data.frame() %>%
              filter(max_slide_id == "slide_4E") %>%
              dplyr::select(top_spot,max_slide_id,Col,Row) %>%
              distinct(),
            by = c("top_spot","max_slide_id")) %>%
  pull(Col) %>%
  as.character()

# Y position of cells
colData(slide_4e_bulked_cds)$Y_pos = 
  colData(slide_4e_bulked_cds) %>%
  as.data.frame() %>% 
  dplyr::select(top_spot,
                max_slide_id) %>%
  left_join(colData(spatial_cds) %>%
              as.data.frame() %>%
              filter(max_slide_id == "slide_4E") %>%
              dplyr::select(top_spot,max_slide_id,Col,Row) %>%
              distinct(),
            by = c("top_spot","max_slide_id")) %>%
  mutate(Row = 85 - as.numeric(as.character(Row)) * -1) %>%
  pull(Row) %>%
  as.character()


# Make Seurat object
sci.space.seq = 
  CreateSeuratObject(counts = counts(slide_4e_bulked_cds), 
                     assay="Spatial")

# Load in spatial coordinates
coord.df = 
  colData(slide_4e_bulked_cds) %>%
  as.data.frame() %>% 
  dplyr::select(x = X_pos,
                y = Y_pos)

identical(colData(slide_4e_bulked_cds)$position %>% as.character(),
          colnames(slide_4e_bulked_cds) %>% as.character())

rownames(coord.df) = colnames(slide_4e_bulked_cds)

sci.space.seq@images$image =  
  new(Class = 'SlideSeq',
      assay = "Spatial",
      key = "image_",
      coordinates = coord.df)


# Run Seurat workflow -----------------------------------------------------

sci.space.seq <- SCTransform(sci.space.seq, assay = "Spatial", verbose = FALSE)
sci.space.seq <- RunPCA(sci.space.seq, assay = "SCT", verbose = FALSE)

sci.space.seq <- FindNeighbors(sci.space.seq, reduction = "pca", dims = 1:30)
sci.space.seq <- FindClusters(sci.space.seq, verbose = FALSE)
sci.space.seq <- RunUMAP(sci.space.seq, reduction = "pca", dims = 1:30)

p1 <- DimPlot(sci.space.seq, 
              reduction = "umap", 
              label = TRUE,pt.size = 0.5)


random_colors = randomcoloR::randomColor(count = 25)
spatial_df =
  Embeddings(sci.space.seq, reduction = "umap") %>%
  as.data.frame()


spatial_df = 
  spatial_df %>%
  mutate(position = rownames(spatial_df) %>% as.character()) %>%
  left_join(data.frame(position = names(Idents(sci.space.seq)), 
                       cluster = Idents(sci.space.seq) %>% as.character()),
            by = "position")

# Supplemental fig S21D
ggplot(spatial_df) +
  geom_point(aes(x = UMAP_1,
                 y = UMAP_2),
             color = "black",
             size = 0.95,
             stroke = 0) +
  geom_point(aes(x = UMAP_1,
                 y = UMAP_2,
                 color = cluster),
             size = 0.8,
             stroke = 0) +
  theme_void() +
  theme(legend.position = "none") +
  scale_color_manual(values = random_colors) + 
  ggsave("Figures/Figure_Components/Supplement_bulk_aggregation/seurat_spatial_umap.pdf",
         height =1.75,
         width = 1.75)

# Supplemental fig S21E
SpatialDimPlot(sci.space.seq,pt.size.factor = 6.5,stroke = 0.1) +
  theme(legend.position = "none") +
  scale_fill_manual(values = random_colors) + 
  ggsave("Figures/Figure_Components/Supplement_bulk_aggregation/seurat_spatial_clusters.png",
         height =2,
         width = 2,
         bg = "transparent")



# Single cell analysis ----------------------------------------------------

# Get single sci-Space nuclei from the same slide 

slide_4E_single_cells = 
  spatial_cds %>%
  colData() %>%
  as.data.frame() %>%
  filter(max_slide_id == "slide_4E") %>%
  pull(Cell)


slide_4E_sc = 
  CreateSeuratObject(counts = spatial_cds[,slide_4E_single_cells] %>%
                       counts(),
                     project = "sciSpace",
                     assay = "RNA",
                     meta.data = spatial_cds[,slide_4E_single_cells] %>% 
                       colData() %>%
                       as.data.frame())


slide_4E_sc =
  SCTransform(slide_4E_sc, 
              ncells = 3000, 
              verbose = FALSE) %>% 
  RunPCA(verbose = FALSE) %>% 
  RunUMAP(dims = 1:30)

slide_4E_sc <- FindNeighbors(slide_4E_sc, dims = 1:10)
slide_4E_sc= FindClusters(slide_4E_sc)
slide_4E_sc$seurat_cluster = Idents(slide_4E_sc)

DimPlot(slide_4E_sc, reduction = "umap")

# the annotation is stored in the 'subclass' column of object metadata
DimPlot(slide_4E_sc, 
        group.by = "seurat_cluster", 
        label = F)

# get anchors for mapping against aggregated map
anchors <- 
  FindTransferAnchors(reference = sci.space.seq, 
                      query = slide_4E_sc, 
                      normalization.method = "SCT")

sci.space.seq$position = rownames(sci.space.seq@meta.data)


# Integrate single cell and spatial aggregated data ----------------------------

# Get the spatial predictions for these nuclei
predictions.assay <- 
  TransferData(anchorset = anchors, 
               refdata = sci.space.seq$position, 
               prediction.assay = TRUE,
               weight.reduction = slide_4E_sc[["pca"]])

slide_4E_sc[["predictions"]] <- predictions.assay


# Recover the top transferred position for each cell
top_transferred_position =
  sapply(X = seq(1,dim(slide_4E_sc[["predictions"]])[2]),
       FUN = function(x,
                      pred_matrix){
         names(which.max(pred_matrix[,x]))
       }, slide_4E_sc[["predictions"]][1:1391,])

# Recover the top transferred position probability for each cell
transfer_probability = 
  sapply(X = seq(1,dim(slide_4E_sc[["predictions"]])[2]),
         FUN = function(x,
                        pred_matrix){
           max(pred_matrix[,x])
         }, slide_4E_sc[["predictions"]][1:1391,])

# Make a dataframe for the comparision
transfer_df = 
  data.frame(Cell = colnames(slide_4E_sc[["predictions"]]),
             top_transferred_position = stringr::str_replace_all(top_transferred_position,pattern = "-","_"),
             transfer_probability = transfer_probability)

transfer_df$top_transferred_position =
  paste(stringr::str_split_fixed(transfer_df$top_transferred_position,
                                 pattern = "_",
                                 n = 3)[,1],
        stringr::str_split_fixed(transfer_df$top_transferred_position,
                                 pattern = "_",
                                 n = 3)[,2],
        sep = "_")

transfer_df$slide = "slide_4E"

transfer_df = 
  transfer_df %>%
  left_join(spatial_cds %>%
              colData() %>%
              as.data.frame() %>%
              dplyr::select(Cell, 
                            top_spot, 
                            cluster,
                            final_cluster_label),
            by = "Cell")

all_image_data = 
  readRDS("Submission_Data/E14_slides/RDS_intermediates/Notebook0_1_images_and_transformations.RDS")

# Get the information for this slide (where all the positions are)
oligo_layout_slide_4e = all_image_data$oligo_plate_layout[[14]]

transfer_df = 
  transfer_df %>%
  left_join(oligo_layout_slide_4e %>%
              dplyr::select(top_transferred_position = Oligo,
                            Row_transferred = Row,
                            Col_transferred = Col),
            by = "top_transferred_position")
  
transfer_df = 
  transfer_df %>%
  left_join(oligo_layout_slide_4e %>%
              dplyr::select(top_spot = Oligo,
                            Row,
                            Col),
            by = "top_spot")

# Calculate the distance between the top integrated spot and the measured position 
transfer_df$distance = 
  sqrt((transfer_df$Row_transferred - transfer_df$Row)^2 +
         (transfer_df$Col_transferred - transfer_df$Col)^2)


transfer_df %>%
  mutate(same_position = (top_spot == top_transferred_position)) %>%
  group_by(cluster) %>%
  summarise(percent_correct_position = sum(same_position)/n()) %>%
  arrange(desc(percent_correct_position))


transfer_df$distance %>%
  mean()

transfer_df %>%
  mutate(same_position = (top_spot == top_transferred_position)) %>%
  summarise(sum(same_position)/n())


# Supplemental Fig S21 F --------------------------------------------------

colors =c("Cardiac muscle lineages" = "#FF3333",
          "#B1E6E7",
          "#748FE2",
          "#D4E981",
          "#62B1D8",
          "#CD60C7",
          "#E76C43",
          "#DCE4C9",
          "#E744DC",
          "#E3538D",
          "#E7E041",
          "#96EA43",
          "#DBB5A0",
          "#8538E8",
          "#E9DEA2",
          "#776AE7",
          "#60E7B2",
          "#6AE6DF",
          "#E1AF52",
          "#D3CFE2",
          "#DE8785",
          "#E2AFCE",
          "#AEEBB9",
          "#D590DD",
          "#AFADE2",
          "#699B9A",
          "#92B56F",
          "#805B92",
          "#876F50",
          "#66DF6C",
          "black",
          "grey",
          "green")

transfer_df %>%
  mutate(same_position = (top_spot == top_transferred_position)) %>%
  group_by(cluster) %>%
  mutate(percent_correct_position = sum(same_position)/n()) %>%
  ggplot() +
  geom_boxplot(aes(x = reorder(final_cluster_label,percent_correct_position),
                   y = distance,
                   fill = final_cluster_label),
               outlier.size = 0.5,
               outlier.stroke = 0,
               size = 0.5,
               show.legend = F) +
  scale_fill_manual(values = colors) +
  theme_classic() +
  scale_y_continuous(breaks = seq(0,60,5),
                     limits = c(0,35)) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6)) +
  ylab(expression(paste(Delta,"Spot Distance (Transfer - Measurement)"))) +
  coord_flip() +
  ggsave("Figures/Figure_Components/Supplement_bulk_aggregation/seurat_spatial_mapping.pdf",
         height =2.25,
         width = 4)
    


