# Align sci-Space data containing E14 embryos to 
# mouse organogenesis cell atlas (MOCA) dataset using Seurat.
# Record nearest neighbors between the datasets to automate cell type 
# annotation

# All computations were run on a computing clustering 

# Load startup packages ---------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyverse)
  library(ggrepel)
  library(FNN)
  library(devtools)
  library(monocle3)
  DelayedArray:::set_verbose_block_processing(TRUE)
  options(DelayedArray.block.size=1000e7)
  
})

setwd(dir="~/projects/Space/published_datasets/Cao_et_al")


append_umap_coordinates = function(cds){
  colData(cds)$umap1 = reducedDim(x = cds,
                                  type = "UMAP")[,1]
  colData(cds)$umap2 = reducedDim(x = cds,
                                  type = "UMAP")[,2]
  return(cds)
}

# Read in Cao -- monocle3 object ------------------------------------------
cds_cao = readRDS("Cao_2019_cds_cleaned_monocle3.RDS")

new_annotations =
  read.table(file = "cell_annotate_20200119.csv",
             sep = ",",
             header = T,
             stringsAsFactors = F)

colData(cds_cao)$Main_trajectory =
  colData(cds_cao) %>%
  as.data.frame() %>%
  left_join(new_annotations %>% 
              dplyr::select(sample,
                            Main_trajectory),
            by = "sample") %>%
  pull(Main_trajectory.y)


colData(cds_cao)$Sub_trajectory_name =
  colData(cds_cao) %>%
  as.data.frame() %>%
  left_join(new_annotations %>% 
              dplyr::select(sample,
                            Sub_trajectory_name),
            by = "sample") %>%
  pull(Sub_trajectory_name.y)

colData(cds_cao)$Main_cell_type =
  colData(cds_cao) %>%
  as.data.frame() %>%
  left_join(new_annotations %>% 
              dplyr::select(sample,
                            Main_cell_type),
            by = "sample") %>%
  pull(Main_cell_type)

# Read in sci-Space data ------------------------------------------
cds_space = readRDS("/net/trapnell/vol1/home/sanjays/projects/Space/Submission/E14_slides/RDS_intermediates/Notebook3_E14_spatial_CDS.RDS")

## Isolate the E13.5 time_point from the Cao data
cds_cao_e13.5 = 
  cds_cao[,colData(cds_cao)$day == "13.5"]
  
# Get the intersecting set of genes between the two datasets
intersecting_genes = 
  intersect(rowData(cds_cao_e13.5)$gene_id %>%
              as.character(),
            rowData(cds_space)$id %>%
              as.character())

# Subset E13.5 MOCA CDS for intersecting genes
cds_cao_e13.5 = cds_cao_e13.5[intersecting_genes,]

# Subset spatial CDS with intersecting genes
cds_space = cds_space[intersecting_genes,]

coldata_cds_space = 
  colData(cds_space) %>%
  as.data.frame() %>%
  dplyr::select(Cell, sample) %>%
  mutate(stage = "14")

rownames(coldata_cds_space) = coldata_cds_space$Cell 

# downsample Cao 2019 E13.5 time-point with the same number of cells as the spatial data

coldata_cao_subsample = 
  colData(cds_cao_e13.5) %>%
  as.data.frame() %>%
  dplyr::select(Cell = sample,
                stage = day) %>%
  dplyr::mutate(sample = "Cao")

rownames(coldata_cao_subsample) = 
  coldata_cao_subsample$Cell

joint_cds =
  new_cell_data_set(expression_data = cbind(counts(cds_cao_e13.5), 
                                            counts(cds_space)),
                    cell_metadata = rbind(coldata_cao_subsample, 
                                          coldata_cds_space),
                    gene_metadata = rowData(cds_space) %>% 
                      as.data.frame())


# Use Seurat Integration ----------------------------------------------------------

# Load libraries
library(scater)
library(Seurat)

# Extract features from the monocle3 CDS to make a seurat object
count_mat = assay(joint_cds)
coldata_df = colData(joint_cds) %>% as.data.frame()
coldata_df$Cell = rownames(coldata_df)

cds_seurat = 
  CreateSeuratObject(counts = count_mat,
                     project = "sciSpace",
                     assay = "RNA",
                     meta.data = coldata_df)

# Break up the CDS by sample, which is either MOCA or one 
# of the experimental samples

cds.list <- 
  SplitObject(cds_seurat, 
              split.by = "sample")

features <- SelectIntegrationFeatures(object.list = cds.list)

cds.list <- lapply(X = cds.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = T)
  x <- RunPCA(x, features = features, verbose = T, npcs = 100)
})

# Anchor selection using reciprocal PCA. 
# We set the MOCA cells to be the reference
anchors <- FindIntegrationAnchors(object.list = cds.list, 
                                  reference = c(1), 
                                  reduction = "rpca", 
                                  dims = 1:100,
                                  verbose = T)


cds.integrated <- IntegrateData(anchorset = anchors, 
                                  dims = 1:100,
                                  verbose = T)

cds.integrated <- ScaleData(cds.integrated, verbose = T)
cds.integrated <- RunPCA(cds.integrated, verbose = T, npcs = 100)
cds.integrated <- RunUMAP(cds.integrated, dims = 1:100)

# Pull out the UMAP embeddings from the Seurat CDS object
seurat_umap_embeddings = 
  as.data.frame(Embeddings(cds.integrated, reduction = "umap"))

seurat_umap_embeddings$Cell = rownames(seurat_umap_embeddings)


# Label transfer ----------------------------------------------------------

# Make a dataframe containing the label information for the MOCA and sci-Space cells 
coldata_cao = 
  colData(cds_cao) %>%
  as.data.frame() %>%
  dplyr::select(Cell = sample,
                stage = day) %>%
  dplyr::mutate(sample = "Cao")

coldata_cds_space = 
  colData(cds_space) %>%
  as.data.frame() %>%
  dplyr::select(Cell, sample) %>%
  mutate(stage = "14")

seurat_umap_embeddings =
  seurat_umap_embeddings %>%
  left_join(rbind(coldata_cao,
                  coldata_cds_space),
            by = "Cell")

cao_seurat_umap = 
  seurat_umap_embeddings %>%
  filter(sample == "Cao") %>%
  left_join(colData(cds_cao) %>%
              as.data.frame() %>% 
              dplyr::select(Cell = sample,
                            Main_trajectory,
                            Sub_trajectory_name,
                            Main_cell_type),
            by = "Cell")


space_seurat_umap = 
  seurat_umap_embeddings %>%
  filter(sample != "Cao") 


# Get nearest neighbors using the KNN package
nn = 
  get.knnx(data = cao_seurat_umap %>%
             dplyr::select(UMAP_1, UMAP_2) %>%
             as.matrix(),
           query = space_seurat_umap %>%
             dplyr::select(UMAP_1, UMAP_2) %>%
             as.matrix(),
           k = 10)


## Use nearest neighbor indices to get top 10 labels
nn = 
  sapply(X = seq(1,dim(nn$nn.index)[2]),
         FUN = function(X){
           cao_seurat_umap$Main_cell_type[nn$nn.index[,X]]
         }
  ) %>%
  as.data.frame()

rownames(nn) =
  space_seurat_umap$Cell

# Collect the results -----------------------------------------------------
# Take the majority celltype annotation for each cell across all the iterations 
nn_majority = 
  sapply(X = seq(1,dim(nn)[1]),
         FUN = function(X){
           tabulated_labels = 
             nn[X,] %>%
             unlist() %>%
             table() %>%
             sort(decreasing = T)
           
           tabulated_labels[1] %>%
             names()
         })

# Calculate the proportion that of nearest neighbors in the majority label across all iterations 
majority_fraction =
  sapply(X = seq(1,dim(nn)[1]), 
                  FUN = function(Y){
                    sum(nn[Y,] == nn_majority[Y])/10
                  
                    })
estimates =
  data.frame(Cell = space_seurat_umap$Cell,
           majority_fraction = majority_fraction,
           nn_majority = nn_majority)

write.table(x = estimates,
            file = "alignment_coldata_seurat.tsv",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)



