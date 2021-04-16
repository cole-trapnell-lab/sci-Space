# Align sci-Space data containing E14 embryos to 
# developing mouse brain atlas (DMBA) dataset.

# This dataset can be downloaded at -- http://mousebrain.org/downloads.html

# G. La Manno, K. Siletti, A. Furlan, D. Gyllborg, E. Vinsland, 
# C. M. Langseth, I. Khven, A. Johnsson, M. Nilsson, P. Lönnerberg, 
# S. Linnarsson, Molecular architecture of the developing mouse brain. 
# Cold Spring Harbor Laboratory (2020), p. 2020.07.02.184051.

# https://www.biorxiv.org/content/10.1101/2020.07.02.184051v1

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
  library(scater)
  library(Seurat)
  
  DelayedArray:::set_verbose_block_processing(TRUE)
  options(DelayedArray.block.size=1000e7)
  
})

setwd(dir="/Volumes/GoogleDrive/Space/published_datasets/")

dev_brain_cds = 
  readRDS(dev_brain_cds,"La_Manno_developing_brain/LaManno_dev_brain_cds.RDS")

# Read in sci-Space data ------------------------------------------
cds_space = readRDS("/net/trapnell/vol1/home/sanjays/projects/Space/Submission/E14_slides/RDS_intermediates/Notebook6_spatial_cds_anatomy.RDS")

cm_space = counts(cds_space)

rowdata_space = 
  rowData(cds_space) %>%
  as.data.frame() %>%
  mutate(Acession_id = stringr::str_split_fixed(string = id,pattern = "\\.",n = 2)[,1])

rownames(rowdata_space) = rowdata_space$Acession_id
rownames(cm_space) = rowdata_space$Acession_id

coldata_space = 
  colData(cds_space) %>%
  as.data.frame() 

rownames(coldata_space) = colnames(cds_space)

cds_space = 
  new_cell_data_set(expression_data = cm_space,
                    cell_metadata = coldata_space,
                    gene_metadata = rowdata_space)

# Get the intersecting set of genes between the two datasets
intersecting_genes = 
  intersect(rowdata$Accession %>%
              as.character(),
            rowdata_space$Acession_id)

age_range = c("e14.5","e14.0","e13.5")

dev_brain_E14_cells =
  colData(dev_brain_cds) %>%
  as.data.frame() %>%
  filter(Age %in% age_range) %>%
  pull(CellID) %>%
  as.character()


## Isolate the cells around the E14 time_point and the genes that intersect the two datasets
dev_brain_cds_E14 = 
  dev_brain_cds[intersecting_genes,dev_brain_E14_cells]

# dim(dev_brain_cds_E14) [1] 30386 38703

# Isolate only the neural cells (partition 2),
# mapping the brain (anatomical_annotation == "Cortex")
neural_cells = 
  colData(cds_space) %>%
  as.data.frame() %>%
  filter(anatomical_annotation == "Cortex",
         partition == "2") %>%
  pull(Cell)


cds_space = cds_space[intersecting_genes,neural_cells]


coldata_space = 
  colData(cds_space) %>%
  as.data.frame() %>%
  dplyr::select(Cell, sample) %>%
  mutate(Age = "e14.0")

rownames(coldata_space) = coldata_space$Cell 

coldata_dev_brain = 
  colData(dev_brain_cds_E14) %>%
  as.data.frame() %>%
  dplyr::select(Cell = CellID,
                Age) %>%
  dplyr::mutate(sample = "LaManno")

rownames(coldata_dev_brain) = 
  colData(dev_brain_cds_E14)$CellID

joint_cds =
  new_cell_data_set(expression_data = 
                      cbind(counts(dev_brain_cds_E14), 
                            counts(cds_space)),
                    cell_metadata = 
                      rbind(coldata_dev_brain, 
                            coldata_space),
                    gene_metadata = 
                      rowData(dev_brain_cds_E14) %>%
                      as.data.frame())


# Use Seurat Integration to align the datasets ----------------------------

# https://satijalab.org/seurat/archive/v3.0/integration.html

count_mat = assay(joint_cds)
coldata_df = colData(joint_cds) %>% as.data.frame()
coldata_df$Cell = rownames(coldata_df)


cds_seurat = 
  CreateSeuratObject(counts = count_mat,
                     project = "sciSpace",
                     assay = "RNA",
                     meta.data = coldata_df)

cds.list <- 
  SplitObject(cds_seurat, 
              split.by = "sample")

features <- SelectIntegrationFeatures(object.list = cds.list)

cds.list <- lapply(X = cds.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = T)
  x <- RunPCA(x, features = features, verbose = T, npcs = 100)
})

anchors <- FindIntegrationAnchors(object.list = cds.list, 
                                  reduction = "rpca", 
                                  dims = 1:100,
                                  verbose = T)


cds.integrated <- IntegrateData(anchorset = anchors, 
                                dims = 1:100,
                                verbose = T)

cds.integrated <- ScaleData(cds.integrated, verbose = T)
cds.integrated <- RunPCA(cds.integrated, verbose = T, npcs = 100)
cds.integrated <- RunUMAP(cds.integrated, dims = 1:100)


# Extract the PCA co-embedding ----------------------------------------------

joint_cds =
  joint_cds %>%
  estimate_size_factors() %>%
  preprocess_cds(num_dim = 100)

seurat_embeddings_pcs =   
  Embeddings(cds.integrated, reduction = "pca") %>%
  as.matrix()

colnames(seurat_embeddings_pcs) = 
  colnames(reducedDim(joint_cds, 
                      type = "PCA"))


reducedDim(x = joint_cds, type = "PCA") = seurat_embeddings_pcs

joint_cds =
  joint_cds %>%
  reduce_dimension() %>%
  cluster_cells()

joint_cds$umap1 = reducedDim(joint_cds, type = "UMAP")[,1]
joint_cds$umap2 = reducedDim(joint_cds, type = "UMAP")[,2]


saveRDS(joint_cds,
        file = "La_Manno_developing_brain/sci_space_cortex_neuron_e13p5_e14p5_joint_embedding.RDS")






