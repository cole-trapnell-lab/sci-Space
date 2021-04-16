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

# All computations in this notebook were run on a computing clustering 

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

setwd(dir="~/projects/Space/published_datasets/")


# read in DMBA count matrix 
count_matrix = Matrix::readMM("La_Manno_developing_brain/count_matrix.mtx")

# read in DMBA cell metadata 
coldata = 
  read.table(file = "La_Manno_developing_brain/LaManno_coldata.tsv",
             header = T,
             row.names = 1,
             sep = "\t")

# read in DMBA cell genedata 
rowdata = 
  read.table(file = "La_Manno_developing_brain/LaManno_rowdata.tsv",
             header = T,
             row.names = 1,
             sep = "\t")

# For CDS objects the rownames of the gene metadata must be the rownames of 
# the count matrix and rownames of the cell metadata must be the colnames 
# of the count matrix

rownames(coldata) = coldata$CellID
rownames(rowdata) = rowdata$Accession
rowdata$gene_short_name = rowdata$Gene

# The count matrix was cells * genes and it needs to be transposed
cm = t(count_matrix %>% as.matrix())

dimnames(cm) = list(rownames(rowdata),coldata$CellID)

# Make CDS object
dev_brain_cds = 
  new_cell_data_set(expression_data = cm,
                    cell_metadata = coldata,
                    gene_metadata = rowdata)

coldata$Age %>%
  unique()

saveRDS(dev_brain_cds,
        "La_Manno_developing_brain/LaManno_dev_brain_cds.RDS")

# Read in sci-Space data ------------------------------------------
cds_space = readRDS("/net/trapnell/vol1/home/sanjays/projects/Space/Submission/E14_slides/RDS_intermediates/Notebook6_spatial_cds_anatomy.RDS")

cm_space = counts(cds_space)

# Genes have isoform tag that needs to be removed to match up with DMBA data
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

## Isolate the cells around the E14 time_point and the genes that intersect the two datasets
dev_brain_cds_E14 = 
  dev_brain_cds[intersecting_genes,]

# dim(dev_brain_cds_E14) [1] 30386 38703

cds_space = cds_space[intersecting_genes,]

# Make matching column data for DMBA and sci-Space cells
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

# Make a joint CDS object
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

# Using UMAP in Monocle 3 -- 
# I like the default UMAP hyperparameters in monocle3

# Do PCA on the joint CDS
joint_cds =
  joint_cds %>%
  estimate_size_factors() %>%
  preprocess_cds(num_dim = 100)

# Replace PCA with the integrated PCA coordinates from Seurat
seurat_embeddings_pcs =   
  Embeddings(cds.integrated, reduction = "pca") %>%
  as.matrix()

colnames(seurat_embeddings_pcs) = 
  colnames(reducedDim(joint_cds, 
                      type = "PCA"))

reducedDim(x = joint_cds, type = "PCA") = seurat_embeddings_pcs

# Run UMAP and cluster the cells
joint_cds =
  joint_cds %>%
  reduce_dimension() %>%
  cluster_cells()

joint_cds$umap1 = reducedDim(joint_cds, type = "UMAP")[,1]
joint_cds$umap2 = reducedDim(joint_cds, type = "UMAP")[,2]

# Isolate the neural cells -----------------------------------------------------

colData(joint_cds)$joint_clusters = 
  clusters(joint_cds)

colData(joint_cds)$joint_partition = 
  partitions(joint_cds)

# Isolate the populations that are actually present in the DMBA
# For example muscle, liver and cardiomyocytes aren't in there

metadata_for_transfer =
  colData(joint_cds)  %>%
  as.data.frame() %>%
  filter(joint_partition %in% c("2","3","8","13","7","4","6") | 
           joint_clusters %in% c("31")) 


metadata_for_transfer =
  metadata_for_transfer %>%
  left_join(coldata,
            by = c("Cell" = "CellID"))

metadata_for_transfer$Clusters = 
  as.character(metadata_for_transfer$Clusters)

ggplot(metadata_for_transfer %>% filter(sample == "LaManno")) + 
  geom_point(aes(x = umap1, 
                 y = umap2, 
                 color = as.factor(Clusters)), 
             size = 1, 
             stroke = 0) + theme(legend.position = "none") + 
  ggsave("La_Manno_developing_brain/joint_clusters.pdf")
  

#  Transfer labels  using knn -----------------------------------------------------

lamanno_coldata =
  metadata_for_transfer %>%
  filter(sample == "LaManno")

spatial_coldata = 
  metadata_for_transfer %>%
  filter(sample != "LaManno")

# Get nearest neighbors
nn = 
  get.knnx(data =  
             lamanno_coldata %>%
             dplyr::select(umap1, umap2) %>%
             as.matrix(),
           query = spatial_coldata %>%
             dplyr::select(umap1, umap2) %>%
             as.matrix(),
           k = 5)


# Class label transfer ----------------------------------------------------
## Use nearest neighbor indices to get top 5 labels
nn_class = 
  sapply(X = seq(1,dim(nn$nn.index)[2]),
         FUN = function(X){
           lamanno_coldata$Class[nn$nn.index[,X]]
         }
  ) %>%
  as.data.frame()

rownames(nn_class) =
  spatial_coldata$Cell

nn_class_majority = 
  sapply(X = seq(1,dim(nn_class)[1]),
         FUN = function(X){
           tabulated_labels = 
             nn_class[X,] %>%
             unlist() %>%
             table() %>%
             sort(decreasing = T)
           
           tabulated_labels[1] %>%
             names()
         })

nn_class_fraction =
  sapply(X = seq(1,dim(nn_class)[1]), 
         FUN = function(Y){
           sum(nn_class[Y,] == nn_class_majority[Y])/5
         })

# Tissue Label Transfer ---------------------------------------------------
## Use nearest neighbor indices to get top 5 labels
nn_tissue = 
  sapply(X = seq(1,dim(nn$nn.index)[2]),
         FUN = function(X){
           lamanno_coldata$Tissue[nn$nn.index[,X]]
         }
  ) %>%
  as.data.frame()

rownames(nn_tissue) =
  spatial_coldata$Cell

nn_tissue_majority = 
  sapply(X = seq(1,dim(nn_tissue)[1]),
         FUN = function(X){
           tabulated_labels = 
             nn_tissue[X,] %>%
             unlist() %>%
             table() %>%
             sort(decreasing = T)
           
           tabulated_labels[1] %>%
             names()
         })

nn_tissue_fraction =
  sapply(X = seq(1,dim(nn_tissue)[1]), 
         FUN = function(Y){
           sum(nn_tissue[Y,] == nn_tissue_majority[Y])/5
         })


# Punchcard Label Transfer ---------------------------------------------------

## Use nearest neighbor indices to get top 5 labels
nn_punchcard = 
  sapply(X = seq(1,dim(nn$nn.index)[2]),
         FUN = function(X){
           lamanno_coldata$Punchcard[nn$nn.index[,X]]
         }
  ) %>%
  as.data.frame()

rownames(nn_punchcard) =
  spatial_coldata$Cell

nn_punchcard_majority = 
  sapply(X = seq(1,dim(nn_punchcard)[1]),
         FUN = function(X){
           tabulated_labels = 
             nn_punchcard[X,] %>%
             unlist() %>%
             table() %>%
             sort(decreasing = T)
           
           tabulated_labels[1] %>%
             names()
         })

nn_punchcard_fraction =
  sapply(X = seq(1,dim(nn_punchcard)[1]), 
         FUN = function(Y){
           sum(nn_punchcard[Y,] == nn_punchcard_majority[Y])/5
         })


# UMAP Cluster Transfer ---------------------------------------------------

## Use nearest neighbor indices to get top 5 labels
nn_Cluster = 
  sapply(X = seq(1,dim(nn$nn.index)[2]),
         FUN = function(X){
           lamanno_coldata$Clusters[nn$nn.index[,X]]
         }
  ) %>%
  as.data.frame()

rownames(nn_Cluster) =
  spatial_coldata$Cell

nn_Cluster_majority = 
  sapply(X = seq(1,dim(nn_Cluster)[1]),
         FUN = function(X){
           tabulated_labels = 
             nn_Cluster[X,] %>%
             unlist() %>%
             table() %>%
             sort(decreasing = T)
           
           tabulated_labels[1] %>%
             names()
         })

nn_Cluster_fraction =
  sapply(X = seq(1,dim(nn_Cluster)[1]), 
         FUN = function(Y){
           sum(nn_Cluster[Y,] == nn_punchcard_majority[Y])/5
         })


# UMAP1 Label Transfer ---------------------------------------------------

## Use nearest neighbor indices to get top 5 labels
nn_umap1 = 
  sapply(X = seq(1,dim(nn$nn.index)[2]),
         FUN = function(X){
           lamanno_coldata$umap1[nn$nn.index[,X]] 
         }
  ) %>%
  as.data.frame()

rownames(nn_umap1) =
  spatial_coldata$Cell

nn_umap1_avg = 
    apply(nn_umap1,MARGIN = 1,FUN = mean)

# UMAP2 Label Transfer ---------------------------------------------------
## Use nearest neighbor indices to get top 5 labels
nn_umap2 = 
  sapply(X = seq(1,dim(nn$nn.index)[2]),
         FUN = function(X){
           lamanno_coldata$umap2[nn$nn.index[,X]] 
         }
  ) %>%
  as.data.frame()

rownames(nn_umap2) =
  spatial_coldata$Cell

nn_umap2_avg = 
  apply(nn_umap2,MARGIN = 1,FUN = mean)


# Collect the results -----------------------------------------------------
# Take the majority celltype annotation for each cell across all the iterations 


# Calculate the proportion that of nearest neighbors in the majority label across all iterations 

estimates =
  data.frame(Cell = spatial_coldata$Cell,
             nn_class_majority = nn_class_majority,
             nn_class_fraction = nn_class_fraction,
             nn_tissue_majority = nn_tissue_majority,
             nn_tissue_fraction = nn_tissue_fraction,
             nn_punchcard_majority = nn_punchcard_majority,
             nn_punchcard_fraction = nn_punchcard_fraction,
             nn_Cluster_majority = nn_Cluster_majority,
             nn_umap1 = nn_umap1_avg,
             nn_umap2 = nn_umap2_avg)

write.table(x = estimates,
            file = "La_Manno_developing_brain/alignment_coldata_lamanno.tsv",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)

