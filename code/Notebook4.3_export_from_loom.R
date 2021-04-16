# Intermediate script used to extract the developing mouse brain atlas data
# from a Loom object 

# This dataset can be downloaded at -- http://mousebrain.org/downloads.html

# G. La Manno, K. Siletti, A. Furlan, D. Gyllborg, E. Vinsland, 
# C. M. Langseth, I. Khven, A. Johnsson, M. Nilsson, P. Lönnerberg, 
# S. Linnarsson, Molecular architecture of the developing mouse brain. 
# Cold Spring Harbor Laboratory (2020), p. 2020.07.02.184051.

# https://www.biorxiv.org/content/10.1101/2020.07.02.184051v1


# Load startup packages ---------------------------------------------------
suppressPackageStartupMessages({ 
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(purrr)
  library(loomR)
  library(monocle3)
})

lfile <- connect(filename = "~/Downloads/dev_all.loom", 
                 mode = "r+",
                 skip.validate = TRUE)

# Get the count matrix and write it out as in Matrix Marker format
full.matrix <- 
  Matrix::Matrix(sparse = T,) %>%
  as(Class = "dgCMatrix")

Matrix::writeMM( lfile[["matrix"]][, ] %>% Matrix::Matrix(sparse = T),
                  file = "Published_Data/Developing_Brain_Atlas/count_matrix.mtx")
  
# Get the appropriate column data
coldata_colnames = lfile[["col_attrs"]]$names

coldata_colnames = paste("col_attrs/",
                         coldata_colnames,
                         sep = "")

# Select columns that aren't tensors (depth of 1)
meta_data_single_column = 
  sapply(seq(1,length(coldata_colnames)),
       function(x){
         (lfile[[coldata_colnames[x]]]$dims %>% length()) == 1
       })

coldata =
  lapply(coldata_colnames[meta_data_single_column],
       function(x){
         lfile[[x]][]
       })

# Bind these columns into a dataframe
coldata = 
  do.call(cbind,
          coldata) %>%
  as.data.frame()

colnames(coldata) = 
  lfile[["col_attrs"]]$names[meta_data_single_column]

# Get the UMAP embeddings from the original study
UMAP_coordinates = t(lfile[["col_attrs/UMAP"]][,])
coldata$UMAP1 = UMAP_coordinates[,1]
coldata$UMAP2 = UMAP_coordinates[,2]

lfile[["col_attrs"]]$names[!meta_data_single_column]

# Get gene metadata
rowdata_colnames = lfile[["row_attrs"]]$names
rowdata_colnames = paste("row_attrs/",
                         rowdata_colnames,
                         sep = "")

# Select columns that aren't tensors (depth of 1)
row_data_single_column = 
  sapply(seq(1,length(rowdata_colnames)),
         function(x){
           (lfile[[rowdata_colnames[x]]]$dims %>% length()) == 1
         })


rowdata =
  lapply(rowdata_colnames[row_data_single_column],
         function(x){
           lfile[[x]][]
         })

# Bind these columns into a dataframe
rowdata = 
  do.call(cbind,
          rowdata) %>%
  as.data.frame()

colnames(rowdata) = 
  lfile[["row_attrs"]]$names[row_data_single_column]

# Set the rownames of the gene metadata data frame to the Ensmbl IDs
rownames(rowdata) = rowdata$Accession


write.table(coldata,
            file = "Published_Data/Developing_Brain_Atlas/coldata.tsv",
            sep = "\t", 
            col.names = T, 
            row.names = T,
            quote = F)

write.table(rowdata,
            file = "Published_Data/Developing_Brain_Atlas/rowdata.tsv",
            sep = "\t", 
            col.names = T, 
            row.names = T,
            quote = F)
