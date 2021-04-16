# Comparisons between sci-Space, 10x Visium and slide-seq v2 data
# Code used to make Supplemental Figure 21 G - I 

suppressPackageStartupMessages(
  {library(Seurat)
    library(SeuratData)
    library(ggplot2)
    library(patchwork)
    library(dplyr)
    library(tidyverse)
    library(ggplot2)
    library(tidyr)
    library(viridis)
    library(sp)
    library(sf)
    library(monocle3)

    space_directory = "/Volumes/GoogleDrive/My Drive/sciSpace/"
    setwd(dir=space_directory)}
)

spatial_cds = 
  readRDS("Submission_Data/E14_slides/RDS_intermediates/Notebook6_spatial_cds_anatomy.RDS")

bulked_spatial_cds= 
  readRDS("Submission_Data/E14_slides/RDS_intermediates/Notebook6.5_joint_bulk_CDS.RDS")

colData(bulked_spatial_cds)$experiment = 
  ifelse(colData(bulked_spatial_cds)$max_slide_id %in% c("slide_1G",
                                                         "slide_2G",
                                                         "slide_1D",
                                                         "slide_1F",
                                                         "slide_1E",
                                                         "slide_2H"),
         "sciSpace v1",
         "sciSpace v1.1")

umis_per_pos_space = 
  data.frame(spot_name = colData(bulked_spatial_cds)$position,
             n.umi =  colData(bulked_spatial_cds)$n.umi,
             n.genes = colData(bulked_spatial_cds)$n.genes,
             name = colData(bulked_spatial_cds)$experiment,
             radius = 0.073)

# 10x Visium Brian dataset
tenx.visium <- LoadData("stxBrain", type = "anterior1")

umis_per_pos_10x = 
  data.frame(spot_name = rownames(tenx.visium@meta.data),
             n.umi = tenx.visium@meta.data$nCount_Spatial,
             n.genes = tenx.visium@meta.data$nFeature_Spatial,
             name = "10x Visium",
             radius  = 0.055/2)

#slideseq datase  
slide.seq <- LoadData("ssHippo")

umis_per_pos_slide_seq = 
  data.frame(spot_name = rownames(slide.seq@meta.data),
             n.umi = slide.seq@meta.data$nCount_Spatial,
             n.genes = slide.seq@meta.data$nFeature_Spatial,
             name = "SlideSeq v2",
             radius  = 0.0055)

comparison_df = 
  rbind(umis_per_pos_space,
      umis_per_pos_10x,
      umis_per_pos_slide_seq)

comparison_df$umi_per_area = comparison_df$n.umi/(comparison_df$radius^2 * pi)


comparison_df %>%
  group_by(name) %>%
  summarise(median(umi_per_area))

ggplot(comparison_df) +
  geom_boxplot(aes(x = name,
                   y = log10(umi_per_area)),
               outlier.stroke = 0,
               outlier.size = 0.75) +
  monocle3:::monocle_theme_opts() +
  scale_y_continuous(breaks = seq(5,9,1),
                     labels = c("10,000","100,000","1,000,000",
                                "10,000,000","100,000,000")) +
  ylab(expression(paste("UMIs per ", mm^2,sep = ""))) +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 6),
        axis.title.y = element_text(size = 8)) +
  ggsave("Figures/Figure_Components/Supplement_bulk_aggregation/comparison.pdf",
         height =1.5,
         width = 3)
  

ggplot(comparison_df) +
  geom_boxplot(aes(x = name,
                   y = log10(n.genes)),
               outlier.stroke = 0,
               outlier.size = 0.75) +
  monocle3:::monocle_theme_opts() +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 6),
        axis.title.y = element_text(size = 8)) +
  ylab("Genes per feature") +
  scale_y_continuous(breaks = seq(1,4,1),labels = c("10","100","1,000","10,000")) +
  ggsave("Figures/Figure_Components/Supplement_bulk_aggregation/comparison_genes.pdf",
         height =1.5,
         width = 3)


# Deeper Sequencing -------------------------------------------------------

load.count.matrix = function(mat.path, gene.annotation.path, cell.annotation.path) {
  
  df = read.table(
    mat.path,
    col.names = c("gene.idx", "cell.idx", "count"),
    colClasses = c("integer", "integer", "integer"))
  print("Expression Matrix Read")
  
  gene.annotations = read.table(
    gene.annotation.path,
    col.names = c("id", "gene_short_name"),
    colClasses = c("character", "character"))
  print("Gene Annotations Read")
  
  cell.annotations = read.table(
    cell.annotation.path,
    col.names = c("Cell", "sample"),
    colClasses = c("character", "factor"))
  print("Cell Annotations Read")
  
  rownames(gene.annotations) = gene.annotations$id
  rownames(cell.annotations) = cell.annotations$Cell
  
  # add a dummy cell to ensure that all genes are included in the matrix
  # even if a gene isn't expressed in any cell
  df = rbind(df, data.frame(
    gene.idx = c(1, nrow(gene.annotations)),
    cell.idx = rep(nrow(cell.annotations)+1, 2),
    count = c(1, 1)))
  
  print("Bulding sparseMatrix")
  mat = Matrix::sparseMatrix(i = df$gene.idx, j = df$cell.idx, x = df$count)
  mat = mat[, 1:(ncol(mat)-1)]
  
  
  print(dim(mat))
  rownames(mat) = gene.annotations$id
  colnames(mat) = cell.annotations$Cell
  
  cds = new_cell_data_set(mat,
                          cell_metadata =cell.annotations,
                          gene_metadata = gene.annotations)
  
  colData(cds)$n.umi = Matrix::colSums(exprs(cds))
  
  return(cds)
}

# Read in transcriptome for E14 deeper sequencing -----------------------

deeper_cds = 
  load.count.matrix(mat.path = "Submission_Data/E14_deeper_sequencing/Sequencing/UMI.count.matrix",
                    gene.annotation.path = "Submission_Data/E14_deeper_sequencing/Sequencing/gene.annotations",
                    cell.annotation.path = "Submission_Data/E14_deeper_sequencing/Sequencing/cell.annotations")

# Calculate the amount of sequencing that went towards deeper sequencing of each cell
# This number corresponds to the number of non-duplicated 
reads.per.cell.deeper = 
  read.table("Submission_Data/E14_deeper_sequencing/Sequencing/reads.per.cell.out",
             stringsAsFactors = F,
             col.names = c("Cell","reads"))

reads.per.cell.deeper %>%
  filter(Cell %in% colnames(deeper_cds)) %>%
  pull(reads) %>%
  mean()
#[1] 43535.67

reads.per.cell.shallow = 
  read.table("Submission_Data/E14_slides/Sequencing/experiment_1/reads.per.cell.out",
             stringsAsFactors = F,
             col.names = c("Cell","reads"))

reads.per.cell.shallow %>%
  filter(Cell %in% colnames(spatial_cds)) %>%
  pull(reads) %>%
  mean()
# [1] 4831.538


spatial_cells_for_comp =
  spatial_cds %>%
  colData() %>%
  as.data.frame() %>%
  filter(sample == "1") %>%
  pull(Cell)
  
rbind(deeper_cds %>% 
        colData() %>%
        as.data.frame() %>% 
        dplyr::select(Cell, n.umi)  %>%
        mutate(Sample = "Deeper"),
      spatial_cds[,spatial_cells_for_comp] %>%
        colData() %>%
        as.data.frame() %>% 
        dplyr::select(Cell, n.umi) %>%
        mutate(Sample = "Shallow")) %>%
  group_by(Sample) %>%
  summarise(median_umis = median(n.umi))

rbind(deeper_cds %>% 
        colData() %>%
        as.data.frame() %>% 
        dplyr::select(Cell, n.umi)  %>%
        mutate(Sample = "Deeper"),
      spatial_cds[,spatial_cells_for_comp] %>%
        colData() %>%
        as.data.frame() %>% 
        dplyr::select(Cell, n.umi) %>%
        mutate(Sample = "Shallow")) %>%
  mutate(Sample = factor(Sample,levels = c("Shallow",
                                           "Deeper"))) %>%
  ggplot() +
  geom_boxplot(aes(x = Sample,
                   y = n.umi),
               outlier.stroke = 0,
               outlier.size = 0.75) +
  scale_y_log10(breaks = c(1000,2500,5000,10000,10000,50000)) +
  ylab("RNA UMIs") +
  monocle3:::monocle_theme_opts() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x  = element_text(size = 6),
        axis.text.y  = element_text(size = 6)) +
  ggsave("Figures/Figure_Components/Supplement_bulk_aggregation/deeper_sequencing.pdf",
         height = 1.5,
         width = 1.5)

# Read in hashing data that corresponds to slide and oligo concentration used 
hashTable = 
  read.table("Submission_Data/E14_deeper_sequencing/Sequencing/hashTable.out",
             sep = "\t",
             header = F,
             col.names = c("sample", "Cell", "oligo", "axis", "count"))

# Read in data for UMIs per sci-RNA-seq barcode to determine background cell distribution
rna_umis = 
  read.table("Submission_Data/E14_deeper_sequencing/Sequencing/UMIs.per.cell.barcode",
             col.names = c("sample","Cell","n.umi"))

hashTable = 
  hashTable %>%
  dplyr::group_by(axis, Cell) %>%
  dplyr::mutate(total_within_axis = sum(count)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Cell) %>%
  dplyr::mutate(total_hash = sum(count)) %>%
  dplyr::ungroup() %>%
  left_join(rna_umis,
            by = "Cell") %>%
  mutate(in_cds = Cell %in% colnames(barnyard_cds))

breaks_for_plot = c(1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000, 500000)

# Knee plot with red line showing cutoff chosen and blue line showing background (ambient) molecule cutoff
rna_umis %>%
  dplyr::select(n.umi) %>%
  arrange(desc(n.umi)) %>%
  mutate(rank = dplyr::row_number()) %>%
  ggplot() +
  geom_line(aes(x = rank,
                y = n.umi)) +
  scale_x_log10(breaks = breaks_for_plot) +
  scale_y_log10(breaks = breaks_for_plot) +
  theme_bw() +
  geom_hline(yintercept = 500, 
             color = "red", 
             size = 0.5) +
  geom_hline(yintercept = 10, 
             color = "dodgerblue", 
             size = 0.5)+
  theme(axis.text.x = element_text(angle = 90, size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 8))+ 
  xlab("Number of Cells") +
  ylab("RNA UMIs")  +
  ggsave("~/Desktop/test.pdf")


# Preprocess replicate1 CDS -----------------------------------------------

# Correct the spatial hash table by removing background estimated from debris
hashTable = 
  hashTable %>%
  # perfor background correction separately for each axis
  dplyr::rename(Oligo = oligo,
                Count = count)

hashTable$Oligo = 
  factor(hashTable$Oligo, 
         levels = unique(hashTable$Oligo))

# Background Cells fall below RNA treshold
background_cells = 
  hashTable %>%
  filter(n.umi < 10) %>%
  pull(Cell) %>%
  as.character()

# Test Cells that are in the CDS object
test_cells =
  hashTable %>%
  filter(in_cds) %>%
  pull(Cell) %>%
  as.character()

background_sample_hashes = 
  hashTable %>% 
  filter(Cell %in% background_cells) 

test_cell_hashes = 
  hashTable %>% 
  filter(Cell %in% test_cells) 

hash_df = assign_hash_labels_return_all(test_cell_hashes, 
                                        background_sample_hashes, 
                                        downsample_rate=1)

hashTable = left_join(hashTable,
                      hash_df,
                      by = c("Cell","Oligo"))

hashTable %>% filter(in_cds) %>% as.data.frame()


# Make a matrix representation of the data
oligo_df = 
  hashTable %>%
  ungroup() %>%
  filter(in_cds) %>%
  dplyr::select(Cell, 
                Oligo, 
                adjusted_count)%>%
  spread(key = Oligo, 
         value = adjusted_count,
         fill = 0)

oligo_cell = oligo_df$Cell
oligo_df = oligo_df[,2:ncol(oligo_df)]

# Get the Maximum slide oligo per Cell
max_oligo = 
  colnames(oligo_df)[apply(X = oligo_df, 
                           MARGIN = 1, 
                           FUN = (which.max))]

# Get the Maximum slide oligo UMIs per Cell
max_oligo_val = 
  apply(X = oligo_df, 
        MARGIN = 1, 
        FUN = function(x){
          maximum_col = which.max(x)
          
          x[maximum_col]
        })

# Define the ratio top 2 metric as the ratio of the top hash oligo in a cell
# divided by the second most abundant hash oligo in that cell
ratio_top_2 = 
  apply(X = oligo_df, 
        MARGIN = 1, 
        FUN = function(x){
          maximum_col = which.max(x)
          max_val = x[maximum_col]
          x[maximum_col] = 0
          seoncd_max = which.max(x)
          second_max_val = x[seoncd_max]
          max_val/second_max_val
        })

oligo_total= 
  rowSums(oligo_df)

oligo_df = data.frame(
  Cell = oligo_cell,
  max_oligo_val,
  max_slide_id = max_oligo,
  ratio_top_2 = ratio_top_2,
  oligo_total)

a = (counts(spatial_cds) > 1)

mean(colSums(a)/colSums(counts(spatial_cds)))
