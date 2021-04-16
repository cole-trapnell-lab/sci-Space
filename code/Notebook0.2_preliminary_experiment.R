# Control experiment: Hash oligos were were loaded onto agarose backed slides at multiple concentrations
# Each slide contained two embryos sections and recieved a single hash oligo (in solution) to mark the slide 
# Each section then recieved oligo from the agarose backed slide. Each concentration oligo has a specific sequence

# Load startup packages ---------------------------------------------------
suppressPackageStartupMessages({ 
  library(tidyverse)
  library(ggplot2)
  library(viridis)
  library(ggridges)
  library(ggpubr)
  library(Matrix)
  library(monocle3)
  
  space_directory = "/Volumes/GoogleDrive/My Drive/sciSpace/"
  setwd(dir=space_directory)
  
  source("data_from_experiments/bin/chiSq_test_functions.R")
  
  # Pass TRUE if you want to see progress output on some of Monocle 3's operations
  DelayedArray:::set_verbose_block_processing(TRUE)
  # Passing a higher value will make some computations faster but use more memory. 
  # Adjust with caution!
  options(DelayedArray.block.size=1000e6)
  set.seed(42)
})

# Function to make cds object from sparse matrix format
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

# Read in transcriptome for E13-E16 CDS replicate 1 -----------------------
# Cutoff in knee plot was 500 RNA UMIs per cell
cds_1 = 
  load.count.matrix(mat.path = "Submission_Data/E13-E16/Experiment1/UMI.count.matrix",
                    gene.annotation.path = "Submission_Data/E13-E16/Experiment1/gene.annotations",
                    cell.annotation.path = "Submission_Data/E13-E16/Experiment1/cell.annotations")

# Read in hashing data that corresponds to slide and oligo concentration used 
hashTable_1 = 
  read.table("Submission_Data/E13-E16/Experiment1/hashTable.out",
             sep = "\t",
             header = F,
             col.names = c("sample", "Cell", "oligo", "axis", "count"))

# Read in data for UMIs per sci-RNA-seq barcode to determine background cell distribution
rna_umis_1 = 
  read.table("Submission_Data/E13-E16/Experiment1/UMIs.per.cell.barcode",
             col.names = c("sample","Cell","n.umi"))

hashTable_1 = 
  hashTable_1 %>%
  dplyr::group_by(axis, Cell) %>%
  dplyr::mutate(total_within_axis = sum(count)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Cell) %>%
  dplyr::mutate(total_hash = sum(count)) %>%
  dplyr::ungroup() %>%
  left_join(rna_umis_1,
            by = "Cell") %>%
  mutate(in_cds = Cell %in% colnames(cds_1))

breaks_for_plot = c(1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000, 500000)

# Knee plot with red line showing cutoff chosen and blue line showing background (ambient) molecule cutoff
rna_umis_1 %>%
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
  ylab("RNA UMIs") 

# Preprocess replicate1 CDS -----------------------------------------------

# Correct the spatial hash table by removing background estimated from debris
corrected_hash_table_1 = 
  hashTable_1 %>%
  # perfor background correction separately for each axis
  group_by(axis) %>% 
  dplyr::rename(Oligo = oligo,
                Count = count) %>%
  nest() %>% 
  mutate(corrected_df = purrr::map(data, .f = function(hash_table_subset){
    hash_table_subset$Oligo = 
      factor(hash_table_subset$Oligo, 
             levels = unique(hash_table_subset$Oligo))
    
    # Background Cells fall below RNA treshold
    background_cells = 
      hash_table_subset %>%
      filter(n.umi < 15) %>%
      pull(Cell) %>%
      as.character()
    
    # Test Cells that are in the CDS object
    test_cells =
      hash_table_subset %>%
      filter(in_cds) %>%
      pull(Cell) %>%
      as.character()
    
    background_sample_hashes = 
      hash_table_subset %>% 
      filter(Cell %in% background_cells) 
    
    test_cell_hashes = 
      hash_table_subset %>% 
      filter(Cell %in% test_cells) 
    
    hash_df = assign_hash_labels_return_all(test_cell_hashes, 
                                            background_sample_hashes, 
                                            downsample_rate=1)
    
    left_join(hash_table_subset,
              hash_df,
              by = c("Cell","Oligo")) 
  })) %>%
  dplyr::select(-data) %>%
  unnest()

corrected_hash_table_1 %>% head()


# Make a matrix representation of the data
slide_oligo_1 = 
  corrected_hash_table_1 %>%
  ungroup() %>%
  filter(axis == 2,
         in_cds) %>%
  dplyr::select(Cell, 
                Oligo, 
                adjusted_count)%>%
  spread(key = Oligo, 
         value = adjusted_count,
         fill = 0)

slide_oligo_cell_1 = slide_oligo_1$Cell
slide_oligo_1 = slide_oligo_1[,2:ncol(slide_oligo_1)]

# Get the Maximum slide oligo per Cell
max_slide_oligo_1 = 
  colnames(slide_oligo_1)[apply(X = slide_oligo_1, 
                              MARGIN = 1, 
                              FUN = (which.max))]

# Get the Maximum slide oligo UMIs per Cell
max_slide_val_1 = 
  apply(X = slide_oligo_1, 
        MARGIN = 1, 
        FUN = function(x){
          maximum_col = which.max(x)
          
          x[maximum_col]
        })

# Define the ratio top 2 metric as the ratio of the top hash oligo in a cell
# divided by the second most abundant hash oligo in that cell
ratio_top_2_1 = 
  apply(X = slide_oligo_1, 
        MARGIN = 1, 
        FUN = function(x){
          maximum_col = which.max(x)
          max_val = x[maximum_col]
          x[maximum_col] = 0
          seoncd_max = which.max(x)
          second_max_val = x[seoncd_max]
          max_val/second_max_val
        })

slide_oligo_total_1 = 
  rowSums(slide_oligo_1)

slide_df_1 = data.frame(
  Cell = slide_oligo_cell_1,
  max_slide_val = max_slide_val_1,
  max_slide_id = max_slide_oligo_1,
  ratio_top_2_slide = ratio_top_2_1,
  slide_oligo_total_1)

rownames(slide_df_1) =
  slide_df_1$Cell

# Isolate those cells that are well marked by a singular slide oligo
marked_cells_1 =
  hashTable_1 %>%
  dplyr::select(Cell,
                total_hash) %>% 
  distinct() %>%
  right_join(slide_df_1,
             by = "Cell") %>% 
  filter(slide_oligo_total_1 >= 10, 
         ratio_top_2_slide > 5, 
         total_hash > 300) %>%
  as.data.frame() %>%
  pull(Cell)

# Check that all the cells are in the cds
((colnames(cds_1) %in% marked_cells_1) %>% sum) / length(marked_cells_1)
cds_1 = cds_1[,colnames(cds_1) %in% marked_cells_1]

coldata_cds_1 = 
  colData(cds_1) %>%
  as.data.frame() %>%
  left_join(slide_df_1,
            by = "Cell")


colData(cds_1)$max_slide_id = 
  coldata_cds_1$max_slide_id

colData(cds_1)$slide_id =
  stringr::str_split_fixed(colData(cds_1)$max_slide_id,
                           "_",
                           2)[,1]

colData(cds_1)$stage =
  stringr::str_split_fixed(colData(cds_1)$max_slide_id,
                           "_",
                           2)[,2]

colData(cds_1)$slide_set =
  ifelse(colData(cds_1)$slide_id %in% c("slide1","slide3","slide5","slide7"), "Set 1","Set 2")
  
# Filter on slides that recieved either the 50uM and 25uM oligo -- Set 1
# or the slides that recieved either the 20uM and 10M oligo -- Set 2

hashTable_1_set1 =
  hashTable_1 %>%
  filter(axis == 1,
         Cell %in% marked_cells_1) %>%
  dplyr::select(Cell,oligo,count) %>%
  right_join(colData(cds_1) %>%
               as.data.frame(),
             by = "Cell") %>% 
  filter(slide_set == "Set 1") %>%
  filter(oligo %in% c("50uM","25uM")) 

hashTable_1_set2 =
  hashTable_1 %>%
  filter(axis == 1,
         Cell %in% marked_cells_1) %>%
  dplyr::select(Cell,oligo,count) %>%
  right_join(colData(cds_1) %>%
               as.data.frame(),
             by = "Cell") %>% 
  filter(slide_set == "Set 2") %>%
  filter(oligo %in% c("10uM","20uM"))


# Supplementary Figure 4 — Panel G ----------------------------------------
rbind(hashTable_1_set2, hashTable_1_set1) %>% 
  ggplot() +
  geom_boxplot(aes(x = oligo,
                   y = log10(count)),
               fill = "grey90",
               color = "black",
               outlier.size = 1,
               outlier.stroke = 0,
               size = 0.25) +
  monocle3:::monocle_theme_opts() + 
  xlab("Concentration of Spotted Oligo") +
  ylab("log10(hash UMIs)") +
  theme(axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text.x  = element_text(size = 6),
        axis.text.y = element_text(size = 6)) +
  ggsave("Figures/Figure_Components/Supplement_test_transfer/spotted_oligo_concentration.pdf",
         height = 2,
         width = 2)
  


# Read in transcriptome for E13-E16 CDS replicate 2 -----------------------

# Cutoff of was 500 RNA UMIs per cell
cds_2 = 
  load.count.matrix(mat.path = "Submission_Data/E13-E16/Experiment2/UMI.count.matrix",
                    gene.annotation.path = "Submission_Data/E13-E16/Experiment2/gene.annotations",
                    cell.annotation.path = "Submission_Data/E13-E16/Experiment2/cell.annotations")

# Read in hashing data that corresponds to slide and oligo concentration used 
hashTable_2 = 
  read.table("Submission_Data/E13-E16/Experiment2/hashTable.out",
             sep = "\t",
             header = F,
             col.names = c("sample", "Cell", "oligo", "axis", "count"))

# Read in data for UMIs per sci-RNA-seq barcode to determine background cell distribution
rna_umis_2 = 
  read.table("Submission_Data/E13-E16/Experiment2/UMIs.per.cell.barcode",
             col.names = c("sample","Cell","n.umi"))


hashTable_2 = 
  hashTable_2 %>%
  dplyr::group_by(axis, Cell) %>%
  dplyr::mutate(total_within_axis = sum(count)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Cell) %>%
  dplyr::mutate(total_hash = sum(count)) %>%
  dplyr::ungroup() %>%
  left_join(rna_umis_2,
            by = "Cell") %>%
  mutate(in_cds = Cell %in% colnames(cds_2))

breaks_for_plot = c(1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000, 500000)

# Knee plot with red line showing cutoff chosen and blue line showing background (ambient) molecule cutoff
rna_umis_2 %>%
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
  geom_hline(yintercept = 15, 
             color = "dodgerblue", 
             size = 0.5)+
  theme(axis.text.x = element_text(angle = 90, size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 8))+ 
  xlab("Number of Cells") +
  ylab("RNA UMIs") 



# Preprocess replicate2 CDS -------------------------------------------------------


# Correct the spatial hash table by removing background estimated from debris
corrected_hash_table_2 = 
  hashTable_2 %>%
  # perfor background correction separately for each axis
  group_by(axis) %>% 
  dplyr::rename(Oligo = oligo,
                Count = count) %>%
  nest() %>% 
  mutate(corrected_df = purrr::map(data, .f = function(hash_table_subset){
    hash_table_subset$Oligo = 
      factor(hash_table_subset$Oligo, 
             levels = unique(hash_table_subset$Oligo))
    
    # Background Cells fall below RNA treshold
    background_cells = 
      hash_table_subset %>%
      filter(n.umi < 15) %>%
      pull(Cell) %>%
      as.character()
    
    # Test Cells that are in the CDS object
    test_cells =
      hash_table_subset %>%
      filter(in_cds) %>%
      pull(Cell) %>%
      as.character()
    
    background_sample_hashes = 
      hash_table_subset %>% 
      filter(Cell %in% background_cells) 
    
    test_cell_hashes = 
      hash_table_subset %>% 
      filter(Cell %in% test_cells) 
    
    hash_df = assign_hash_labels_return_all(test_cell_hashes, 
                                            background_sample_hashes, 
                                            downsample_rate=1)
    
    left_join(hash_table_subset,
              hash_df,
              by = c("Cell","Oligo")) 
  })) %>%
  dplyr::select(-data) %>%
  unnest()

corrected_hash_table_2 %>% head()


# Make a matrix representation of the data
slide_oligo_2 = 
  corrected_hash_table_2 %>%
  ungroup() %>%
  filter(axis == 2,
         in_cds) %>%
  dplyr::select(Cell, 
                Oligo, 
                adjusted_count)%>%
  spread(key = Oligo, 
         value = adjusted_count,
         fill = 0)

slide_oligo_cell_2 = slide_oligo_2$Cell
slide_oligo_2 = slide_oligo_2[,2:ncol(slide_oligo_2)]

# Get the Maximum slide oligo per Cell
max_slide_oligo_2 = 
  colnames(slide_oligo_2)[apply(X = slide_oligo_2, 
                                MARGIN = 1, 
                                FUN = (which.max))]

# Get the Maximum slide oligo UMIs per Cell
max_slide_val_2 = 
  apply(X = slide_oligo_2, 
        MARGIN = 1, 
        FUN = function(x){
          maximum_col = which.max(x)
          
          x[maximum_col]
        })

# Define the ratio top 2 metric as the ratio of the top hash oligo in a cell
# divided by the second most abundant hash oligo in that cell
ratio_top_2_2 = 
  apply(X = slide_oligo_2, 
        MARGIN = 1, 
        FUN = function(x){
          maximum_col = which.max(x)
          max_val = x[maximum_col]
          x[maximum_col] = 0
          seoncd_max = which.max(x)
          second_max_val = x[seoncd_max]
          max_val/second_max_val
        })

slide_oligo_total_2 = 
  rowSums(slide_oligo_2)

slide_df_2 = data.frame(
  Cell = slide_oligo_cell_2,
  max_slide_val = max_slide_val_2,
  max_slide_id = max_slide_oligo_2,
  ratio_top_2_slide = ratio_top_2_2,
  slide_oligo_total_2)

rownames(slide_df_2) =
  slide_df_2$Cell

# Isolate those cells that are well marked by a singular slide oligo
marked_cells_2 =
  hashTable_2 %>%
  dplyr::select(Cell,
                total_hash) %>% 
  distinct() %>%
  right_join(slide_df_2,
             by = "Cell") %>% 
  filter(slide_oligo_total_2 >= 3, 
         ratio_top_2_slide > 5, 
         total_hash > 30) %>%
  as.data.frame() %>%
  pull(Cell)

# Check that all the cells are in the cds
((colnames(cds_2) %in% marked_cells_2) %>% sum) / length(marked_cells_2)
cds_2 = cds_2[,colnames(cds_2) %in% marked_cells_2]

coldata_cds_2 = 
  colData(cds_2) %>%
  as.data.frame() %>%
  left_join(slide_df_2,
            by = "Cell")


colData(cds_2)$max_slide_id = 
  coldata_cds_2$max_slide_id


colData(cds_2)$slide_id =
  stringr::str_split_fixed(colData(cds_2)$max_slide_id,
                           "_",
                           2)[,1]

colData(cds_2)$stage =
  stringr::str_split_fixed(colData(cds_2)$max_slide_id,
                           "_",
                           2)[,2]



# Combine the two biological replicates into one CDS ----------------------

# joint count matrix
joint_count_mat = cbind(counts(cds_1),
                        counts(cds_2))

# joint column data
joint_coldata = rbind(colData(cds_1) %>%
                        as.data.frame() %>%
                        dplyr::select(-slide_set),
                      colData(cds_2))

# Create a joint CDS
joint_cds = new_cell_data_set(expression_data = joint_count_mat,
                              cell_metadata = joint_coldata,
                              gene_metadata = rowData(cds_2) %>% as.data.frame())


colData(joint_cds)$stage =
  factor(colData(joint_cds)$stage, 
         levels = c("E16","E15","E14","E13"))

colData(joint_cds)$sample = 
  ifelse(colData(joint_cds)$sample == "space60",
         "Replicate 1",
         "Replicate 2")



# Process Joint CDS -------------------------------------------------------

joint_cds = 
  joint_cds %>%
  estimate_size_factors() %>%
  preprocess_cds(method = "PCA",
                 num_dim = 50) %>%
  align_cds(preprocess_method = "PCA",
            alignment_group = "sample") %>%
  reduce_dimension(umap.fast_sgd = F)

colData(joint_cds)$umap1 = 
  reducedDim(joint_cds,
             type = "UMAP")[,1]

colData(joint_cds)$umap2 = 
  reducedDim(joint_cds,
             type = "UMAP")[,2]

# Supplementary Figure 4 — Panel A ----------------------------------------

colData(joint_cds) %>% 
  as.data.frame() %>%
  ggplot() +
  geom_boxplot(aes(x = stage,
                   y = log10(n.umi),
                   fill = stage),
               outlier.size = 1,
               outlier.stroke = 0,
               size = 0.25) +
  facet_wrap(~sample) +
  scale_fill_viridis_d() +
  coord_flip() +
  xlab("Embryonic Stage") +
  ylab("Log10(RNA UMIs)") +
  monocle3:::monocle_theme_opts() +
  theme(axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.text.x  = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        strip.text.x = element_text(size = 8),
        legend.position = "none") +
  ggsave("Figures/Figure_Components/Supplement_test_transfer/umis_per_stage.pdf", 
         height = 2., 
         width = 2.5)

  

# Supplementary Figure 4 — Panel B ----------------------------------------
colData(joint_cds) %>% 
  as.data.frame() %>%
  filter(sample == "Replicate 1") %>%
  group_by(slide_id,
           sample,
           stage) %>%
  summarise(n = n()) %>%
  ggplot() +
  geom_bar(aes(x = reorder(slide_id,n),
               y = n,
               fill = stage),
           color = "black",
           size = 0.25,
           width = .9,
           stat = "identity") +
  facet_wrap(~sample) +
  scale_fill_viridis_d() +
  coord_flip() +
  ylab("Cells Recovered") +
  monocle3:::monocle_theme_opts() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.text.x  = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        strip.text.x = element_text(size = 8),
        legend.position = "none") +
  ggsave("Figures/Figure_Components/Supplement_test_transfer/cells_per_slide_rep1.pdf", 
         height = 2.25, 
         width = 1.5)

colData(joint_cds) %>% 
  as.data.frame() %>%
  filter(sample == "Replicate 2") %>%
  group_by(slide_id,
           sample,
           stage) %>%
  summarise(n = n()) %>%
  ggplot() +
  geom_bar(aes(x = reorder(slide_id,n),
               y = n,
               fill = stage),
           color = "black",
           size = 0.25,
           width = .9,
           stat = "identity") +
  facet_wrap(~sample) +
  scale_fill_viridis_d() +
  coord_flip() +
  ylab("Cells Recovered") +
  monocle3:::monocle_theme_opts() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.text.x  = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        strip.text.x = element_text(size = 8),
        legend.position = "none") +
  ggsave("Figures/Figure_Components/Supplement_test_transfer/cells_per_slide_rep2.pdf", 
         height = 2.25, 
         width = 1.5)


# Supplementary Figure 4 — Panel C ----------------------------------------

ggplot(colData(joint_cds) %>% 
         as.data.frame()) +
  geom_point(aes(x = umap1, 
                 y = umap2, 
                 color = stage), 
             stroke = 0, 
             size = .25) +
  monocle3:::monocle_theme_opts() +
  scale_color_viridis_d(name = "Stage") +
  xlab("umap 1") +
  ylab("umap 2") +
  guides(guides(colour = guide_legend(override.aes = list(size=1)))) +
  theme_void()+
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Supplement_test_transfer/stage_umap.pdf", 
         height = 2, 
         width = 2)


# Supplementary Figure 4 — Panel D ----------------------------------------

colData(joint_cds) %>% 
  as.data.frame() %>%
  ggplot() +
  geom_point(aes(x = umap1, 
                 y = umap2, 
                 color = sample), 
             stroke = 0, 
             size = .25) +
  monocle3:::monocle_theme_opts() +
  scale_color_brewer(palette = "Set1",
                     name = "Replicate",
                     labels = c("Rep1",
                                "Rep2")) +
  xlab("umap 1") +
  ylab("umap 2") +
  guides(guides(colour = guide_legend(override.aes = list(size=1)))) +
  theme_void()+
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Supplement_test_transfer/replicate_umap.pdf", 
         height = 2, 
         width = 2)


# Supplementary Figure 4 — Panel E ----------------------------------------


colData(joint_cds)$ttn_expression = 
  counts(joint_cds)[which(rowData(joint_cds)$gene_short_name == "Ttn"),]

colData(joint_cds) %>% 
  as.data.frame() %>%
  ggplot() +
  geom_point(aes(x = umap1, 
                 y = umap2, 
                 color = log10(ttn_expression)), 
             stroke = 0, 
             size = .25) +
  monocle3:::monocle_theme_opts() +
  theme_void()+
  scale_color_viridis_c(na.value = "grey90") +
  theme(legend.position = "none",
        strip.text.x = element_blank()) +
  ggsave("Figures/Figure_Components/Supplement_test_transfer/ttn_umap.pdf", 
         height = 2, 
         width = 2)

# Supplementary Figure 4 — Panel F ----------------------------------------

colData(joint_cds)$sample = 
  ifelse(colData(joint_cds)$sample == "Replicate 1",
         "Rep1",
         "Rep2")

counts_per_gene = 
  colData(joint_cds) %>% 
  as.data.frame() %>%
  group_by(sample,stage) %>%
  nest() %>%
  mutate(counts_per_gene = purrr::map(data, .f = function(pdata_subset, cds_filtered) {
    cds_subset = joint_cds[,pdata_subset$Cell]
    cds_subset = estimate_size_factors(cds_subset)
    norm_exprs_mat = 
      Matrix::t(Matrix::t(counts(cds_subset))/size_factors(cds_subset))
    totals = Matrix::rowSums(norm_exprs_mat)
    tibble(gene_id = rownames(exprs(cds_subset)),
           totals = totals )
  },cds)) %>%
  dplyr::select(-data)

counts_per_gene = 
  counts_per_gene %>%
  unnest(counts_per_gene)

counts_per_gene %>%
  spread(key = sample, value = totals,fill = 0)  %>% 
  ggplot() +
  geom_point(aes(x = log10(Rep1 + 1), 
                 y = log10(Rep2 + 1)),
             stroke = 0, 
             size = .25, 
             color = "grey80") +
  geom_abline(intercept = 0, 
              slope = 1, 
              size = .05) +
  geom_smooth(aes(x = log10(Rep1 + 1), 
                  y = log10(Rep2 + 1)),
              method = "lm", 
              color = "red", 
              size = .25 ) +
  stat_cor(aes(x = log10(Rep1 + 1), 
               y = log10(Rep2 + 1)),
           size=1.5, 
           label.y.npc = 0.9, 
           color="grey31", 
           method = "pearson") +
  monocle:::monocle_theme_opts() +
  facet_wrap(~stage) +
  theme(axis.title.y = element_text(size = 6),
        axis.title.x = element_text(size = 6),
        axis.text.x  = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        strip.text.x = element_text(size = 6),
        legend.position = "none") +
  xlab("Replicate1 log10(UMIs + 1)") +
  ylab("Replicate2 log10(UMIs + 1)") +
  ggsave("Figures/Figure_Components/Supplement_test_transfer/replicate_correlation.png",
         height = 2, 
         width = 2,
         units = "in",
         dpi = 300)
 
