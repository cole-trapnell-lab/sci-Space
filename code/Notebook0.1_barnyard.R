# R analysis file to analyse control experiment. 
# HEK293T (human) and NIH3T3 (mouse) cells were grown on a glass cover slip 
# and then hashed with a single oligo transferred from a thin layer of agarose.  

# Load startup packages ---------------------------------------------------

suppressPackageStartupMessages({ 
  library(tidyverse)
  library(ggplot2)
  library(viridis)
  library(ggridges)
  library(ggpubr)
  library(Matrix)
  library(monocle3)
  source("~/Google Drive File Stream/My Drive/sciSpace/data_from_experiments/bin/chiSq_test_functions.R")
  
  space_directory = "~/Google Drive File Stream/My Drive/sciSpace/"
  setwd(dir=space_directory)
  
  # Pass TRUE if you want to see progress output on some of Monocle 3's operations
  DelayedArray:::set_verbose_block_processing(TRUE)
  # Passing a higher value will make some computations faster but use more memory. 
  # Adjust with caution!
  options(DelayedArray.block.size=1000e6)
  set.seed(42)
})


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
barnyard_cds = 
  load.count.matrix(mat.path = "Submission_Data/control_experiment/barnyard/UMI.count.matrix",
                    gene.annotation.path = "Submission_Data/control_experiment/barnyard/gene.annotations",
                    cell.annotation.path = "Submission_Data/control_experiment/barnyard/cell.annotations")

# Read in hashing data that corresponds to slide and oligo concentration used 
hashTable = 
  read.table("Submission_Data/control_experiment/barnyard/hashTable.out",
             sep = "\t",
             header = F,
             col.names = c("sample", "Cell", "oligo", "axis", "count"))

# Read in data for UMIs per sci-RNA-seq barcode to determine background cell distribution
rna_umis = 
  read.table("Submission_Data/control_experiment/barnyard/UMIs.per.cell.barcode",
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
  ylab("RNA UMIs")


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

# Background cells are those that have a very low number of RNA molecules
background_cells = 
  hashTable %>%
  filter(n.umi < 10) %>%
  pull(Cell) %>%
  as.character()

# Test cells are those cells in the CDS object
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

# Make a dataframe with all of these computed values
oligo_df = data.frame(
  Cell = oligo_cell,
  max_oligo_val,
  max_slide_id = max_oligo,
  ratio_top_2 = ratio_top_2,
  oligo_total)


# Barnyard ----------------------------------------------------------------

# Add relevant gene statistics to CDS -- Count human and mouse genes per cell
human.genes = rownames(barnyard_cds)[grepl("^ENSG", rowData(barnyard_cds)$id)]
mouse.genes = rownames(barnyard_cds)[grepl("^ENSMUSG", rowData(barnyard_cds)$id)]
colData(barnyard_cds)$n_human_umi = Matrix::colSums(exprs(barnyard_cds)[rownames(barnyard_cds) %in% human.genes,])
colData(barnyard_cds)$n_mouse_umi = Matrix::colSums(exprs(barnyard_cds)[rownames(barnyard_cds) %in% mouse.genes,])

# Make dataframe for plotting
barnyard_coldata =
  colData(barnyard_cds) %>%
  as.data.frame() %>%
  left_join(oligo_df,
            by = "Cell") 


barnyard_coldata %>%
  mutate(ratio = n_mouse_umi/n.umi) %>% 
  ggplot() +
  geom_histogram(aes(x = ratio,
                     fill = max_slide_id),
                      bins = 100,
                 color = "black",
                 size = 0.15) +
  theme_classic() +
  theme(legend.position = "none") +
  xlab("Ratio (Mouse UMIs / Total UMIs)") +
  scale_fill_brewer(palette = "Set1") +
  ggsave("~/Google Drive File Stream/My Drive/sciSpace/Figures/Figure_Components/Supplement_experiment_QC/barnyard_purity.pdf",
         height = 2,
         width = 2.5)



# Supplemental Figure 5E - Barnyard scatterplot
ggplot(barnyard_coldata %>%
         filter(ratio_top_2 > 5)) +
  geom_point(aes(x = n_human_umi,
                 y = n_mouse_umi),
             size = 0.125) +
  monocle3:::monocle_theme_opts() +
  theme(legend.position = "none") +
  xlab("Human UMIs") +
  ylab("Mouse UMIs") +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) +
  ggsave("~/Google Drive File Stream/My Drive/sciSpace/Figures/Figure_Components/Supplement_experiment_QC/barnyard.pdf",
         height = 2.75,
         width = 2.75)

# Supplemental Figure 5F - Histogram
barnyard_coldata %>%
  mutate(ratio = n_mouse_umi/n.umi,
         call = ifelse(ratio <.1, 
                       "Human",
                       ifelse(ratio > 0.9,
                              "Mouse",
                              "Doublet"))) %>%
  mutate(ratio_species = ifelse(call == "Human",
                                100* (1 - n_human_umi / n.umi),
                                100* (1 - n_mouse_umi/ n.umi))) %>%
  filter(call !="Doublet" ) %>%
  ggplot() +
  geom_histogram(aes(x = ratio_species,
                     fill = call),
                 bins = 100,
                 color = "black",
                 size = 0.15) +
  scale_fill_manual(values = c("Human" = "#af8dc3",
                               "Mouse" = "#7fbf7b")) +
  monocle3:::monocle_theme_opts() +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = seq(0,100,10),
                     limits = c(0,50),
                     labels  = paste(seq(0,100,10),"%",sep = "")) +
  xlab("Contaminating UMIs (%)")  +
  ylab("Cells") +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) +
  ggsave("~/Google Drive File Stream/My Drive/sciSpace/Figures/Figure_Components/Supplement_experiment_QC/barnyard_purity.pdf",
         height = 2.75,
         width = 2.75)

# Calculate the percentage of contaminating transcripts
barnyard_coldata %>%
  mutate(ratio = n_mouse_umi/n.umi,
         call = ifelse(ratio <.1, 
                       "Human",
                       ifelse(ratio > 0.95,
                              "Mouse",
                              "Doublet"))) %>%
  mutate(ratio_species = ifelse(call == "Human",
                                (1 - n_human_umi / n.umi),
                                (1 - n_mouse_umi/ n.umi))) %>%
  filter(call != "Doublet") %>%
  group_by(call) %>%
  summarise(mean_contamination = mean(ratio_species))
  
  