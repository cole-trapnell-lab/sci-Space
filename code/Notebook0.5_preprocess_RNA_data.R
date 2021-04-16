# Preprocess single-cell transcriptomes from E14 embryos:Single cell transcriptomes
# were collected from 14 E14 slides. This notebook takes all the sequenced reads mapping 
# to the transcriptome and assigns them to a single slide based on the predominant 
# slide specific hash oligo associated with that nucleus

# Load startup packages ---------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(Matrix)
  library(monocle3)
  
  space_directory = "/Volumes/GoogleDrive/My Drive/sciSpace/"
  setwd(dir=space_directory)
  
  source("data_from_experiments/bin/chiSq_test_functions.R")
  DelayedArray:::set_verbose_block_processing(TRUE)
  options(DelayedArray.block.size=1000e7)
})

sessionInfo()

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


# Read in transcriptome for spatial CDS -----------------------------------
# These are the main transcriptomic data that constitute the dataset
# Slides 1-6 are experiment 1 and Slides 8-14 are experiment 2
# Slides 1-6 were prepared over 4 PCR plates and sequenced as one batch over 
# 2 Nextseq 75high kits
# Slides 8-14 were prepared over 6 PCR plates and sequenced as three batches of 2 PCR plates 
# with 6 Nextseq 75high kits in total  


# # Slides 1-6 Cutoff of 500 UMIs per cell
spatial_cds_1 = 
  load.count.matrix(mat.path = "Submission_Data/E14_slides/Sequencing/experiment_1/UMI.count.matrix",
                    gene.annotation.path = "Submission_Data/E14_slides/Sequencing/experiment_1/gene.annotations",
                    cell.annotation.path = "Submission_Data/E14_slides/Sequencing/experiment_1/cell.annotations")
# Slides 8-14 (1) Cutoff of 400 UMIs per cell
spatial_cds_2 = 
  load.count.matrix(mat.path = "Submission_Data/E14_slides/Sequencing/experiment_2/A_plates/UMI.count.matrix",
                    gene.annotation.path = "Submission_Data/E14_slides/Sequencing/experiment_2/A_plates/gene.annotations",
                    cell.annotation.path = "Submission_Data/E14_slides/Sequencing/experiment_2/A_plates/cell.annotations")
# Slides 8-14 (2) Cutoff of 400 UMIs per cell
spatial_cds_3 = 
  load.count.matrix(mat.path = "Submission_Data/E14_slides/Sequencing/experiment_2/D_plates/UMI.count.matrix",
                    gene.annotation.path = "Submission_Data/E14_slides/Sequencing/experiment_2/D_plates/gene.annotations",
                    cell.annotation.path = "Submission_Data/E14_slides/Sequencing/experiment_2/D_plates/cell.annotations")
# Slides 8-14 (3) Cutoff of 400 UMIs per cell
spatial_cds_4 = 
  load.count.matrix(mat.path = "Submission_Data/E14_slides/Sequencing/experiment_2/E_plates/UMI.count.matrix",
                    gene.annotation.path = "Submission_Data/E14_slides/Sequencing/experiment_2/E_plates/gene.annotations",
                    cell.annotation.path = "Submission_Data/E14_slides/Sequencing/experiment_2/E_plates/cell.annotations")


# Read in hashing data that corresponds to spatial coordinates
hashTable_1 = 
  read.table("Submission_Data/E14_slides/Sequencing/experiment_1/hashTable.out",
            sep = "\t",
            header = F,
            col.names = c("sample", "Cell", "oligo", "axis", "count"))

hashTable_2 = 
  read.table("Submission_Data/E14_slides/Sequencing/experiment_2/A_plates/hashTable.out",
             sep = "\t",
             header = F,
             col.names = c("sample", "Cell", "oligo", "axis", "count"))

hashTable_3 = 
  read.table("Submission_Data/E14_slides/Sequencing/experiment_2/D_plates/hashTable.out",
             sep = "\t",
             header = F,
             col.names = c("sample", "Cell", "oligo", "axis", "count"))

hashTable_4 = 
  read.table("Submission_Data/E14_slides/Sequencing/experiment_2/E_plates/hashTable.out",
             sep = "\t",
             header = F,
             col.names = c("sample", "Cell", "oligo", "axis", "count"))



# Read in the number of RNA UMI molecules per cell 
rna_umis_1 = 
  read.table("Submission_Data/E14_slides/Sequencing/experiment_1/UMIs.per.cell.barcode",
             col.names = c("sample","Cell","n.umi")) %>%
  dplyr::select(-"sample")

rna_umis_2 = 
  read.table("Submission_Data/E14_slides/Sequencing/experiment_2/A_plates/UMIs.per.cell.barcode",
             col.names = c("sample","Cell","n.umi")) %>%
  dplyr::select(-"sample")

rna_umis_3 = 
  read.table("Submission_Data/E14_slides/Sequencing/experiment_2/D_plates/UMIs.per.cell.barcode",
             col.names = c("sample","Cell","n.umi")) %>%
  dplyr::select(-"sample")

rna_umis_4 = 
  read.table("Submission_Data/E14_slides/Sequencing/experiment_2/E_plates/UMIs.per.cell.barcode",
             col.names = c("sample","Cell","n.umi")) %>%
  dplyr::select(-"sample")

# Make lists for all experiments/batches of sciSpace preps 
cds.list = 
  list(spatial_cds_1,
       spatial_cds_2,
       spatial_cds_3,
       spatial_cds_4)

hashTable.list = 
  list(hashTable_1,
       hashTable_2,
       hashTable_3,
       hashTable_4)

rnas.per.cell.list = 
  list(rna_umis_1,
       rna_umis_2,
       rna_umis_3,
       rna_umis_4)

# Tabulate per cell information and add to the spatial hash table
for(i in seq(1,4)){
   hashTable.list[[i]] =
     hashTable.list[[i]] %>%
     dplyr::group_by(axis, Cell) %>%
     dplyr::mutate(total_within_axis = sum(count)) %>%
     dplyr::ungroup() %>%
     dplyr::group_by(Cell) %>%
     dplyr::mutate(total_hash = sum(count)) %>%
     dplyr::ungroup() %>%
     left_join(rnas.per.cell.list[[i]],
               by = "Cell") %>%
     mutate(in_cds = Cell %in% colnames(cds.list[[i]]))
  }


# Make RNA knee plot to set background UMI treshold for estimating background hashes

# The RNA cutoff was set at 400 or 500 
# The background count treshold was set at 10

breaks_for_plot = c(1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000, 250000, 500000)

do.call(rbind,
        rnas.per.cell.list) %>%
  dplyr::select(n.umi) %>%
  arrange(desc(n.umi)) %>%
  mutate(rank = dplyr::row_number()) %>%
  ggplot() +
  geom_line(aes(x = rank,
                y = n.umi)) +
  scale_x_log10(breaks = breaks_for_plot) +
  scale_y_log10(breaks = breaks_for_plot) +
  theme_bw() +
  geom_hline(yintercept = 400, 
             color = "red", 
             size = 0.5) +
  theme(axis.text.x = element_text(angle = 90, size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 8))+ 
  xlab("Number of Cells") +
  ylab("RNA UMIs") +
  ggsave(filename = "Submission_Data/E14_slides/Auxillary_Figures/Notebook0_RNA_kneeplot.pdf",
         height = 2,
         width = 4) 

# Correct the spatial hash table by removing background estimated from debris --------

# The hashTable is in a "long" form (as opposed to wide) where every row is a Cell/hashOligo
# combination with the number of UMIs recovered for that hash Oligo. Additionally, based on the
# design of the experiment hashes can denote labeling of different things.

# In this exeriment the "Axis" denotes each independent labeling step 
# Axis 1 -- Slide hash oligo -- in solution when the nuclei are permabilized on a slide
# Axis 2 -- Sector hash oligo -- found in contiguous square pathches on the space-grid
# Axis 3 -- Spot hash oligo -- sequence for each spot. 

# The spot in combination with the sector defines a specific position 
corrected_hash_table.list = 
  lapply(hashTable.list, 
         FUN = function(x){
           x %>%
             # perfor background correction separately for each axis
             group_by(axis) %>% 
             dplyr::rename(Oligo = oligo,
                           Count = count) %>%
             nest() %>% 
             mutate(corrected_df = 
                      purrr::map(.x = data, 
                                 .f = function(hash_table_subset){
               hash_table_subset$Oligo = 
                 factor(hash_table_subset$Oligo, 
                        levels = unique(hash_table_subset$Oligo))
               
               # Background Cells fall below RNA treshold
               background_cells = 
                 hash_table_subset %>%
                 filter(n.umi < 10) %>%
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
             unnest(cols = c("corrected_df"))
         })

# Recover the slide information from the spatial hashing data -------------
# Slide Oligo was labeled as Axis 1 

slide_oligo.list = 
  lapply(corrected_hash_table.list,
         function(x){
           # Make a matrix representation of the data
           slide_oligo = 
             x %>%
             ungroup() %>%
             filter(axis == 1,
                    in_cds) %>%
             dplyr::select(Cell, 
                           Oligo, 
                           adjusted_count)%>%
             spread(key = Oligo, 
                    value = adjusted_count,
                    fill = 0)
           
           slide_oligo_cell = slide_oligo$Cell
           slide_oligo = slide_oligo[,2:ncol(slide_oligo)]
           
           # Get the Maximum slide oligo per Cell
           max_slide_oligo = 
             colnames(slide_oligo)[apply(X = slide_oligo, 
                                         MARGIN = 1, 
                                         FUN = (which.max))]
           
           max_slide_val = 
             apply(X = slide_oligo, 
                   MARGIN = 1, 
                   FUN = function(x){
                     maximum_col = which.max(x)
                     
                     x[maximum_col]
                   })
           
           ratio_top_2 = 
             apply(X = slide_oligo, 
                   MARGIN = 1, 
                   FUN = function(x){
                     maximum_col = which.max(x)
                     max_val = x[maximum_col]
                     x[maximum_col] = 0
                     seoncd_max = which.max(x)
                     second_max_val = x[seoncd_max]
                     max_val/second_max_val
                   })
           
           slide_oligo_total = 
             rowSums(slide_oligo)
           
           data.frame(Cell = slide_oligo_cell,
             max_slide_val = max_slide_val,
             max_slide_id = max_slide_oligo,
             ratio_top_2_slide = ratio_top_2,
             slide_oligo_total)
         })

# Isolate those cells that are well marked by a singular slide oligo --------

# Marked by a singular slide oligo is defined here as exceeding 5 fold enrichment 
# over the next most abundant oligo and containing more than 300 hash UMIs 
marked_cells = 
  do.call(rbind,slide_oligo.list) %>% 
  filter(slide_oligo_total >= 5, 
         ratio_top_2_slide > 5) %>%
  as.data.frame() %>%
  pull(Cell)


# Make Supplemental Figure 8A,B -------------------------------------------
# Supplemental Figure 8B 

ggplot(data = 
         do.call(rbind,
                 rnas.per.cell.list) %>%
         mutate(hashed = (Cell %in% marked_cells)) %>%
         dplyr::select(n.umi,hashed) %>%
         arrange(desc(n.umi)) %>%
         mutate(rank = dplyr::row_number()),
       aes(x = rank,
           y = n.umi),
       size = 0.15) +
  geom_line() +
  geom_hline(yintercept = 400,
             color = "red") +
  scale_x_log10(breaks = breaks_for_plot) +
  scale_y_log10(breaks = breaks_for_plot) +
  theme_bw() +
  scale_color_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(hjust = 1,
                                   angle = 90)) +
  xlab("Number of Cells") +
  ylab("RNA UMIs") +
  ggsave(filename = "Figures/Figure_Components/Supplement_experiment_QC/RNA_kneeplot.png",
         dpi = 450,
         height = 2.5,
         width = 4)

# Supplemental Figure 8B 
do.call(rbind,
        rnas.per.cell.list) %>%
  mutate(hashed = (Cell %in% marked_cells)) %>%
  filter(n.umi > 400) %>%
  ggplot() +
  geom_violin(aes(x = hashed,
                  y = n.umi),
              fill = "grey80",
              color = "black",
              size = 0.5) +
  geom_hline(yintercept = 400,
             color = "red") +
  scale_y_log10(breaks = breaks_for_plot) +
  
  monocle3:::monocle_theme_opts() +
  xlab("Passed Slide Hash Filter") +
  ylab("RNA UMIs") +
  ggsave("Figures/Figure_Components/Supplement_experiment_QC/hashed_violin.pdf",
         height = 2,
         width = 3)

rna.v.hash.df = 
  inner_join(do.call(rbind,
                    rnas.per.cell.list),
            do.call(rbind,
                    hash_table.list) %>%
              dplyr::select(total_hash,Cell) %>%
              distinct(),
            by = "Cell")

ggplot(rna.v.hash.df) +
  geom_point(aes(x = log10(n.umi),
               y = log10(total_hash)),
             size = 0.25,
             stroke = 0) +
  geom_density2d(aes(x = log10(n.umi),
                     y = log10(total_hash))) +
  monocle3:::monocle_theme_opts() +
  scale_y_continuous(breaks = seq(1,10)) +
  scale_x_continuous(breaks = seq(1,10)) +
  ylab("Log10(Hash UMIs)") + 
  xlab("Log10(RNA UMIs)") +
  ggsave("Figures/Figure_Components/Reviewer_Figures/hash_v_rna.png",
         height = 4,
         width = 4)
  
# Join the experiments into a single cds
spatial_cds =
  combine_cds(cds_list = cds.list)

# > dim(spatial_cds)
# [1]  52636 232415

# Subset the CDS by these singularly marked cells
spatial_cds = spatial_cds[,colData(spatial_cds)$Cell %in% marked_cells]

# > dim(spatial_cds)
# [1]  52636 151490

# Add a new slide id name for plotting and representation purposes

slide_id_key = 
  data_frame(max_slide_id = 
               c("slide_1D",
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
                 "slide_4C",
                 "slide_4D",
                 "slide_4E"),
             slide_id = 
               factor(x = 
                        c("Slide 1",
                          "Slide 2",
                          "Slide 3",
                          "Slide 4",
                          "Slide 5",
                          "Slide 6",
                          "Slide 7",
                          "Slide 8",
                          "Slide 9",
                          "Slide 10",
                          "Slide 11",
                          "Slide 12",
                          "Slide 13",
                          "Slide 14"),
                      levels = c("Slide 1",
                                 "Slide 2",
                                 "Slide 3",
                                 "Slide 4",
                                 "Slide 5",
                                 "Slide 6",
                                 "Slide 7",
                                 "Slide 8",
                                 "Slide 9",
                                 "Slide 10",
                                 "Slide 11",
                                 "Slide 12",
                                 "Slide 13",
                                 "Slide 14")))


# Add which slide wach cell originate from into each cell's metadata in the CDS
coldata_spatial_cds =
  colData(spatial_cds) %>%
  as.data.frame() %>%
  left_join(do.call(rbind,slide_oligo.list),
            by = "Cell") %>%
  left_join(slide_id_key,
            by = "max_slide_id")

# Slides 1-6 are experiment 1 and Slides 8-14 are experiment 2
coldata_spatial_cds$experiment = 
  ifelse(coldata_spatial_cds$sample == "1",
         "experiment_1",
         "experiment_2")

rownames(coldata_spatial_cds) =
  coldata_spatial_cds$Cell

count_matrix_spatial_cds = counts(spatial_cds)
colnames(count_matrix_spatial_cds) = 
  stringr::str_replace(string = colnames(count_matrix_spatial_cds),
                       pattern = "_[1-4]$", 
                       replacement = "")

# Check to make sure the colData and counts match up by cell name
identical(colnames(count_matrix_spatial_cds),
          rownames(coldata_spatial_cds))

spatial_cds = 
  new_cell_data_set(expression_data = count_matrix_spatial_cds,
                    cell_metadata = coldata_spatial_cds,
                    gene_metadata = rowData(spatial_cds))
  

# Write out CDS that contains the new slide information -------------------
saveRDS(object =spatial_cds,
        file = "Submission_Data/E14_slides/RDS_intermediates/Notebook0_E14_spatial_CDS.RDS")

# Write out mtx file that contains the new slide information
writeMM(obj = counts(spatial_cds),
        file = "Submission_Data/E14_slides/RDS_intermediates/Notebook0_E14_count_matrix.mtx")

# Write out an RDS of the hashTable
saveRDS(object = hashTable.list,
        file = "Submission_Data/E14_slides/RDS_intermediates/Notebook0_hash_table.RDS")


