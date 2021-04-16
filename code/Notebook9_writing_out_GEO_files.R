# Write out the count matrix, cell metadata and gene metadat for GEO 
# and select desired columns for output

# Load startup packages ---------------------------------------------------
suppressPackageStartupMessages({ 
  library(monocle3)

  space_directory = "/Volumes/GoogleDrive/My Drive/sciSpace/"
  setwd(dir=space_directory)

  # Pass TRUE if you want to see progress output on some of Monocle 3's operations
  DelayedArray:::set_verbose_block_processing(TRUE)
  # Passing a higher value will make some computations faster but use more memory. 
  # Adjust with caution!
  options(DelayedArray.block.size=1000e6)
  set.seed(42)
})

spatial_cds = 
  readRDS("Submission_Data/E14_slides/RDS_intermediates/Notebook6.01_spatial_cds_anatomy.RDS")

saveRDS(object =spatial_cds,
        file = "Submission_Data/E14_slides/RDS_intermediates/Notebook0_E14_spatial_CDS.RDS")

column_data = 
  colData(spatial_cds) %>%
  as.data.frame()

columns_to_keep = 
  c("Cell",
    "sample",
    "experiment",
    "max_slide_id",
    "slide_id",
    "top_spot",
    "Row",
    "Col",
    "coords.x1",
    "coords.x2",
    "n.umi",
    "log.n.umi",
    "Size_Factor",
    "umap1",
    "umap2",
    "cluster",
    "partition",
    "sub_cluster",
    "final_cluster_label",
    "manual_annotation_2",
    "anatomical_annotation",
    "brain_region",
    "g1s_score",
    "g2m_score",
    "proliferation_index",
    "lamanno_Punchcard",
    "lamanno_Tissue")

column_data = column_data[,columns_to_keep]

rowdata = rowData(spatial_cds) %>% as.data.frame()

# Write out mtx file that contains the new slide information
writeMM(obj = counts(spatial_cds),
        file = "abb9536_Science_Submission/GEO_submission/sciSpace_count_matrix.mtx")

write.table(x = column_data,
            file = "abb9536_Science_Submission/GEO_submission/sciSpace_cell_metadata.tsv",
            sep = "\t",
            quote = F,
            col.names = T,
            row.names = T)

write.table(x = rowdata,
            file = "abb9536_Science_Submission/GEO_submission/sciSpace_gene_metadata.tsv",
            sep = "\t",
            quote = F,
            col.names = T,
            row.names = T)

