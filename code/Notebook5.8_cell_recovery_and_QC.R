# Using nuclei estimate counts from 2 layer neural network run on each image
# display the number of nuclei that are recovered by sci-space as sequenced 
# transcriptomes versus the number of estimated nuclei in the slide

# Also quantifying and plotting the number of cells recovered from each slide

# Load startup packages ---------------------------------------------------
suppressPackageStartupMessages({ 
  library(tidyverse)
  library(ggplot2)
  library(tidyr)
  library(viridis)
  library(ggridges)
  library(monocle3)
  
  space_directory = "/Volumes/GoogleDrive/My Drive/sciSpace/"
  setwd(dir=space_directory)
  set.seed(42)
})

spatial_cds = 
  readRDS(file = "Submission_Data/E14_slides/RDS_intermediates/Notebook5.5_E14_spatial_CDS.RDS")

all_image_data = 
  readRDS("Submission_Data/E14_slides/RDS_intermediates/Notebook0_1_images_and_transformations.RDS")


# Set slide labels and colors ---------------------------------------------

slide_colors = 
  c("slide_1D" = "#8DC9DA",
    "slide_1E" = "#4FB961",
    "slide_1F" = "#C01E67",
    "slide_1G" = "grey",
    "slide_2G" = "#C0A630",
    "slide_2H" = "#CD3029",
    "slide_3D" = "#E2AFCE",
    "slide_3F" = "#876F50",
    "slide_3G" = "#8538E8",
    "slide_3H" = "#60E7B2",
    "slide_4A" = "#E7E041",
    "slide_4C" = "#EE7733",
    "slide_4D" =  "#0077BB", 
    "slide_4E" = "#E744DC"
  )

# Supplementary Figure 9 — Panel A ----------------------------------------

# Calculate the median number of genes per cell
colData(spatial_cds)$genes_per_cell = 
  Matrix::colSums(counts(spatial_cds) > 0)
  
print(paste(median(colData(spatial_cds)$genes_per_cell),
            " median genes detected per Cell",
            sep = ""))

colData(spatial_cds) %>%
  as.data.frame() %>%
  dplyr::select(slide_id,
                genes_per_cell,
                n.umi) %>%
  mutate(slide_id = 
           factor(slide_id,
                  levels =  c("Slide 1",
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
                              "Slide 14"))) %>%
  gather("measurement",
         "value",
         c(genes_per_cell,
           n.umi)) %>%
  ggplot() +
  geom_boxplot(aes(x = slide_id,
                   y = log10(value),
                   fill = measurement),
               size = 0.25,
               outlier.stroke = 0,
               outlier.size = 0.75) +
  monocle3:::monocle_theme_opts() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 6,
                                   hjust = 1,
                                   angle = 45),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6)) +
  scale_fill_brewer(palette = "Set1" ,
                    labels = c("Genes", 
                               "UMIs"),
                    name = "") +
  ylab("Log10(Observations)") +
  ggsave("Figures/Figure_Components/Supplement_cell_recovery/umis_genes_per_cell.pdf",
         height = 2,
         width = 2.5)


# Read in estimated nuclei counts from image ------------------------------

read_nuclei_counts = function(file_path,
                              slide_number){
  read.csv(file = file_path,
           row.names = 1) %>%
    t() %>%
    as.data.frame() %>%
    mutate(slide_id = paste("Slide ",slide_number, sep = ""))
}

nuclei_counts_1 = 
  read_nuclei_counts(file_path = "Submission_Data/E14_slides/Images/Image_nuclei_count/csvs/nuclei_locations_1_summary.csv",1)

nuclei_counts_2 = 
  read_nuclei_counts(file_path = "Submission_Data/E14_slides/Images/Image_nuclei_count/csvs/nuclei_locations_2_summary.csv",2)

nuclei_counts_3 = 
  read_nuclei_counts(file_path = "Submission_Data/E14_slides/Images/Image_nuclei_count/csvs/nuclei_locations_3_summary.csv",3)

nuclei_counts_4 = 
  read_nuclei_counts(file_path = "Submission_Data/E14_slides/Images/Image_nuclei_count/csvs/nuclei_locations_4_summary.csv",4)

nuclei_counts_6 = 
  read_nuclei_counts(file_path = "Submission_Data/E14_slides/Images/Image_nuclei_count/csvs/nuclei_locations_6_summary.csv",6)

nuclei_counts_7 = 
  read_nuclei_counts(file_path = "Submission_Data/E14_slides/Images/Image_nuclei_count/csvs/nuclei_locations_7_summary.csv",7)

nuclei_counts_8 = 
  read_nuclei_counts(file_path = "Submission_Data/E14_slides/Images/Image_nuclei_count/csvs/nuclei_locations_8_summary.csv",8)

nuclei_counts_9 = 
  read_nuclei_counts(file_path = "Submission_Data/E14_slides/Images/Image_nuclei_count/csvs/nuclei_locations_9_summary.csv",9)

nuclei_counts_10 = 
  read_nuclei_counts(file_path = "Submission_Data/E14_slides/Images/Image_nuclei_count/csvs/nuclei_locations_10_summary.csv",10)

nuclei_counts_11 = 
  read_nuclei_counts(file_path = "Submission_Data/E14_slides/Images/Image_nuclei_count/csvs/nuclei_locations_11_summary.csv",11)

nuclei_counts_13 = 
  read_nuclei_counts(file_path = "Submission_Data/E14_slides/Images/Image_nuclei_count/csvs/nuclei_locations_13_summary.csv",13)

nuclei_counts_14 = 
  read_nuclei_counts(file_path = "Submission_Data/E14_slides/Images/Image_nuclei_count/csvs/nuclei_locations_14_summary.csv",14)

nuclei_counts_5 =
  nuclei_counts_4

nuclei_counts_12 =
  nuclei_counts_4

# Nuclei counts for slide 5 and slide 12 could not be calculated because the high definition CZI file for the
# image was corrupted or accidentally deleted, respectively
nuclei_counts_5[!is.na(nuclei_counts_5)] <- NA
nuclei_counts_12[!is.na(nuclei_counts_12)] <- NA

nuclei_counts_5 = 
  nuclei_counts_5 %>%
  mutate(slide_id = "Slide 5")

nuclei_counts_12 = 
  nuclei_counts_12 %>%
  mutate(slide_id = "Slide 12")

slide_nuclei_totals =
  rbind(nuclei_counts_1,
      nuclei_counts_2,
      nuclei_counts_3,
      nuclei_counts_4,
      nuclei_counts_5,
      nuclei_counts_6,
      nuclei_counts_7,
      nuclei_counts_8,
      nuclei_counts_9,
      nuclei_counts_10,
      nuclei_counts_11,
      nuclei_counts_12,
      nuclei_counts_13,
      nuclei_counts_14)

slide_nuclei_totals %>%
  pull(Total) %>%
  sum(na.rm = T)/12 * 14

# Supplementary Figure 10 — Panel B ---------------------------------------

colData(spatial_cds) %>%
  as.data.frame() %>%
  group_by(slide_id) %>% 
  summarise(cells_sequenced = n()) %>%
  right_join(slide_nuclei_totals,
             by = "slide_id") %>%
  dplyr::select(slide_id,
                cells_sequenced,
                Total) %>%
  dplyr::mutate(slide_id = 
                  factor(slide_id, 
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
                                    "Slide 14"))) %>%
  gather(key = Total, 
         value = "value",
         -slide_id) %>%
  ggplot() +
  geom_point(aes(x = slide_id,
                 y = log10(value),
                 color = Total),
             shape = 21,
             stroke = 2,
             size = 2) +
  scale_color_brewer(palette = "Set1",
                     name = "Method",
                     labels = c("Sequencing",
                                "Imaging Estimate")) +
  scale_y_continuous(limits = c(3,6),
                     breaks = seq(1,6,1),
                     labels = c("10",
                                "100",
                                "1,000",
                                "10,000",
                                "100,000",
                                "1,000,000")) +
  monocle3:::monocle_theme_opts() +
  ylab("Number of Cells") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 8,
                                   angle = 45,
                                   hjust = 1),
        axis.text.y = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)) +
  ggsave("Figures/Figure_Components/Supplement_cell_recovery/imaging_vs_sequencing.pdf",
         height = 2.5,
         width = 5)


# Calculate the percent of cells sequenced per slide
colData(spatial_cds) %>%
  as.data.frame() %>%
  group_by(slide_id) %>% 
  summarise(cells_sequenced = n()) %>%
  right_join(slide_nuclei_totals,
             by = "slide_id") %>%
  dplyr::select(slide_id,
                cells_sequenced,
                Total) %>%
  dplyr::mutate(slide_id = 
                  factor(slide_id, 
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
                                    "Slide 14"))) %>%
  mutate(percent_captured = cells_sequenced/Total)

# Supplementary Figure 14 — Histograms ---------------------------------------

colData(spatial_cds) %>%
  as.data.frame() %>% 
  group_by(top_spot,
           slide_id) %>%
  summarise(n = n()) %>% 
  filter(slide_id %in%  c("Slide 1",
                          "Slide 2",
                          "Slide 3",
                          "Slide 4",
                          "Slide 5",
                          "Slide 6",
                          "Slide 7")) %>%
  ggplot() +
  geom_histogram(aes(x = n),
                 fill = "grey90",
                      bins = 75,
                 color = "black",
                 size = 0.2) +
  facet_wrap(~slide_id,
             scales = "free_y",
             ncol =1) +
  monocle3:::monocle_theme_opts() +
  theme(legend.position = "none",
        axis.title.y = element_text(size= 8),
        axis.title.x = element_text(size= 8),
        strip.text.x = element_text(size= 8, color = "white"),
        axis.text.x = element_text(size= 6),
        axis.text.y = element_text(size= 6)) +
  xlab("Cells Recovered Per Spot") +
  ylab("Number of Spots") +
  scale_x_log10(breaks = c(1,2,3,4,5,10,25,50,100)) +
  scale_y_continuous(breaks = c(50,125,250,500)) +
  ggsave("Figures/Figure_Components/Supplement_cell_recovery/cells_per_spot_1.pdf",
         height = 8,
         width = 3)


colData(spatial_cds) %>%
  as.data.frame() %>% 
  group_by(top_spot,
           slide_id) %>%
  summarise(n = n()) %>%
  group_by(slide_id) %>%
  mutate(total = sum(n)) %>%
  ungroup() %>%
  filter(!(slide_id %in%  c("Slide 1",
                          "Slide 2",
                          "Slide 3",
                          "Slide 4",
                          "Slide 5",
                          "Slide 6",
                          "Slide 7"))) %>%
  ggplot() +
  geom_histogram(aes(x = n),
                 fill = "grey90",
                 bins = 75,
                 color = "black",
                 size = 0.2) +
  facet_wrap(~slide_id,
             scales = "free_y",
             ncol =1) +
  monocle3:::monocle_theme_opts() +
  theme(legend.position = "none",
        axis.title.y = element_text(size= 8),
        axis.title.x = element_text(size= 8),
        strip.text.x = element_text(size= 8, color = "white"),
        axis.text.x = element_text(size= 6),
        axis.text.y = element_text(size= 6)) +
  xlab("Cells Recovered Per Spot") +
  ylab("Number of Spots") +
  scale_x_log10(breaks = c(1,2,3,4,5,10,25,50,100)) +
  scale_y_continuous(breaks = c(50,100,200,300,400,500)) +
  ggsave("Figures/Figure_Components/Supplement_cell_recovery/cells_per_spot_2.pdf",
         height = 8,
         width = 3)

# Supplementary Figure 14 — Embryo heatmaps ---------------------------------------

# Get the maximum cells in a position to set color scale
max_n =
  colData(spatial_cds) %>%
  as.data.frame() %>%
  group_by(top_spot,
           top_spot,
           coords.x1,
           coords.x2) %>% 
  summarise(n = n()) %>%
  arrange(-n) %>%
  pull(n) %>%
  max()


for( i in seq(1,14)){
colData(spatial_cds) %>%
  as.data.frame() %>%
  filter(slide_id == paste("Slide ",i,sep = "")) %>%
  group_by(slide_id,
           top_spot,
           coords.x1,
           coords.x2) %>% 
  summarise(n = n()) %>%
  ggplot() +
    geom_sf(data = all_image_data$slide_polygon[i][[1]],
            fill = NA,
            color = "black",
            size = 0.25) +
  geom_tile(aes(x = coords.x1,
                y = coords.x2,
                fill = (log10(n))),
            height = 25,
            width = 25,
            alpha = 0.85,
            show.legend = F) +
  scale_fill_viridis(breaks =c(0,1,log10(50),2,log10(200)),
                     labels = c(0, 10,50,100,200),
                     name = "",
                     limits = c(0,log10(max_n))) +
  theme_void() +
  ggsave(paste("Figures/Figure_Components/Supplement_cell_recovery/cells_per_spot_slide",i,".png",sep = ""),
         dpi = 300,
         height = 1,
         width = 1)
}

# Fix Slide 12 -- not in pixel space like the others
colData(spatial_cds) %>%
  as.data.frame() %>%
  filter(slide_id == "Slide 12") %>%
  group_by(slide_id,
           top_spot,
           coords.x1,
           coords.x2) %>% 
  summarise(n = n()) %>%
  ggplot() +
  geom_sf(data = all_image_data$slide_polygon[12][[1]],
          fill = NA,         
          color = "black",
          size = 0.25) +
  geom_tile(aes(x = coords.x1,
                y = coords.x2,
                fill = (log10(n))),
            height = 14,
            width = 14,
            alpha = 0.85,
            show.legend = F) +
  scale_fill_viridis(breaks =c(0,1,log10(50),2,log10(200)),
                     labels = c(0, 10,50,100,200),
                     name = "",
                     limits = c(0,log10(max_n))) +
  theme_void() +  
  ggsave("Figures/Figure_Components/Supplement_cell_recovery/cells_per_spot_slide12.png",
         dpi = 300,
         height = 3,
         width = 3)

# Version for legend
colData(spatial_cds) %>%
  as.data.frame() %>%
  filter(slide_id == "Slide 12") %>%
  group_by(slide_id,
           top_spot,
           coords.x1,
           coords.x2) %>% 
  summarise(n = n()) %>%
  ggplot() +
  geom_sf(data = all_image_data$slide_polygon[12][[1]],
          fill = NA,
          color = "black",
          size = 0.25) +
  geom_tile(aes(x = coords.x1,
                y = coords.x2,
                fill = (log10(n))),
            height = 14,
            width = 14,
            alpha = 0.85,
            show.legend = T) +
  scale_fill_viridis(breaks =c(0,1,log10(50),2,log10(200)),
                     labels = c(0, 10,50,100,200),
                     name = "",
                     limits = c(0,log10(max_n))) +
  theme_void() +
  ggsave("Figures/Figure_Components/Supplement_cell_recovery/cells_per_spot_slide_12_legend.png",
         dpi = 300,
         height = 3,
         width = 3)


# Calculate statistics related to cell recovery ---------------------------
colData(spatial_cds) %>%
  as.data.frame() %>% 
  group_by(top_spot,
           max_slide_id) %>%
  summarise(n = n()) %>% 
  group_by(max_slide_id) %>%
  summarise(med = median(n))


colData(spatial_cds) %>%
  as.data.frame() %>% 
  dplyr::select(top_spot,
                max_slide_id) %>%
  distinct() %>%
  group_by(max_slide_id) %>%
  summarise(n = n()) %>% 
  pull(n) %>%
  mean()




colData(spatial_cds) %>%
  as.data.frame() %>%
  group_by(slide_id) %>%
  summarise(cells_sequenced = n()) %>%
  right_join(slide_nuclei_totals,
             by = "slide_id") %>%
  mutate(percent = 100* cells_sequenced/Total) %>%
  pull(percent) %>%
  mean(na.rm = T)
 
           