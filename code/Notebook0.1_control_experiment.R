# Control experiment: Hash oligos from 3 slides were dissolved and used to label HEK293T cells
# A subset of cells were also stained with a uniform mix of sector oligos
# Only hash oligos were sequenced in this experiment

# Load startup packages ---------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggridges)
  library(ggpubr)
  library(devtools)
  library(monocle3)
  
  space_directory = "/Volumes/GoogleDrive/My Drive/sciSpace/"
  setwd(dir=space_directory)
  set.seed(42)
})


# Read in hashTable and append experiment specific metadata ---------------

hashTable = read.table("Submission_Data/control_experiment/Sequencing/hashTable.out", sep = "\t", header = F)
colnames(hashTable) = c("sample", "cell", "oligo", "axis", "count")

cell_meta_data = 
  str_split_fixed(hashTable$cell, 
                  pattern = "_",
                  n = 3)

hashTable$rt_well = cell_meta_data[,3] 


labels = data.frame(rt_well = paste("RT_",(seq(0:95) + 96),sep = ""),
           experiment = c(rep(1,24),
                         rep(2,24),
                         rep(3,24),
                         rep(4,24)))

experiment_pretty = 
  data_frame(experiment = seq(1,4),
             experiment_pretty = 
               c("Uniform Conc.",
                 "Slide Rep1",
                 "Slide Rep2",
                 "Slide Rep3"))

hashTable  =
  left_join(hashTable, labels) %>%
  left_join(experiment_pretty)

hashTable$experiment = as.character(hashTable$experiment)

# Get cell names for cells that passed hash treshold
cells_passing_treshold =
  hashTable %>%
  dplyr::group_by(cell, experiment) %>%
  dplyr::summarise(total = sum(count)) %>%
  tidyr::drop_na() %>% 
  dplyr::filter(total > 1000) %>%
  dplyr::pull(cell) %>%
  unique() %>%
  as.character()

hashTable = 
  hashTable %>%
  dplyr::filter(cell %in% cells_passing_treshold)


# Join hash information with oligo layout ---------------------------------
oligo_plate_layout = 
  read.table(file = "Submission_Data/Oligo_layout/oligo_plate_layout.tsv", 
             sep = '\t', 
             header = TRUE)

# Colors to label sectors 
colors = c("#41e05e", 
           "#fc9cf6",
           "#f9e0a9", 
           "#ed6f99",
           "#1956d1",
           "#eaaf85",
           "#dd7f7a",
           "#27e8d7",
           "#8e68ff",
           "#7dd7d8",
           "#3a63a8",
           "#f2b529",
           "#bbfc9f",
           "#20188e",
           "#9b08cc",
           "#e0503a")

sectors = oligo_plate_layout$sector %>% unique()

colors_df = data.frame(color = colors, 
                       sector = sectors)

temp = 
  hashTable %>% 
  filter(axis == 2,
         experiment == 1)  %>%
  mutate(sector_ordering = 
           str_split_fixed(oligo,"sector",2)[,2] %>%
           as.numeric()) %>%
  # color df from quadrants analysis.R
  left_join(colors_df, by = c("oligo" = "sector"))  


col <- as.character(temp$color)
names(col) <- as.character(temp$oligo)

# Number of Hash Oligos per cell for when each 
# hash oligo is mixed at an equal concentration
hashTable %>% 
  filter(axis == 2,
         experiment == 1) %>%
  mutate(sector_ordering = 
           str_split_fixed(oligo,"sector",2)[,2] %>%
           as.numeric()) %>%
  ggplot() +
  geom_density_ridges(aes(y = reorder(oligo,sector_ordering),
                          x = log10(count),
                          fill = oligo)) + 
  monocle:::monocle_theme_opts()  +
  scale_fill_manual(values = col) +
  theme(legend.position = "none") +
  facet_wrap(~experiment_pretty) +
  ylab("Sector Oligo") +
  xlab("Log10(Hashes per Cell)")


# Supplemental Figure 5 — Panel B -----------------------------------------

hashTable %>% 
  filter(axis == 3,
         experiment %in% c(2,4)) %>%
  mutate(experiment = ifelse(experiment  == 2, "slide1","slide2" )) %>%
  group_by(oligo, experiment) %>%
  dplyr::summarise(median_spot = median(count)) %>% 
  spread(key = experiment, value = median_spot) %>%
  ggplot(aes(y = slide2,
             x = slide1)) +
  geom_jitter(stroke = 0,
              size = 0.5) +
  scale_x_continuous(breaks = seq(0,100,2)) +
  scale_y_continuous(breaks = seq(0,100,2)) +
  geom_smooth(method = "lm")+
  stat_cor(size=2, color="grey31",label.x.npc = "left", label.y.npc = "top") +
  monocle3:::monocle_theme_opts()  +
  theme(axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        axis.title.y  = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "none") +
  ylab("UMIs per Spot (Rep2)") +
  xlab("UMIs per Spot (Rep1)") +
  ggsave("Figures/Figure_Components/Supplement_grid_QC/correlation_between_slide_per_spot.pdf",
         height = 2.5,
         width = 2.5)

# Supplementary Figure 5 — Panel C ----------------------------------------
hashTable %>% 
  filter(axis == 2,
         experiment != 1) %>%
  mutate(sector_ordering = 
           str_split_fixed(oligo,"sector",2)[,2] %>%
           as.numeric()) %>% 
  ggplot() +
  geom_density_ridges(aes(y = reorder(oligo,sector_ordering),
                          x = log10(count),
                          fill = oligo)) + 
  monocle:::monocle_theme_opts()  +
  facet_wrap(~experiment_pretty) +
  scale_fill_manual(values = col) +
  theme(axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        axis.title.y  = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "none") +  
  ylab("Sector Oligo") +
  xlab("Log10(Hash UMIs per Cell)") +
  xlim(1,3.2) +
  ggsave("Figures/Figure_Components/Supplement_grid_QC/hashes_from_slide.pdf",
         height = 3, width = 3.5)
