# Calculate the angular distance between cells that mapped to various positions
# on the grid. If a pair of cells (of the same celltyoe ) are further from
# each other in spatial position, are the transcriptomes of the cells more disimilar?

# Load startup packages ---------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(tidyr)
  library(viridis)
  library(purrr)
  library(ggridges)
  library(monocle3)
  library(broom)
  space_directory = "/Volumes/GoogleDrive/My Drive/sciSpace/"
  setwd(dir=space_directory)
  
  set.seed(42)
})


spatial_cds =
  readRDS(file = "Submission_Data/E14_slides/RDS_intermediates/Notebook6_spatial_cds_anatomy.RDS")


slides_to_analyse =
  c("slide_1D",
    "slide_2H",
    "slide_3F",
    "slide_4A",
    "slide_4E")

angular.dist_df =
  colData(spatial_cds) %>%
  as.data.frame() %>%
  filter(max_slide_id %in% slides_to_analyse) %>%
  group_by(max_slide_id,
           final_cluster_label) %>%
  filter(n() > 100) %>%
  nest() %>%
  mutate(ang.dist =
           purrr::map(.x = data,
                      .f = function(coldata_subset){

                        count_mat =
                          spatial_cds[,coldata_subset$Cell] %>%
                          monocle3:::normalize_expr_data(norm_method = "log",
                                                         pseudo_count = 1) %>%
                          as.matrix()


                        temp = lsa::cosine(count_mat)
                        temp =
                          (2/pi*acos(temp))^2 %>%
                          as.data.frame() %>%
                          rownames_to_column(var = "Cell1")

                        gather(temp,
                               key = "Cell2",
                               value = "angular.distance",
                               -Cell1)

                      })) %>%
  dplyr::select(-data) %>%
  unnest(cols = c("ang.dist"))


physical_distance = 
  colData(spatial_cds) %>%
  as.data.frame() %>%
  filter(max_slide_id %in% slides_to_analyse) %>%
  group_by(max_slide_id, 
           final_cluster_label) %>%
  filter(n() > 100) %>%
  nest() %>%
  mutate(physical.distance =
           purrr::map(.x = data,
                      .f = function(coldata_subset){
                        temp = 
                          combn(x = coldata_subset$Cell,m = 2) %>% 
                          t() %>%
                          as.data.frame()
                          
                        colnames(temp) = c("Cell1", "Cell2")
                        
                        temp = 
                          left_join(temp,
                                    coldata_subset %>%
                                      dplyr::select(Cell1 = Cell,
                                                    Row_1 = Row,
                                                    Col_1 = Col),
                                    by = "Cell1")
                        
                        temp = 
                          left_join(temp,
                                    coldata_subset %>%
                                      dplyr::select(Cell2 = Cell,
                                                    Row_2 = Row,
                                                    Col_2 = Col),
                                    by = "Cell2") %>%
                          mutate(Row_1 = Row_1 %>% as.character() %>% as.numeric(),
                                 Row_2 = Row_2 %>% as.character() %>% as.numeric(),
                                 Col_1 = Col_1 %>% as.character() %>% as.numeric(),
                                 Col_2 = Col_2 %>% as.character() %>% as.numeric())
                        

                        temp$distance = sqrt((temp$Row_1 - temp$Row_2)^2 +
                                               (temp$Col_1 - temp$Col_2)^2)
                        
                        temp
                      })
         ) %>%
  dplyr::select(-data) %>%
  unnest(cols = c("physical.distance"))


dist_df = inner_join(physical_distance,angular.dist_df,by = c("Cell1","Cell2","final_cluster_label","max_slide_id"))


dist_df$distance_breaks =
  cut(dist_df$distance,
      breaks = c(-1,0.75,1.01,2.01,3.01,4.01,84),
      labels = seq(0,5))  %>%
  as.character()


# Fit a linear model asking whether the model angular distance ~ physical_distance has
# a significant Beta
distance_linear_model=
  dist_df %>%
  group_by(max_slide_id,
           final_cluster_label) %>%
  nest() %>%
  mutate(distance_model =
           purrr::map(.x = data,
                      .f = function(distance_df){
                        model = lm(formula = angular.distance ~ distance,
                                   data = distance_df)
                        tidy(model)
                      })) %>%
  dplyr::select(-data) %>%
  unnest(distance_model) %>%
  filter(term == "distance")

distance_linear_model$q.value = 
  p.adjust(p = distance_linear_model$p.value)



# For plotting purposes have stars denoting level of significance
significance_df =
  distance_linear_model %>%
  filter(q.value < 1e-2) %>%
  mutate(significance =
           ifelse(q.value < 1e-10,
                  "***",
                  ifelse(q.value < 1e-5,
                         "**",
                         "*")))


# Supplemental Figure 27 — Panel A — Slide 1 --------------------------------


dist_df %>%
  filter(max_slide_id == "slide_1D") %>%
  ggplot() +
  geom_boxplot(aes(x = final_cluster_label,
                   y = angular.distance,
                   fill = distance_breaks),
               color = "black",
               size = 0.25,
               outlier.stroke = 0,
               outlier.size = 0.75) +
  geom_text(data =
              dist_df %>%
              left_join(significance_df,
                        by = c("max_slide_id",
                               "final_cluster_label")) %>%
              filter(max_slide_id == "slide_1D") %>%
              dplyr::select(final_cluster_label,
                            significance) %>%
              distinct() %>%
              drop_na(),
            aes(x = final_cluster_label,
                label = significance,
                y = 0.9)) +
  monocle3:::monocle_theme_opts() +
  scale_fill_manual(values =
                      viridis(option = "magma",
                              n = 8)[3:8] %>%
                      rev()) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   size = 8),
        axis.title.x = element_blank(),
        plot.margin = unit(c(1,1,1,1), "cm")) +
  scale_y_continuous(limits = c(0.2,0.9))+
  ylab("Angular Distance")+
  ggsave("Figures/Figure_Components/Supplemental_Angular_Distance/ang_dist_spatial_bin_1D.png",
         height = 5,
         width = 8,
         limitsize = F,
         dpi = 300,
         bg = "transparent")

# Supplemental Figure 27 — Panel B — Slide 14 --------------------------------

dist_df %>%
  filter(max_slide_id == "slide_4E") %>%
  ggplot() +
  geom_boxplot(aes(x = final_cluster_label,
                   y = angular.distance,
                   fill = distance_breaks),
               color = "black",
               size = 0.25,
               outlier.stroke = 0,
               outlier.size = 0.75) +
  geom_text(data =
              dist_df %>%
              left_join(significance_df,
                        by = c("max_slide_id",
                               "final_cluster_label")) %>%
              filter(max_slide_id == "slide_4E") %>%
              dplyr::select(final_cluster_label,
                            significance) %>%
              distinct() %>%
              drop_na(),
            aes(x = final_cluster_label,
                label = significance,
                y = 0.9)) +
  monocle3:::monocle_theme_opts() +
  scale_fill_manual(values =
                      viridis(option = "magma",
                              n = 8)[3:8] %>%
                      rev()) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   size = 8),
        axis.title.x = element_blank(),
        plot.margin = unit(c(1,1,1,1), "cm")) +
  scale_y_continuous(limits = c(0.2,0.9))+
  ylab("Angular Distance")+
  ggsave("Figures/Figure_Components/Supplemental_Angular_Distance/ang_dist_spatial_bin_4E.png",
         height = 5,
         width = 8,
         limitsize = F,
         dpi = 300,
         bg = "transparent")



# Figure 4A ----------------------------------------------------------------


ggplot(dist_df %>%
         filter(max_slide_id == "slide_4E",
                final_cluster_label %in% c("Neuron",
                             "Endothelial Cells",
                             "Radial glia"))) +
  geom_boxplot(aes(x = final_cluster_label,
                   y = angular.distance,
                   fill = distance_breaks),
               color = "black",
               size = 0.25,
               outlier.stroke = 0,
               outlier.size = 0.75) +
  monocle3:::monocle_theme_opts() +
  scale_fill_manual(values =
                      viridis(option = "magma",
                              n = 8)[3:8] %>%
                      rev()) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        plot.background = element_blank()) +
  scale_y_continuous(limits = c(0.4,0.8))+
  ylab("Angular Distance")+
  ggsave("Figures/Figure_Components/Figure5/ang_dist.pdf",
         height = 1.75,
         width = 2.25 ,
         limitsize = F)



# Gaussian Process --------------------------------------------------------
# Reviewer suggestion to fit gaussian process to the data

library(GauPro)

gp = 
  dist_df %>%
  group_by(max_slide_id,
           final_cluster_label,
           distance) %>%
  summarise(mean_ang_dist = mean(angular.distance)) %>%
  group_by(max_slide_id,
           final_cluster_label) %>%
  nest() %>%
  mutate(gp = purrr::map(.x = data,
                         .f = function(x){
                           gp = GauPro(x$distance, x$mean_ang_dist, 
                                       parallel=FALSE, 
                                       se.fit = T)
                           
                           gp_pred = gp$predict(seq(1,50,0.1),
                                                    se.fit = T)
                           
                           data.frame(gp_p =gp_pred[,1],
                                      se_lower = gp_pred[,2],
                                      se_upper = gp_pred[,3],
                                      dist = seq(1,50,0.1))
                         })) %>%
  dplyr::select(-data) %>%
  unnest()



# Filter for Gaussian Process Fits that didn't fail
ggplot(gp %>%
         filter(max_slide_id == "slide_4E") %>%
         filter(!final_cluster_label %in% c("Schwann Cells",
                                           "Lateral Plate Mesoderm"))) + 
  geom_line(aes(x = dist * 0.22,
                 y = gp_p),
             color = "blue") +
  geom_line(aes(x = dist * 0.22,
                y = gp_p + 2*se_upper),
            linetype = "dashed",
            color = "grey") +
  geom_line(aes(x = dist * 0.22,
                y = gp_p - 2*se_upper),
            color = "grey",
            linetype = "dashed") +
  geom_vline(xintercept = 1.1,
             size = 0.2,
             color = "red") +
  theme_classic() +
  facet_wrap(~final_cluster_label,scales = "free_y") +
  xlab("Physical distance (mm)") +
  ylab("Angular Distance") +
  scale_x_continuous(breaks = seq(1,to = 10,1)) +
  theme(strip.text = element_text(size = 4),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title  = element_text(size = 8))+
  ggsave("Figures/Reviewer_Figures/gaussian_process_slide_4E.png",
         width = 6,
         height = 4)

