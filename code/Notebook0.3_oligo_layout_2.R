# Oligo layout: Hash oligos were arrayed in a predetermined format
# This notebook transforms the 19 384-well plates which contained aliquotted hash oligos 
# into their grid coordinates on the space-grid

# This notebook constructs a second oligo layout grid optimized such that a spot oligo never appears in
# an adjacent section. Used to label Slides 8 - 14

# Load startup packages ---------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

setwd("/Volumes/GoogleDrive/My Drive/sciSpace/Submission_Data/Oligo_layout_2/")

# Constructing the oligo array ---------------------------------------------------
# Layout of how the printer arrays the spatial grid
oligo_plate_layout =read.delim("plate_layout.tsv", header = F,sep = "\t")

rownames(oligo_plate_layout) = as.character(1:84)
colnames(oligo_plate_layout) = as.character(1:84)


oligo_plate_layout = as.data.frame.table(as.matrix(oligo_plate_layout))
colnames(oligo_plate_layout) = c("Row","Col","Oligo")

rownames(oligo_plate_layout) = oligo_plate_layout$Oligo

oligo_plate_layout$plate = stringr::str_split_fixed(string = oligo_plate_layout$Oligo,pattern = "_",n = 2)[,1]
oligo_plate_layout$well = stringr::str_split_fixed(string = oligo_plate_layout$Oligo,pattern = "_",n = 2)[,2]
oligo_plate_layout$row_well = (substr(x = oligo_plate_layout$well,start = 1,stop = 1))
oligo_plate_layout$col_well = (substr(x = oligo_plate_layout$well,start = 2,stop = 3))




# Read in all the hash oligos available
all_hashes = read.table(file = "all_hashes.tsv", sep= "\t",header = F)
colnames(all_hashes) = c("plate","position","RT_seq","seq","barcode")
all_hashes = 
  all_hashes %>% 
  mutate(plate = stringr::str_split_fixed(plate, "_",2)[,1]) %>%
  dplyr::select(-seq)


# Plate 19 is a set of hash oligos that is just plate12 in a new orientation
new19 = 
  read.csv(file = "plate12_to_19_SS.csv", header = T) %>%
  left_join( all_hashes %>% filter(plate == "plate12"),
             by = c("plate12" = "position")) %>%
  mutate(plate = "plate19") %>%
  dplyr::select(plate,
                position = plate19,
                RT_seq,
                barcode) %>%
  arrange(position)

all_hashes = 
  all_hashes %>% 
  filter(plate != "plate19") %>%
  rbind(new19)

# This function builds a 384 well plate from 4 given 96 well plates in a specific order
build_384_well_plate = function(plate1, plate2, plate3, plate4, hashes, label = ""){
  
  plate1_table = 
    data.frame(
      hash =  hashes %>% filter(plate == plate1) %>% pull(RT_seq),
      row = seq(1:8) * 2  - 1,
      col = (seq(0,95) %/% 8 + 1)*2 - 1,
      position = "1",
      well_plate_96 = plate1)
  
  plate2_table = 
    data.frame(
      hash =  rev(hashes %>% filter(plate == plate2) %>% pull(RT_seq)),
      row = seq(1:8) * 2 - 1,
      col = (seq(0,95) %/% 8 + 1)*2,
      position = "2",
      well_plate_96 = plate2)
  
  plate3_table = 
    data.frame(
      hash =  rev(hashes %>% filter(plate == plate3) %>% pull(RT_seq)),
      row = seq(1:8) * 2,
      col = (seq(0,95) %/% 8 + 1)*2 - 1,
      position = "3",
      well_plate_96 = plate3)
  
  plate4_table =
    data.frame(
      hash =  hashes %>% filter(plate == plate4) %>% pull(RT_seq),
      row = seq(1:8) * 2,
      col = (seq(0,95) %/% 8 + 1)*2,
      position = "4",
      well_plate_96 = plate4)
  
  plate_384_well_plate = 
    rbind(plate1_table,
          plate2_table,
          plate3_table,
          plate4_table)
  
  plate_384_well_plate$plate = label
  
  rows_to_letters = 
    data.frame(row_letters = LETTERS[1:16],
               row = seq(1,16))
  
  
    
    plate_384_well_plate = 
      left_join(plate_384_well_plate,rows_to_letters, by = "row")
    
    plate_384_well_plate$well_position = paste0(label, "_",
                                                plate_384_well_plate$row_letters,
                                                plate_384_well_plate$col)
    
    return(plate_384_well_plate)
  
}


# Build the 19 384 well plates that constitute the spatial grid -----------

plate1 = build_384_well_plate("plate1","plate15","plate9","plate10",hashes = all_hashes, label = "Plate1")
plate2 = build_384_well_plate("plate2","plate5","plate2","plate6",hashes = all_hashes, label = "Plate2")
plate3 = build_384_well_plate("plate3","plate18","plate16","plate15",hashes = all_hashes, label = "Plate3")
plate4 = build_384_well_plate("plate4","plate13","plate5","plate9",hashes = all_hashes, label = "Plate4")
plate5 = build_384_well_plate("plate5","plate7","plate19","plate11",hashes = all_hashes, label = "Plate5")
plate6 = build_384_well_plate("plate6","plate16","plate15","plate5",hashes = all_hashes, label = "Plate6")
plate7 = build_384_well_plate("plate7","plate10","plate11","plate1",hashes = all_hashes, label = "Plate7")
plate8 = build_384_well_plate("plate8","plate14","plate8","plate18",hashes = all_hashes, label = "Plate8")
plate9 = build_384_well_plate("plate9","plate12","plate4","plate2",hashes = all_hashes, label = "Plate9")
plate10 = build_384_well_plate("plate10","plate2","plate12","plate4",hashes = all_hashes, label = "Plate10")
plate11 = build_384_well_plate("plate12","plate6","plate18","plate16",hashes = all_hashes, label = "Plate11")
plate12 = build_384_well_plate("plate11","plate9","plate7","plate13",hashes = all_hashes, label = "Plate12")
plate13 = build_384_well_plate("plate13","plate19","plate3","plate8",hashes = all_hashes, label = "Plate13")
plate14 = build_384_well_plate("plate14","plate3","plate14","plate3",hashes = all_hashes, label = "Plate14")
plate15 = build_384_well_plate("plate15","plate11","plate6","plate17",hashes = all_hashes, label = "Plate15")
plate16 = build_384_well_plate("plate16","plate17","plate1","plate12",hashes = all_hashes, label = "Plate16")
plate17 = build_384_well_plate("plate17","plate8","plate17","plate7",hashes = all_hashes, label = "Plate17")
plate18 = build_384_well_plate("plate18","plate1","plate13","plate14",hashes = all_hashes, label = "Plate18")
plate19 = build_384_well_plate("plate19","plate4","plate10","plate19",hashes = all_hashes, label = "Plate19")

# 
duplicated_positions = 
  c("Plate11_M21","Plate11_M22","Plate11_M23","Plate11_M24",
  "Plate11_N21","Plate11_N22","Plate11_N23","Plate11_N24",
  "Plate11_O21","Plate11_O22","Plate11_O23","Plate11_O24",
  "Plate11_P21","Plate11_P22","Plate11_P23","Plate11_P24")

new_positions = 
  c("Plate11_D4","Plate11_D3","Plate11_D2","Plate11_D1",
    "Plate11_C4","Plate11_C3","Plate11_C2","Plate11_C1",
    "Plate11_B4","Plate11_B3","Plate11_B2","Plate11_B1",
    "Plate11_A4","Plate11_A3","Plate11_A2","Plate11_A1")


flipped_well_position_df = 
  data.frame(old = duplicated_positions,
             new = new_positions)

plate11 = 
  plate11 %>%
  filter(!(well_position %in% new_positions)) %>%
  rbind(inner_join(plate11, 
             flipped_well_position_df,
             by = c("well_position" = "old")) %>%
  dplyr::select(everything(),
                -well_position,
                well_position = new))

all_plates = rbind(plate1,plate2,plate3,plate4,
                   plate5,plate6,plate7,plate8,
                   plate9,plate10,plate11,plate12,
                   plate13,plate14,plate15,plate16,
                   plate17,plate18,plate19)

  

all_plates = all_plates %>%
  dplyr::rename(Oligo = well_position) %>%
  dplyr::select(-plate)


oligo_plate_layout =
  left_join(oligo_plate_layout,all_plates, by = "Oligo" )



# Add SYBR Wells ----------------------------------------------------------

SYBR_wells = 
  read.table(file = "SYBR_wells.tsv",
           header = T,
           sep = "\t")

oligo_plate_layout = 
  oligo_plate_layout %>%
  left_join(SYBR_wells,
            by = "Oligo")


# Plot of how 96 well plates end up
ggplot(oligo_plate_layout) + 
  geom_tile(data=oligo_plate_layout %>% dplyr::select(-plate96), aes(x = Row, y = Col)) +
  geom_tile(aes(x = Row, y = Col, fill = plate96)) +
  facet_wrap(~plate96, ncol = 5) +
  theme_void()  +
  theme(legend.position = "none", 
        text = element_text(size = 6),
        strip.text.x = element_text(size = 6)) 


# Revalue oligo combos corresponding to rows and combos for sectors --------

oligo_plate_layout$rowCombo <- plyr::revalue((oligo_plate_layout$row_well),c("A" = 'A',
                                                                             "E" = 'A',
                                                                             "I" = 'A',
                                                                             "M" = 'A',
                                                                             "B" = 'B',
                                                                             "F" = 'B',
                                                                             "J" = 'B',
                                                                             "N" = 'B',
                                                                             "C" = 'C',
                                                                             "G" = 'C',
                                                                             "K" = 'C',
                                                                             "O" = 'C',
                                                                             "D" = 'D',
                                                                             "H" = 'D',
                                                                             "L" = 'D',
                                                                             "P" = 'D'))

oligo_plate_layout$colCombo <- plyr::revalue((oligo_plate_layout$col_well),c("1" = '1',
                                                                             "10" = '2',
                                                                             "11" = '3',
                                                                             "12" = '4',
                                                                             "13" = '1',
                                                                             "14" = '2',
                                                                             "15" = '3',
                                                                             "16" = '4',
                                                                             "17" = '1',
                                                                             "18" = '2',
                                                                             "19" = '3',
                                                                             "20" = '4',
                                                                             "2" = '2',
                                                                             "21" = '1',
                                                                             "22" = '2',
                                                                             "23" = '3',
                                                                             "24" = '4',
                                                                             "3" = '3',
                                                                             "4" = '4',
                                                                             "5" = '1',
                                                                             "6" = '2',
                                                                             "7" = '3',
                                                                             "8" = '4',
                                                                             "9" = '1'))



oligo_plate_layout$combo = paste(oligo_plate_layout$rowCombo, oligo_plate_layout$colCombo)

oligo_plate_layout$sector<- plyr::revalue((oligo_plate_layout$combo),c(
  "D 4"= "sector16", 
  "C 4" = "sector12", 
  "B 4" = "sector8", 
  "A 4" = "sector4", 
  "D 3"= "sector15",
  "C 3" = "sector11", 
  "B 3" = "sector7", 
  "A 3" = "sector3", 
  "D 2"= "sector14", 
  "C 2" = "sector10", 
  "B 2"= "sector6", 
  "A 2" = "sector2", 
  "D 1" = "sector13", 
  "C 1" = "sector9", 
  "B 1" = "sector5", 
  "A 1" = "sector1"))



ggplot(oligo_plate_layout) + 
  geom_tile(data=oligo_plate_layout %>% dplyr::select(-combo), aes(x = Row, y = Col)) +
  geom_tile(aes(x = Row, y = Col, fill = sector)) + 
  facet_wrap(~combo) +
  theme_void()  +
  theme(legend.position="none") 



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

sectors = 
  oligo_plate_layout$sector %>% 
  unique()

colors_df = data.frame(color = colors, 
                       sector = sectors)




# Code used to initially create the randomized SYBR wells  ----------------

# sybr_wells =
#   oligo_plate_layout %>%
#   sample_n(size = 192) %>%
#   pull(Oligo)
# 
# oligo_plate_layout =
#   oligo_plate_layout %>%
#   mutate(SYBR = (Oligo %in% sybr_wells))
# 

# oligo_plate_layout %>%
#   dplyr::select(Oligo,
#                 SYBR) %>%
#   mutate(Oligo= as.character(Oligo)) %>%
#   filter(SYBR) %>%
#   arrange(Oligo) %>% 
#   write.table(file = "SYBR_wells_for_print.tsv",
#               sep = "\t",
#               quote = F,
#               row.names = T)
# 
# oligo_plate_layout %>% 
#   dplyr::select(Oligo,
#                 SYBR) %>%
#   arrange(Oligo) %>%
#   write.table(file = "SYBR_wells.tsv",
#               sep = "\t",
#               quote = F,
#               row.names = F)



ggplot(oligo_plate_layout) + 
  geom_tile(data=oligo_plate_layout, 
            aes(x = Row, 
                y = Col,
                fill = SYBR)) +
  theme_void()  +
  theme(legend.position = "none",
        text = element_text(size = 6)) 

ggplot(oligo_plate_layout) + 
  geom_tile(data=oligo_plate_layout, 
            aes(x = Row, 
                y = Col,
                fill = combo)) +
  geom_point(data=oligo_plate_layout %>%
               filter(SYBR),
             aes(x = Row, 
                 y = Col),
             color = "white",
             size = 0.25) +
  theme_void()  +
  scale_fill_manual(values = colors) +
  theme(legend.position = "none",
        text = element_text(size = 6)) +
  ggsave("barcode_grid_2.png",
         height = 3,
         width = 3)

oligo_plate_layout = 
  left_join(oligo_plate_layout, 
            colors_df, by = "sector")

# Write out the result 

write.table(oligo_plate_layout,
            file = "oligo_plate_layout_2.tsv",
            sep = "\t",
            quote = F)


