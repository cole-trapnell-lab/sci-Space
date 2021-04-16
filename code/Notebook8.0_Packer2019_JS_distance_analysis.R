# Repeat analysis from Packer et. al. 2019 -- 10.1126/science.aax1971
# This notebook explores how the angular distance between related cells
# (cells that share the same parent) correlates with the Jensen-Shannon 
# distance and how angular distance can be used to calculate variance explained
# statistics using the Law of Total Variance

# Load startup packages ---------------------------------------------------
suppressPackageStartupMessages({ 
  library(monocle3)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  library(viridis)
  library(ggpubr)
  library(tidyr)
  library(stringr)
  
  space_directory = "~/Google Drive File Stream/My Drive/sciSpace/Submission_Data/"
  setwd(dir=space_directory)
  set.seed(42)
})

load(file = "Packer2019/lineage.tpm.bg.corrected.RData")
ls()


# Define functions for calculating distance measurements --------
kld = function(prob.1, prob.2) {
  if (any(!(prob.2 > 0)))
    warning("Zero values in prob.2")
  
  tmp = ifelse(prob.1 > 0, log2(prob.1/prob.2), 0)
  return(sum(prob.1 * tmp))
}

bin.js.distance = function(lin.1, lin.2, tpm.mat) {
  tpm.1 = tpm.mat[, lin.1]
  tpm.2 = tpm.mat[, lin.2]
  
  prob.1 = tpm.1 / sum(tpm.1)
  prob.2 = tpm.2 / sum(tpm.2)
  
  genes.1 = which(prob.1 > 0)
  genes.2 = which(prob.2 > 0)
  common.genes = unique(c(genes.1, genes.2))
  
  prob.1 = prob.1[common.genes]
  prob.2 = prob.2[common.genes]
  center = (prob.1 + prob.2) / 2.0
  
  js.div = 0.5 * kld(prob.1, center) + 0.5 * kld(prob.2, center)
  return(sqrt(js.div))
}


bin.ang.distance = function(lin.1, lin.2, tpm.mat) {
  tpm.1 = tpm.mat[, lin.1]
  tpm.1 = 
    log(tpm.1/1e6 + 1)
  
  tpm.2 = tpm.mat[, lin.2]
  tpm.2 = 
    log(tpm.2/1e6 + 1)
  
  ang.dist =
    angular_distance(tpm.1 %>% 
                     unit_vector(),
                     tpm.2 %>%
                     unit_vector())
  return(ang.dist)
}

# Calculate the norm of a vector
norm_vec  =
  function(x) {
    sqrt(sum(x^2))
  }

# Compute the angular distance between two cells
angular_distance = 
  function(cell_1, cell_2){
    if(identical(cell_1,cell_2)){
      0
    }
    else{
    (2/pi * acos(cell_1 %*% cell_2 /(norm_vec(cell_1) * norm_vec(cell_2))))^2
    }
  }

# Rescale so that the magnitude of the vector is unit
unit_vector = 
  function(vector){
    vector / norm_vec(vector)
  }

# Take a normalized expression count matrix
# Return the unit mean vector for that matrix
unit_mean =
  function(count_matrix){
    Matrix::rowSums(count_matrix) %>%
      unit_vector()
  }


# Load dataframes related to lineage information from Packer 2019 ---------
class(lineage.tpm.bg.corrected)
class(lineage.tpm.robust.bg.corrected)

valid.cells = readRDS("Packer2019/valid.cells.rds")
length(valid.cells)
valid.cells

all.lineage.df = readRDS("Packer2019/all.lineage.df.rds")
head(all.lineage.df)

load("Packer2019/lineage.tpm.bg.corrected.RData")

all.genes = rownames(lineage.tpm.robust.bg.corrected)


tmp = sort(unique(all.lineage.df$canonical.lineage))
tmp[!(tmp %in% colnames(lineage.tpm.bg.corrected))]

all.lineage.df = 
  all.lineage.df %>% 
  filter(canonical.lineage %in% colnames(lineage.tpm.bg.corrected))

is.terminal = read.table(
  "Packer2019/is.terminal",
  col.names = c("canonical.lineage", "is.terminal"),
  colClasses = c("character", "integer")
)

is.terminal$is.terminal = is.terminal$is.terminal == 1


all.lineage.df = 
  inner_join(all.lineage.df, 
             is.terminal, 
             by = "canonical.lineage")


# Isolate in the AB lineage -----------------------------------------------

ab.lineage.df = 
  all.lineage.df %>% 
  filter(clade == "AB")

ab.lineages = 
  sort(unique(ab.lineage.df$canonical.lineage))

ab.symmetric.lineages = 
  ab.lineage.df %>%
  filter(str_count(canonical.lineage, "/") + str_count(canonical.lineage, "x") < 2) %>%
  inner_join(ab.lineage.df, by = "canonical.lineage") %>%
  filter(cell.x != cell.y) %>%
  select(cell = cell.x, symmetric.cell = cell.y)

# remember to include reverse
tmp.df = data.frame(
  cell = c(
    "ABalapppaap", "ABalaaaarlp", # RME
    "ABalappppaaa", "ABalapaappaa", "ABalppapppaa", # IL1
    "ABalappppap", "ABalapaappp", "ABalppapppp", # IL2
    "ABalapppppp", "ABalppaappp", # CEPso
    "ABalpapapaa", "ABalpppapad", "ABalppapaaa", # RMD
    "ABalpapapapp", "ABalpapappp", # SMB
    "ABalappppaa", "ABalapaappa", "ABalppapppa", # IL1 parent
    "ABalappppa", "ABalapaapp", "ABalppappp", # IL1/2 neuroblast
    "ABplaaaaappa", "ABplpaappppa", # CEP
    "ABplaapaaap", "ABalaaapall", "ABalppapapp", # ILso
    "ABalapppapaa", "ABplpaaappaa", # OLQ
    "ABalapppapa", "ABplpaaappa", # OLQ parent
    "ABalapppap", "ABplpaaapp", # OLQ-URY neuroblast
    "ABplpapaapa", "ABplpapappa", # SIA
    "ABplppaaaaa", "ABplpapaapp", # SIB
    "ABplpapaaaa", "ABalppappaa", # SMD
    "ABarpapapa", "ABarpapapp", # hyp4/6
    "ABarpapppa", "ABarpapppp", # H0/H1,
    "ABalpaapaaa", "ABalpaappap", "ABaraaappaa", # posterior arcades
    "ABalpaapppa", "ABalpaapppp", # hyp1/2
    "ABplpaaapap", "ABarpaaaapp", # CEPsh
    "ABaraapaapp", # mc2D
    "ABarpaappa", # hyp7
    "ABplppppaap" # rect_V
  ),
  
  symmetric.cell = c(
    "ABplpappaaa", "ABalaaaarrp", # RME
    "ABalapappaaa", "ABalaappppaa", "ABarapppppaa", # IL1
    "ABalapappap", "ABalaappppp", "ABarapppppp", # IL2
    "ABalapapppp", "ABalaapappp", # CEPso
    "ABarappapaa", "ABpraaaapad", "ABarapppaaa", # RMD
    "ABarappapapp", "ABarappappp", # SMB
    "ABalapappaa", "ABalaappppa", "ABarapppppa", # IL1 parent
    "ABalapappa", "ABalaapppp", "ABarappppp", # IL1/2 neuroblast
    "ABarpapaappa", "ABprpaappppa", # CEP
    "ABpraapaaap", "ABalaaapprr", "ABarapppapp", # ILso
    "ABalapapapaa", "ABprpaaappaa", # OLQ
    "ABalapapapa", "ABprpaaappa", # OLQ parent
    "ABalapapap", "ABprpaaapp", # OLQ-URY neuroblast
    "ABprpapaapa", "ABprpapappa", # SIA
    "ABprppaaaaa", "ABprpapaapp", # SIB
    "ABprpapaaaa", "ABarappppaa", # SMD
    "ABplaaaapa", "ABplaaaapp", # hyp4/6
    "ABplaaappa", "ABplaaappp", # H0/H1
    "ABaraaapaaa", "ABaraaappap", "ABarapapapa", # posterior arcades
    "ABaraaapppp", "ABalpapaaap", # hyp1/2
    "ABprpaaapap", "ABarpaaapap", # CEPsh
    "ABaraappapp", # mc2V
    "ABarpaappp", # hyp7
    "ABprppppaap" # rect_V
  ),
  
  stringsAsFactors = F)

length(intersect(tmp.df$cell, tmp.df$symmetric.cell)) == 0

tmp.df = rbind(tmp.df, tmp.df[seq(nrow(tmp.df), 1, -1),])

# Define lineage relationships --------------------------------------------
ab.symmetric.lineages = rbind(ab.symmetric.lineages, tmp.df)
ab.lineage.nonsymmetric.nonredundant.sisters = ab.lineage.df %>%
  select(canonical.lineage, tmp.lineage, cell, generation, parent) %>%
  filter(parent != "NA") %>%
  inner_join(ab.lineage.df %>% select(
    sister.canonical.lineage = canonical.lineage,
    sister.tmp.lineage = tmp.lineage,
    sister.cell = cell,
    parent), by = "parent") %>%
  filter(cell != sister.cell,
         tmp.lineage != sister.tmp.lineage,
         sister.cell %in% valid.cells,
         sister.canonical.lineage %in% ab.lineages) %>%
  left_join(ab.symmetric.lineages,  by = "cell") %>%
  filter(is.na(symmetric.cell) | sister.cell != symmetric.cell) %>%
  select(canonical.lineage, sister.canonical.lineage, generation) %>%
  unique() %>%
  filter(as.character(canonical.lineage) < as.character(sister.canonical.lineage)) %>%
  arrange(canonical.lineage)

ab.lineage.nonsymmetric.nonredundant.sisters$canonical.lineage =
  as.character(ab.lineage.nonsymmetric.nonredundant.sisters$canonical.lineage)

ab.lineage.nonsymmetric.nonredundant.sisters$sister.canonical.lineage =
  as.character(ab.lineage.nonsymmetric.nonredundant.sisters$sister.canonical.lineage)

head(ab.lineage.nonsymmetric.nonredundant.sisters)

ab.lineage.nonsymmetric.nonredundant.cousins = ab.lineage.df %>%
  select(canonical.lineage, tmp.lineage, cell, generation, parent, grandparent) %>%
  filter(grandparent != "NA") %>%
  inner_join(ab.lineage.df %>% select(
    cousin.canonical.lineage = canonical.lineage,
    cousin.tmp.lineage = tmp.lineage,
    cousin.cell = cell,
    cousin.parent = parent,
    grandparent), by = "grandparent") %>%
  filter(cell != cousin.cell,
         parent != cousin.parent,
         tmp.lineage != cousin.tmp.lineage,
         cousin.cell %in% valid.cells,
         cousin.canonical.lineage %in% ab.lineages) %>%
  left_join(ab.symmetric.lineages,  by = "cell") %>%
  filter(is.na(symmetric.cell) | cousin.cell != symmetric.cell) %>%
  select(canonical.lineage, cousin.canonical.lineage, generation) %>%
  unique() %>%
  filter(as.character(canonical.lineage) < as.character(cousin.canonical.lineage)) %>%
  arrange(canonical.lineage)

ab.lineage.nonsymmetric.nonredundant.cousins$canonical.lineage =
  as.character(ab.lineage.nonsymmetric.nonredundant.cousins$canonical.lineage)

ab.lineage.nonsymmetric.nonredundant.cousins$cousin.canonical.lineage =
  as.character(ab.lineage.nonsymmetric.nonredundant.cousins$cousin.canonical.lineage)

head(ab.lineage.nonsymmetric.nonredundant.cousins)

ab.lineage.nonsymmetric.nonredundant.2nd.cousins = ab.lineage.df %>%
  select(canonical.lineage, tmp.lineage, cell, generation, parent, grandparent, g.grandparent) %>%
  filter(g.grandparent != "NA") %>%
  inner_join(ab.lineage.df %>% select(
    cousin.canonical.lineage = canonical.lineage,
    cousin.tmp.lineage = tmp.lineage,
    cousin.cell = cell,
    cousin.parent = parent,
    cousin.grandparent = grandparent,
    g.grandparent), by = "g.grandparent") %>%
  filter(cell != cousin.cell,
         parent != cousin.parent,
         grandparent != cousin.grandparent,
         tmp.lineage != cousin.tmp.lineage,
         cousin.cell %in% valid.cells,
         cousin.canonical.lineage %in% ab.lineages) %>%
  left_join(ab.symmetric.lineages,  by = "cell") %>%
  filter(is.na(symmetric.cell) | cousin.cell != symmetric.cell) %>%
  select(canonical.lineage, cousin.canonical.lineage, generation) %>%
  unique() %>%
  filter(as.character(canonical.lineage) < as.character(cousin.canonical.lineage)) %>%
  arrange(canonical.lineage)

ab.lineage.nonsymmetric.nonredundant.2nd.cousins$canonical.lineage =
  as.character(ab.lineage.nonsymmetric.nonredundant.2nd.cousins$canonical.lineage)

ab.lineage.nonsymmetric.nonredundant.2nd.cousins$cousin.canonical.lineage =
  as.character(ab.lineage.nonsymmetric.nonredundant.2nd.cousins$cousin.canonical.lineage)

head(ab.lineage.nonsymmetric.nonredundant.2nd.cousins)

ab.lineage.nonsymmetric.nonredundant.3rd.cousins = ab.lineage.df %>%
  select(canonical.lineage, tmp.lineage, cell, generation, parent, grandparent, g.grandparent, gg.grandparent) %>%
  filter(gg.grandparent != "NA") %>%
  inner_join(ab.lineage.df %>% select(
    cousin.canonical.lineage = canonical.lineage,
    cousin.tmp.lineage = tmp.lineage,
    cousin.cell = cell,
    cousin.parent = parent,
    cousin.grandparent = grandparent,
    cousin.g.grandparent = g.grandparent,
    gg.grandparent), by = "gg.grandparent") %>%
  filter(cell != cousin.cell,
         parent != cousin.parent,
         grandparent != cousin.grandparent,
         g.grandparent != cousin.g.grandparent,
         tmp.lineage != cousin.tmp.lineage,
         cousin.cell %in% valid.cells,
         cousin.canonical.lineage %in% ab.lineages) %>%
  left_join(ab.symmetric.lineages,  by = "cell") %>%
  filter(is.na(symmetric.cell) | cousin.cell != symmetric.cell) %>%
  select(canonical.lineage, cousin.canonical.lineage, generation) %>%
  unique() %>%
  filter(as.character(canonical.lineage) < as.character(cousin.canonical.lineage)) %>%
  arrange(canonical.lineage)

ab.lineage.nonsymmetric.nonredundant.3rd.cousins$canonical.lineage =
  as.character(ab.lineage.nonsymmetric.nonredundant.3rd.cousins$canonical.lineage)

ab.lineage.nonsymmetric.nonredundant.3rd.cousins$cousin.canonical.lineage =
  as.character(ab.lineage.nonsymmetric.nonredundant.3rd.cousins$cousin.canonical.lineage)

head(ab.lineage.nonsymmetric.nonredundant.3rd.cousins)

ab.lineage.nonsymmetric.nonredundant.4th.cousins = ab.lineage.df %>%
  select(canonical.lineage, tmp.lineage, cell, generation, parent, grandparent, g.grandparent, gg.grandparent, ggg.grandparent) %>%
  filter(ggg.grandparent != "NA") %>%
  inner_join(ab.lineage.df %>% select(
    cousin.canonical.lineage = canonical.lineage,
    cousin.tmp.lineage = tmp.lineage,
    cousin.cell = cell,
    cousin.parent = parent,
    cousin.grandparent = grandparent,
    cousin.g.grandparent = g.grandparent,
    cousin.gg.grandparent = gg.grandparent,
    ggg.grandparent), by = "ggg.grandparent") %>%
  filter(cell != cousin.cell,
         parent != cousin.parent,
         grandparent != cousin.grandparent,
         g.grandparent != cousin.g.grandparent,
         gg.grandparent != cousin.gg.grandparent,
         tmp.lineage != cousin.tmp.lineage,
         cousin.cell %in% valid.cells,
         cousin.canonical.lineage %in% ab.lineages) %>%
  left_join(ab.symmetric.lineages,  by = "cell") %>%
  filter(is.na(symmetric.cell) | cousin.cell != symmetric.cell) %>%
  select(canonical.lineage, cousin.canonical.lineage, generation) %>%
  unique() %>%
  filter(as.character(canonical.lineage) < as.character(cousin.canonical.lineage)) %>%
  arrange(canonical.lineage)

ab.lineage.nonsymmetric.nonredundant.4th.cousins$canonical.lineage =
  as.character(ab.lineage.nonsymmetric.nonredundant.4th.cousins$canonical.lineage)

ab.lineage.nonsymmetric.nonredundant.4th.cousins$cousin.canonical.lineage =
  as.character(ab.lineage.nonsymmetric.nonredundant.4th.cousins$cousin.canonical.lineage)

head(ab.lineage.nonsymmetric.nonredundant.4th.cousins)

ab.lineage.nonsymmetric.nonredundant.5th.cousins = ab.lineage.df %>%
  select(canonical.lineage, tmp.lineage, cell, generation, parent, grandparent, g.grandparent, gg.grandparent, ggg.grandparent, gggg.grandparent) %>%
  filter(gggg.grandparent != "NA") %>%
  inner_join(ab.lineage.df %>% select(
    cousin.canonical.lineage = canonical.lineage,
    cousin.tmp.lineage = tmp.lineage,
    cousin.cell = cell,
    cousin.parent = parent,
    cousin.grandparent = grandparent,
    cousin.g.grandparent = g.grandparent,
    cousin.gg.grandparent = gg.grandparent,
    cousin.ggg.grandparent = ggg.grandparent,
    gggg.grandparent), by = "gggg.grandparent") %>%
  filter(cell != cousin.cell,
         parent != cousin.parent,
         grandparent != cousin.grandparent,
         g.grandparent != cousin.g.grandparent,
         gg.grandparent != cousin.gg.grandparent,
         ggg.grandparent != cousin.ggg.grandparent,
         tmp.lineage != cousin.tmp.lineage,
         cousin.cell %in% valid.cells,
         cousin.canonical.lineage %in% ab.lineages) %>%
  left_join(ab.symmetric.lineages,  by = "cell") %>%
  filter(is.na(symmetric.cell) | cousin.cell != symmetric.cell) %>%
  select(canonical.lineage, cousin.canonical.lineage, generation) %>%
  unique() %>%
  filter(as.character(canonical.lineage) < as.character(cousin.canonical.lineage)) %>%
  arrange(canonical.lineage)

ab.lineage.nonsymmetric.nonredundant.5th.cousins$canonical.lineage =
  as.character(ab.lineage.nonsymmetric.nonredundant.5th.cousins$canonical.lineage)

ab.lineage.nonsymmetric.nonredundant.5th.cousins$cousin.canonical.lineage =
  as.character(ab.lineage.nonsymmetric.nonredundant.5th.cousins$cousin.canonical.lineage)

head(ab.lineage.nonsymmetric.nonredundant.5th.cousins)

ab.lineage.nonsymmetric.nonredundant.6th.cousins = ab.lineage.df %>%
  select(canonical.lineage, tmp.lineage, cell, generation, parent, grandparent, g.grandparent, gg.grandparent, ggg.grandparent, gggg.grandparent, ggggg.grandparent) %>%
  filter(ggggg.grandparent != "NA") %>%
  inner_join(ab.lineage.df %>% select(
    cousin.canonical.lineage = canonical.lineage,
    cousin.tmp.lineage = tmp.lineage,
    cousin.cell = cell,
    cousin.parent = parent,
    cousin.grandparent = grandparent,
    cousin.g.grandparent = g.grandparent,
    cousin.gg.grandparent = gg.grandparent,
    cousin.ggg.grandparent = ggg.grandparent,
    cousin.gggg.grandparent = gggg.grandparent,
    ggggg.grandparent), by = "ggggg.grandparent") %>%
  filter(cell != cousin.cell,
         parent != cousin.parent,
         grandparent != cousin.grandparent,
         g.grandparent != cousin.g.grandparent,
         gg.grandparent != cousin.gg.grandparent,
         ggg.grandparent != cousin.ggg.grandparent,
         gggg.grandparent != cousin.gggg.grandparent,
         tmp.lineage != cousin.tmp.lineage,
         cousin.cell %in% valid.cells,
         cousin.canonical.lineage %in% ab.lineages) %>%
  left_join(ab.symmetric.lineages,  by = "cell") %>%
  filter(is.na(symmetric.cell) | cousin.cell != symmetric.cell) %>%
  select(canonical.lineage, cousin.canonical.lineage, generation) %>%
  unique() %>%
  filter(as.character(canonical.lineage) < as.character(cousin.canonical.lineage)) %>%
  arrange(canonical.lineage)

ab.lineage.nonsymmetric.nonredundant.6th.cousins$canonical.lineage =
  as.character(ab.lineage.nonsymmetric.nonredundant.6th.cousins$canonical.lineage)

ab.lineage.nonsymmetric.nonredundant.6th.cousins$cousin.canonical.lineage =
  as.character(ab.lineage.nonsymmetric.nonredundant.6th.cousins$cousin.canonical.lineage)

head(ab.lineage.nonsymmetric.nonredundant.6th.cousins)

ab.lineage.nonsymmetric.nonredundant.7th.cousins = ab.lineage.df %>%
  select(canonical.lineage, tmp.lineage, cell, generation, parent, grandparent, g.grandparent, gg.grandparent, ggg.grandparent, gggg.grandparent, ggggg.grandparent, gggggg.grandparent) %>%
  filter(gggggg.grandparent != "NA") %>%
  inner_join(ab.lineage.df %>% select(
    cousin.canonical.lineage = canonical.lineage,
    cousin.tmp.lineage = tmp.lineage,
    cousin.cell = cell,
    cousin.parent = parent,
    cousin.grandparent = grandparent,
    cousin.g.grandparent = g.grandparent,
    cousin.gg.grandparent = gg.grandparent,
    cousin.ggg.grandparent = ggg.grandparent,
    cousin.gggg.grandparent = gggg.grandparent,
    cousin.ggggg.grandparent = ggggg.grandparent,
    gggggg.grandparent), by = "gggggg.grandparent") %>%
  filter(cell != cousin.cell,
         parent != cousin.parent,
         grandparent != cousin.grandparent,
         g.grandparent != cousin.g.grandparent,
         gg.grandparent != cousin.gg.grandparent,
         ggg.grandparent != cousin.ggg.grandparent,
         gggg.grandparent != cousin.gggg.grandparent,
         ggggg.grandparent != cousin.ggggg.grandparent,
         tmp.lineage != cousin.tmp.lineage,
         cousin.cell %in% valid.cells,
         cousin.canonical.lineage %in% ab.lineages) %>%
  left_join(ab.symmetric.lineages,  by = "cell") %>%
  filter(is.na(symmetric.cell) | cousin.cell != symmetric.cell) %>%
  select(canonical.lineage, cousin.canonical.lineage, generation) %>%
  unique() %>%
  filter(as.character(canonical.lineage) < as.character(cousin.canonical.lineage)) %>%
  arrange(canonical.lineage)

ab.lineage.nonsymmetric.nonredundant.7th.cousins$canonical.lineage =
  as.character(ab.lineage.nonsymmetric.nonredundant.7th.cousins$canonical.lineage)

ab.lineage.nonsymmetric.nonredundant.7th.cousins$cousin.canonical.lineage =
  as.character(ab.lineage.nonsymmetric.nonredundant.7th.cousins$cousin.canonical.lineage)

head(ab.lineage.nonsymmetric.nonredundant.7th.cousins)

ab.lineage.nonsymmetric.nonredundant.8th.cousins = ab.lineage.df %>%
  select(canonical.lineage, tmp.lineage, cell, generation, parent, grandparent, g.grandparent, gg.grandparent, ggg.grandparent, gggg.grandparent, ggggg.grandparent, gggggg.grandparent, ggggggg.grandparent) %>%
  filter(ggggggg.grandparent != "NA") %>%
  inner_join(ab.lineage.df %>% select(
    cousin.canonical.lineage = canonical.lineage,
    cousin.tmp.lineage = tmp.lineage,
    cousin.cell = cell,
    cousin.parent = parent,
    cousin.grandparent = grandparent,
    cousin.g.grandparent = g.grandparent,
    cousin.gg.grandparent = gg.grandparent,
    cousin.ggg.grandparent = ggg.grandparent,
    cousin.gggg.grandparent = gggg.grandparent,
    cousin.ggggg.grandparent = ggggg.grandparent,
    cousin.gggggg.grandparent = gggggg.grandparent,
    ggggggg.grandparent), by = "ggggggg.grandparent") %>%
  filter(cell != cousin.cell,
         parent != cousin.parent,
         grandparent != cousin.grandparent,
         g.grandparent != cousin.g.grandparent,
         gg.grandparent != cousin.gg.grandparent,
         ggg.grandparent != cousin.ggg.grandparent,
         gggg.grandparent != cousin.gggg.grandparent,
         ggggg.grandparent != cousin.ggggg.grandparent,
         gggggg.grandparent != cousin.gggggg.grandparent,
         tmp.lineage != cousin.tmp.lineage,
         cousin.cell %in% valid.cells,
         cousin.canonical.lineage %in% ab.lineages) %>%
  left_join(ab.symmetric.lineages,  by = "cell") %>%
  filter(is.na(symmetric.cell) | cousin.cell != symmetric.cell) %>%
  select(canonical.lineage, cousin.canonical.lineage, generation) %>%
  unique() %>%
  filter(as.character(canonical.lineage) < as.character(cousin.canonical.lineage)) %>%
  arrange(canonical.lineage)

ab.lineage.nonsymmetric.nonredundant.8th.cousins$canonical.lineage =
  as.character(ab.lineage.nonsymmetric.nonredundant.8th.cousins$canonical.lineage)

ab.lineage.nonsymmetric.nonredundant.8th.cousins$cousin.canonical.lineage =
  as.character(ab.lineage.nonsymmetric.nonredundant.8th.cousins$cousin.canonical.lineage)

head(ab.lineage.nonsymmetric.nonredundant.8th.cousins)


# Calculate JS distance between cells based on lineage distance -----------
ab.lineage.nonsymmetric.nonredundant.sisters$js.distance =
  with(ab.lineage.nonsymmetric.nonredundant.sisters, mapply(function(x, y) {
    bin.js.distance(x, y, lineage.tpm.robust.bg.corrected)
  }, canonical.lineage, sister.canonical.lineage))

ab.lineage.nonsymmetric.nonredundant.cousins$js.distance =
  with(ab.lineage.nonsymmetric.nonredundant.cousins, mapply(function(x, y) {
    bin.js.distance(x, y, lineage.tpm.robust.bg.corrected)
  }, canonical.lineage, cousin.canonical.lineage))

ab.lineage.nonsymmetric.nonredundant.2nd.cousins$js.distance =
  with(ab.lineage.nonsymmetric.nonredundant.2nd.cousins, mapply(function(x, y) {
    bin.js.distance(x, y, lineage.tpm.robust.bg.corrected)
  }, canonical.lineage, cousin.canonical.lineage))

ab.lineage.nonsymmetric.nonredundant.3rd.cousins$js.distance =
  with(ab.lineage.nonsymmetric.nonredundant.3rd.cousins, mapply(function(x, y) {
    bin.js.distance(x, y, lineage.tpm.robust.bg.corrected)
  }, canonical.lineage, cousin.canonical.lineage))

ab.lineage.nonsymmetric.nonredundant.4th.cousins$js.distance =
  with(ab.lineage.nonsymmetric.nonredundant.4th.cousins, mapply(function(x, y) {
    bin.js.distance(x, y, lineage.tpm.robust.bg.corrected)
  }, canonical.lineage, cousin.canonical.lineage))

ab.lineage.nonsymmetric.nonredundant.5th.cousins$js.distance =
  with(ab.lineage.nonsymmetric.nonredundant.5th.cousins, mapply(function(x, y) {
    bin.js.distance(x, y, lineage.tpm.robust.bg.corrected)
  }, canonical.lineage, cousin.canonical.lineage))

ab.lineage.nonsymmetric.nonredundant.6th.cousins$js.distance =
  with(ab.lineage.nonsymmetric.nonredundant.6th.cousins, mapply(function(x, y) {
    bin.js.distance(x, y, lineage.tpm.robust.bg.corrected)
  }, canonical.lineage, cousin.canonical.lineage))

ab.lineage.nonsymmetric.nonredundant.7th.cousins$js.distance =
  with(ab.lineage.nonsymmetric.nonredundant.7th.cousins, mapply(function(x, y) {
    bin.js.distance(x, y, lineage.tpm.robust.bg.corrected)
  }, canonical.lineage, cousin.canonical.lineage))

ab.lineage.nonsymmetric.nonredundant.8th.cousins$js.distance =
  with(ab.lineage.nonsymmetric.nonredundant.8th.cousins, mapply(function(x, y) {
    bin.js.distance(x, y, lineage.tpm.robust.bg.corrected)
  }, canonical.lineage, cousin.canonical.lineage))

# Calculate angular distance between cells based on lineage distance -----------

ab.lineage.nonsymmetric.nonredundant.sisters$ang.distance =
  with(ab.lineage.nonsymmetric.nonredundant.sisters, mapply(function(x, y) {
    bin.ang.distance(x, y, lineage.tpm.robust.bg.corrected)
  }, canonical.lineage, sister.canonical.lineage))

ab.lineage.nonsymmetric.nonredundant.cousins$ang.distance =
  with(ab.lineage.nonsymmetric.nonredundant.cousins, mapply(function(x, y) {
    bin.ang.distance(x, y, lineage.tpm.robust.bg.corrected)
  }, canonical.lineage, cousin.canonical.lineage))

ab.lineage.nonsymmetric.nonredundant.2nd.cousins$ang.distance =
  with(ab.lineage.nonsymmetric.nonredundant.2nd.cousins, mapply(function(x, y) {
    bin.ang.distance(x, y, lineage.tpm.robust.bg.corrected)
  }, canonical.lineage, cousin.canonical.lineage))

ab.lineage.nonsymmetric.nonredundant.3rd.cousins$ang.distance =
  with(ab.lineage.nonsymmetric.nonredundant.3rd.cousins, mapply(function(x, y) {
    bin.ang.distance(x, y, lineage.tpm.robust.bg.corrected)
  }, canonical.lineage, cousin.canonical.lineage))

ab.lineage.nonsymmetric.nonredundant.4th.cousins$ang.distance =
  with(ab.lineage.nonsymmetric.nonredundant.4th.cousins, mapply(function(x, y) {
    bin.ang.distance(x, y, lineage.tpm.robust.bg.corrected)
  }, canonical.lineage, cousin.canonical.lineage))


ab.lineage.nonsymmetric.nonredundant.5th.cousins$ang.distance =
  with(ab.lineage.nonsymmetric.nonredundant.5th.cousins, mapply(function(x, y) {
    bin.ang.distance(x, y, lineage.tpm.robust.bg.corrected)
  }, canonical.lineage, cousin.canonical.lineage))

ab.lineage.nonsymmetric.nonredundant.6th.cousins$ang.distance =
  with(ab.lineage.nonsymmetric.nonredundant.6th.cousins, mapply(function(x, y) {
    bin.ang.distance(x, y, lineage.tpm.robust.bg.corrected)
  }, canonical.lineage, cousin.canonical.lineage))

ab.lineage.nonsymmetric.nonredundant.7th.cousins$ang.distance =
  with(ab.lineage.nonsymmetric.nonredundant.7th.cousins, mapply(function(x, y) {
    bin.ang.distance(x, y, lineage.tpm.robust.bg.corrected)
  }, canonical.lineage, cousin.canonical.lineage))

ab.lineage.nonsymmetric.nonredundant.8th.cousins$ang.distance =
  with(ab.lineage.nonsymmetric.nonredundant.8th.cousins, mapply(function(x, y) {
    bin.ang.distance(x, y, lineage.tpm.robust.bg.corrected)
  }, canonical.lineage, cousin.canonical.lineage))


# Collate results into a single dataframe ---------------------------------

ab.lineage.nonsymmetric.nonredundant.sisters$relationship = "sisters"
ab.lineage.nonsymmetric.nonredundant.cousins$relationship = "cousins"
ab.lineage.nonsymmetric.nonredundant.2nd.cousins$relationship = "2nd cousins"
ab.lineage.nonsymmetric.nonredundant.3rd.cousins$relationship = "3rd cousins"
ab.lineage.nonsymmetric.nonredundant.4th.cousins$relationship = "4th cousins"
ab.lineage.nonsymmetric.nonredundant.5th.cousins$relationship = "5th cousins"
ab.lineage.nonsymmetric.nonredundant.6th.cousins$relationship = "6th cousins"
ab.lineage.nonsymmetric.nonredundant.7th.cousins$relationship = "7th cousins"
ab.lineage.nonsymmetric.nonredundant.8th.cousins$relationship = "8th cousins"

ab.lineage.nonsymmetric.nonredundant.sisters$lineage.distance = 1
ab.lineage.nonsymmetric.nonredundant.cousins$lineage.distance = 2
ab.lineage.nonsymmetric.nonredundant.2nd.cousins$lineage.distance = 3
ab.lineage.nonsymmetric.nonredundant.3rd.cousins$lineage.distance = 4
ab.lineage.nonsymmetric.nonredundant.4th.cousins$lineage.distance = 5
ab.lineage.nonsymmetric.nonredundant.5th.cousins$lineage.distance = 6
ab.lineage.nonsymmetric.nonredundant.6th.cousins$lineage.distance = 7
ab.lineage.nonsymmetric.nonredundant.7th.cousins$lineage.distance = 8
ab.lineage.nonsymmetric.nonredundant.8th.cousins$lineage.distance = 9


ab.plot.df = rbind(
  ab.lineage.nonsymmetric.nonredundant.sisters %>%
    select(canonical.lineage, relative = sister.canonical.lineage, generation, js.distance, relationship, lineage.distance, ang.distance),
  ab.lineage.nonsymmetric.nonredundant.cousins %>%
    select(canonical.lineage, relative = cousin.canonical.lineage, generation, js.distance, relationship, lineage.distance, ang.distance),
  ab.lineage.nonsymmetric.nonredundant.2nd.cousins %>%
    select(canonical.lineage, relative = cousin.canonical.lineage, generation, js.distance, relationship, lineage.distance, ang.distance),
  ab.lineage.nonsymmetric.nonredundant.3rd.cousins %>%
    select(canonical.lineage, relative = cousin.canonical.lineage, generation, js.distance, relationship, lineage.distance, ang.distance),
  ab.lineage.nonsymmetric.nonredundant.4th.cousins %>%
    select(canonical.lineage, relative = cousin.canonical.lineage, generation, js.distance, relationship, lineage.distance, ang.distance),
  ab.lineage.nonsymmetric.nonredundant.5th.cousins %>%
    select(canonical.lineage, relative = cousin.canonical.lineage, generation, js.distance, relationship, lineage.distance, ang.distance),
  ab.lineage.nonsymmetric.nonredundant.6th.cousins %>%
    select(canonical.lineage, relative = cousin.canonical.lineage, generation, js.distance, relationship, lineage.distance, ang.distance),
  ab.lineage.nonsymmetric.nonredundant.7th.cousins %>%
    select(canonical.lineage, relative = cousin.canonical.lineage, generation, js.distance, relationship, lineage.distance, ang.distance),
  ab.lineage.nonsymmetric.nonredundant.8th.cousins %>%
    select(canonical.lineage, relative = cousin.canonical.lineage, generation, js.distance, relationship, lineage.distance, ang.distance))

ab.plot.df$relationship = factor(ab.plot.df$relationship, levels = c(
  "sisters", "cousins", "2nd cousins", "3rd cousins", "4th cousins",
  "5th cousins", "6th cousins", "7th cousins", "8th cousins"))

ab.plot.df$lineage.distance = as.integer(ab.plot.df$lineage.distance)

ab.plot.df = inner_join(ab.plot.df, is.terminal, by = "canonical.lineage")

ab.plot.df = inner_join(
  ab.plot.df,
  is.terminal %>% select(relative = canonical.lineage, relative.is.terminal = is.terminal),
  by = "relative"
)

ab.plot.df$both.preterminal = with(ab.plot.df, !is.terminal & !relative.is.terminal)
ab.plot.df$both.terminal = with(ab.plot.df, is.terminal & relative.is.terminal)

ab.plot.df$facet = rep(NA, nrow(ab.plot.df))
ab.plot.df$facet = with(ab.plot.df, ifelse(both.preterminal, "Pre-terminal", facet))
ab.plot.df$facet = with(ab.plot.df, ifelse(both.terminal, "Terminal", facet))
ab.plot.df$facet = factor(ab.plot.df$facet, levels = c("Pre-terminal", "Terminal"))


head(ab.plot.df)
nrow(ab.plot.df)


# Supplementary Figure 28 — Panel A ---------------------------------------

plot = ggplot(
  ab.plot.df %>% filter(
    (generation %in% c(5, 6, 7, 8) & both.preterminal) |
      (generation %in% c(9) & both.terminal),
    relationship %in% c("sisters", "cousins", "2nd cousins", "3rd cousins", "4th cousins",
                        "5th cousins")),
  aes(
    x = factor(generation),
    y = ang.distance,
    fill = relationship
  )
) +
  facet_grid(~ facet, scales = "free_x", space = "free_x") +
  geom_boxplot(outlier.size = 0.75,
               outlier.stroke = 0) +
  scale_y_continuous(breaks = seq(0.0, 1.0, 0.1)) +
  scale_fill_discrete(labels = c("Sisters", "Cousins", "2nd cousins", "3rd cousins", "4th cousins", "5th cousins")) +
  xlab("AB lineage generation") +
  ylab("Angular distance") +
  guides(fill = guide_legend(title = "Lineage\nrelationship")) +
  monocle:::monocle_theme_opts() +
  theme() +
  ggsave("Figures/Figure_Components/Supplement_Packer2019/angular_distance_by_generation.pdf",
         height = 3,
         width = 5)

# Supplementary Figure 28 — Panel B ---------------------------------------
ggplot(ab.plot.df,
       aes(x = js.distance,
           y = ang.distance)) +
  geom_point(stroke = 0,
             size = 0.25) +
  stat_cor(size = 3,
    label.y.npc = 0.99,
    color = "grey31",
    method = "pearson") +
  monocle3:::monocle_theme_opts() +
  ylim(-0.05,1) +
  xlab("Jenson-Shannon Distance")+
  ylab("Angular Distance") +
  ggsave("Figures/Figure_Components/Supplement_Packer2019/ang_js_dist_cor.pdf",
         height = 3,
         width = 3)


# Perform variance decomposition ------------------------------------------
# now for the ab-lineage we have the a data frame for a pair of cells and their relationship 
# each set of relationships


# Get the unit mean cell vector for every generation
global_means =
  lapply(X = seq(5,10),
       FUN = function(gen,
                      count_mat,
                      lineage.df){
         
         # Isolate the cells that are going to be looked at
         cells.in.this.generation = 
           lineage.df %>%
           filter(generation == gen) %>%
           dplyr::select(canonical.lineage,
                         parent,
                         generation) %>%
           distinct() %>%
           group_by(parent) %>%
           filter(n() > 1) %>% 
           filter(canonical.lineage %in% colnames(lineage.tpm.robust.bg.corrected)) %>%
           pull(canonical.lineage) %>%
           unique() %>%
           as.character()
         
         # Get the global mean for these cells
         global_unit_mean =
           count_mat[,cells.in.this.generation] %>%
           unit_mean() 
         
         tibble(generation = gen,
                global_mean =  list(global_unit_mean))
           
       },lineage.tpm.robust.bg.corrected,
       ab.lineage.df)

global_means = 
  do.call(rbind,
          global_means)

# Calculate the average distance between every cell and the global mean
# Understanding the variance as the second moment of the dataset
total_variance_per_generation =
  ab.lineage.df %>%
  filter(canonical.lineage %in% colnames(lineage.tpm.robust.bg.corrected)) %>% 
  filter(generation >= 5) %>%
  group_by(generation,parent) %>%
  filter(n() > 1) %>% 
  ungroup() %>% 
  left_join(global_means,
            by = "generation") %>%
  group_by(generation) %>%
  nest() %>%
  mutate(total_var = 
           purrr:::map(.x = data,
                       .f = function(generation_group,
                                     count_mat){
                         
                         cells_in_this_generation = 
                           generation_group$canonical.lineage
                         
                         generation_global_mean =
                           generation_group$global_mean[[1]]
                         
                         sapply(X = cells_in_this_generation,
                                  function(cell){
                                    cell_1 = 
                                      count_mat[,cell] %>%
                                      unit_vector()
                                    
                                    angular_distance(cell_1 = cell_1,
                                                     cell_2 = generation_global_mean)
                                  }) %>%
                           mean()
                       },lineage.tpm.robust.bg.corrected))

  
total_variance_per_generation = 
  total_variance_per_generation %>%
  dplyr::select(-data) %>%
  unnest(cols = c("total_var"))

# Calculate the average variance using the law of total variance
# Breaking down the variance within a group (between two daughter cells) and
# between that group and the global mean
conditional_var_explained =
  lapply(X = seq(5,10),
         FUN = function(gen,
                        lineage.df){
          tmp =
            lineage.df %>%
             filter(generation == gen) %>% 
             dplyr::select(canonical.lineage,
                           parent,
                           generation) %>%
             distinct() %>%
             group_by(parent) %>%
             filter(canonical.lineage %in% colnames(lineage.tpm.robust.bg.corrected)) %>%
             filter(n() > 1) %>% 
             left_join(global_means,
                       by = "generation") %>%
             nest()  %>%
             mutate(var_explained =
                      purrr::map(.x = data,
                                 .f = function(lineage_group,
                                               count_mat){
                                   cells_in_group = 
                                     lineage_group$canonical.lineage %>%
                                     as.character()
                                   
                                   group_mean =
                                     count_mat[,cells_in_group] %>%
                                     unit_mean()
                                   
                                   gg_ang_dist =
                                     angular_distance(group_mean,
                                                      lineage_group$global_mean[[1]])
                                   
                                   num_cells_in_group = dim(lineage_group)[1]
                                   
                                   cell_group_angular_dist = 
                                     sapply(X = cells_in_group,
                                            function(cell,
                                                     group_angular_mean){
                                              cell_1 =
                                                count_mat[,cell] %>%
                                                unit_vector()
                                              
                                              angular_distance(cell_1 = cell_1,
                                                               cell_2 = group_angular_mean)
                                            },group_mean)
                                   
                                   cell_group_angular_dist = 
                                     mean(cell_group_angular_dist)
                                   
                                   
                                   tibble(cell_group_angular_dist = cell_group_angular_dist,
                                          global_group_ang_dist = gg_ang_dist,
                                          generation = lineage_group$generation,
                                          num_cells = num_cells_in_group,
                                          group_mean = list(group_mean))
                                 },lineage.tpm.robust.bg.corrected))
        tmp %>%
          dplyr::select(-data) %>%
          unnest(cols = c("var_explained"))
         },ab.lineage.df)

conditional_var_explained =
  do.call(rbind,
          conditional_var_explained)
      
conditional_var_explained =
  conditional_var_explained %>%
  dplyr::select(parent,
                cell_group_angular_dist,
                global_group_ang_dist,
                generation,
                num_cells) %>%
  distinct()

# get only cells for which both daughter cells are identified in the dataset
conditional_var_explained[conditional_var_explained$num_cells != 2,]

# Sum the within group and between group average variances
conditional_var_explained = 
  conditional_var_explained %>%
  group_by(generation) %>%
  summarise(mean_conditional_ang_dist = mean(cell_group_angular_dist),
            mean_global_group_ang_dist =mean(global_group_ang_dist)) %>%
  group_by(generation) %>%
  mutate(sum_of_conditionals = sum(mean_conditional_ang_dist,mean_global_group_ang_dist))

conditional_var_explained =
  conditional_var_explained %>%
  left_join(total_variance_per_generation,
            by = "generation")


# Supplemental Figure 28 — Panel C ----------------------------------------

conditional_var_explained %>%
  dplyr::select(generation,
                sum_of_conditionals,
                total_var) %>%
  gather(key =  "variance_calculation",
         value = "Variance",
         -generation) %>%
  ggplot() +
  geom_bar(aes(x = as.factor(generation),
               y = Variance,
               fill = variance_calculation),
           stat = "identity",
           position = "dodge",
           size = 0.4,
           color = "black",
           width = 0.9,) +
  scale_fill_brewer(palette = "Set1",
                    name = "Calculation",
                    label = c("Variance \nDecomposition",
                              "Total Variance")) +
  monocle3:::monocle_theme_opts() +
  xlab("AB lineage generation") +
  ylab("Average Variance") +
  ggsave("Figures/Figure_Components/Supplement_Packer2019/law_of_total_variance_pseudobulk.pdf",
         height = 2, 
         width = 4)


# Calculate a variance explained null model for each generation ------------
# How much of the variance can be explained by knowing that two cells are sisters
# tabulated per generation.


# To calculate a variance explained statistic divide the within group variance 
# by the variance of a model with the same degrees of freedom

# Calculate the variance of a model with the same degrees of freedom by permuting
# the labels.
conditional_var_permuted =
  lapply(X = seq(5,10),
         FUN = function(gen,
                        lineage.df){
           tmp =
             lineage.df %>%
             filter(generation == gen) %>% 
             filter(canonical.lineage %in% colnames(lineage.tpm.robust.bg.corrected)) %>%
             dplyr::select(canonical.lineage,
                           parent,
                           generation) %>%
             distinct() %>%
             left_join(global_means,
                       by = "generation") %>%
             group_by(parent) %>%
             filter(n() > 1) %>% 
             ungroup() %>%
             # Labels permuted
             mutate(shuffled.parent = sample(parent,
                                             replace = F)) %>%
             # Group by the permuted parent cell
             group_by(shuffled.parent) %>%
             nest()  %>%
             mutate(var_explained =
                      purrr::map(.x = data,
                                 .f = function(lineage_group,
                                               count_mat){
                                   cells_in_group = 
                                     lineage_group$canonical.lineage %>%
                                     as.character()
                                   
                                   group_mean =
                                     count_mat[,cells_in_group] %>%
                                     unit_mean()
                                   
                                   gg_ang_dist =
                                     angular_distance(group_mean,
                                                      lineage_group$global_mean[[1]]) %>%
                                     as.numeric()
                                   
                                   num_cells_in_group = dim(lineage_group)[1]
                                   
                                   cell_group_angular_dist = 
                                     sapply(X = cells_in_group,
                                            function(cell,
                                                     group_angular_mean){
                                              cell_1 =
                                                count_mat[,cell] %>%
                                                unit_vector()
                                              
                                              angular_distance(cell_1 = cell_1,
                                                               cell_2 = group_angular_mean)
                                              
                                              
                                            },group_mean)
                                   
                                   cell_group_angular_dist = 
                                     mean(cell_group_angular_dist)
                                   
                                   
                                   tibble(cell_group_angular_dist = cell_group_angular_dist,
                                          global_group_ang_dist = gg_ang_dist,
                                          generation = lineage_group$generation,
                                          num_cells = num_cells_in_group,
                                          group_mean = list(group_mean))
                                 },lineage.tpm.robust.bg.corrected))
           tmp %>%
             dplyr::select(-data) %>%
             unnest(cols = c("var_explained"))
         },ab.lineage.df)

# Collate the results for all the generations
conditional_var_permuted =
  do.call(rbind,
          conditional_var_permuted)

conditional_var_permuted =
  conditional_var_permuted %>%
  dplyr::select(shuffled.parent,
                cell_group_angular_dist,
                global_group_ang_dist,
                generation,
                num_cells) %>%
  distinct()

conditional_var_permuted[conditional_var_permuted$num_cells != 2,]

# Calculate the mean within group and between group variances
conditional_var_permuted = 
  conditional_var_permuted %>%
  group_by(generation) %>%
  summarise(mean_conditional_ang_dist = mean(cell_group_angular_dist),
            mean_global_group_ang_dist =mean(global_group_ang_dist)) %>%
  group_by(generation) %>%
  mutate(sum_of_conditionals = sum(mean_conditional_ang_dist,mean_global_group_ang_dist))


colnames(conditional_var_permuted)[colnames(conditional_var_permuted) != "generation"] = 
  paste0("permuted.",
         colnames(conditional_var_permuted)[colnames(conditional_var_permuted) != "generation"],
         sep = "")

conditional_var_explained =
  conditional_var_explained %>%
  left_join(conditional_var_permuted,
            by = "generation")

# Calculate a variance explained statistic for each generation ------------

conditional_var_explained$var_explained = 
  1 - (conditional_var_explained$mean_conditional_ang_dist/conditional_var_explained$permuted.mean_conditional_ang_dist)


# Supplemental Figure 28 — Panel D ----------------------------------------
ggplot(conditional_var_explained %>%
         filter(generation != 10)) +
  geom_bar(aes(x = as.factor(generation),
               y = var_explained *100,
               fill = as.factor(generation)),
           stat = "identity",
           color = "black",
           size = 0.4) +
  monocle3:::monocle_theme_opts() +
  scale_fill_manual(values = c("#003f5c", "#58508d", "#bc5090", "#ff6361", "#ffa600")) +
  xlab("AB lineage generation") +
  ylab("Variance Explained") +
  scale_y_continuous(breaks = seq(0,100,10),
                     labels = c(0,paste(seq(10,100,10),
                                        "%",
                                        sep =""))) +
  theme(legend.position = "none") +
  ggsave("Figures/Figure_Components/Supplement_Packer2019/r2_pseudobulk.pdf",
         height = 2,
         width = 3)

  