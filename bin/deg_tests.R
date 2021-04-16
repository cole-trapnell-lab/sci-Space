
# Load startup packages ---------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(monocle3)
  
  DelayedArray:::set_verbose_block_processing(TRUE)
  options(DelayedArray.block.size=1000e7)
})

setwd(dir="~/projects/Space/")

spatial_cds = readRDS("/net/trapnell/vol1/home/sanjays/projects/Space/Submission/E14_slides/RDS_intermediates/Notebook6_spatial_cds_anatomy.RDS")

two.set.differential.gene.test = 
  function(cds, set.1.filter, set.2.filter, formal = F, cores = 1,
           covariates = NULL, sample.n.max = NULL) {
    s1.cds = cds[, set.1.filter]
    s2.cds = cds[, set.2.filter]
    
    if (!is.null(sample.n.max)) {
      set.seed(42)
      if (ncol(s1.cds) > sample.n.max)
        s1.cds = s1.cds[, sample(1:ncol(s1.cds), sample.n.max)]
      if (ncol(s2.cds) > sample.n.max)
        s2.cds = s2.cds[, sample(1:ncol(s2.cds), sample.n.max)]
    }
    
    message(paste("# of cells sampled in set 1:", ncol(s1.cds)))
    message(paste("# of cells sampled in set 2:", ncol(s2.cds)))
    
    s1.norm.expr = monocle3:::normalize_expr_data(s1.cds,
                                                  norm_method = "size_only",
                                                  pseudo_count = 0)
    s2.norm.expr = monocle3:::normalize_expr_data(s2.cds,
                                                  norm_method = "size_only",
                                                  pseudo_count = 0)
    
    s1.tpm = Matrix::rowSums(s1.norm.expr)
    s1.tpm = s1.tpm / sum(s1.tpm) * 1000000
    s2.tpm = Matrix::rowSums(s2.norm.expr)
    s2.tpm = s2.tpm / sum(s2.tpm) * 1000000
    
    s1.n.umi = Matrix::rowSums(counts(s1.cds))
    s2.n.umi = Matrix::rowSums(counts(s2.cds))
    
    higher.expr = ifelse(s1.tpm > s2.tpm, "Set 1", "Set 2")
    higher.tpm = ifelse(higher.expr == "Set 1", s1.tpm, s2.tpm)
    
    s1.ratio = s1.tpm / (s2.tpm + 1)
    s2.ratio = s2.tpm / (s1.tpm + 1)
    log2.ratio = ifelse(
      s1.tpm == 0 & s2.tpm == 0, 0, ifelse(
        higher.expr == "Set 1", log2(s1.ratio), log2(s2.ratio)))
    
    s1.n.expr = Matrix::rowSums(counts(s1.cds))
    s2.n.expr = Matrix::rowSums(counts(s2.cds))
    
    res = data.frame(
      gene = fData(cds)$gene_short_name,
      set.1.n.umi = s1.n.umi,
      set.2.n.umi = s2.n.umi,
      set.1.tpm = s1.tpm,
      set.2.tpm = s2.tpm,
      higher.expr = higher.expr,
      log2.ratio = log2.ratio,
      heuristic.score = log2.ratio * log2(higher.tpm+1)
    ) %>% arrange(-heuristic.score)
    
    if (formal) {
      pData(cds)$tmp = ifelse(set.1.filter, 1, ifelse(set.2.filter, 2, NA))
      
      cds.subset = cds[, set.1.filter | set.2.filter]
      cds.subset = estimate_size_factors(cds.subset)
      #cds.subset = estimateDispersions(cds.subset)
      cds.subset = detect_genes(cds.subset,min_expr =  0.1)
      
      expressed.genes = subset(fData(cds.subset), num_cells_expressed >= 5)[, 1]
      message(paste(length(expressed.genes), "genes expressed in at least 5 cells across both sets"))
      message("Computing differential expression p-values")
      
      DEG = fit_models(cds = cds.subset[expressed.genes,], model_formula_str = "~ tmp")
      
      res = inner_join(res,
                       DEG %>% 
                         coefficient_table() %>%
                         select(everything(),
                                gene = gene_short_name),
                       by = "gene")
      
      pData(cds)$tmp = NULL
    }
    
    return(res)
  }




DEG.results = 
  lapply(X = sort((unique(colData(spatial_cds)$cluster %>% as.character()))),
         FUN = function(this.cluster) { 
           message("Finding markers for cluster ", this.cluster)
           tmp = two.set.differential.gene.test(cds = spatial_cds,
                                          set.1.filter = as.character(colData(spatial_cds)$cluster) == this.cluster,
                                          set.2.filter = as.character(colData(spatial_cds)$cluster) != this.cluster,
                                          sample.n.max = 10000)
           tmp$cluster = this.cluster
           
           tmp})
          
DEG.results = 
  do.call(rbind,DEG.results) %>%
  filter(higher.expr == "Set 1")

DEG.results  %>%
  mutate(cluster = as.numeric(cluster)) %>%
  group_by(cluster) %>%
  top_n(n = 10,wt = heuristic.score) %>%
  arrange(cluster,-heuristic.score)%>%
  as.data.frame() %>%
  dplyr::select(-higher.expr) %>%
  write.table("Submission_Data/E14_slides/RDS_intermediates/manual_cell_annotations.tsv",
              sep = "\t",
              quote = F,
              col.names = T,
              row.names = F)

DEG.results  %>%
  mutate(cluster = as.numeric(cluster)) %>%
  filter(cluster == "39") %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = heuristic.score) %>%
  arrange(cluster,-heuristic.score)%>%
  as.data.frame() 
