#' calculate spatial weights
#' functions take from graph_test
#' @k is how many neighbors in knn graph 
#' @reduction_method default is UMAP
#' 
calculateSpatialWeights <- function(cds, k, reduction_method = "UMAP"){
  lw <- calculateLW(cds = cds, k = k, 
                    verbose = FALSE,
                    neighbor_graph = "knn",
                    reduction_method = reduction_method)
  wc <- spdep::spweights.constants(lw, zero.policy = TRUE, adjust.n = TRUE)
}


#' @var is a specific label value
#' @column_name of labels of interest
#' @lw is output of calculateLW 
calculateLocalG <- function(df, lw, wc, column_name, var) {
  
  # absence or presence of specified variable
  # df$LG <- ifelse(df[,column_name] == var, 1,0)
  # z = df[,"LG"]
  # names(z) <- df$Cell
  
  df$LG <- ifelse(df[,column_name] == var, 1,0)
  z = df[,"LG"]
  names(z) <- df$Cell
  
  # calculate local g
  localG = spdep::localG(x=z, listw = lw, wc$n, wc$S0, zero.policy = TRUE)
  localG.df = localG[1:length(localG)] %>% as.data.frame()
  colnames(localG.df)[1] <- "local_g"
  localG.df$Cell = rownames(localG.df)
  return(localG.df)
}

#' select local g values of interest + merge with full df
#' only values 1 mean anything in the local G calculation
mergeLocalG <- function(localG_result, df, column) {
  
  localG_result$Cell = row.names(localG_result)
  this.label = colnames(localG_result)[1]
  
  localG_result = 
    data_frame(Cell = localG_result$Cell,
               localG = localG_result[,1])
  
  sel = (df[,column] == this.label)
  
  wt = 
    df[sel,] %>%
    left_join(localG_result,
              by = "Cell")
  wt
}

#' to turn the zscore into a pvalue
#' method : c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
getPval <- function(df, method, tail) {
  
  # I think this is right, but should check it...
  if (tail=="onesided"){
    df = df %>% mutate(pval = pnorm(-abs(local_g)))
  } else if (tail=="twosided"){
    df = df %>% mutate(pval = 2*pnorm(-abs(local_g)))
  }
  adjusted_pval = p.adjust( p = df$pval, method = method)
  df %>% mutate(qval = adjusted_pval)
}


