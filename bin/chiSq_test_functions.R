suppressPackageStartupMessages({
  library(ggplot2)
  library(dirmult)
  library(tidyr)
  library(dplyr)
  library(pbapply)
})


prepare_hash_matrix <- function(hashes){
  hash_matrix =
    hashes %>%
    dplyr::select(Cell, Oligo, Count) %>%
    distinct() %>%
    spread(key="Oligo", value="Count", fill = 0, drop = F)
  rownames(hash_matrix) = hash_matrix$Cell
  cell_ids = hash_matrix$Cell
  hash_matrix = hash_matrix[,-1]
  hash_matrix = as.matrix(hash_matrix)
  
  return(list(hash_matrix,cell_ids))
}

generate_pval_matrix <- function(test_hash_matrix, hash_frequencies, overdispersion, cores=1){
  # Note: use outfile argument to makeCluster for debugging
  platform <- Sys.info()[['sysname']]
  if (platform == "Windows")
    cl <- makeCluster(cores)
  if (platform %in% c("Linux", "Darwin"))
    cl <- makeCluster(cores)
  
  cleanup <- function(){
    stopCluster(cl)
  }
  on.exit(cleanup)
  required_packages = c("dirmult")
  if (is.null(required_packages) == FALSE){
    clusterCall(cl, function(pkgs) {
      for (req in pkgs) {
        library(req, character.only=TRUE)
      }
    }, required_packages)
  }
  pval_matrix = t(pbapply(test_hash_matrix, 1, function(x, hash_frequencies, overdispersion) {
    total_counts = sum(x)
    sim_sample = simPop(J=1000,
                        n=total_counts,
                        pi=hash_frequencies,
                        theta=overdispersion)
    cdf_list = apply(sim_sample$data, 2, function(y) { ecdf(y) })
    pvals = list()
    for (i in 1:length(cdf_list)){
      pvals = append(pvals, cdf_list[[i]] (x[i]) )
    }
    pvals = 1 - unlist(pvals)
    qvals = p.adjust(pvals)
    qvals
  },
  hash_frequencies=background_hash_frequencies,
  overdispersion=overdispersion,
  cl=cl))
  colnames(pval_matrix) = colnames(test_hash_matrix)
  return(pval_matrix)
}


fit_background_dist <- function(hash_matrix, training_samples=5000, return_training_matrix=FALSE){
  #hash_mtrix = hash_matrix[sample(nrow(hash_matrix), 100),]
  hash_mtrix_sample = hash_matrix[sample(nrow(hash_matrix), training_samples),]
  background_fit = dirmult(hash_mtrix_sample, trace=TRUE)
  #dirmult.summary(hash_mtrix_sample, background_fit)
  if (return_training_matrix)
    return(list(background_fit, hash_mtrix_sample))
  else
    return(background_fit)
}


chisq_vs_background <- function(test_hash_matrix, hash_frequencies){
  hash_frequencies_nz = which(hash_frequencies > 0)
  hash_frequencies = hash_frequencies[hash_frequencies_nz]
  pvals= pbapply(test_hash_matrix[,hash_frequencies_nz], 1, function(x) {
    tryCatch({
      res = chisq.test(x, p=hash_frequencies,  simulate.p.value = FALSE)
      unlist(res[["p.value"]])
    }, error = function(e) { 1.0 }
    )
  })
  return(pvals)
}



assign_hash_labels <- function(test_cell_hashes, background_cell_hashes, min_best_vs_second_best_ratio=2, qval_thresh = 0.05, downsample_rate=NULL){
  background_hash_matrix = prepare_hash_matrix(background_cell_hashes)[[1]]
  # if (is.null(downsample_rate) == FALSE){
  #     background_hash_matrix = floor(background_hash_matrix * downsample_rate)
  # }
  background_hash_matrix = background_hash_matrix[rowSums(background_hash_matrix)>0,]
  
  background_hash_frequencies = colSums(background_hash_matrix)/sum(colSums(background_hash_matrix))
  
  test_hash_list = prepare_hash_matrix(test_cell_hashes)
  test_hash_matrix = test_hash_list[[1]]
  if (is.null(downsample_rate) == FALSE){
    test_hash_matrix = floor(test_hash_matrix * downsample_rate)
  }
  pvals = chisq_vs_background(test_hash_matrix, hash_frequencies=background_hash_frequencies)
  qvals = p.adjust(pvals)
  
  expected_background_hashes = outer(rowSums(test_hash_matrix), background_hash_frequencies)
  background_subtracted_test_hashes = test_hash_matrix - expected_background_hashes
  background_subtracted_test_hashes[background_subtracted_test_hashes < 0] = 0
  #hash_hits = background_subtracted_test_hashes * (qvals < qval_thresh)
  top_to_second_best_ratios = apply(background_subtracted_test_hashes, 1, function(x) { y = sort(x, decreasing=TRUE); y[1] / y[2]})
  top_hash = apply(background_subtracted_test_hashes, 1, function(x) {
    m = which(x == max(x));
    if (length(m) > 1)
      return (NA)
    colnames(background_subtracted_test_hashes)[which(x == max(x))]
  })
  #unambiguous_hits = which(top_to_second_best_ratios > min_best_vs_second_best_ratio)
  hash_label_df = data.frame(Cell = test_hash_list[[2]],
                             hash_umis = rowSums(test_hash_matrix),
                             pval = pvals,
                             qval = qvals,
                             #hit = qvals < qval_thresh,
                             top_to_second_best_ratio = top_to_second_best_ratios,
                             #unambiguous_hit = top_to_second_best_ratios > min_best_vs_second_best_ratio,
                             top_oligo = top_hash)
}



assign_hash_labels_return_all = function(test_cell_hashes, background_cell_hashes, min_best_vs_second_best_ratio=2, qval_thresh = 0.05, downsample_rate=NULL){
  background_hash_matrix = prepare_hash_matrix(background_cell_hashes)[[1]]
  
  background_hash_matrix = background_hash_matrix[rowSums(background_hash_matrix)>0,]
  
  background_hash_frequencies = colSums(background_hash_matrix)/sum(colSums(background_hash_matrix))
  
  test_hash_list = prepare_hash_matrix(test_cell_hashes)
  test_hash_matrix = test_hash_list[[1]]
  if (is.null(downsample_rate) == FALSE){
    test_hash_matrix = floor(test_hash_matrix * downsample_rate)
  }
  pvals = chisq_vs_background(test_hash_matrix, hash_frequencies=background_hash_frequencies)
  qvals = p.adjust(pvals)
  
  # if we saw X reads for a single cell multiply the background frequencies by the distribution
  expected_background_hashes = outer(rowSums(test_hash_matrix), background_hash_frequencies)
  
  # subtract those background_hashes from the test hashes present
  background_subtracted_test_hashes = test_hash_matrix - expected_background_hashes
  background_subtracted_test_hashes[background_subtracted_test_hashes < 0] = 0
  background_subtracted_test_hashes  =
    background_subtracted_test_hashes %>%
    as.data.frame
  
  background_subtracted_test_hashes$Cell = 
    test_hash_list[[2]]
  
  background_subtracted_test_hashes = 
    gather(data = background_subtracted_test_hashes,
           key = "Oligo", 
           value = adjusted_count, 
           -Cell )
  
  
  return(background_subtracted_test_hashes)
}
