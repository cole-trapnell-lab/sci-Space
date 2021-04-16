# Krigging is a method used in geostatistics where values are interpolated and modeled  
# by a gaussian process. Here we apply this method to interpolate gene expression from 
# each celltype

# Load startup packages ---------------------------------------------------
suppressPackageStartupMessages({ 
  library(tidyverse)
  library(ggplot2)
  library(tictoc)
  library(monocle3)
  library(tidymodels)
  library(furrr)
  library(ggplot2)
  library(automap)
  library(spatstat)
  library(imager)
  library(vec2dtransf)
  library(sp)
  library(sf)
  
  
  space_directory = "~/Google Drive File Stream/My Drive/sciSpace/"
  setwd(dir=space_directory)
  cc.genes <- readRDS("Submission_Data/bin/cc.genes.mouse.RDS")
  
  # Pass TRUE if you want to see progress output on some of Monocle 3's operations
  DelayedArray:::set_verbose_block_processing(TRUE)
  # Passing a higher value will make some computations faster but use more memory. 
  # Adjust with caution!
  options(DelayedArray.block.size=1000e6)
  options(future.globals.maxSize= 891289600)
  set.seed(42)
})



# Functions to interpolate and plot points --------------------------------
nscore <- function(x) {
   # Takes a vector of values x and calculates their normal scores. Returns 
   # a list with the scores and an ordered table of original values and
   # scores, which is useful as a back-transform table. See backtr().
   nscore <- qqnorm(x, plot.it = FALSE)$x  # normal score 
   trn.table <- data.frame(x=sort(x),nscore=sort(nscore))
   
   return (list(nscore=nscore, trn.table=trn.table))
}

backtr <- function(scores, nscore, tails='none', draw=TRUE) {
   # Given a vector of normal scores and a normal score object 
   # (from nscore), the function returns a vector of back-transformed 
   # values. One major issue is how to extrapolate to the tails. Options 
   # other than none may result in dramatically incorrect tail estimates!
   # tails options:
   # 'none' : No extrapolation; more extreme score values will revert 
   # to the original min and max values. 
   # 'equal' : Calculate magnitude in std deviations of the scores about 
   # initial data mean. Extrapolation is linear to these deviations. 
   # will be based upon deviations from the mean of the original 
   # hard data - possibly quite dangerous!
   # 'separate' :  This calculates a separate sd for values 
   # above and below the mean.
   
   if(tails=='separate') { 
      mean.x <- mean(nscore$trn.table$x)
      small.x <- nscore$trn.table$x < mean.x
      large.x <- nscore$trn.table$x > mean.x
      small.sd <- sqrt(sum((nscore$trn.table$x[small.x]-mean.x)^2)/
                       (length(nscore$trn.table$x[small.x])-1))
      large.sd <- sqrt(sum((nscore$trn.table$x[large.x]-mean.x)^2)/
                       (length(nscore$trn.table$x[large.x])-1))
      min.x <- mean(nscore$trn.table$x) + (min(scores) * small.sd)
      max.x <- mean(nscore$trn.table$x) + (max(scores) * large.sd)
      # check to see if these values are LESS extreme than the
      # initial data - if so, use the initial data.
      #print(paste('lg.sd is:',large.sd,'max.x is:',max.x,'max nsc.x is:',max(nscore$trn.table$x)))
      if(min.x > min(nscore$trn.table$x)) {min.x <- min(nscore$trn.table$x)}
      if(max.x < max(nscore$trn.table$x)) {max.x <- max(nscore$trn.table$x)}
   }
   if(tails=='equal') { # assumes symmetric distribution around the mean
      mean.x <- mean(nscore$trn.table$x)
      sd.x <- sd(nscore$trn.table$x)
      min.x <- mean(nscore$trn.table$x) + (min(scores) * sd.x)
      max.x <- mean(nscore$trn.table$x) + (max(scores) * sd.x)
      # check to see if these values are LESS extreme than the
      # initial data - if so, use the initial data.
      if(min.x > min(nscore$trn.table$x)) {min.x <- min(nscore$trn.table$x)}
      if(max.x < max(nscore$trn.table$x)) {max.x <- max(nscore$trn.table$x)}
   }
   if(tails=='none') {   # No extrapolation
      min.x <- min(nscore$trn.table$x)
      max.x <- max(nscore$trn.table$x)
   }
   min.sc <- min(scores)
   max.sc <- max(scores)
   x <- c(min.x, nscore$trn.table$x, max.x)
   nsc <- c(min.sc, nscore$trn.table$nscore, max.sc)
   
   if(draw) {plot(nsc,x, main='Transform Function')}
   back.xf <- approxfun(nsc,x) # Develop the back transform function
   val <- back.xf(scores)
   
   return(val)
}


jitter.points <- function(pts,jitter_amount=1e-6) {
	x = coordinates(pts)
	x =  x + rnorm(length(x),0,jitter_amount) 
	res = SpatialPoints(x)
	proj4string(res)=CRS(proj4string(pts))
	if (class(pts)=="SpatialPointsDataFrame") {
		res = SpatialPointsDataFrame(res,data.frame(pts))}
	return(res) 
}

# Main kriging function
krige_gene_in_cell_type = function(agg_count_cds, 
                                   genes, 
                                   cell_types, 
                                   border_file,
                                   cell_type_column="cell_cluster", 
                                   log_transform=FALSE, 
                                   normal_scores=TRUE, 
                                   pseudocount=1e-3){
    coord_vars = c("x","y")
    data_vars = setdiff(colnames(colData(agg_count_cds)), coord_vars)

    #fit_points = SpatialPoints(as.data.frame(grid_points))

    if (!is.null(length(genes)) && length(genes) >= 2){
      gene_ids = rownames(rowData(agg_count_cds)$gene_short_name %in% genes)
        resid_mat = aggregate_gene_expression(agg_count_cds, gene_ids, norm_method="size_only")
    }else{
        resid_mat = counts(agg_count_cds)[rowData(agg_count_cds)$gene_short_name %in% genes,]
        resid_mat = t(t(resid_mat) / size_factors(agg_count_cds))
    }

    voxels = colData(agg_count_cds)$voxel_id
    expr_from_cell_type = resid_mat[,colData(agg_count_cds)[,cell_type_column] %in% cell_types]
    expr_from_cell_type = as.matrix(t(monocle3:::my.aggregate.Matrix(t(expr_from_cell_type), 
        voxels[colData(agg_count_cds)[,cell_type_column] %in% cell_types])))
    expr_mat_voxel_totals = as.matrix(t(monocle3:::my.aggregate.Matrix(t(resid_mat), voxels)))
    prop_expr_from_cell_type = matrix(0, 
        ncol=ncol(expr_mat_voxel_totals), nrow=nrow(expr_mat_voxel_totals))
    row.names(prop_expr_from_cell_type) = row.names(expr_mat_voxel_totals)
    colnames(prop_expr_from_cell_type) = colnames(expr_mat_voxel_totals)
    expr_from_cell_type = expr_from_cell_type / expr_mat_voxel_totals[,colnames(expr_from_cell_type)]
    prop_expr_from_cell_type[,colnames(expr_from_cell_type)] = expr_from_cell_type
    prop_expr_from_cell_type[is.nan(prop_expr_from_cell_type)] = 0

    coord_df = colData(agg_count_cds) %>% as.data.frame %>% dplyr::select(voxel_id, x, y) %>%
         group_by(voxel_id) %>% summarize(x=mean(x), y=mean(y))  %>% as.data.frame
    coord_df$voxel_id = as.character(coord_df$voxel_id)
    row.names(coord_df) = coord_df$voxel_id
    coord_df = coord_df[colnames(prop_expr_from_cell_type),]

    sp_points = jitter.points(SpatialPoints(as.matrix(coord_df[,coord_vars])))
    
    # Create a fine grid
    pixels_per_side = 200
    bottom.left = apply(sp_points@coords,2,min)
    top.right = apply(sp_points@coords,2,max)
    margin = abs((top.right-bottom.left))/10
    bottom.left = bottom.left-margin
    top.right = top.right+margin
    pixel.size = abs(top.right-bottom.left)/pixels_per_side
    g = GridTopology(cellcentre.offset=bottom.left,
                    cellsize=pixel.size,
                    cells.dim=c(pixels_per_side,pixels_per_side))

    embryo_pg = as(sf::st_geometry(border_file), "Spatial")
    grid_points = SpatialPoints(g)
    in_points = !is.na(over(grid_points,embryo_pg))
    fit_points = SpatialPoints(as.data.frame(grid_points)[in_points,])

    sp_df = SpatialPointsDataFrame(sp_points, as.data.frame(coord_df))

    krig_res = lapply(row.names(resid_mat), function(x){
        resids = prop_expr_from_cell_type[x,] 
        if (log_transform){
            resids = log10(resids + pseudocount)
        }
        if (normal_scores){
            resids = nscore(resids)$nscore
        }

        krig = autoKrige(resids~1, sp_df, new_data=fit_points, remove_duplicates=FALSE)
        return(krig)
    })
    names(krig_res) = row.names(resid_mat)
    return(krig_res)
}


interpolate_genes = function(cds, krig_res){
  krig_res = lapply(names(krig_res), function(x){
      krig = krig_res[[x]]
    interp_data = as.data.frame(krig$krige_output)
    colnames(interp_data) = c("x","y","expr_pred","expr_var","expr_stdev")
    interp_data$feature_id = x
    interp_data$scaled_pred = 100 * (interp_data$expr_pred - min(interp_data$expr_pred)) / 
      (max(interp_data$expr_pred) - min(interp_data$expr_pred))
    return(interp_data)
  })
  krig_res = do.call(rbind, krig_res)
  return(krig_res)
}



# Prepare CDS by loading spatial position into "tSNE" slot in reduced dims ---------------

all_image_data = 
  readRDS(file = "Submission_Data/E14_slides/RDS_intermediates/Notebook0_1_images_and_transformations.RDS")

slide_1D_border_file = all_image_data$slide_polygon[[1]]

cds = 
  readRDS(file = "Submission_Data/E14_slides/RDS_intermediates/Notebook6_spatial_cds_anatomy.RDS")

# Load the spatial coordinates of every cell into the "tSNE" dimensionality reduction slot 
reducedDims(cds)[["tSNE"]] = 
  matrix(cbind(colData(cds)$coords.x1, 
               colData(cds)$coords.x2), 
         ncol=2)

rownames(reducedDims(cds)[["tSNE"]]) = 
  rownames(reducedDims(cds)[["UMAP"]])

cds@clusters[["tSNE"]] = cds@clusters[["UMAP"]]

colData(cds)$cell_cluster = clusters(cds)
colData(cds)$cell_type = colData(cds)$final_cluster_label

min_fraction_cells = 0.01
cds = detect_genes(cds)
expressed_genes = 
  row.names(subset(rowData(cds), 
                   num_cells_expressed > ncol(cds)*min_fraction_cells))
print(length(expressed_genes))
cds = cds[expressed_genes,]

# Collapse spatial bins to aggregate more counts per celltype -------------

bin_size_um = 100
colData(cds)$trans_x_bin = floor((colData(cds)$coords.x1 - min(colData(cds)$coords.x1))  / bin_size_um)
colData(cds)$trans_y_bin = floor((colData(cds)$coords.x2 - min(colData(cds)$coords.x2))  / bin_size_um)
colData(cds)$voxel_id = interaction(colData(cds)$max_slide_id, colData(cds)$trans_x_bin, colData(cds)$trans_y_bin, drop=TRUE)
colData(cds)$voxel_id.cell_cluster = interaction(colData(cds)$cell_cluster, colData(cds)$voxel_id, drop=TRUE)

# Do not collapse cells of different types:
spot_summary = colData(cds) %>% 
  as.data.frame %>% 
  filter(is.na(voxel_id) == FALSE) %>%
  group_by(max_slide_id, voxel_id, cell_cluster) %>% 
  summarize(x = mean(coords.x1, na.rm=TRUE), 
            y = mean(coords.x2, na.rm=TRUE), 
            total_cells = n())


slice_cds = cds[,colData(cds)$max_slide_id == "slide_1D"]
colData(slice_cds)$cell_cluster = droplevels(colData(slice_cds)$cell_cluster)

count_helper <- function(gene_id, cds, agg_grouping, pseudocount){
tryCatch({
    orig_data = as.vector(counts(cds)[gene_id,]) / size_factors(cds)
    gene_counts = orig_data

    agg_residuals = lapply(split(gene_counts, agg_grouping), sum) 
    group_names = names(agg_residuals)
    agg_residuals = as(as.numeric(agg_residuals), "sparseMatrix")
    row.names(agg_residuals) = group_names

    return(agg_residuals)
}, error = function(e){
    print (e)
    retval = rep_len(NA, length(levels(agg_grouping)))
    names(retval) = levels(agg_grouping)
    return(retval)
})
}

# Function to aggregate cells of the same type in the same spatial bin
aggregate_counts <- function(cds, aggregation_col, pseudocount=0.1) {
    require(speedglm)
    agg_grouping = factor(colData(cds)[,aggregation_col])
    model_tbl <- rowData(cds) %>% as.data.frame %>%
        dplyr::mutate(predictions = furrr::future_map(id, count_helper, cds, agg_grouping, pseudocount))
    pred_matrix = t(do.call(cbind, model_tbl$predictions))
    return(pred_matrix)
}

# Aggregate counts at each position for each celltype ---------------------
plan(multicore, workers=4)

# Temporarily disable OpenMP threading in functions to be run in parallel
old_omp_num_threads = as.numeric(Sys.getenv("OMP_NUM_THREADS"))

if (is.na(old_omp_num_threads)){old_omp_num_threads = 1}
RhpcBLASctl::omp_set_num_threads(1)

# Temporarily set the number of threads the BLAS library can use to be 1
old_blas_num_threads = as.numeric(Sys.getenv("OPENBLAS_NUM_THREADS"))

if (is.na(old_omp_num_threads)){old_blas_num_threads = 1}
RhpcBLASctl::blas_set_num_threads(1)

# Aggregate counts at each position for each celltype
agg_counts = aggregate_counts(slice_cds, aggregation_col = "voxel_id.cell_cluster")

Sys.setenv(OMP_NUM_THREADS = 1)
RhpcBLASctl::omp_set_num_threads(old_omp_num_threads)
RhpcBLASctl::blas_set_num_threads(old_blas_num_threads)

spatal_bin_metadata = as.data.frame(spot_summary)
spatal_bin_metadata$voxel_id.cell_cluster = paste(spatal_bin_metadata$cell_cluster, spatal_bin_metadata$voxel_id, sep=".")
row.names(spatal_bin_metadata) = spatal_bin_metadata$voxel_id.cell_cluster
spatal_bin_metadata = spatal_bin_metadata[colnames(agg_counts),]

# Make a new cell data set that holds voxels instead of cells
# with appropriate metadata
agg_count_cds = new_cell_data_set(agg_counts, 
                                     gene_metadata = rowData(slice_cds)%>% as.data.frame,
                                     cell_metadata = spatal_bin_metadata)

reducedDims(agg_count_cds)[["tSNE"]] = 
  matrix(cbind(colData(agg_count_cds)$x, 
               colData(agg_count_cds)$y), 
         ncol=2)
rownames(reducedDims(agg_count_cds)[["tSNE"]]) = colnames(agg_count_cds)

count_summary = rowData(agg_count_cds)
count_summary$total_res = rowMeans(counts(agg_count_cds))
count_summary %>% as.data.frame %>% arrange(total_res) %>% tail(n=25)

colData(agg_count_cds)$cell_cluster = droplevels(colData(agg_count_cds)$cell_cluster)



# Krig gene expression for genes in two celltypes -------------------------
nbin = 10

figure2_genes = c("Fgfr1","Fgfr2","Fgf9","Fgfr3")


cardiomyocyte_clusters = 
  which(prop.table(table(colData(cds)$cell_cluster, 
                         colData(cds)$cell_type), margin=1)[,"Cardiac muscle lineages"] > 0.15) %>%
  as.character()

fig2_cardio_gene_krig_res = 
  krige_gene_in_cell_type(agg_count_cds = agg_count_cds,
                          genes = figure2_genes, 
                          border_file = slide_1D_border_file,
                          cell_types= cardiomyocyte_clusters, 
                          log_transform=FALSE, 
                          normal_score=FALSE)

fig2_cardio_interpolated_res = 
  interpolate_genes(agg_count_cds, 
                    fig2_cardio_gene_krig_res)

fig2_cardio_interpolated_res = 
  inner_join(fig2_cardio_interpolated_res, 
             as.data.frame(rowData(agg_count_cds)), 
             by=c("feature_id"="id"))

fig2_cardio_interpolated_res$feature_id = 
  fig2_cardio_interpolated_res$gene_short_name

fig2_cardio_interpolated_res$expr_pred[fig2_cardio_interpolated_res$expr_pred < 0] = 0

endothelial_clusters = 
  which(prop.table(table(colData(cds)$cell_cluster, 
                         colData(cds)$cell_type), 
                   margin=1)[,"Endothelial Cells"] > 0.15)

fig2_endothelial_gene_krig_res = 
  krige_gene_in_cell_type(agg_count_cds = agg_count_cds, 
                          genes = figure2_genes, border_file = slide_1D_border_file,
                          cell_types = endothelial_clusters, 
                          log_transform=FALSE, 
                          normal_score=FALSE)

fig2_endothelial_interpolated_res = 
  interpolate_genes(agg_count_cds, 
                    fig2_endothelial_gene_krig_res)

fig2_endothelial_interpolated_res = 
  inner_join(fig2_endothelial_interpolated_res, 
             as.data.frame(rowData(agg_count_cds)), 
             by=c("feature_id"="id"))

fig2_endothelial_interpolated_res$feature_id = 
  fig2_endothelial_interpolated_res$gene_short_name

fig2_endothelial_interpolated_res$expr_pred[fig2_endothelial_interpolated_res$expr_pred < 0] = 0

# Figure 2 — Panel F ------------------------------------------------------

max_Fgfr2 = max(c(fig2_cardio_interpolated_res %>%
                    filter(feature_id == "Fgfr2") %>% 
                    pull(expr_pred),
                  fig2_endothelial_interpolated_res %>% 
                    filter(feature_id == "Fgfr2") %>% 
                    pull(expr_pred))) * 100

max_Fgfr1 = max(c(fig2_cardio_interpolated_res %>%
                    filter(feature_id == "Fgfr1") %>% 
                    pull(expr_pred),
                  fig2_endothelial_interpolated_res %>% 
                    filter(feature_id == "Fgfr1") %>% 
                    pull(expr_pred))) * 100

max_Fgfr3 = max(c(fig2_cardio_interpolated_res %>%
                    filter(feature_id == "Fgfr3") %>% 
                    pull(expr_pred),
                  fig2_endothelial_interpolated_res %>% 
                    filter(feature_id == "Fgfr3") %>% 
                    pull(expr_pred))) * 100

max_Fgf9 = max(c(fig2_cardio_interpolated_res %>%
                   filter(feature_id == "Fgf9") %>% 
                   pull(expr_pred),
                 fig2_endothelial_interpolated_res %>% 
                   filter(feature_id == "Fgf9") %>% 
                   pull(expr_pred))) * 100


png("Figures/Figure_Components/Figure2/embryo_cardiac_Fgfr2.png", width=1.5, height=2, units="in", res = 300)
ggplot() +
  geom_sf(data = all_image_data$slide_polygon[[1]] * rbind(c(0, -1),c(1,0)),
          fill = NA,
          size = .25) +
  geom_tile(data= fig2_cardio_interpolated_res %>% 
              filter(feature_id == "Fgfr2"),
            aes(x=y, 
                y=-x, 
                fill=expr_pred * 100),
            color=NA) +
  stat_contour(data= fig2_cardio_interpolated_res %>% 
                 filter(feature_id == "Fgfr2"),
               aes(x=y, 
                   y=-x, 
                   z=expr_pred * 100), 
               bins=nbin, size=I(0.1), 
               color="white") +
  viridis::scale_fill_viridis(name="% contributed",
                              option="B",
                              limits = c(0, max_Fgfr2)) + 
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "none",
        legend.title = element_text(size =8),
        legend.text = element_text(size = 6),
        strip.text = element_blank()) 
dev.off()



png("Figures/Figure_Components/Figure2/embryo_endo_Fgfr2.png", width=1.5, height=2, units="in", res = 300)
ggplot() +
  geom_sf(data = all_image_data$slide_polygon[[1]] * rbind(c(0, -1),c(1,0)),
          fill = NA,
          size = .25) +
  geom_tile(data= fig2_endothelial_interpolated_res %>% 
              filter(feature_id == "Fgfr2"),
            aes(x=y, 
                y=-x, 
                fill=expr_pred * 100),
            color=NA) +
  stat_contour(data= fig2_endothelial_interpolated_res %>% 
                 filter(feature_id == "Fgfr2"),
               aes(x=y, 
                   y=-x, 
                   z=expr_pred * 100), 
               bins=nbin, size=I(0.1), 
               color="white") +
  viridis::scale_fill_viridis(name="% contributed",
                              option="B",
                              limits = c(0, max_Fgfr2)) + 
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "none",
        legend.title = element_text(size =8),
        legend.text = element_text(size = 6),
        strip.text = element_blank()) 
dev.off()


png("Figures/Figure_Components/Figure2/embryo_cardiac_Fgfr3.png", width=1.5, height=2, units="in", res = 300)
ggplot() +
  geom_sf(data = all_image_data$slide_polygon[[1]] * rbind(c(0, -1),c(1,0)),
          fill = NA,
          size = .25) +
  geom_tile(data= fig2_cardio_interpolated_res %>% 
              filter(feature_id == "Fgfr3"),
            aes(x=y, 
                y=-x, 
                fill=expr_pred * 100),
            color=NA) +
  stat_contour(data= fig2_cardio_interpolated_res %>% 
                 filter(feature_id == "Fgfr3"),
               aes(x=y, 
                   y=-x, 
                   z=expr_pred * 100), 
               bins=nbin, size=I(0.1), 
               color="white") +
  viridis::scale_fill_viridis(name="% contributed",
                              option="B",
                              limits = c(0, max_Fgfr3)) + 
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "none",
        legend.title = element_text(size =8),
        legend.text = element_text(size = 6),
        strip.text = element_blank()) 
dev.off()



png("Figures/Figure_Components/Figure2/embryo_endo_Fgfr3.png", width=1.5, height=2, units="in", res = 300)
ggplot() +
  geom_sf(data = all_image_data$slide_polygon[[1]] * rbind(c(0, -1),c(1,0)),
          fill = NA,
          size = .25) +
  geom_tile(data= fig2_endothelial_interpolated_res %>% 
              filter(feature_id == "Fgfr3"),
            aes(x=y, 
                y=-x, 
                fill=expr_pred * 100),
            color=NA) +
  stat_contour(data= fig2_endothelial_interpolated_res %>% 
                 filter(feature_id == "Fgfr3"),
               aes(x=y, 
                   y=-x, 
                   z=expr_pred * 100), 
               bins=nbin, size=I(0.1), 
               color="white") +
  viridis::scale_fill_viridis(name="% contributed",
                              option="B",
                              limits = c(0, max_Fgfr3)) + 
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "none",
        legend.title = element_text(size =8),
        legend.text = element_text(size = 6),
        strip.text = element_blank()) 
dev.off()

png("Figures/Figure_Components/Figure2/embryo_cardiac_Fgfr1.png", width=1.5, height=2, units="in", res = 300)
ggplot() +
  geom_sf(data = all_image_data$slide_polygon[[1]] * rbind(c(0, -1),c(1,0)),
          fill = NA,
          size = .25) +
  geom_tile(data= fig2_cardio_interpolated_res %>% 
              filter(feature_id == "Fgfr1"),
            aes(x=y, 
                y=-x, 
                fill=expr_pred * 100),
            color=NA) +
  stat_contour(data= fig2_cardio_interpolated_res %>% 
                 filter(feature_id == "Fgfr1"),
               aes(x=y, 
                   y=-x, 
                   z=expr_pred * 100), 
               bins=nbin, size=I(0.1), 
               color="white") +
  viridis::scale_fill_viridis(name="% contributed",
                              option="B",
                              limits = c(0, max_Fgfr1)) + 
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "none",
        legend.title = element_text(size =8),
        legend.text = element_text(size = 6),
        strip.text = element_blank()) 
dev.off()


png("Figures/Figure_Components/Figure2/embryo_endo_Fgfr1.png", width=1.5, height=2, units="in", res = 300)
ggplot() +
  geom_sf(data = all_image_data$slide_polygon[[1]] * rbind(c(0, -1),c(1,0)),
          fill = NA,
          size = .25) +
  geom_tile(data= fig2_endothelial_interpolated_res %>% 
              filter(feature_id == "Fgfr1"),
            aes(x=y, 
                y=-x, 
                fill=expr_pred * 100),
            color=NA) +
  stat_contour(data= fig2_endothelial_interpolated_res %>% 
                 filter(feature_id == "Fgfr1"),
               aes(x=y, 
                   y=-x, 
                   z=expr_pred * 100), 
               bins=nbin, size=I(0.1), 
               color="white") +
  viridis::scale_fill_viridis(name="% contributed",
                              option="B",
                              limits = c(0, max_Fgfr1)) + 
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "none",
        legend.title = element_text(size =8),
        legend.text = element_text(size = 6),
        strip.text = element_blank()) 
dev.off()




png("Figures/Figure_Components/Figure2/embryo_cardiac_Fgf9.png", width=1.5, height=2, units="in", res = 300)
ggplot() +
  geom_sf(data = all_image_data$slide_polygon[[1]] * rbind(c(0, -1),c(1,0)),
          fill = NA,
          size = .25) +
  geom_tile(data= fig2_cardio_interpolated_res %>% 
              filter(feature_id == "Fgf9"),
            aes(x=y, 
                y=-x, 
                fill=expr_pred * 100),
            color=NA) +
  stat_contour(data= fig2_cardio_interpolated_res %>% 
                 filter(feature_id == "Fgf9"),
               aes(x=y, 
                   y=-x, 
                   z=expr_pred * 100), 
               bins=nbin, size=I(0.1), 
               color="white") +
  viridis::scale_fill_viridis(name="% contributed",
                              option="B",
                              limits = c(0, max_Fgf9)) + 
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "none",
        legend.title = element_text(size =8),
        legend.text = element_text(size = 6),
        strip.text = element_blank())
dev.off()


png("Figures/Figure_Components/Figure2/embryo_endo_Fgf9.png", width=1.5, height=2, units="in", res = 300)
ggplot() +
  geom_sf(data = all_image_data$slide_polygon[[1]] * rbind(c(0, -1),c(1,0)),
          fill = NA,
          size = .55) +
  geom_tile(data= fig2_endothelial_interpolated_res %>% 
              filter(feature_id == "Fgf9"),
            aes(x=y, 
                y=-x, 
                fill=expr_pred * 100),
            color=NA) +
  stat_contour(data= fig2_endothelial_interpolated_res %>% 
                 filter(feature_id == "Fgf9"),
               aes(x=y, 
                   y=-x, 
                   z=expr_pred * 100), 
               bins=nbin, size=I(0.1), 
               color="white") +
  viridis::scale_fill_viridis(name="% contributed",
                              option="B",
                              limits = c(0, max_Fgf9)) + 
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "none",
        legend.title = element_text(size =8),
        legend.text = element_text(size = 6),
        strip.text = element_blank())
dev.off()


png("Figures/Figure_Components/Figure2/embryo_endothelial_genes_counts.png", width=2, height=1, units="in", res = 300)
ggplot() +
  geom_sf(data = all_image_data$slide_polygon[[1]] * rbind(c(0, -1),c(1,0)),
          fill = NA,
          size = .25) +
  geom_tile(aes(x=y, y=-x, fill=expr_pred * 100),color=NA, data=fig2_endothelial_interpolated_res) +
  stat_contour(aes(x=y, y=-x, z=expr_pred * 100), bins=nbin, size=I(0.1), color="white", data=fig2_endothelial_interpolated_res) +
  viridis::scale_fill_viridis(name="% contributed", option="B") + 
  facet_wrap(~feature_id, 
             nrow = 1) + 
  monocle3:::monocle_theme_opts() +
  theme_void() +
  theme(legend.position = "none",
        legend.title = element_text(size =8),
        legend.text = element_text(size = 6),
        strip.text = element_blank())
dev.off()

