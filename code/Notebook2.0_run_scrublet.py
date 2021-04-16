# Run Scrublet from AllonKleinLab
# https://github.com/AllonKleinLab/scrublet

import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import scrublet as scr


data_path = "/Volumes/GoogleDrive/My Drive/sciSpace/Submission_Data/E14_slides/RDS_intermediates/"

os.chdir(data_path)

counts_matrix = scipy.io.mmread("Notebook1_E14_count_matrix.mtx").T

scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.05)


doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=100)

aux_figures_path = "/Volumes/GoogleDrive/My Drive/sciSpace/Submission_Data/E14_slides/Auxillary_Figures/"
os.chdir(aux_figures_path)

scrub.plot_histogram()[0].savefig("scrublet_histogram.png");

print('Running UMAP...')
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
    
print('Done.')


scrub.plot_embedding('UMAP', order_points=True)[0].savefig("UMAP_doublets.png");

# Write out predicted doublet scores 
os.chdir(data_path)
np.savetxt(fname = "scrublet_doublet_scores.txt", X = doublet_scores)
np.savetxt(fname = "scrublet_predicted_doublets.txt", X = predicted_doublets)