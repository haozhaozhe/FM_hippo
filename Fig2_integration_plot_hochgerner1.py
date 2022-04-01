import matplotlib.pyplot as plt
import scanpy as sc

import sys
sys.path.append('scctools/scctools_1.0')
from scctools import *

#load data
hoch1_500= sc.read('data/Data_integration/admk_Hoch1_500_finish.h5')
hoch1_500.obs_names = [x[:-2] for x in list(hoch1_500.obs_names)]

#plot heatmpa
plot_matrix_Hoch1_500 = two_species_heatmap(hoch1_500, species_1 = 'monkey', species_2 = 'mouse',species_1_key = 'label_14_merge', species_2_key = 'characteristics: cell cluster',                        louvain = 1.5,figure_path = 'test_heatmap.png')

sns_plot = sns.heatmap(plot_matrix_Hoch1_500_plot, cmap='Greys', cbar=True, xticklabels=1,yticklabels=1, linewidth = 1, linecolor = 'gray')
plot_matrix_Hoch1_500.to_csv('Hoch1_500_plot_matrix.csv')

#plot umap
admk = hoch1_500[hoch1_500.obs['species'] =='monkey',:]
admm = hoch1_500[hoch1_500.obs['species'] =='mouse',:]

fig = plt.figure(constrained_layout=True,figsize=(5,13))
gs = GridSpec(3, 1, figure=fig)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[1,0])
ax3 = fig.add_subplot(gs[2, 0])
ax1 = sc.pl.umap(hoch1_500, color=['species'],ax=ax1, legend_fontsize=10,title=' ', show=False, size=10, frameon = False)

ax2 = sc.pl.umap(admk, color=['label'],ax=ax2,legend_loc = 'on data',legend_fontoutline = 2,legend_fontsize=10,title=' ', show=False, size=10, frameon = False,)

ax3 = sc.pl.umap(admm, color=['label'],ax=ax3,legend_loc = 'on data',legend_fontoutline = 2,legend_fontsize=10,title=' ', show=False, size=10, frameon = False)

