import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

adMerge= sc.read('./admk_sestan2021_500Cell.h5')

plot_matrix = two_species_heatmap(adMerge, species_1 = 'macaque_Hao', species_2 = 'macaque_Sestan',species_1_key = 'label', species_2_key = 'cluster',                        louvain = 2.7,figure_path = 'test_heatmap_glia.png')


df_norm_row = plot_matrix.apply(lambda x: (x-x.mean())/x.std(), axis = 1)
df_norm_row = df_norm_row/(df_norm_row.max().max())
sns.heatmap(df_norm_row, cmap='Greys', cbar=True, xticklabels=1,yticklabels=1, linewidth = 1,
            linecolor = 'gray', vmax = 1, vmin = 0)

plt.savefig('./Franjic_heat.pdf',dpi = 600, )
plot_matrix_plot.to_csv('./Franjic_100.csv')


admk = adMerge[adMerge.obs['species'] =='macaque_Hao',:]


matplotlib.rcParams.update({'font.size': 10})
fig = plt.figure(constrained_layout=True,figsize=(6,12))

gs = GridSpec(3, 1, figure=fig)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[1, 0])
ax3 = fig.add_subplot(gs[2, 0])
ax1 = sc.pl.umap(adMerge[plot_names_sestan], color=['species'],ax=ax1,           legend_fontsize=10,title=' ', show=False, size=50, frameon = False)

ax2 = sc.pl.umap(admk, color=['label_'],ax=ax2,legend_fontoutline = 2,           legend_fontsize=10,title=' ', show=False, size=15, frameon = False,)

ax3 = sc.pl.umap(adMerge[plot_names_sestan], color = 'Sestan_label',ax=ax3,legend_fontoutline = 1,           legend_fontsize=10,title=' ', show=False, size=110, frameon = False)


fig.tight_layout(w_pad=0.3)
plt.savefig("./Fig_Franjic_umap_label.pdf", dpi = 600)

