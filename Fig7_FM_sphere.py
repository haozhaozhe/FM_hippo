import pickle
import scanpy as sc

#load data
FM_sphere=sc.read('FM_sphere_6988_wk.h5')

#plot SCCAF

matplotlib.rcParams.update({'font.size': 10})
fig = plt.figure(constrained_layout=True,figsize=(10,3.5))
gs = GridSpec(2,7, figure=fig)
ax1 = fig.add_subplot(gs[0:2, 0:2])
ax2 = fig.add_subplot(gs[0, 2])
ax3 = fig.add_subplot(gs[0,3 ])
ax4 = fig.add_subplot(gs[0,4 ])
ax5 = fig.add_subplot(gs[0,5 ])
ax6 = fig.add_subplot(gs[0,6 ])
ax7 = fig.add_subplot(gs[1,2 ])
ax8 = fig.add_subplot(gs[1,3 ])
ax9 = fig.add_subplot(gs[1,4 ])
ax10 = fig.add_subplot(gs[1,5 ])
ax11 = fig.add_subplot(gs[1,6 ])

ax1 = sc.pl.umap(FM_sphere, color = 'marker1_predict_harmony', title = 'SCCAF Projection', ax = ax1,legend_loc = None, frameon = False, show = False)
ax2 = sc.pl.umap(FM_sphere, color = 'IPC_1', title = 'IPC_1',  ax = ax2,legend_loc = None, frameon = False, show = False)

ax3 = sc.pl.umap(FM_sphere, color = 'Neuroblast', title = 'Neuroblast',ax = ax3,legend_loc = None, frameon = False, show = False)
ax4 = sc.pl.umap(FM_sphere, color = 'Granule_immature', title = 'GC_immature',  ax = ax4,legend_loc = None, frameon = False, show = False)
ax5 = sc.pl.umap(FM_sphere, color = 'Granule_mature', title = 'GC_mature', ax = ax5,legend_loc = None, frameon = False, show = False)
ax6 = sc.pl.umap(FM_sphere, color = 'Granule_mature_1', title = 'GC_mature_1', ax = ax6,legend_loc = None, frameon = False, show = False)
ax7 = sc.pl.umap(FM_sphere, color = 'RGL_2', title = 'RGL_2',  ax = ax7,legend_loc = None, frameon = False, show = False)
ax8 = sc.pl.umap(FM_sphere, color = 'Astro_im1', title = 'Astro_im1', ax = ax8,legend_loc = None, frameon = False, show = False)
ax9 = sc.pl.umap(FM_sphere, color = 'Astro_im2', title = 'Astro_im2', ax = ax9,legend_loc = None, frameon = False, show = False)
ax10 = sc.pl.umap(FM_sphere, color = 'Astro_1', title = 'Astro_1', ax = ax10,legend_loc = None, frameon = False, show = False)
ax11 = sc.pl.umap(FM_sphere, color = 'Astro_2', title = 'Astro_2', ax = ax11,legend_loc = None, frameon = False, show = False)


#plot gene expression
matplotlib.rcParams.update({'font.size': 14})
fig = plt.figure(figsize=(13,5))

gs = GridSpec(2,4, figure=fig)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[0, 2])
ax4 = fig.add_subplot(gs[0, 3])

ax5 = fig.add_subplot(gs[1, 0])
ax6 = fig.add_subplot(gs[1, 1])
ax7 = fig.add_subplot(gs[1, 2])
ax8 = fig.add_subplot(gs[1, 3])

gl = ['HMGB2','MKI67','NNAT','SOX4', 'SLC1A3','SDC4','RSPO3','PAX6']

ax1 = sc.pl.umap(FM_sphere, color=gl[0],ax=ax1,legend_loc = 'none',cmap = 'RdGy_r',  legend_fontsize=10,title=gl[0], show=False, frameon = False, s = 12)

ax2 = sc.pl.umap(FM_sphere, color=gl[1],ax=ax2,legend_loc = 'none',cmap = 'RdGy_r',  legend_fontsize=10,title=gl[1], show=False, frameon = False,s = 12)

ax3 = sc.pl.umap(FM_sphere, color=gl[2],ax=ax3,legend_loc = 'none',cmap = 'RdGy_r', legend_fontsize=10,title=gl[2], show=False,  frameon = False,s = 12)

ax4 = sc.pl.umap(FM_sphere, color=gl[3],ax=ax4,legend_loc = 'none',cmap = 'RdGy_r', legend_fontsize=10,title=gl[3], show=False,  frameon = False,s = 12)

ax5 = sc.pl.umap(FM_sphere, color=gl[4],ax=ax5,legend_loc = 'none',cmap = 'RdGy_r',legend_fontsize=10,title=gl[4], show=False,  frameon = False,s = 12)

ax6 = sc.pl.umap(FM_sphere, color=gl[5],ax=ax6,legend_loc = 'none',cmap = 'RdGy_r', legend_fontsize=10,title=gl[5], show=False,  frameon = False,s = 12)

ax7 = sc.pl.umap(FM_sphere, color=gl[6],ax=ax7,legend_loc = 'none',cmap = 'RdGy_r', legend_fontsize=10,title=gl[6], show=False,  frameon = False,s = 12)

ax8 = sc.pl.umap(FM_sphere, color=gl[7],ax=ax8,legend_loc = 'none',cmap = 'RdGy_r', legend_fontsize=10,title=gl[7], show=False, frameon = False,s = 12)




