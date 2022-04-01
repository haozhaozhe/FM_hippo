import scanpy as sc
import matplotlib.pyplot as plt
import scvelo as scv

#load data
ad = sc.read('FM_hippo_207785.h5')

figsize(5,5)
sc.pl.umap(ad, color = 'label')

#plot umap

matplotlib.rcParams.update({'font.size': 10})
fig = plt.figure(constrained_layout=True,figsize=(5,5))
#sns.despine()

gs = GridSpec(1, 1, figure=fig)
ax2 = fig.add_subplot(gs[0, 0])

ax2 = sc.pl.umap(ad, color=['show_GC'],ax=ax2,legend_loc = None, legend_fontsize=10,title=' ', show=False, size=5, frameon = False)

#plot velocity
ad_wk = sc.read('toGC1.h5')

figsize(4,4)
scv.pl.velocity_embedding_stream(ad_wk, basis='umap', color = 'label',legend_loc = None,                                 title = ' ', alpha = 0.1, s = 50, legend_fontoutline = 3,                                  save  = "Fig_GC_velo.png",dpi = 600)

scv.tl.rank_velocity_genes(ad_zc, groupby='label', min_corr=.3)

df = scv.DataFrame(ad_zc.uns['rank_velocity_genes']['names'])
df[0:10]

scv.tl.rank_velocity_genes(ad_wk, groupby='label', min_corr=.3)

df = scv.DataFrame(ad_wk.uns['rank_velocity_genes']['names'])
df.to_csv('NB_GC_rank_velocity_genes.csv')

top_genes = ad_wk.var['fit_likelihood'].sort_values(ascending=False).index[:1000]

top_genes = ad_wk.var['fit_likelihood'].sort_values(ascending=False).index[:1000]
scv.pl.heatmap(ad_wk, var_names=top_genes[0:1000], sortby='latent_time', col_color='label', n_convolve=100,                figsize = (10,5), save = 'GC_velo_heat.png')

matplotlib.rcParams.update({'font.size': 10})
fig = plt.figure(figsize=(2.5,11))
#sns.despine()

gs = GridSpec(5,1, figure=fig)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[1, 0])
ax3 = fig.add_subplot(gs[2, 0])
ax4 = fig.add_subplot(gs[3, 0])
ax5 = fig.add_subplot(gs[4, 0])

cm = 'RdGy_r'

ax1 = sc.pl.umap(ad_zc, color = 'ASAP2', title = 'ASAP2',
                    color_map=cm,ax = ax1, frameon = False, show = False, vmin = -1, vmax = 2)

ax2 = sc.pl.umap(ad_zc, color = 'TMEM159',title = 'TMEM159', 
                    color_map=cm,ax = ax2, frameon = False, show = False, vmin = -0.7, vmax = 2, s = 20)

ax3 = sc.pl.umap(ad_zc, color = 'COTL1',title = 'COTL1',
                    color_map=cm,ax = ax3, frameon = False, show = False, vmin = -0.5, vmax = 1.5)

ax4 = sc.pl.umap(ad_zc, color = 'NOV',title = 'NOV',
                    color_map=cm,ax = ax4, frameon = False, show = False, vmin = -0.7, vmax = 1.5)

ax5 = sc.pl.umap(ad_zc, color = 'BCL6',title = 'BCL6',
                    color_map=cm,ax = ax5, frameon = False, show = False, vmin = -2, vmax = 3.5, s = 10)

plot_matrix = pd.read_csv('plot_matrix_GC.csv', index_col = 0)

plot_gene = list(plot_matrix.index)

GC = ad[ad.obs['label'].isin(['Granule_1','Granule_immature']),:]

sns.heatmap(plot_matrix, vmax = 1,vmin = -1, cmap = 'seismic',robust = True)

plot_matrix_norm1 = (plot_matrix.T/plot_matrix.T.max()).T
plot_matrix_norm = plot_matrix_norm1/plot_matrix_norm1.max()

plt.figure()
figsize(2,4)
sns_fig = sns.heatmap(plot_matrix, vmax = 1,vmin = -1, cmap = 'seismic',robust = True)

