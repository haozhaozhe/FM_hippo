import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

import scvelo as scv

#load data
ad = sc.read('data/FM_hippo_207785.h5')

#plot umap
fig = plt.figure(constrained_layout=True,figsize=(5,5))
gs = GridSpec(1, 1, figure=fig)
ax2 = fig.add_subplot(gs[0, 0])
ax2 = sc.pl.umap(ad, color=['show_Astro'],ax=ax2,legend_loc = None,           legend_fontsize=10,title=' ', show=False, size=5, frameon = False)

#plot velocity
ad_zc = sc.read('toAstro_v.h5')

scv.tl.latent_time(ad_zc)

fig = plt.figure(figsize=(14.2,7.2))

gs = GridSpec(2,3, figure=fig)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[0, 2])
ax4 = fig.add_subplot(gs[1, 0])
ax5 = fig.add_subplot(gs[1, 1])
ax6 = fig.add_subplot(gs[1, 2])

cm = 'RdGy_r'

ax1 = sc.pl.umap(ad_zc, color = 'ID3',title = 'ID3',
                    color_map=cm,ax = ax1, frameon = False, show = False,  s = 20)
ax2 = sc.pl.umap(ad_zc, color = 'CDK4',title = 'CDK4', 
                    color_map=cm,ax = ax2, frameon = False, show = False, s = 40)
ax3 = sc.pl.umap(ad_zc, color = 'VIM',title = 'VIM',
                    color_map=cm,ax = ax3, frameon = False, show = False,  s = 20)
ax4 = sc.pl.umap(ad_zc, color = 'HTRA1',title = 'HTRA1',
                    color_map=cm,ax = ax4, frameon = False, show = False,  s = 20)
ax5 = sc.pl.umap(ad_zc, color = 'ACSL6',title = 'ACSL6',
                    color_map=cm,ax = ax5, frameon = False, show = False,  s = 10)
ax6 = sc.pl.umap(ad_zc, color = 'NOTCH2',title = 'NOTCH2',
                    color_map=cm,ax = ax6, frameon = False, show = False, s = 10)

scv.pl.velocity_embedding_stream(ad_zc, basis='umap', color = 'label',                                 title = ' ', alpha = 0.1, s = 50, legend_fontoutline = 3,                               dpi = 600, legend_fontsize = 16)


scv.tl.latent_time(ad_zc)

fig = plt.figure(figsize=(5.5,5))
gs = GridSpec(1, 1, figure=fig)
ax2 = fig.add_subplot(gs[0, 0])

ax2 = scv.pl.scatter(ad_zc, color='latent_time', color_map='coolwarm',                     title = ' ', size=20, ax = ax2, show = False, vmin =-0, vmax = 0.9)

scv.tl.rank_velocity_genes(ad_zc, groupby='label', min_corr=.8)

df = scv.DataFrame(ad_zc.uns['rank_velocity_genes']['names'])
top_genes = ad_zc.var['fit_likelihood'].sort_values(ascending=False).index[:1000]
scv.pl.heatmap(ad_zc, var_names=top_genes[1:1000], sortby='latent_time', col_color='label', n_convolve=100,                figsize = (10,5))

