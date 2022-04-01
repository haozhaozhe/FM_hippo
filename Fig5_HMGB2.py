import scanpy as sc


#load data

ad = sc.read("data/FM_hippo_207785.h5")
sc.pl.umap(ad, color = 'label')

#plot genes
matplotlib.rcParams.update({'font.size': 14})
fig = plt.figure(figsize=(9,2.6))
gs = GridSpec(1,3, figure=fig)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[0, 2])
gl =[FM_key,'EZH2','MKI67',]

ax1 = sc.pl.umap(ad, color=gl[0],ax=ax1,legend_loc = 'none',show=False,  frameon = False, s = 2,legend_fontweight = 'bold')

ax2 = sc.pl.umap(ad, color='EZH2',ax=ax2,legend_loc = 'none',cmap = 'RdGy_r',            legend_fontsize=10,title='EZH2', show=False, frameon = False,s = 1.7)

ax3 = sc.pl.umap(ad, color='MKI67',ax=ax3,legend_loc = 'none',cmap = 'RdGy_r',            legend_fontsize=10,title='MKI67', show=False,  frameon = False,s = 2)


#plot violin
figsize(5,3)

plt.figure()
ad_sub = ad[ad.obs['label'].isin(['RGL','IPC_1', 'IPC_2']),:]
color_chart = dict(zip(ad_sub.obs['label'].cat.categories, ad_sub.uns['label_colors']))
sc.pl.violin(ad_sub, keys = ['GFAP','HMGB2'],jitter = False,                      groupby = 'label', swap_axes = True, palette= ad_sub.uns['label_colors'], show = False)
