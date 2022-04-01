import scanpy as sc
#load data
ad = sc.read('FM_hippo_207785.h5')

#plot umap
sc.pl.umap(ad, color = 'label', s= 3, frameon = False)

#plot violin
RGL_marker = ['SLC1A3','SOX2','GFAP']
RGL_young_marker = ['DBI','FABP7']
aNSC_marker = ['HMGB2']
NB_marker = ['SOX4','NNAT']
MicroGlia_marker = ['C1QA',]
CA3_Pyr_marker = ['RGS4','ELAVL2']
Cajal_Retzius_marker = ['RELN']
Endothelial_marker = ['FLT1']
Ependymal_marker = ['FOXJ1']
GABA_marker = ['GAD2']
GC_marker = ['PROX1','NPY1R','ERC2']
Granule_1_marker = ['CPLX2']
GC_immature_marker = ['BHLHE22']
OPC_marker = ['CSPG4',]
NFOL_marker = ['BCAS1']
MOL_marker = ['MOG']
PVM_marker = ['LYVE1']
VLMC_marker = ['DCN']
nIPC_perin_marker = ['ASCL1','RFC4','MKI67']
nIPC_marker = ['']

good_gene_list = RGL_marker +  aNSC_marker+ nIPC_perin_marker + NB_marker + GC_marker+Granule_1_marker                 + CA3_Pyr_marker + GABA_marker                   + OPC_marker + NFOL_marker   + MOL_marker + MicroGlia_marker + Cajal_Retzius_marker                 + Endothelial_marker + Ependymal_marker + PVM_marker + VLMC_marker 

plot_gene = good_gene_list
sc.pl.stacked_violin(ad, groupby='label',                     var_names= good_gene_list, jitter=False,layer = 'raw',
                     swap_axes = True,dendrogram = False,figsize=(10,6.5),\
                    linewidth = 0.15, palette=ad.uns['label_colors'])

#plot gene expression
matplotlib.rcParams.update({'font.size': 14})
fig = plt.figure(figsize=(15,5))

gs = GridSpec(2,5, figure=fig)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[0, 2])
ax4 = fig.add_subplot(gs[0, 3])
ax5 = fig.add_subplot(gs[0, 4])

ax6 = fig.add_subplot(gs[1, 0])
ax7 = fig.add_subplot(gs[1, 1])
ax8 = fig.add_subplot(gs[1, 2])
ax9 = fig.add_subplot(gs[1, 3])
ax10 = fig.add_subplot(gs[1, 4])


gl = ['GFAP','MKI67','SOX4','NNAT','PROX1', 'RGS4','GAD2','AQP4','CSPG4','C1QA']


ax1 = sc.pl.umap(ad, color=gl[0],ax=ax1,legend_loc = 'none',cmap = 'RdGy_r',           legend_fontsize=10,title=gl[0], show=False,  frameon = False, s = 1)

ax2 = sc.pl.umap(ad, color=gl[1],ax=ax2,legend_loc = 'none',cmap = 'RdGy_r',            legend_fontsize=10,title=gl[1], show=False, frameon = False,s = 1)

ax3 = sc.pl.umap(ad, color=gl[2],ax=ax3,legend_loc = 'none',cmap = 'RdGy_r',            legend_fontsize=10,title=gl[2], show=False,  frameon = False,s = 0.2)

ax4 = sc.pl.umap(ad, color=gl[3],ax=ax4,legend_loc = 'none',cmap = 'RdGy_r',           legend_fontsize=10,title=gl[3], show=False,  frameon = False)

ax5 = sc.pl.umap(ad, color=gl[4],ax=ax5,legend_loc = 'none',cmap = 'RdGy_r',            legend_fontsize=10,title=gl[4], show=False,  frameon = False)

ax6 = sc.pl.umap(ad, color=gl[5],ax=ax6,legend_loc = 'none',cmap = 'RdGy_r',           legend_fontsize=10,title=gl[5], show=False,  frameon = False)

ax7 = sc.pl.umap(ad, color=gl[6],ax=ax7,legend_loc = 'none',cmap = 'RdGy_r',           legend_fontsize=10,title=gl[6], show=False,  frameon = False)

ax8 = sc.pl.umap(ad, color=gl[7],ax=ax8,legend_loc = 'none',cmap = 'RdGy_r',           legend_fontsize=10,title=gl[7], show=False, frameon = False)

ax9 = sc.pl.umap(ad, color=gl[8],ax=ax9,legend_loc = 'none',cmap = 'RdGy_r',           legend_fontsize=10,title=gl[8], show=False,  frameon = False)

ax10 = sc.pl.umap(ad, color=gl[9],ax=ax10,legend_loc = 'none',cmap = 'RdGy_r',           legend_fontsize=10,title=gl[9], show=False,  frameon = False)


#plot heatmap
plot_matrix_norm=pd.read_csv('Label_plot_matrix_sub_norm.csv', index_col = 0)

figsize(5,10)
sns_plot = sns.heatmap(plot_matrix_norm_plot, cmap = 'RdGy_r',robust = False,xticklabels = 1,
yticklabels = 0)

plot_matrix_norm_plot.to_csv('data/DE/Hao_SourceData_Fig1_plot_matrix.csv')

