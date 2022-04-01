import pickle
from SCCAF import *
from NaiveDE import *
import scanpy as sc
import harmonypy as hm
from sklearn.preprocessing import scale

#load data
ad = sc.read('data/FM_hippo_207785.h5')

#plot umap
matplotlib.rcParams.update({'font.size': 10})
fig = plt.figure(constrained_layout=True,figsize=(5,5))
gs = GridSpec(1, 1, figure=fig)
ax2 = fig.add_subplot(gs[0, 0])
ax2 = sc.pl.umap(ad_sub, color=['show_Neuron'],ax=ax2,legend_loc = None,           legend_fontsize=10,title=' ', show=False, size=5, frameon = False)

#plot mKNN

ad_zc = sc.read('toPyr1_f.h5')

matplotlib.rcParams.update({'font.size': 14})
fig = plt.figure(constrained_layout=True,figsize=(8,8))

gs = GridSpec(1, 1, figure=fig)
ax2 = fig.add_subplot(gs[0, 0])

ax2 = sc.pl.draw_graph(ad_sub, color = ['label'],legend_loc = None, legend_fontoutline = 3,frameon = False, edges = True,
                      ax = ax2, show = False, title = '', edges_width = 0.2)


matplotlib.rcParams.update({'font.size': 8})
fig = plt.figure(figsize=(4.5,6))

gs = GridSpec(3,2, figure=fig)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[1, 0])
ax4 = fig.add_subplot(gs[1, 1])
ax5 = fig.add_subplot(gs[2, 0])
ax6 = fig.add_subplot(gs[2, 1])

cm = 'RdGy_r'

gene_list = ['MKI67','ASCL1','EOMES','NEUROD2','PROX1','ELAVL2']
ax1 = sc.pl.draw_graph(ad_sub, color = gene_list[0],title = gene_list[0], legend_loc = None, legend_fontoutline = 3,frameon = False, edges = True,
     ax = ax1, show = False, cmap = cm )

ax2 = sc.pl.draw_graph(ad_sub, color = gene_list[1],title = gene_list[1], legend_loc = None,  legend_fontoutline = 3,frameon = False, edges = True, ax = ax2, show = False, cmap = cm, s = 30)

ax3 = sc.pl.draw_graph(ad_sub, color = gene_list[2],title = gene_list[2], legend_loc = None, legend_fontoutline = 3,frameon = False, edges = True, ax = ax3, show = False, cmap = cm )
ax4 = sc.pl.draw_graph(ad_sub, color = gene_list[3],title = gene_list[3], legend_loc = None, legend_fontoutline = 3,frameon = False, edges = True, ax = ax4, show = False, cmap = cm)
ax5 = sc.pl.draw_graph(ad_sub, color = gene_list[4],title = gene_list[4], legend_loc = None, legend_fontoutline = 3,frameon = False, edges = True, ax = ax5, show = False, cmap = cm)
ax6 = sc.pl.draw_graph(ad_sub, color = gene_list[5],title = gene_list[5], legend_loc = None, legend_fontoutline = 3,frameon = False, edges = True, ax = ax6, show = False, cmap = cm, s = 30)

marker_gene_list = ['BCAN',
 'ANXA5',
 'PTN',
 'CCND1',
 'RGCC',
 'NNAT','TUBB2B',
 'MARCKSL1',
 'RPS13',
 'TMSB10',
 'NRN1',
 'CALB1',
 'LAMP5',
 'TNNT2',
 'C1QL3',
 'SNAP25',
 'SYT1',
 'STXBP1',
 'LMO4',
 'RGS4','NUAK1',
 'CPLX1',]

matplotlib.rcParams['figure.dpi'] = 600

sc.pl.heatmap(ad_sub, marker_gene_list, groupby='label', figsize=(3.3, 6), standard_scale = 'var',
              cmap='RdGy_r',var_group_rotation=90)

