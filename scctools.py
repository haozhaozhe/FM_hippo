
#useful tools for single cell RNA-seq analysis#
#mostly scanpy based#


import warnings
warnings.filterwarnings("ignore")
import pickle
from SCCAF import *
import seaborn as sns
from NaiveDE import *
import scanpy as sc
import harmonypy as hm
from sklearn.preprocessing import scale
from numpy import unique
import os


import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from mpl_toolkits.mplot3d import Axes3D 

    
##########################
## cluster tools        ##
########################## 


def show_one_cluster(ad,key_level, key,basis):
    #high light selected cluster
    ad.obs['class_show'] = (ad.obs[key_level]== key)
    ad.obs['class_show'] = ad.obs['class_show'].astype(str)
    ad.uns['class_show_colors'] = ['#d3d3d3','#0000FF']
    sc.pl.scatter(ad, color = [key_level,'class_show'], title = [key_level,key], basis = basis)


###########################
## plot alignment heatmap##
###########################
def two_species_heatmap(ad, species_1 = 'FM', species_2 = 'HS',species_1_key = 'subclass', species_2_key = 'marker2',\
                        louvain = 0,figure_path = 'test_heatmap.png'):
#generate hodge fig 5d heatmap
#input: ad: pre_merged, harmony corrected two species data, with .obs['species'] mark the species
#input : species_key, which .observe to generate the heatmap
#input: louvain, the louvain resolution for cluster the merged data
#input: figure_path, the path to save the heatplot, if figure_path = 0, do not save figure

    #prepare data
    if not louvain ==0:
        sc.pp.neighbors(ad,  metric='euclidean',use_rep = 'X_harmony' )
        sc.tl.louvain(ad, resolution = louvain,key_added = 'louvain')

    sc.pl.umap(ad, color = ['species','louvain'])
    
    ad_1 = ad[ad.obs['species'].isin([species_1]),:] 
    sc.pl.umap(ad_1, color = [species_1_key,'louvain'], legend_loc = 'on data')
    
    ad_2 = ad[ad.obs['species'].isin([species_2]),:]
    sc.pl.umap(ad_2, color = [species_2_key,'louvain'], legend_loc = 'on data')
    
    df_1 = pd.crosstab(ad_1.obs[species_1_key], ad_1.obs['louvain'], normalize ='index')
    df_2 = pd.crosstab(ad_2.obs[species_2_key], ad_2.obs['louvain'], normalize ='index')

    df_1.columns = df_1.columns.tolist()
    df_2.columns = df_2.columns.tolist()

    #add missing louvain
    dif_list = list(set(list(df_2.columns)) ^set(list(df_1.columns)))
    for dif in dif_list:
        if dif not in list(df_1.columns):
            df_1[dif] = 0
        if dif not in list(df_2.columns):
            df_2[dif] = 0  
    df_2 = df_2[df_1.columns] #reorder column sequence, important!

    #generate heatmap matrix
    mk_cluster_all = df_1.index
    hu_cluster_all = df_2.index 
    low_sum_matrix = pd.DataFrame(index=mk_cluster_all, columns = hu_cluster_all)

    low_sum_matrix.index.name = species_1 +'_cluster'
    low_sum_matrix.columns.name = species_2 +'_cluster'
    
    
    low_sum_matrix = pd.DataFrame(index=mk_cluster_all, columns = hu_cluster_all)
    for mk_cluster in mk_cluster_all:
        for hu_cluster in hu_cluster_all:
        
            two_row = np.column_stack([df_1.loc[mk_cluster],df_2.loc[hu_cluster]])
            low_sum = sum(two_row.min(axis=1))
            low_sum_matrix.loc[mk_cluster, hu_cluster] = low_sum
            
    low_sum_matrix.to_csv('data/python/Homology_map/temp.csv')
    del(low_sum_matrix)
    low_sum_matrix =pd.read_csv('data/python/Homology_map/temp.csv')
    low_sum_matrix.set_index((species_1+'_cluster'), inplace = True)
    # for some reason, direct use matrix does not work.... save and read did the trick...
    
    low_sum_matrix_sort = low_sum_matrix.copy()

    ss = sns.clustermap(low_sum_matrix_sort, cmap = 'Greys',  cbar = True, col_cluster = False,  xticklabels=1, yticklabels=1,  method =  'single')
    
    hu_cluster_all = low_sum_matrix_sort.columns
    low_sum_matrix_sort1 = low_sum_matrix_sort
    low_sum_matrix_sort2 = low_sum_matrix_sort1.iloc[ss.dendrogram_row.reordered_ind]  
    
    low_sum_matrix1 =  low_sum_matrix_sort2
    for hu_cluster in hu_cluster_all:
        line = list(low_sum_matrix1.loc[:,hu_cluster])
        ss = np.array(line)
        tmp = []
        for kk in  range(len(ss)):
            if kk < len(ss) - 5:
                tmp.append( ss[kk]+ss[kk+1]/2+ss[kk+2]/3+ss[kk+3]/4+ss[kk+4]/5)
            elif kk < len(ss) - 4: 
                tmp.append( ss[kk]+ss[kk+1]/2+ss[kk+2]/3+ss[kk+3]/4)
            elif kk < len(ss) - 3: 
                tmp.append( ss[kk]+ss[kk+1]/2+ss[kk+2]/3 )
            elif kk < len(ss) - 2: 
                tmp.append( ss[kk]+ss[kk+1]/2 )
            elif kk < len(ss) - 1: 
                tmp.append( ss[kk] )
            elif kk < len(ss): 
                tmp.append(ss[kk]+ss[kk]/2+ss[kk]/3+ss[kk]/4+ss[kk]/5+ss[kk]/6)
        low_sum_matrix_sort2.loc['max_order',hu_cluster]= float(tmp.index(max(tmp)))
        del tmp
    low_sum_matrix_sort2 = low_sum_matrix_sort2.sort_values(by = ['max_order'],axis=1)
    plt.figure(figsize = (15,15))
    graph = sns.heatmap(low_sum_matrix_sort2[0:-1], cmap="Greys", cbar=True, xticklabels=1,yticklabels=1, linewidth = 0.01, linecolor = 'gray')
    
    plot_matrix = low_sum_matrix_sort2[0:-1]
    
    if not figure_path == 0:
        fig = graph.get_figure()
        fig.savefig(figure_path)
        
    return(plot_matrix)    


######################
## plot stacked bar###
######################

def sc_plot_stacke_bar(ad, \
                       y_key="predict_plot_direct",\
                       x_key="cell_type",\
                       cluster_palette=None,xlabel_rotation=0,\
                       savePath = 'Figure/Fig_NG_SCCAF/stackedbar_',
                       rotation = 0, fig_size = (8,5)
                      ):
    
    #plot stacked bar plot for SCCAF projection results
    #implement from  gist.github.com/wflynny 
    #adata : AnnData object
    #cluster_key : original cell type key
    #predict_key : SCCAF projected key
    # cluster_palette: list of color, same sequence as category
    
    sizes = ad.obs.groupby([y_key, x_key]).size()
    props = sizes.groupby(level=1).apply(lambda x: 100 * x / x.sum()).reset_index() 
    props = props.pivot(columns=x_key, index=y_key).T
    props.index = props.index.droplevel(0)
    props.fillna(0, inplace=True)
    
    #plot figure
    
    fig, ax = plt.subplots(dpi=600, figsize = fig_size)
    fig.patch.set_facecolor("white")
    
    cmap = None
    
    if cluster_palette is not None:
        cmap = sns.palettes.blend_palette(
            cluster_palette, 
            n_colors=len(cluster_palette), 
            as_cmap=True)
   
    cluster_props = props.copy()
    cluster_props.plot(
        kind="bar", 
        stacked=True, 
        ax=ax, 
        legend=None, 
        colormap=cmap
    )
    
    ax.legend(bbox_to_anchor=(1.01, 1), frameon=False, title=" ")
    sns.despine(fig, ax)
    ax.tick_params(axis="x", rotation=rotation)
    ax.set_xlabel(cluster_props.index.name.capitalize())
    ax.set_ylabel("Proportion")
    fig.tight_layout() 
    
    plt.savefig(savePath+y_key+'.pdf', dpi = 600)
    plt.savefig(savePath+y_key+'.png', dpi = 600)
   
    
    return fig


def plot_enrich_square(data, n_terms=20, cmap = 'Reds_r', vmax_scale = 2, vmin_scale = 2, save=False,):
# Plotting GO enrichment terms, only keep p value, remove intersection
# data,  Pandas Dataframe output by gprofiler
# n_terms, number of terms to plot
# save, the path to save the file, if empty, not save anything
# vmax_scale, vmin_scale, factors to adjust display color

    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sb
    from matplotlib import colors
    from matplotlib import rcParams

    
    def scale_data_5_75(data):
        mind = np.min(data)
        maxd = np.max(data)
    
        if maxd == mind:
            maxd=maxd+1
            mind=mind-1
        
        drange = maxd - mind
        return ((((data - mind)/drange*0.70)+0.05)*100)

    

    # Test data input
    if not isinstance(data, pd.DataFrame):
        raise ValueError('Please input a Pandas Dataframe output by gprofiler.')
        
    if not np.all([term in data.columns for term in ['p_value', 'name', 'intersection_size']]):
        raise TypeError('The data frame {} does not contain enrichment results from gprofiler.'.format(data))
    
    data_to_plot = data.iloc[:n_terms,:].copy()
    data_to_plot['go.id'] = data_to_plot.index
    data_to_plot = data_to_plot.sort_values(by=['p_value'])
    
    min_pval = data_to_plot['p_value'].min()
    max_pval = data_to_plot['p_value'].max()
    
    # Scale intersection_size to be between 5 and 75 for plotting
    #Note: this is done as calibration was done for values between 5 and 75
    #data_to_plot['scaled.overlap'] = scale_data_5_75(data_to_plot['intersection_size'])
    
    norm = colors.LogNorm(min_pval, max_pval)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    rcParams.update({'font.size': 14, 'font.weight': 'bold'})

    #sb.set(style="whitegrid")

    path = plt.scatter(x=[1] *n_terms, y="name", c='p_value', cmap=cmap, 
                       norm=colors.LogNorm(min_pval, max_pval),
                       vmax = max_pval*vmax_scale, vmin = min_pval/vmin_scale,
                       data=data_to_plot, linewidth=1, edgecolor="k",marker="s", s = 300)
    ax = plt.gca()
    ax.invert_yaxis()

    if save: ax.set_title(save.split('/')[-1])
    ax.set_ylabel('')
    #ax.set_xlabel('Gene ratio', fontsize=14, fontweight='bold')
    ax.xaxis.grid(False)
    ax.yaxis.grid(False)
    
    ax.spines['bottom'].set_color('#FFFFFF')
    ax.spines['left'].set_color('#FFFFFF')
    ax.spines['top'].set_color('#FFFFFF')
    ax.spines['right'].set_color('#FFFFFF') 

    ax.set_xticks([])
    ax.tick_params(axis='y', length=0)

    
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Get tick marks for this plot
    #Note: 6 ticks maximum
    min_tick = np.floor(np.log10(min_pval)).astype(int)
    max_tick = np.ceil(np.log10(max_pval)).astype(int)
    tick_step = np.ceil((max_tick - min_tick)/6).astype(int)
    
    # Ensure no 0 values
    if tick_step == 0:
        tick_step = 1
        min_tick = max_tick-1
    
    ticks_vals = [10**i for i in range(max_tick, min_tick-1, -tick_step)]
    ticks_labs = ['$10^{'+str(i)+'}$' for i in range(max_tick, min_tick-1, -tick_step)]

    #Colorbar
    fig = plt.gcf()
    cbaxes = fig.add_axes([1.0, 0.15, 0.1, 0.2])
    cbar = ax.figure.colorbar(sm, ticks=ticks_vals, shrink=0.5, anchor=(0,0.1), cax=cbaxes)
    cbar.ax.set_yticklabels(ticks_labs,fontsize=12,fontweight='normal')
    cbar.set_label("Adjusted p-value", fontsize=12,fontweight='normal')
    cbar.ax.invert_yaxis() 
    if save:
        plt.savefig(save+'.pdf', dpi=600, format='pdf')
        plt.savefig(save+'.png', dpi=600, format='png')
    plt.show()
    
    

