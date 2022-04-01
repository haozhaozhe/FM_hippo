import matplotlib.pyplot as plt
import scanpy as sc
import scvi


def ksample(x, n=1000):
    if len(x) <= n:
        return x
    else:
        return x.sample(n=n)


#load data
ad = sc.read('FM_hippo_207785.h5')


df = ad.obs.groupby('label').apply(lambda x: ksample(x, n=500))
ds = df.index.droplevel(level=0)
ad = ad[ds,:]

ad.X = ad.layers['counts']

#load mice data
ad1 = sc.read('admm_hoch.h5')

df = ad1.obs.groupby('characteristics: cell cluster').apply(lambda x: ksample(x,n=500))
ds = df.index.droplevel(level=0)
ad1 = ad1[ds,:]
ad1


#merge two dataset
df = pd.read_csv("data/CrossSpecies/mouse2monkey.csv", sep = '\t')
Ln_col_name = 'Crab-eating macaque gene name'
df = df[~df[Ln_col_name].isna()]
genes = df[['Gene name',Ln_col_name]]
genes = genes.drop_duplicates('Gene name')
genes = genes.drop_duplicates(Ln_col_name)
genes.index = genes[Ln_col_name]

genes = genes[genes[Ln_col_name].isin(ad.var_names)]
genes = genes[genes['Gene name'].isin(ad1.var_names)]

ad = ad[:,ad.var_names.isin(genes[Ln_col_name])]
ad.var['mouse'] = genes.reindex(ad.var_names)['Gene name']

ad1 = ad1[:, ad.var['mouse'].tolist()]

ad1.var['mouse'] = ad1.var_names
ad1.var_names = ad.var_names

ad.obs['org'] = 'monkey'
ad1.obs['org'] = 'mouse'

ad2 = ad.concatenate(ad1,index_unique = None)

#run integration
ad2.layers['counts']=ad2.X

sc.pp.normalize_per_cell(ad2, counts_per_cell_after=1e4)
sc.pp.log1p(ad2)
sc.pp.highly_variable_genes(
    ad2,
    batch_key="org",
    flavor="seurat",
    n_top_genes=4000,
    subset=True
)

x = sc.pp.calculate_qc_metrics(ad2, percent_top=None, inplace=False)
ad2.obs['n_counts'] = x[0].iloc[:,1]
ad2.obs['n_genes'] = x[0].iloc[:,0]

scvi.data.setup_anndata(
    ad2,
    batch_key = 'batch1',
    categorical_covariate_keys=['org'],
    continuous_covariate_keys = ['n_genes','n_counts']
)

model = scvi.model.SCVI(ad2)

model.train(600)

ad2.obsm['X_scvi'] = model.get_latent_representation()

sc.pp.neighbors(ad, use_rep="X_scvi")
sc.tl.umap(ad, min_dist=0.3)

import harmonypy as hm
ho = hm.run_harmony(ad2.obsm['X_scvi'], ad2.obs, ['batch','org'])

ad2.obsm['X_harmonypca'] = ho.Z_corr.T
sc.pp.neighbors(ad2, use_rep='X_harmonypca')

ad2.obsm['X_umapraw'] = ad2.obsm['X_umap']
sc.tl.umap(ad2)
ad2.obsm['X_umapharmony'] = ad2.obsm['X_umap']

ad2.write('admk_Hoch1_500_finish.h5')

