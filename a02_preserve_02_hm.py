import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import json
from scipy.stats import zscore
import numpy as np

#----------------------------------------------------------------
f_gl='./out/a02_preserve_01_hm-pp/gene.json'
f_ada='./raw/count/h5ad/concat_merged.h5ad'
f_meta='./raw/meta_info.csv'
fd_out='./out/a02_preserve_02_hm'

#--------------------------------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)

#color map
df_meta=pd.read_csv(f_meta, index_col=0)
dic_cmap=df_meta.to_dict()['color']

def plot_topbar(df, f_out, l_cell=None, sz=(6,6), dic_cmap=dic_cmap):
	cmap=[dic_cmap[i] for i in l_cell]
	#1. make df
	df=df.loc[df['anno'].isin(l_cell), :]
	#2. sort
	df['anno']=pd.Categorical(df['anno'], categories=l_cell, ordered=True)
	df=df.sort_values('anno', ascending=True)
	df['anno']=df.anno.cat.codes
	#pp
	df=df.loc[:, ['anno']].T
	#3. plot
	fig, ax=plt.subplots(figsize=sz)
	ax=sns.heatmap(df, cmap=cmap, cbar=False)	
	#3. adjust
	ax.xaxis.label.set_visible(False)
	ax.yaxis.label.set_visible(False)
	plt.xticks([])
	plt.yticks([])
	#4. save
	plt.tight_layout()
	plt.savefig(f_out, dpi=300)
	plt.close()	
	return


def plot_hm(df, f_out, l_cell=None, size=(10,15), vmax=1, vmin=-0.3, y=11.5):
	#2. heatmap
	fig, ax=plt.subplots(figsize=size)
	ax=sns.heatmap(df, cmap='Purples', vmax=vmax, vmin=vmin, cbar=False)
	#3. adjust
	ax.xaxis.label.set_visible(False)
	ax.yaxis.label.set_visible(False)
	plt.xticks([])
	plt.yticks(np.arange(0.5, df.shape[0]+0.5, 1), df.index.tolist(), fontsize=y, rotation=0, weight='medium')
	#4. save
	plt.tight_layout()
	plt.savefig(f_out, dpi=300)
	plt.close()
	return

################################################################################
#load df
ada=sc.read(f_ada)
df=pd.DataFrame(ada.raw.X.toarray(), index=ada.obs.index, columns=ada.raw.var.index)
df['anno']=ada.obs['anno']
df['anno']=df['anno'].replace(['Spindle-Root-1', 'Spindle-Root-2'], ['Root', 'Spindle'])

l_rss_main=['Marginal', 'Intermediate', 'Basal', 'Spindle', 'Root', 'Fibrocyte', 'Reissner', 'Macrophage']

#load dic_cell
with open(f_gl, 'r') as f:
	dic_cell=json.load(f)

#get gene list
l_gl_main=[]       #only required cell genes
for cell in l_rss_main:
	l_gl_main.extend(dic_cell[cell])

#plot topbar
f_out=f'{fd_out}/top.png'
plot_topbar(df, f_out, l_cell=l_rss_main)

##prepare df
#df=df.loc[df['anno'].isin(l_rss_main), :].copy()
#df['anno']=pd.Categorical(df['anno'], categories=l_rss_main, ordered=True)
#df=df.sort_values('anno', ascending=True).drop('anno', axis=1)

#df_main=df.loc[:, l_gl_main].copy()
#df_main=df_main.loc[:, (df_main!=0).any(axis=0)]
#df_main=df_main.T

#print(df_main.shape)  #155

##heatmap filter (high rss genes)
#df1=df_main.iloc[0:77, :].copy()
#df2=df_main.iloc[77:, :].copy()

#f_out=f'{fd_out}/main_1.png'
#plot_hm(df1, f_out, l_cell=l_rss_main, y=9.5)

#f_out=f'{fd_out}/main_2.png'
#plot_hm(df2, f_out, l_cell=l_rss_main, y=9.5)



























