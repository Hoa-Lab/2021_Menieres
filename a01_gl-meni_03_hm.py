import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import json
from scipy.stats import zscore
import numpy as np

#----------------------------------------------------------------
fd_gl='./out/a01_gl-meni_02_hm-pp'
fd_ada='./out/a00_pp_00_load'
f_meta='./raw/meta_info.csv'
fd_out='./out/a01_gl-meni_03_hm'

#---------------------------------------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)

#color map
df_meta=pd.read_csv(f_meta, index_col=0)
dic_cmap=df_meta.to_dict()['color']


def plot_topbar(ada, f_out, l_cell=None, sz=(6,6), dic_cmap=dic_cmap):
	cmap=[dic_cmap[i] for i in l_cell]
	#1. make df
	df=ada.obs.loc[:, ['anno']].copy()
	df=df.loc[df['anno'].isin(l_cell), :]
	#2. sort
	df['anno']=pd.Categorical(df['anno'], categories=l_cell, ordered=True)
	df=df.sort_values('anno', ascending=True)
	df['anno']=df.anno.cat.codes
	#3. plot
	fig, ax=plt.subplots(figsize=sz)
	ax=sns.heatmap(df.T, cmap=cmap, cbar=False)	
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


###########################################################################
sample='sc_P30'
ada=sc.read(f'{fd_ada}/{sample}.h5ad')

l_rss_main=['Marginal', 'Intermediate', 'Basal', 'Spindle-Root', 'Macrophage', 'Fibrocyte']
l_rss_other=['B Cell', 'Neutrophil', 'Erythrocyte', 'Erythroblast', 'Unknown', 'not_spec']

#load dic_cell
with open(f'{fd_gl}/gl_{sample}.json', 'r') as f:
	dic_cell=json.load(f)

#2. get gene list
l_gl_main=[]       #only required cell genes
l_gl_other=[]          #other genes
for cell in l_rss_main:
	l_gl_main.extend(dic_cell[cell])
for cell in l_rss_other:
	l_gl_other.extend(dic_cell[cell])

#plot topbar
f_out=f'{fd_out}/{sample}_top.png'
plot_topbar(ada, f_out, l_cell=l_rss_main)

#prepare df
df=ada.to_df()
df['anno']=ada.obs['anno']
df=df.loc[df['anno'].isin(l_rss_main), :].copy()
df['anno']=pd.Categorical(df['anno'], categories=l_rss_main, ordered=True)
df=df.sort_values('anno', ascending=True).drop('anno', axis=1)

df_main=df.loc[:, l_gl_main].copy()
df_main=df_main.loc[:, (df_main!=0).any(axis=0)]
df_main=df_main.T

df_other=df.loc[:, l_gl_other].copy()
df_other=df_other.loc[:, (df_other!=0).any(axis=0)]
df_other=df_other.T

#print(df_main.shape)  #218
#print(df_other.shape)  #286

#heatmap filter (high rss genes)
df1=df_main.iloc[0:109, :].copy()
df2=df_main.iloc[109:, :].copy()

f_out=f'{fd_out}/{sample}_main_1.png'
plot_hm(df1, f_out, l_cell=l_rss_main, y=9.5)

f_out=f'{fd_out}/{sample}_main_2.png'
plot_hm(df2, f_out, l_cell=l_rss_main, y=9.5)

#heatmap extra
df1=df_other.iloc[0:95, :].copy()
df2=df_other.iloc[95:190, :].copy()
df3=df_other.iloc[190:, :].copy()

f_out=f'{fd_out}/{sample}_other1.png'
plot_hm(df1, f_out, l_cell=l_rss_main, y=10)

f_out=f'{fd_out}/{sample}_other2.png'
plot_hm(df2, f_out, l_cell=l_rss_main, y=10)

f_out=f'{fd_out}/{sample}_other3.png'
plot_hm(df3, f_out, l_cell=l_rss_main, y=10)



##-----------------------------------------------------------------------
#sample='sNuc_P30'
#ada=sc.read(f'{fd_ada}/{sample}.h5ad')

#l_rss_main=['Marginal', 'Intermediate', 'Basal', 'Spindle', 'Root', 'Reissner', 'Macrophage']
#l_rss_other=['B Cell', 'Neutrophil', 'Erythrocyte', 'Erythroblast', 'Unknown', 'not_spec']

##load dic_cell
#with open(f'{fd_gl}/gl_{sample}.json', 'r') as f:
#	dic_cell=json.load(f)

##2. get gene list
#l_gl_main=[]       #only required cell genes
#l_gl_other=[]          #other genes
#for cell in l_rss_main:
#	l_gl_main.extend(dic_cell[cell])
#for cell in l_rss_other:
#	l_gl_other.extend(dic_cell[cell])

##plot topbar
#f_out=f'{fd_out}/{sample}_top.png'
#plot_topbar(ada, f_out, l_cell=l_rss_main)

##prepare df
#df=ada.to_df()
#df['anno']=ada.obs['anno']
#df=df.loc[df['anno'].isin(l_rss_main), :].copy()
#df['anno']=pd.Categorical(df['anno'], categories=l_rss_main, ordered=True)
#df=df.sort_values('anno', ascending=True).drop('anno', axis=1)

#df_main=df.loc[:, l_gl_main].copy()
#df_main=df_main.loc[:, (df_main!=0).any(axis=0)]
#df_main=df_main.T

#df_other=df.loc[:, l_gl_other].copy()
#df_other=df_other.loc[:, (df_other!=0).any(axis=0)]
#df_other=df_other.T

##print(df_main.shape)  #172
##print(df_other.shape)  #456

##heatmap filter (high rss genes)
#df1=df_main.iloc[0:86, :].copy()
#df2=df_main.iloc[86:, :].copy()

#f_out=f'{fd_out}/{sample}_main_1.png'
#plot_hm(df1, f_out, l_cell=l_rss_main, y=10)

#f_out=f'{fd_out}/{sample}_main_2.png'
#plot_hm(df2, f_out, l_cell=l_rss_main, y=10)


##heatmap extra
#df1=df_other.iloc[0:91, :].copy()
#df2=df_other.iloc[91:192, :].copy()
#df3=df_other.iloc[192:283, :].copy()
#df4=df_other.iloc[283:374, :].copy()
#df5=df_other.iloc[374:, :].copy()

#f_out=f'{fd_out}/{sample}_other1.png'
#plot_hm(df1, f_out, l_cell=l_rss_main, y=10)

#f_out=f'{fd_out}/{sample}_other2.png'
#plot_hm(df2, f_out, l_cell=l_rss_main, y=10)

#f_out=f'{fd_out}/{sample}_other3.png'
#plot_hm(df3, f_out, l_cell=l_rss_main, y=10)

#f_out=f'{fd_out}/{sample}_other4.png'
#plot_hm(df4, f_out, l_cell=l_rss_main, y=10)

#f_out=f'{fd_out}/{sample}_other5.png'
#plot_hm(df5, f_out, l_cell=l_rss_main, y=10)












