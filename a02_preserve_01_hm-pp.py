import pandas as pd
import scanpy as sc
from pathlib import Path
from scipy.stats import zscore
import json

#---------------------------------------------------------
z=2.5
min_avg=0.1

f_rss='./out/a02_preserve_00_pp/rss.csv'
f_ada='./raw/count/h5ad/concat_merged.h5ad'
fd_out='./out/a02_preserve_01_hm-pp'


#--------------------------------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)
ada=sc.read(f_ada)
df_rss=pd.read_csv(f_rss, index_col=0)

#--------------------------------------------------------
def zscore_rss(df):
	df=df.dropna().T
	df=df.apply(zscore).T
	return df
	
def get_dic_cell(df, z, main=[], sup=[]):
	'''main & sup are used for customized gene list'''
	#pp
	dic_cell={}
	l_cell=df.columns.tolist()
	#get non-specific genes
	l_nonspec=df.loc[df.max(axis=1) < z, :].index.tolist()
	l_nonspec=list(set(l_nonspec+sup))
	l_nonspec=[i for i in l_nonspec if i not in main]
	dic_cell['not_spec']=l_nonspec
	#get specific gene
	all_gene=df.index.tolist()
	l_spec=[i for i in all_gene if i not in l_nonspec]
	#get spec genes
	df=df.reindex(l_spec)
	df['cell']=df.idxmax(axis=1)
	for cell in l_cell:
		l_gene=df.loc[df['cell']==cell, :].index.tolist()
		dic_cell[cell]=l_gene
	return dic_cell


def sort_gene(dic_cell, ada, avg, main=[], sup=[]):
	'''sort gene and move low exp gene to "not_spec"
	main & sup are used for customized gene list'''
	#convert to df
	df=pd.DataFrame(ada.raw.X.toarray(), index=ada.obs.index, columns=ada.raw.var.index)
	df['anno']=ada.obs['anno']
	l_cell=df['anno'].unique().tolist()
	l_cell=['Reissner', 'Marginal', 'Intermediate', 'Basal', 'Fibrocyte', 'Unknown', 'Macrophage', 'Spindle-Root-2', 'Spindle-Root-1']
	#loop on cell
	d_cell={}
	d_cell['not_spec']=dic_cell['not_spec']
	for cell in l_cell:
		#pp
		l_gene=dic_cell[cell]
		dfi=df.loc[df['anno']==cell, l_gene].copy().T
		#sort by avg exp
		dfi['avg']=dfi.mean(axis=1)
		dfi=dfi.sort_values('avg', ascending=False)
		l_gene=dfi.index.tolist() #keep order
		#move low avg exp gene to not_spec
		l_lowexp=dfi.loc[dfi['avg']<avg, :].index.tolist()
		l_lowexp=[i for i in l_lowexp if i not in main]
		d_cell['not_spec']=d_cell['not_spec']+l_lowexp
		#update spec genes
		d_cell[cell]=[i for i in l_gene if i not in l_lowexp]
	return d_cell

#######################################################################
#clean rss
df_rss=df_rss.iloc[:, 0:-1]
df_rss=df_rss.fillna(0)
df_rss=zscore_rss(df_rss)

#get dic_cell {cell: l_gene}
dic_cell=get_dic_cell(df_rss, z)
dic_cell=sort_gene(dic_cell, ada, min_avg)

#rename
dic_cell['Spindle']=dic_cell['Spindle-Root-2']
dic_cell['Root']=dic_cell['Spindle-Root-1']
del dic_cell['Spindle-Root-1']
del dic_cell['Spindle-Root-2']


#calculate gene number
n=0
l_cell=['Marginal', 'Intermediate', 'Basal', 'Spindle', 'Root', 'Fibrocyte', 'Reissner', 'Macrophage']
for cell in l_cell:
	n+=len(dic_cell[cell])
print(n)  #155

#save
f_out=f'{fd_out}/gene.json'
with open(f_out, 'w') as f:
	json.dump(dic_cell, f)
















