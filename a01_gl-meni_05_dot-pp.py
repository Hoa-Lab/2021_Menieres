import scanpy as sc
import pandas as pd
from pathlib import Path

#-------------------------------------------
fd_ada='./out/a00_pp_00_load'
fd_out='./out/a01_gl-meni_05_dot-pp'

l_sample=['sc_P30', 'sNuc_P30']
l_cell=['Marginal', 'Intermediate', 'Basal']

dic_cell={'sc_P30': ['Marginal', 'Intermediate', 'Basal', 'Spindle-Root', 'Fibrocyte', 'Macrophage'], 'sNuc_P30': ['Marginal', 'Intermediate', 'Basal', 'Spindle', 'Root', 'Macrophage']}

dic_gl={'Marginal': ['Kcne1', 'Atp1b2', 'Esrrb', 'Add2', 'Sgk1', 'Atp13a5', 'Hmx2', 'Cacna2d1', 'Car12', 'Eya4', 'Tnfrsf12a', 'Shroom3', 'Wnk2', 'Dtna'],
		'Intermediate': ['Met', 'Cdc14a', 'Ednrb', 'Tmem176a', 'Tmem176b', 'Gsta1'],
		'Basal': ['Wnk4', 'Col11a2', 'Slc44a2', 'Cldn11']}

#---------------------------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)


####################################################################
for sample in l_sample:
	#load
	l_dotcell=dic_cell[sample]
	ada=sc.read(f'{fd_ada}/{sample}.h5ad')
	df=ada.to_df()
	df_anno=ada.obs.loc[:, ['anno']]
	df=df.merge(df_anno, left_index=True, right_index=True)
	
	#get gene list
	l_gene=[]
	for cell in l_cell:
		l_gene=l_gene+dic_gl[cell]
	
	#truncate df
	dfi=df.loc[df['anno'].isin(l_dotcell), l_gene+['anno']].copy()
	
	#get mean
	dfi_mean=dfi.groupby('anno').mean().reindex(l_dotcell)
	
	#get non-zero perc df 
	dfi_perc=(dfi.loc[:, l_gene]>0)
	dfi_perc['anno']=dfi['anno']
	dfi_perc=dfi_perc.groupby('anno').mean().reindex(l_dotcell)
	
	#save
	dfi_mean.to_csv(f'{fd_out}/{sample}_mean.csv')
	dfi_perc.to_csv(f'{fd_out}/{sample}_perc.csv')	
		

