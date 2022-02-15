#scanpy dotplot, for the legend plot

import scanpy as sc
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

#------------------------------------------------
fd_ada='./out/a00_pp_00_load'
fd_out='./out/a01_gl-meni_04_dotplot'

dic_cell={'sc_P30': ['Marginal', 'Intermediate', 'Basal', 'Spindle-Root', 'Fibrocyte', 'Macrophage'], 'sNuc_P30': ['Marginal', 'Intermediate', 'Basal', 'Spindle', 'Root', 'Macrophage']}

dic_gl={'Marginal': ['Kcne1', 'Atp1b2', 'Esrrb', 'Add2', 'Sgk1', 'Atp13a5', 'Hmx2', 'Cacna2d1', 'Car12', 'Eya4', 'Tnfrsf12a', 'Shroom3', 'Wnk2', 'Dtna'],
		'Intermediate': ['Met', 'Cdc14a', 'Ednrb', 'Tmem176a', 'Tmem176b', 'Gsta1'],
		'Basal': ['Wnk4', 'Slc44a2', 'Slc9a2','Coch']}

#-------------------------------------------------
Path(fd_out).mkdir(parents=True, exist_ok=True)


################################################################
sample='sc_P30'
sz_title='Fraction of cells\nin cell types (%)'
cbar_title='Mean expression\nin cell types'

l_cell=dic_cell[sample]
ada=sc.read(f'{fd_ada}/{sample}.h5ad')

#clean ada
ada=ada[ada.obs['anno'].isin(l_cell), :].copy()
ada.obs['anno']=pd.Categorical(ada.obs['anno'], categories=l_cell, ordered=True)

#plot
f_out=f'{fd_out}/{sample}.png'
dic_ax=sc.pl.dotplot(ada, dic_gl, 'anno', show=False, var_group_rotation=0, swap_axes=False, cmap='Oranges', colorbar_title=cbar_title, size_title=sz_title, dot_max=1)

#adjust
fig=plt.gcf()
fig.set_figwidth(12)
fig.set_figheight(7)

ax0=dic_ax['mainplot_ax']
ax0.yaxis.set_tick_params(labelsize=14)

ax1=dic_ax['gene_group_ax']
ax1.yaxis.set_tick_params(labelsize=14)


#save
plt.savefig(f_out, dpi=300)


#---------------------------------------------------
sample='sNuc_P30'
sz_title='Fraction of cells\nin cell types (%)'
cbar_title='Mean expression\nin cell types'

l_cell=dic_cell[sample]
ada=sc.read(f'{fd_ada}/{sample}.h5ad')

#clean ada
ada=ada[ada.obs['anno'].isin(l_cell), :].copy()
ada.obs['anno']=pd.Categorical(ada.obs['anno'], categories=l_cell, ordered=True)

#plot
f_out=f'{fd_out}/{sample}.png'
dic_ax=sc.pl.dotplot(ada, dic_gl, 'anno', show=False, var_group_rotation=0, swap_axes=False, cmap='Oranges', colorbar_title=cbar_title, size_title=sz_title)

#adjust
fig=plt.gcf()
fig.set_figwidth(12)
fig.set_figheight(7)

ax0=dic_ax['mainplot_ax']
ax0.yaxis.set_tick_params(labelsize=14)

ax1=dic_ax['gene_group_ax']
ax1.yaxis.set_tick_params(labelsize=14)


#save
plt.savefig(f_out, dpi=300)

