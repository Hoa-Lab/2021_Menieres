import pandas as pd
import scanpy as sc
from pathlib import Path
from pyscenic.rss import regulon_specificity_scores

#--------------------------------------------------------
f_gl='./out/a01_gl-meni_00_genelist/gl_meni.txt'
fd_ada='./out/a00_pp_00_load'
fd_out='./out/a01_gl-meni_01_rss'

l_sample=['sc_P30', 'sNuc_P30']

#------------------------------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)
l_gene=Path(f_gl).read_text().split('\n')


##########################################################################
#sc P30
sample='sc_P30'

ada=sc.read(f'{fd_ada}/{sample}.h5ad')
l_gl=[i for i in l_gene if i in ada.var.index]
l_anno=ada.obs['anno']

df=pd.DataFrame(ada.X, columns=ada.var.index, index=ada.obs.index)
df=df.loc[:, l_gl]

df_rss=regulon_specificity_scores(df, l_anno).T
df_rss.to_csv(f'{fd_out}/{sample}.csv')


#-------------------------------------------------------
#sNuc P30
sample='sNuc_P30'

ada=sc.read(f'{fd_ada}/{sample}.h5ad')
l_gl=[i for i in l_gene if i in ada.var.index]
l_anno=ada.obs['anno']

df=pd.DataFrame(ada.X.toarray(), columns=ada.var.index, index=ada.obs.index)
df=df.loc[:, l_gl]

df_rss=regulon_specificity_scores(df, l_anno).T
df_rss.to_csv(f'{fd_out}/{sample}.csv')
