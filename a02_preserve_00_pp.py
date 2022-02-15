import pandas as pd
import scanpy as sc
from pathlib import Path
#from pyscenic.rss import regulon_specificity_scores

#--------------------------------------------------------
f_gl='./out/a01_gl-meni_00_genelist/gl_meni.txt'
f_ada='./raw/count/h5ad/concat_merged.h5ad'
fd_out='./out/a02_preserve_01_hm-pp'

#------------------------------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)
l_gene=Path(f_gl).read_text().split('\n')
ada=sc.read(f_ada)

####################################################################
#pp
l_gl=[i for i in l_gene if i in ada.raw.var.index]
l_anno=ada.obs['anno']

#make df
df=pd.DataFrame(ada.raw.X.toarray(), index=ada.obs.index, columns=ada.raw.var.index)
df=df.loc[:, l_gl]

print(df)

##rss
#df_rss=regulon_specificity_scores(df, l_anno).T
#df_rss.to_csv(f'{fd_out}/rss.csv')
