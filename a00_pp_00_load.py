import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
from shutil import copyfile
from pathlib import Path

#------------------------------------------------
fd_in='./raw/count'
fd_out='./out/a00_pp_00_load'

#-------------------------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)

#####################################################################
#sNuc P30
src=f'{fd_in}/sNuc_P30-processed/sNuc_P30.h5ad'
dst=f'{fd_out}/sNuc_P30.h5ad'
copyfile(src, dst)

#-----------------------------------------------------------------
#sc P30
#load
df_cnt=pd.read_csv(f'{fd_in}/sc_P30-Ian/cnts_01_sc_p30_pp_Ian.csv', index_col=0)
df_tsne=pd.read_csv(f'{fd_in}/sc_P30-Ian/tsne_01_sc_p30_pp_Ian.csv', index_col=0)
df_idt=pd.read_csv(f'{fd_in}/sc_P30-Ian/idts_01_sc_p30_pp_Ian.csv', index_col=0)

#load to sc
df_cnt=df_cnt.T
ada=ad.AnnData(df_cnt)

#add anno
df_idt.columns=['anno']
ada.obs=ada.obs.merge(df_idt, left_index=True, right_index=True)
ada.obs['anno']=ada.obs['anno'].replace(['B-Cell', 'Spindle/Root'], ['B Cell', 'Spindle-Root'])

#add tsne
idx=ada.obs.index
df_tsne=df_tsne.reindex(idx)
ada.obsm['X_tsne']=df_tsne.values

ada.write(f'{fd_out}/sc_P30.h5ad')
