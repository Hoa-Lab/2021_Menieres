#get gene list from lab ref

import pandas as pd
from pathlib import Path

#-----------------------varaible------------------------------
f_gl='./raw/genelist/meniere_2020-09-17.txt'
fd_out='./out/a01_gl-meni_00_genelist'

#---------------------setup----------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)


#####################################################################
#get gene
l_gene=Path(f_gl).read_text().strip().split('\n')
l_gene=list(set(l_gene))
l_gene=[i.capitalize() for i in l_gene]


Path(f'{fd_out}/gl_meni.txt').write_text('\n'.join(l_gene))


