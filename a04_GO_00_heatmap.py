import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path

#----------------------------------------------------
n=20

f_in='./raw/misc/GO_Biological_Process_2018_table.txt'
f_gene='./out/a01_gl-meni_00_genelist/gl_meni.txt'
fd_out='./out/a04_GO_00_heatmap'

#-----------------------------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)
l_gene=Path(f_gene).read_text().strip().split('\n')


###############################################################
#clean go
df=pd.read_csv(f_in, sep='\t')
df=df.loc[:, ['Term', 'Genes', 'Combined Score', 'Adjusted P-value']]
df=df.sort_values(['Combined Score', 'Adjusted P-value'], ascending=False)
df=df.iloc[0:n, :].copy()

#prepare genes
go_gene=[]
for g in df['Genes'].tolist():
	go_gene.extend(g.split(';'))

go_gene=list(set(go_gene))
go_gene=[i.capitalize() for i in go_gene]
l_gene=[i for i in l_gene if i in go_gene]

#make df
df_hm=pd.DataFrame(0, index=df['Term'], columns=l_gene)
for _, row in df.iterrows():
	go=row['Term']
	l_temp=row['Genes'].strip().split(';')
	l_temp=[i.capitalize() for i in l_temp]
	for g in l_temp:
		df_hm[g][go]=1

#sort
df_hm=df_hm.T
df_hm=df_hm.sort_values(df['Term'].tolist(), ascending=False)
df_hm=df_hm.T

##plot
#fig, ax=plt.subplots(figsize=(16,7))
#ax=sns.heatmap(df_hm, cmap='Purples',cbar=False, linewidths=.5, vmin=-0.1)

##adjust
#ax.xaxis.label.set_visible(False)
#ax.yaxis.label.set_visible(False)
#ax.xaxis.tick_top()
#ax.yaxis.tick_right()
##plt.xticks(np.arange(0.5, len(l_gene)+0.5, 1), l_gene, fontsize=8, rotation=90)
#plt.xticks([])
#plt.yticks(fontsize=13, rotation=0)

##save
#plt.tight_layout()
#plt.savefig(f'{fd_out}/go_heatmap.png', dpi=300)
#plt.close()


#------------------------------------------------------------------
#split
df=df_hm.iloc[:, 0:110]
l_gene=l_gene[0:110]

#plot
fig, ax=plt.subplots(figsize=(18.5,7))
ax=sns.heatmap(df, cmap='Purples',cbar=False, linewidths=.5, vmin=-0.1)

#adjust
ax.xaxis.label.set_visible(False)
ax.yaxis.label.set_visible(False)
ax.xaxis.tick_top()
ax.yaxis.tick_right()
plt.xticks(np.arange(0.5, len(l_gene)+0.5, 1), l_gene, fontsize=7, rotation=90)
plt.yticks(fontsize=13, rotation=0)

#save
plt.tight_layout()
plt.savefig(f'{fd_out}/go_heatmap_1.png', dpi=300)
plt.close()


	
