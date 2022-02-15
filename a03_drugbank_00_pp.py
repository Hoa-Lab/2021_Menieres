import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

#----------------------------------------------
f_drug='./raw/drugbank.csv'
f_gl='./out/a01_gl-meni_00_genelist/gl_meni.txt'
fd_out='./out/a03_drugbank_00_pp'

#----------------------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)


#------------------function---------------------------
def plot_bar(df, f_out, color='Grey', title=None, n=10, sz=(6,8), y=16):
	dfi=df.iloc[0:n].copy()
	dfi.index.name='drug'
	dfi=dfi.reset_index()
	#1. plot
	sns.set()
	fig, ax=plt.subplots(figsize=sz)
	ax=sns.barplot(x='gene', y='drug', data=dfi, color=color, alpha=0.7)
	#2.adjust
	plt.title('Target Counts', fontsize=20, pad=10, weight='semibold')
	plt.xlabel('Targets', fontsize=18, labelpad=10, weight='semibold')
	plt.ylabel('', fontsize=18, labelpad=10)
	plt.yticks(fontsize=y, rotation=0, weight='semibold')
	#3. save
	plt.tight_layout()
	plt.savefig(f_out, dpi=300)
	plt.close()
	return
	
	
#############################################################
#load
df_drug=pd.read_csv(f_drug)
l_gene=Path(f_gl).read_text().strip().split('\n')
l_gene=[i.upper() for i in l_gene]

#clean
df=df_drug.loc[df_drug['gene'].isin(l_gene), :]

#get df
df=df.groupby('name').count().sort_values('gene', ascending=False)
df=df.loc[:, ['gene']]

df.to_csv(f'{fd_out}/drug.csv')

#plot
f_out=f'{fd_out}/bar.png'

plot_bar(df, f_out)
