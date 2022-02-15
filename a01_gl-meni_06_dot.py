import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

#----------------------------------------------------------
fd_in='./out/a01_gl-meni_05_dot-pp'
fd_out='./out/a01_gl-meni_06_dot'

l_sample=['sc_P30', 'sNuc_P30']
l_cell=['Marginal', 'Intermediate', 'Basal']
dic_title={'sc_P30': 'SC', 'sNuc_P30': 'SN'}

dic_cell={'sc_P30': ['Marginal', 'Intermediate', 'Basal', 'Spindle-Root', 'Fibrocyte', 'Macrophage'], 'sNuc_P30': ['Marginal', 'Intermediate', 'Basal', 'Spindle', 'Root', 'Macrophage']}

dic_gl={'Marginal': ['Kcne1', 'Atp1b2', 'Esrrb', 'Add2', 'Sgk1', 'Atp13a5', 'Hmx2', 'Cacna2d1', 'Car12', 'Eya4', 'Tnfrsf12a', 'Shroom3', 'Wnk2', 'Dtna'],
		'Intermediate': ['Met', 'Cdc14a', 'Ednrb', 'Tmem176a', 'Tmem176b', 'Gsta1'],
		'Basal': ['Wnk4', 'Col11a2', 'Slc44a2', 'Cldn11']}



#-----------------------------------------------------------
Path(fd_out).mkdir(exist_ok=True, parents=True)


#-------------------------------------------------------------
def make_df(df1, df2):
	''' columns: [cell, gene, avg, perc]'''
	#1. transform df1
	df1=df1.stack().reset_index()
	df1.columns=['cell', 'gene', 'Mean']
	df1['idx']=df1['cell']+'_'+df1['gene']
	df1=df1.set_index('idx')
	#2. transform df2
	df2=df2.stack().reset_index()
	df2.columns=['cell', 'gene', 'Percentage']
	df2['idx']=df2['cell']+'_'+df2['gene']
	df2=df2.set_index('idx')
	#3. merge
	df2=df2.loc[:, ['Percentage']]
	df=df1.merge(df2, left_index=True, right_index=True)
	return df
	

def plot_dot(df, f_out, cmap='Oranges', title=None, sz=(8,8), ytick=12, ylim=None, ylen=11, c_size=250):
	#1. plot
	fig, ax=plt.subplots(figsize=sz)
	ax = sns.scatterplot(x='cell', y='gene', hue='Mean', size='Percentage', data=df, sizes=(0, c_size), palette=cmap, vmin=-1)
	
	#2. adjust
	#ax.spines['top'].set_visible(False)
	#ax.spines['right'].set_visible(False)
	plt.title(title, fontsize=22, pad=15, weight='semibold')
	plt.xlabel('')
	plt.ylabel('')
	plt.yticks(fontsize=ytick, rotation=0, weight='semibold')
	plt.xticks(fontsize=18, rotation=90, weight='semibold')
	plt.legend(loc=(1.05, 0), frameon=False, ncol=2, prop={'size': ylen, 'weight': 'semibold'})
	plt.xlim([-0.5, 5.5])
	plt.ylim(ylim)
	
	#4. save
	plt.tight_layout()
	plt.savefig(f_out, dpi=300)
	plt.close()
	return

#########################################################################
for sample in l_sample:
	#load
	df_mean=pd.read_csv(f'{fd_in}/{sample}_mean.csv', index_col=0)
	df_perc=pd.read_csv(f'{fd_in}/{sample}_perc.csv', index_col=0)	
	
	#reverse columns order
	df_mean=df_mean.iloc[:, ::-1]
	df_perc=df_perc.iloc[:, ::-1]		
	
	#make df
	df=make_df(df_mean, df_perc)
	
	#plot
	f_out=f'{fd_out}/{sample}.png'
	title=f'P30 ({dic_title[sample]})'
	plot_dot(df, f_out, title=title)
	





