import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.patches as mpatches


# read the file used to construct the consensusmap
#consensusInfo=pd.read_csv("C:/Users/tim.kuijpers/Desktop/CCLE_result_clusterk8_consensusmatrix.csv",sep=",",index_col=0)
#consensusInfo=pd.read_csv("E:/CCLE manusscript/Data/Cluster results/Feature genes and CpGs 40 runs/ConsensusMatrix_CCLE_K8.csv",sep=",",index_col=0)
consensusInfo=pd.read_csv("E:/CCLE manusscript/Data/Cluster results/featuers for genes and cpgs 40 runs/ConsensusMatrix_k8_CLLE.csv",sep=",",index_col=0)
# read the file for the cell line histology
cellAnnotation=pd.read_csv("E:/CCLE manusscript/Data/RAW data/Cell_lines_annotations_20181226.txt",sep="\t")

print(consensusInfo.shape)
print(cellAnnotation.shape)

# Get a colormap for the colors
CCLEData=cellAnnotation.loc[cellAnnotation['CCLE_ID'].isin(consensusInfo.index),:]
histologySamples=CCLEData.loc[:,['CCLE_ID','Histology']]

# make sure samples in consensusmap and data frame have same order
histologySamples.index=histologySamples['CCLE_ID']
histologySamples=histologySamples.reindex(consensusInfo.index)

colorsCustom=['red','darkorange','lime','darkgreen','blue','cyan','deepskyblue','dodgerblue','slategray','darkviolet','deeppink','grey','olive','yellow','fuchsia','rosybrown','chocolate','darkseagreen','white','teal','coral','black']
colormapCustom=dict(zip(histologySamples['Histology'].unique(), sns.color_palette(colorsCustom,len(set(histologySamples['Histology'].unique())),desat=1)))
colorSamples =histologySamples['Histology'].map(colormapCustom)# create the basic consensusmap
g=sns.clustermap(consensusInfo, cmap='rocket',row_colors=colorSamples,xticklabels=False, yticklabels=False,linewidths=0)
g.ax_heatmap.set_xlabel("")
g.ax_heatmap.set_ylabel("")
g.ax_heatmap.tick_params(axis='both', which='both', length=0)

legendforplot=[mpatches.Patch(color=colormapCustom[l], label=l) for l in colormapCustom.keys()]
g.ax_heatmap.legend(bbox_to_anchor=(0,-11),loc=2, borderaxespad=0.,handles=legendforplot,ncol=3,frameon=True)
plt.show()
