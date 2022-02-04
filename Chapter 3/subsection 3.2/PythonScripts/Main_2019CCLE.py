import pandas as pd
import numpy as np
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
from getMetaDataInfo import getMetaDataInfo
from getClusterWithSample import getclusterwithsample,createdataframeclusters,meltdataframetolongformat,determineClusterMembers,getclustermembersfromarray

def extractFeatures(Wmatrix, methodScore,score_threshold, clustersize):
	# Script to score the basis matrix
	print('Get the basis matrix')
	scores = np.zeros(Wmatrix.shape[0])
	print('score features')
	if (methodScore == 'Kim'):
		# Part to score the matrix
		for x in range(Wmatrix.shape[0]):
			probability = Wmatrix[x, :] / Wmatrix[x, :].sum()
			scores[x] = np.dot(probability, np.log2(probability + 0.0001).T)
		scores = 1. + 1. / np.log2(Wmatrix.shape[1]) * scores
		scoredMatrix = scores
		threshold = np.median(scores) + 3 * np.median(abs(scores - np.median(scores)))
		sel = scores > threshold
		m = np.median(scores)
		sel = np.array([sel[i] and np.max(Wmatrix[i, :]) > m for i in range(Wmatrix.shape[0])])
		features_per_cluster = [[] for i in range(Wmatrix.shape[1])]
		for i in range(Wmatrix.shape[0]):
			if sel[i] == True:
				pos = np.argmax(Wmatrix[i, :])
				features_per_cluster[pos].append(i)

	if (methodScore == 'Max_Kim_Score'):
		scoredMatrix = np.argmax(Wmatrix, axis=1)
		# Part to score the matrix
		for x in range(Wmatrix.shape[0]):
			probability = Wmatrix[x, :] / Wmatrix[x, :].sum()
			scores[x] = np.dot(probability, np.log2(probability + 0.0001).T)
		scores = 1. + 1. / np.log2(Wmatrix.shape[1]) * scores
		scoredMatrix = scores
		sel = scores >=score_threshold
		m = np.median(scores)
		sel = np.array([sel[i] and np.max(Wmatrix[i, :]) > m for i in range(Wmatrix.shape[0])])
		features_per_cluster = [[] for i in range(Wmatrix.shape[1])]
		for i in range(Wmatrix.shape[0]):
			if sel[i] == True:
				pos = np.argmax(Wmatrix[i, :])
				features_per_cluster[pos].append(i)
	return scoredMatrix, sel, features_per_cluster

def determineClusterMembers(Hmatrix, samplenames):
	"For each column in H, maximum value will be determined and right cluster assigned"
	index_maxvalue = np.argmax(Hmatrix, axis=0)
	ClusterMembers = []
	for cluster in range(np.min(index_maxvalue), np.max(index_maxvalue) + 1):
		ClusterMembers.append([i for indx, i in enumerate(samplenames) if index_maxvalue[indx] == cluster])
	return ClusterMembers, index_maxvalue

# First read the file that contains the list of features for each cluster
# Note that for the paper only cluster 6 and cluster 7 are analyzed
path_to_genes_cluster6 = "E:/CCLE manusscript/Data\Results/Results used for paper/NMF clustering results/Genes_cluster6_k8_nrun40.csv"
cluster6_feature_genes = pd.read_csv(path_to_genes_cluster6,sep=',',index_col=False)
path_to_genes_cluster7 = "E:/CCLE manusscript/Data/Results/Results used for paper/NMF clustering results/Genes_cluster7_k8_nrun40.csv"
cluster7_feature_genes = pd.read_csv(path_to_genes_cluster7,sep=",",index_col=False)

# Read the gene expression data
path_to_gene_expression_file = "E:/CCLE manusscript/Data/Cluster results/GeneExpressionSet_CCLE.csv"
gene_expression_data=pd.read_csv(path_to_gene_expression_file,sep=",",index_col=0)
# Load the W matrix to score the features with only features equal to one
pickle_in = open("E:/CCLE manusscript/Data/Results/Results used for paper/NMF clustering results/Wmatrix_CCLE.pickle","rb")
Wmatrix=pickle.load(pickle_in)
pickle_in1=open("E:/CCLE manusscript/Data/Results/Results used for paper/NMF clustering results/Hmatrix_CCLE.pickle","rb")
Hmatrix=pickle.load(pickle_in1)
path_to_gene_expression_file = "E:/CCLE manusscript/Data/Cluster results/GeneExpressionSet_CCLE.csv"
gene_expression_data=pd.read_csv(path_to_gene_expression_file,sep=",",index_col=0)
# Load the methylation data
path_to_methylation_file="E:/CCLE manusscript/Data/Cluster results/MethylationDataSet_CCLE.csv"
methylationdata=pd.read_csv(path_to_methylation_file,sep=',',index_col=0)
# Load the meta data
file_microsatellite="E:/CCLE manusscript/Data/Results/Results used for paper/Microsatellite_Instability_Info.csv"
microsatellite=pd.read_csv(file_microsatellite,sep=";")
meta_data_lineage="E:/CCLE manusscript/Data/RAW data/Cell_lines_annotations_20181226.txt"
meta_data_lineage=pd.read_csv(meta_data_lineage,sep='\t')
# Get the cell line annotation
pathToAnnotation="E:/2019 CCLE analysis/CellLineAnnotation.txt"
cellAnnotation=pd.read_csv(pathToAnnotation,sep="\t")
# The Wmatrix and H matrix with the lowest error are determined from the NMF results earlier
cpg_wmatrix=Wmatrix[0][32][1]
genes_wmatrix=Wmatrix[0][32][0]
scored_matrix_genes,sel_genes,features_per_cluster_genes=extractFeatures(genes_wmatrix,'Kim',score_threshold=0.95,clustersize=8)
scored_matrix_cpg,sel_cpg,features_per_cluster_cpg=extractFeatures(cpg_wmatrix,'Kim',score_threshold=0.8,clustersize=8)
scored_matrix_maximum_genes,sel_genes,features_maximum_per_cluster_genes=extractFeatures(genes_wmatrix,'Max_Kim_Score',score_threshold=0.95,clustersize=8)
scored_matrix_maximum_cpg,sel_cpg,features_maximum_per_cluster_cpg=extractFeatures(cpg_wmatrix,'Max_Kim_Score',score_threshold=0.75,clustersize=8)
hmatrix_minimal_error=Hmatrix[35]
cluster_members,index_maxvalue=determineClusterMembers(hmatrix_minimal_error,samplenames=gene_expression_data.columns.values)

analyze_haematopoietic_lymphoid_tissue=getMetaDataInfo(cluster_number='cluster7')
analyze_haematopoietic_lymphoid_tissue.samples_in_cluster(cluster_members[6])
analyze_haematopoietic_lymphoid_tissue.set_meta_data(meta_data=microsatellite)
analyze_haematopoietic_lymphoid_tissue.extractMeta(column='CCLE_ID')


analyze_cluster_3=getMetaDataInfo(cluster_number='cluster3')
analyze_cluster_3.samples_in_cluster(cluster_members[2])
analyze_cluster_3.set_meta_data(meta_data=microsatellite)
analyze_cluster_3.extractMeta(column='CCLE_ID')
analyze_cluster_4=getMetaDataInfo(cluster_number='cluster4')
analyze_cluster_4.samples_in_cluster(cluster_members[3])
analyze_cluster_4.set_meta_data(meta_data=microsatellite)
analyze_cluster_4.extractMeta(column='CCLE_ID')
analyze_cluster_5=getMetaDataInfo(cluster_number='cluster5')
analyze_cluster_5.samples_in_cluster(cluster_members[4])
analyze_cluster_5.set_meta_data(meta_data=microsatellite)
analyze_cluster_5.extractMeta(column='CCLE_ID')

methylationDataFeaturesCluster7=methylationdata.iloc[features_maximum_per_cluster_cpg[6],:]
geneExpressionDataFeaturesCluster7=gene_expression_data.iloc[features_maximum_per_cluster_genes[6],:]

# Get the gene expression for CTLA4, CD28, CD86
idGeneCTLA4='ENSG00000163599.10'
idGeneCD28='ENSG00000178562.13'
idGeneCD86="ENSG00000114013.11"

idCpGCTLA4=''
idCpGCD28=''
idCpgCD86=''

geneExpressionCTLA4=pd.DataFrame({'CTLA4':gene_expression_data.loc[idGeneCTLA4,cluster_members[6]]})
geneExpressionCD28=pd.DataFrame({'CD28':gene_expression_data.loc[idGeneCD28,cluster_members[6]]})
geneExpressionCD86=pd.DataFrame({'CD86':gene_expression_data.loc[idGeneCD86,cluster_members[6]]})

cd28cd86axis=pd.concat([geneExpressionCTLA4,geneExpressionCD28,geneExpressionCD86],axis=1)
cd28cd86axisLong=pd.melt(cd28cd86axis,value_vars=['CTLA4','CD28','CD86'],var_name='Gene')

# Investigate the expression of CD28 and CD86 and CTLA4 in the other clusters
cluster1CTLA4=pd.DataFrame({'CTLA4':gene_expression_data.loc[idGeneCTLA4,cluster_members[0]]})
cluster1CTLA4=cluster1CTLA4.assign(Cluster='Cluster 1')
cluster2CTLA4=pd.DataFrame({'CTLA4':gene_expression_data.loc[idGeneCTLA4,cluster_members[1]]})
cluster2CTLA4=cluster2CTLA4.assign(Cluster='Cluster 2')
cluster3CTLA4=pd.DataFrame({'CTLA4':gene_expression_data.loc[idGeneCTLA4,cluster_members[2]]})
cluster3CTLA4=cluster3CTLA4.assign(Cluster='Cluster 3')
cluster4CTLA4=pd.DataFrame({'CTLA4':gene_expression_data.loc[idGeneCTLA4,cluster_members[3]]})
cluster4CTLA4=cluster4CTLA4.assign(Cluster='Cluster 4')
cluster5CTLA4=pd.DataFrame({'CTLA4':gene_expression_data.loc[idGeneCTLA4,cluster_members[4]]})
cluster5CTLA4=cluster5CTLA4.assign(Cluster='Cluster 5')
cluster6CTLA4=pd.DataFrame({'CTLA4':gene_expression_data.loc[idGeneCTLA4,cluster_members[5]]})
cluster6CTLA4=cluster6CTLA4.assign(Cluster='Cluster 6')
cluster7CTLA4=pd.DataFrame({'CTLA4':gene_expression_data.loc[idGeneCTLA4,cluster_members[6]]})
cluster7CTLA4=cluster7CTLA4.assign(Cluster='Cluster 7')
cluster8CTLA4=pd.DataFrame({'CTLA4':gene_expression_data.loc[idGeneCTLA4,cluster_members[7]]})
cluster8CTLA4=cluster8CTLA4.assign(Cluster='Cluster 8')
totalCTLA4expression=pd.concat([cluster1CTLA4,cluster2CTLA4,cluster3CTLA4,cluster4CTLA4,cluster5CTLA4,cluster6CTLA4,cluster7CTLA4,cluster8CTLA4],axis=0)

cluster1CD28=pd.DataFrame({'CD28':gene_expression_data.loc[idGeneCD28,cluster_members[0]]})
cluster1CD28=cluster1CD28.assign(Cluster='Cluster 1')
cluster2CD28=pd.DataFrame({'CD28':gene_expression_data.loc[idGeneCD28,cluster_members[1]]})
cluster2CD28=cluster2CD28.assign(Cluster='Cluster 2')
cluster3CD28=pd.DataFrame({'CD28':gene_expression_data.loc[idGeneCD28,cluster_members[2]]})
cluster3CD28=cluster3CD28.assign(Cluster='Cluster 3')
cluster4CD28=pd.DataFrame({'CD28':gene_expression_data.loc[idGeneCD28,cluster_members[3]]})
cluster4CD28=cluster4CD28.assign(Cluster='Cluster 4')
cluster5CD28=pd.DataFrame({'CD28':gene_expression_data.loc[idGeneCD28,cluster_members[4]]})
cluster5CD28=cluster5CD28.assign(Cluster='Cluster 5')
cluster6CD28=pd.DataFrame({'CD28':gene_expression_data.loc[idGeneCD28,cluster_members[5]]})
cluster6CD28=cluster6CD28.assign(Cluster='Cluster 6')
cluster7CD28=pd.DataFrame({'CD28':gene_expression_data.loc[idGeneCD28,cluster_members[6]]})
cluster7CD28=cluster7CD28.assign(Cluster='Cluster 7')
cluster8CD28=pd.DataFrame({'CD28':gene_expression_data.loc[idGeneCD28,cluster_members[7]]})
cluster8CD28=cluster8CD28.assign(Cluster='Cluster 8')
totalCD28expression=pd.concat([cluster1CD28,cluster2CD28,cluster3CD28,cluster4CD28,cluster5CD28,cluster6CD28,cluster7CD28,cluster8CD28],axis=0)

cluster1CD86=pd.DataFrame({'CD86':gene_expression_data.loc[idGeneCD86,cluster_members[0]]})
cluster1CD86=cluster1CD86.assign(Cluster='Cluster 1')
cluster2CD86=pd.DataFrame({'CD86':gene_expression_data.loc[idGeneCD86,cluster_members[1]]})
cluster2CD86=cluster2CD86.assign(Cluster='Cluster 2')
cluster3CD86=pd.DataFrame({'CD86':gene_expression_data.loc[idGeneCD86,cluster_members[2]]})
cluster3CD86=cluster3CD86.assign(Cluster='Cluster 3')
cluster4CD86=pd.DataFrame({'CD86':gene_expression_data.loc[idGeneCD86,cluster_members[3]]})
cluster4CD86=cluster4CD86.assign(Cluster='Cluster 4')
cluster5CD86=pd.DataFrame({'CD86':gene_expression_data.loc[idGeneCD86,cluster_members[4]]})
cluster5CD86=cluster5CD86.assign(Cluster='Cluster 5')
cluster6CD86=pd.DataFrame({'CD86':gene_expression_data.loc[idGeneCD86,cluster_members[5]]})
cluster6CD86=cluster6CD86.assign(Cluster='Cluster 6')
cluster7CD86=pd.DataFrame({'CD86':gene_expression_data.loc[idGeneCD86,cluster_members[6]]})
cluster7CD86=cluster7CD86.assign(Cluster='Cluster 7')
cluster8CD86=pd.DataFrame({'CD86':gene_expression_data.loc[idGeneCD86,cluster_members[7]]})
cluster8CD86=cluster8CD86.assign(Cluster='Cluster 8')
totalCD86expression=pd.concat([cluster1CD86,cluster2CD86,cluster3CD86,cluster4CD86,cluster5CD86,cluster6CD86,cluster7CD86,cluster8CD86],axis=0)


"""
idCpGSitesEGFR="EGFR_7_55085724_55086724"
idCpGsiteCEBPD="CEBPD_8_48650726_48651726"
cluster1CEBPD=pd.DataFrame({'CEBPD':methylationdata.loc[idCpGsiteCEBPD,cluster_members[0]]})
cluster1CEBPD=cluster1CEBPD.assign(Cluster='Cluster 1')
cluster2CEBPD=pd.DataFrame({'CEBPD':methylationdata.loc[idCpGsiteCEBPD,cluster_members[1]]})
cluster2CEBPD=cluster2CEBPD.assign(Cluster='Cluster 2')
cluster3CEBPD=pd.DataFrame({'CEBPD':methylationdata.loc[idCpGsiteCEBPD,cluster_members[2]]})
cluster3CEBPD=cluster3CEBPD.assign(Cluster='Cluster 3')
cluster4CEBPD=pd.DataFrame({'CEBPD':methylationdata.loc[idCpGsiteCEBPD,cluster_members[3]]})
cluster4CEBPD=cluster4CEBPD.assign(Cluster='Cluster 4')
cluster5CEBPD=pd.DataFrame({'CEBPD':methylationdata.loc[idCpGsiteCEBPD,cluster_members[4]]})
cluster5CEBPD=cluster5CEBPD.assign(Cluster='Cluster 5')
cluster6CEBPD=pd.DataFrame({'CEBPD':methylationdata.loc[idCpGsiteCEBPD,cluster_members[5]]})
cluster6CEBPD=cluster6CEBPD.assign(Cluster='Cluster 6')
cluster7CEBPD=pd.DataFrame({'CEBPD':methylationdata.loc[idCpGsiteCEBPD,cluster_members[6]]})
cluster7CEBPD=cluster7CEBPD.assign(Cluster='Cluster 7')
cluster8CEBPD=pd.DataFrame({'CEBPD':methylationdata.loc[idCpGsiteCEBPD,cluster_members[7]]})
cluster8CEBPD=cluster8CEBPD.assign(Cluster='Cluster 8')
totalCEBPDexpression=pd.concat([cluster1CEBPD,cluster2CEBPD,cluster3CEBPD,cluster4CEBPD,cluster5CEBPD,cluster6CEBPD,cluster7CEBPD,cluster8CEBPD],axis=0)

cluster1EGFR=pd.DataFrame({'EGFR':methylationdata.loc[idCpGSitesEGFR,cluster_members[0]]})
cluster1EGFR=cluster1EGFR.assign(Cluster='Cluster 1')
cluster2EGFR=pd.DataFrame({'EGFR':methylationdata.loc[idCpGSitesEGFR,cluster_members[1]]})
cluster2EGFR=cluster2EGFR.assign(Cluster='Cluster 2')
cluster3EGFR=pd.DataFrame({'EGFR':methylationdata.loc[idCpGSitesEGFR,cluster_members[2]]})
cluster3EGFR=cluster3EGFR.assign(Cluster='Cluster 3')
cluster4EGFR=pd.DataFrame({'EGFR':methylationdata.loc[idCpGSitesEGFR,cluster_members[3]]})
cluster4EGFR=cluster4EGFR.assign(Cluster='Cluster 4')
cluster5EGFR=pd.DataFrame({'EGFR':methylationdata.loc[idCpGSitesEGFR,cluster_members[4]]})
cluster5EGFR=cluster5EGFR.assign(Cluster='Cluster 5')
cluster6EGFR=pd.DataFrame({'EGFR':methylationdata.loc[idCpGSitesEGFR,cluster_members[5]]})
cluster6EGFR=cluster6EGFR.assign(Cluster='Cluster 6')
cluster7EGFR=pd.DataFrame({'EGFR':methylationdata.loc[idCpGSitesEGFR,cluster_members[6]]})
cluster7EGFR=cluster7EGFR.assign(Cluster='Cluster 7')
cluster8EGFR=pd.DataFrame({'EGFR':methylationdata.loc[idCpGSitesEGFR,cluster_members[7]]})
cluster8EGFR=cluster8EGFR.assign(Cluster='Cluster 8')
totalEGFRexpression=pd.concat([cluster1EGFR,cluster2EGFR,cluster3EGFR,cluster4EGFR,cluster5EGFR,cluster6EGFR,cluster7EGFR,cluster8EGFR],axis=0)


idCpGSiteCYFIP1=['CYFIP1_15_22891613_22892613','CYFIP1_15_22891169_22892169','CYFIP1_15_22891218_22892218']
cluster1CYFIP1=pd.DataFrame(methylationdata.loc[idCpGSiteCYFIP1,cluster_members[0]])
cluster1CYFIP1=cluster1CYFIP1.T
cluster1CYFIP1=cluster1CYFIP1.assign(Cluster='Cluster 1')
cluster2CYFIP1=pd.DataFrame(methylationdata.loc[idCpGSiteCYFIP1,cluster_members[1]])
cluster2CYFIP1=cluster2CYFIP1.T
cluster2CYFIP1=cluster2CYFIP1.assign(Cluster='Cluster 2')
cluster3CYFIP1=pd.DataFrame(methylationdata.loc[idCpGSiteCYFIP1,cluster_members[2]])
cluster3CYFIP1=cluster3CYFIP1.T
cluster3CYFIP1=cluster3CYFIP1.assign(Cluster='Cluster 3')
cluster4CYFIP1=pd.DataFrame(methylationdata.loc[idCpGSiteCYFIP1,cluster_members[3]])
cluster4CYFIP1=cluster4CYFIP1.T
cluster4CYFIP1=cluster4CYFIP1.assign(Cluster='Cluster 4')
cluster5CYFIP1=pd.DataFrame(methylationdata.loc[idCpGSiteCYFIP1,cluster_members[4]])
cluster5CYFIP1=cluster5CYFIP1.T
cluster5CYFIP1=cluster5CYFIP1.assign(Cluster='Cluster 5')
cluster6CYFIP1=pd.DataFrame(methylationdata.loc[idCpGSiteCYFIP1,cluster_members[5]])
cluster6CYFIP1=cluster6CYFIP1.T
cluster6CYFIP1=cluster6CYFIP1.assign(Cluster='Cluster 6')
cluster7CYFIP1=pd.DataFrame(methylationdata.loc[idCpGSiteCYFIP1,cluster_members[6]])
cluster7CYFIP1=cluster7CYFIP1.T
cluster7CYFIP1=cluster7CYFIP1.assign(Cluster='Cluster 7')
cluster8CYFIP1=pd.DataFrame(methylationdata.loc[idCpGSiteCYFIP1,cluster_members[7]])
cluster8CYFIP1=cluster8CYFIP1.T
cluster8CYFIP1=cluster8CYFIP1.assign(Cluster='Cluster 8')
totalCYFIP1methylation=pd.concat([cluster1CYFIP1,cluster2CYFIP1,cluster3CYFIP1,cluster4CYFIP1,cluster5CYFIP1,cluster6CYFIP1,cluster7CYFIP1,cluster8CYFIP1],axis=0)
combinedDataFrameCYFIP1Long=meltdataframetolongformat(input=totalCYFIP1methylation,idVars=['Cluster'],valueVars=['CYFIP1_15_22891613_22892613','CYFIP1_15_22891169_22892169','CYFIP1_15_22891218_22892218'],varName='CpGsites')

idCpGYAP1="YAP1_11_101982170_101983170"
cluster1YAP1=pd.DataFrame({'YAP1':methylationdata.loc[idCpGYAP1,cluster_members[0]]})
cluster1YAP1=cluster1YAP1.assign(Cluster='Cluster 1')
cluster2YAP1=pd.DataFrame({'YAP1':methylationdata.loc[idCpGYAP1,cluster_members[1]]})
cluster2YAP1=cluster2YAP1.assign(Cluster='Cluster 2')
cluster3YAP1=pd.DataFrame({'YAP1':methylationdata.loc[idCpGYAP1,cluster_members[2]]})
cluster3YAP1=cluster3YAP1.assign(Cluster='Cluster 3')
cluster4YAP1=pd.DataFrame({'YAP1':methylationdata.loc[idCpGYAP1,cluster_members[3]]})
cluster4YAP1=cluster4YAP1.assign(Cluster='Cluster 4')
cluster5YAP1=pd.DataFrame({'YAP1':methylationdata.loc[idCpGYAP1,cluster_members[4]]})
cluster5YAP1=cluster5YAP1.assign(Cluster='Cluster 5')
cluster6YAP1=pd.DataFrame({'YAP1':methylationdata.loc[idCpGYAP1,cluster_members[5]]})
cluster6YAP1=cluster6YAP1.assign(Cluster='Cluster 6')
cluster7YAP1=pd.DataFrame({'YAP1':methylationdata.loc[idCpGYAP1,cluster_members[6]]})
cluster7YAP1=cluster7YAP1.assign(Cluster='Cluster 7')
cluster8YAP1=pd.DataFrame({'YAP1':methylationdata.loc[idCpGYAP1,cluster_members[7]]})
cluster8YAP1=cluster8YAP1.assign(Cluster='Cluster 8')
totalYAP1methylation=pd.concat([cluster1YAP1,cluster2YAP1,cluster3YAP1,cluster4YAP1,cluster5YAP1,cluster6YAP1,cluster7YAP1,cluster8YAP1],axis=0)

idCpGSiteTEAD4=['TEAD4_12_3067477_3068477', 'TEAD4_12_3067747_3068747']
cluster1TEAD4=pd.DataFrame(methylationdata.loc[idCpGSiteTEAD4, cluster_members[0]])
cluster1TEAD4=cluster1TEAD4.T
cluster1TEAD4=cluster1TEAD4.assign(Cluster='Cluster 1')
cluster2TEAD4=pd.DataFrame(methylationdata.loc[idCpGSiteTEAD4, cluster_members[1]])
cluster2TEAD4=cluster2TEAD4.T
cluster2TEAD4=cluster2TEAD4.assign(Cluster='Cluster 2')
cluster3TEAD4=pd.DataFrame(methylationdata.loc[idCpGSiteTEAD4, cluster_members[2]])
cluster3TEAD4=cluster3TEAD4.T
cluster3TEAD4=cluster3TEAD4.assign(Cluster='Cluster 3')
cluster4TEAD4=pd.DataFrame(methylationdata.loc[idCpGSiteTEAD4, cluster_members[3]])
cluster4TEAD4=cluster4TEAD4.T
cluster4TEAD4=cluster4TEAD4.assign(Cluster='Cluster 4')
cluster5TEAD4=pd.DataFrame(methylationdata.loc[idCpGSiteTEAD4, cluster_members[4]])
cluster5TEAD4=cluster5TEAD4.T
cluster5TEAD4=cluster5TEAD4.assign(Cluster='Cluster 5')
cluster6TEAD4=pd.DataFrame(methylationdata.loc[idCpGSiteTEAD4, cluster_members[5]])
cluster6TEAD4=cluster6TEAD4.T
cluster6TEAD4=cluster6TEAD4.assign(Cluster='Cluster 6')
cluster7TEAD4=pd.DataFrame(methylationdata.loc[idCpGSiteTEAD4, cluster_members[6]])
cluster7TEAD4=cluster7TEAD4.T
cluster7TEAD4=cluster7TEAD4.assign(Cluster='Cluster 7')
cluster8TEAD4=pd.DataFrame(methylationdata.loc[idCpGSiteTEAD4, cluster_members[7]])
cluster8TEAD4=cluster8TEAD4.T
cluster8TEAD4=cluster8TEAD4.assign(Cluster='Cluster 8')
totalTEAD4methylation=pd.concat([cluster1TEAD4, cluster2TEAD4, cluster3TEAD4, cluster4TEAD4, cluster5TEAD4, cluster6TEAD4, cluster7TEAD4, cluster8TEAD4], axis=0)
combinedDataFrameTEAD4Long=meltdataframetolongformat(input=totalTEAD4methylation, idVars=['Cluster'], valueVars=['TEAD4_12_3067477_3068477', 'TEAD4_12_3067747_3068747'], varName='CpGsites')

idcpgKLF5="KLF5_13_73631929_73632929"
cluster1KLF5=pd.DataFrame(methylationdata.loc[idcpgKLF5, cluster_members[0]])
cluster1KLF5=cluster1KLF5.assign(Cluster='Cluster 1')
cluster2KLF5=pd.DataFrame(methylationdata.loc[idcpgKLF5, cluster_members[1]])
cluster2KLF5=cluster2KLF5.assign(Cluster='Cluster 2')
cluster3KLF5=pd.DataFrame(methylationdata.loc[idcpgKLF5, cluster_members[2]])
cluster3KLF5=cluster3KLF5.assign(Cluster='Cluster 3')
cluster4KLF5=pd.DataFrame(methylationdata.loc[idcpgKLF5, cluster_members[3]])
cluster4KLF5=cluster4KLF5.assign(Cluster='Cluster 4')
cluster5KLF5=pd.DataFrame(methylationdata.loc[idcpgKLF5, cluster_members[4]])
cluster5KLF5=cluster5KLF5.assign(Cluster='Cluster 5')
cluster6KLF5=pd.DataFrame(methylationdata.loc[idcpgKLF5, cluster_members[5]])
cluster6KLF5=cluster6KLF5.assign(Cluster='Cluster 6')
cluster7KLF5=pd.DataFrame(methylationdata.loc[idcpgKLF5, cluster_members[6]])
cluster7KLF5=cluster7KLF5.assign(Cluster='Cluster 7')
cluster8KLF5=pd.DataFrame(methylationdata.loc[idcpgKLF5, cluster_members[7]])
cluster8KLF5=cluster8KLF5.assign(Cluster='Cluster 8')
totalKLF5methylation=pd.concat([cluster1KLF5, cluster2KLF5, cluster3KLF5, cluster4KLF5, cluster5KLF5, cluster6KLF5, cluster7KLF5, cluster8KLF5], axis=0)
combinedDataFrameKLF5Long=meltdataframetolongformat(input=totalKLF5methylation, idVars=['Cluster'], valueVars=['KLF5_13_73631929_73632929'], varName='KLF5')

idcpgLATS2="LATS2_13_21635722_21636722"
cluster1LATS2=pd.DataFrame(methylationdata.loc[idcpgLATS2, cluster_members[0]])
cluster1LATS2=cluster1LATS2.assign(Cluster='Cluster 1')
cluster2LATS2=pd.DataFrame(methylationdata.loc[idcpgLATS2, cluster_members[1]])
cluster2LATS2=cluster2LATS2.assign(Cluster='Cluster 2')
cluster3LATS2=pd.DataFrame(methylationdata.loc[idcpgLATS2, cluster_members[2]])
cluster3LATS2=cluster3LATS2.assign(Cluster='Cluster 3')
cluster4LATS2=pd.DataFrame(methylationdata.loc[idcpgLATS2, cluster_members[3]])
cluster4LATS2=cluster4LATS2.assign(Cluster='Cluster 4')
cluster5LATS2=pd.DataFrame(methylationdata.loc[idcpgLATS2, cluster_members[4]])
cluster5LATS2=cluster5LATS2.assign(Cluster='Cluster 5')
cluster6LATS2=pd.DataFrame(methylationdata.loc[idcpgLATS2, cluster_members[5]])
cluster6LATS2=cluster6LATS2.assign(Cluster='Cluster 6')
cluster7LATS2=pd.DataFrame(methylationdata.loc[idcpgLATS2, cluster_members[6]])
cluster7LATS2=cluster7LATS2.assign(Cluster='Cluster 7')
cluster8LATS2=pd.DataFrame(methylationdata.loc[idcpgLATS2, cluster_members[7]])
cluster8LATS2=cluster8LATS2.assign(Cluster='Cluster 8')
totalLATS2methylation=pd.concat([cluster1LATS2, cluster2LATS2, cluster3LATS2, cluster4LATS2, cluster5LATS2, cluster6LATS2, cluster7LATS2, cluster8LATS2], axis=0)
combinedDataFrameLATS2Long=meltdataframetolongformat(input=totalLATS2methylation, idVars=['Cluster'], valueVars=['LATS2_13_21635722_21636722'], varName='LATS2')


geneProbeEGFR="ENSG00000146648.11"
cluster1EGFR=pd.DataFrame({'EGFR':gene_expression_data.loc[geneProbeEGFR,cluster_members[0]]})
cluster1EGFR=cluster1EGFR.assign(Cluster='Cluster 1')
cluster2EGFR=pd.DataFrame({'EGFR':gene_expression_data.loc[geneProbeEGFR,cluster_members[1]]})
cluster2EGFR=cluster2EGFR.assign(Cluster='Cluster 2')
cluster3EGFR=pd.DataFrame({'EGFR':gene_expression_data.loc[geneProbeEGFR,cluster_members[2]]})
cluster3EGFR=cluster3EGFR.assign(Cluster='Cluster 3')
cluster4EGFR=pd.DataFrame({'EGFR':gene_expression_data.loc[geneProbeEGFR,cluster_members[3]]})
cluster4EGFR=cluster4EGFR.assign(Cluster='Cluster 4')
cluster5EGFR=pd.DataFrame({'EGFR':gene_expression_data.loc[geneProbeEGFR,cluster_members[4]]})
cluster5EGFR=cluster5EGFR.assign(Cluster='Cluster 5')
cluster6EGFR=pd.DataFrame({'EGFR':gene_expression_data.loc[geneProbeEGFR,cluster_members[5]]})
cluster6EGFR=cluster6EGFR.assign(Cluster='Cluster 6')
cluster7EGFR=pd.DataFrame({'EGFR':gene_expression_data.loc[geneProbeEGFR,cluster_members[6]]})
cluster7EGFR=cluster7EGFR.assign(Cluster='Cluster 7')
cluster8EGFR=pd.DataFrame({'EGFR':gene_expression_data.loc[geneProbeEGFR,cluster_members[7]]})
cluster8EGFR=cluster8EGFR.assign(Cluster='Cluster 8')
totalEGFRexpression=pd.concat([cluster1EGFR,cluster2EGFR,cluster3EGFR,cluster4EGFR,cluster5EGFR,cluster6EGFR,cluster7EGFR,cluster8EGFR],axis=0)

geneProbeCEBPD="ENSG00000221869.4"
cluster1CEBPD=pd.DataFrame({'CEBPD':gene_expression_data.loc[geneProbeCEBPD,cluster_members[0]]})
cluster1CEBPD=cluster1CEBPD.assign(Cluster='Cluster 1')
cluster2CEBPD=pd.DataFrame({'CEBPD':gene_expression_data.loc[geneProbeCEBPD,cluster_members[1]]})
cluster2CEBPD=cluster2CEBPD.assign(Cluster='Cluster 2')
cluster3CEBPD=pd.DataFrame({'CEBPD':gene_expression_data.loc[geneProbeCEBPD,cluster_members[2]]})
cluster3CEBPD=cluster3CEBPD.assign(Cluster='Cluster 3')
cluster4CEBPD=pd.DataFrame({'CEBPD':gene_expression_data.loc[geneProbeCEBPD,cluster_members[3]]})
cluster4CEBPD=cluster4CEBPD.assign(Cluster='Cluster 4')
cluster5CEBPD=pd.DataFrame({'CEBPD':gene_expression_data.loc[geneProbeCEBPD,cluster_members[4]]})
cluster5CEBPD=cluster5CEBPD.assign(Cluster='Cluster 5')
cluster6CEBPD=pd.DataFrame({'CEBPD':gene_expression_data.loc[geneProbeCEBPD,cluster_members[5]]})
cluster6CEBPD=cluster6CEBPD.assign(Cluster='Cluster 6')
cluster7CEBPD=pd.DataFrame({'CEBPD':gene_expression_data.loc[geneProbeCEBPD,cluster_members[6]]})
cluster7CEBPD=cluster7CEBPD.assign(Cluster='Cluster 7')
cluster8CEBPD=pd.DataFrame({'CEBPD':gene_expression_data.loc[geneProbeCEBPD,cluster_members[7]]})
cluster8CEBPD=cluster8CEBPD.assign(Cluster='Cluster 8')
totalCEBPDexpression=pd.concat([cluster1CEBPD,cluster2CEBPD,cluster3CEBPD,cluster4CEBPD,cluster5CEBPD,cluster6CEBPD,cluster7CEBPD,cluster8CEBPD],axis=0)

geneKLF5="ENSG00000102554.9"
cluster1KLF5=pd.DataFrame({'KLF5':gene_expression_data.loc[geneKLF5,cluster_members[0]]})
cluster1KLF5=cluster1KLF5.assign(Cluster='Cluster 1')
cluster2KLF5=pd.DataFrame({'KLF5':gene_expression_data.loc[geneKLF5,cluster_members[1]]})
cluster2KLF5=cluster2KLF5.assign(Cluster='Cluster 2')
cluster3KLF5=pd.DataFrame({'KLF5':gene_expression_data.loc[geneKLF5,cluster_members[2]]})
cluster3KLF5=cluster3KLF5.assign(Cluster='Cluster 3')
cluster4KLF5=pd.DataFrame({'KLF5':gene_expression_data.loc[geneKLF5,cluster_members[3]]})
cluster4KLF5=cluster4KLF5.assign(Cluster='Cluster 4')
cluster5KLF5=pd.DataFrame({'KLF5':gene_expression_data.loc[geneKLF5,cluster_members[4]]})
cluster5KLF5=cluster5KLF5.assign(Cluster='Cluster 5')
cluster6KLF5=pd.DataFrame({'KLF5':gene_expression_data.loc[geneKLF5,cluster_members[5]]})
cluster6KLF5=cluster6KLF5.assign(Cluster='Cluster 6')
cluster7KLF5=pd.DataFrame({'KLF5':gene_expression_data.loc[geneKLF5,cluster_members[6]]})
cluster7KLF5=cluster7KLF5.assign(Cluster='Cluster 7')
cluster8KLF5=pd.DataFrame({'KLF5':gene_expression_data.loc[geneKLF5,cluster_members[7]]})
cluster8KLF5=cluster8KLF5.assign(Cluster='Cluster 8')
totalKLF5expression=pd.concat([cluster1KLF5,cluster2KLF5,cluster3KLF5,cluster4KLF5,cluster5KLF5,cluster6KLF5,cluster7KLF5,cluster8KLF5],axis=0)


geneYAP1="ENSG00000137693.9"
cluster1YAP1=pd.DataFrame({'YAP1':gene_expression_data.loc[geneYAP1,cluster_members[0]]})
cluster1YAP1=cluster1YAP1.assign(Cluster='Cluster 1')
cluster2YAP1=pd.DataFrame({'YAP1':gene_expression_data.loc[geneYAP1,cluster_members[1]]})
cluster2YAP1=cluster2YAP1.assign(Cluster='Cluster 2')
cluster3YAP1=pd.DataFrame({'YAP1':gene_expression_data.loc[geneYAP1,cluster_members[2]]})
cluster3YAP1=cluster3YAP1.assign(Cluster='Cluster 3')
cluster4YAP1=pd.DataFrame({'YAP1':gene_expression_data.loc[geneYAP1,cluster_members[3]]})
cluster4YAP1=cluster4YAP1.assign(Cluster='Cluster 4')
cluster5YAP1=pd.DataFrame({'YAP1':gene_expression_data.loc[geneYAP1,cluster_members[4]]})
cluster5YAP1=cluster5YAP1.assign(Cluster='Cluster 5')
cluster6YAP1=pd.DataFrame({'YAP1':gene_expression_data.loc[geneYAP1,cluster_members[5]]})
cluster6YAP1=cluster6YAP1.assign(Cluster='Cluster 6')
cluster7YAP1=pd.DataFrame({'YAP1':gene_expression_data.loc[geneYAP1,cluster_members[6]]})
cluster7YAP1=cluster7YAP1.assign(Cluster='Cluster 7')
cluster8YAP1=pd.DataFrame({'YAP1':gene_expression_data.loc[geneYAP1,cluster_members[7]]})
cluster8YAP1=cluster8YAP1.assign(Cluster='Cluster 8')
totalYAP1expression=pd.concat([cluster1YAP1,cluster2YAP1,cluster3YAP1,cluster4YAP1,cluster5YAP1,cluster6YAP1,cluster7YAP1,cluster8YAP1],axis=0)

geneTEAD4="ENSG00000197905.4"
cluster1TEAD4=pd.DataFrame({'TEAD4':gene_expression_data.loc[geneTEAD4,cluster_members[0]]})
cluster1TEAD4=cluster1TEAD4.assign(Cluster='Cluster 1')
cluster2TEAD4=pd.DataFrame({'TEAD4':gene_expression_data.loc[geneTEAD4,cluster_members[1]]})
cluster2TEAD4=cluster2TEAD4.assign(Cluster='Cluster 2')
cluster3TEAD4=pd.DataFrame({'TEAD4':gene_expression_data.loc[geneTEAD4,cluster_members[2]]})
cluster3TEAD4=cluster3TEAD4.assign(Cluster='Cluster 3')
cluster4TEAD4=pd.DataFrame({'TEAD4':gene_expression_data.loc[geneTEAD4,cluster_members[3]]})
cluster4TEAD4=cluster4TEAD4.assign(Cluster='Cluster 4')
cluster5TEAD4=pd.DataFrame({'TEAD4':gene_expression_data.loc[geneTEAD4,cluster_members[4]]})
cluster5TEAD4=cluster5TEAD4.assign(Cluster='Cluster 5')
cluster6TEAD4=pd.DataFrame({'TEAD4':gene_expression_data.loc[geneTEAD4,cluster_members[5]]})
cluster6TEAD4=cluster6TEAD4.assign(Cluster='Cluster 6')
cluster7TEAD4=pd.DataFrame({'TEAD4':gene_expression_data.loc[geneTEAD4,cluster_members[6]]})
cluster7TEAD4=cluster7TEAD4.assign(Cluster='Cluster 7')
cluster8TEAD4=pd.DataFrame({'TEAD4':gene_expression_data.loc[geneTEAD4,cluster_members[7]]})
cluster8TEAD4=cluster8TEAD4.assign(Cluster='Cluster 8')
totalTEAD4expression=pd.concat([cluster1TEAD4,cluster2TEAD4,cluster3TEAD4,cluster4TEAD4,cluster5TEAD4,cluster6TEAD4,cluster7TEAD4,cluster8TEAD4],axis=0)

geneJAG1="ENSG00000101384.7"
cluster1JAG1=pd.DataFrame({'JAG1':gene_expression_data.loc[geneJAG1,cluster_members[0]]})
cluster1JAG1=cluster1JAG1.assign(Cluster='Cluster 1')
cluster2JAG1=pd.DataFrame({'JAG1':gene_expression_data.loc[geneJAG1,cluster_members[1]]})
cluster2JAG1=cluster2JAG1.assign(Cluster='Cluster 2')
cluster3JAG1=pd.DataFrame({'JAG1':gene_expression_data.loc[geneJAG1,cluster_members[2]]})
cluster3JAG1=cluster3JAG1.assign(Cluster='Cluster 3')
cluster4JAG1=pd.DataFrame({'JAG1':gene_expression_data.loc[geneJAG1,cluster_members[3]]})
cluster4JAG1=cluster4JAG1.assign(Cluster='Cluster 4')
cluster5JAG1=pd.DataFrame({'JAG1':gene_expression_data.loc[geneJAG1,cluster_members[4]]})
cluster5JAG1=cluster5JAG1.assign(Cluster='Cluster 5')
cluster6JAG1=pd.DataFrame({'JAG1':gene_expression_data.loc[geneJAG1,cluster_members[5]]})
cluster6JAG1=cluster6JAG1.assign(Cluster='Cluster 6')
cluster7JAG1=pd.DataFrame({'JAG1':gene_expression_data.loc[geneJAG1,cluster_members[6]]})
cluster7JAG1=cluster7JAG1.assign(Cluster='Cluster 7')
cluster8JAG1=pd.DataFrame({'JAG1':gene_expression_data.loc[geneJAG1,cluster_members[7]]})
cluster8JAG1=cluster8JAG1.assign(Cluster='Cluster 8')
totalJAG1expression=pd.concat([cluster1JAG1,cluster2JAG1,cluster3JAG1,cluster4JAG1,cluster5JAG1,cluster6JAG1,cluster7JAG1,cluster8JAG1],axis=0)

geneLATS2="ENSG00000150457.7"
cluster1LATS2=pd.DataFrame({'LATS2':gene_expression_data.loc[geneLATS2,cluster_members[0]]})
cluster1LATS2=cluster1LATS2.assign(Cluster='Cluster 1')
cluster2LATS2=pd.DataFrame({'LATS2':gene_expression_data.loc[geneLATS2,cluster_members[1]]})
cluster2LATS2=cluster2LATS2.assign(Cluster='Cluster 2')
cluster3LATS2=pd.DataFrame({'LATS2':gene_expression_data.loc[geneLATS2,cluster_members[2]]})
cluster3LATS2=cluster3LATS2.assign(Cluster='Cluster 3')
cluster4LATS2=pd.DataFrame({'LATS2':gene_expression_data.loc[geneLATS2,cluster_members[3]]})
cluster4LATS2=cluster4LATS2.assign(Cluster='Cluster 4')
cluster5LATS2=pd.DataFrame({'LATS2':gene_expression_data.loc[geneLATS2,cluster_members[4]]})
cluster5LATS2=cluster5LATS2.assign(Cluster='Cluster 5')
cluster6LATS2=pd.DataFrame({'LATS2':gene_expression_data.loc[geneLATS2,cluster_members[5]]})
cluster6LATS2=cluster6LATS2.assign(Cluster='Cluster 6')
cluster7LATS2=pd.DataFrame({'LATS2':gene_expression_data.loc[geneLATS2,cluster_members[6]]})
cluster7LATS2=cluster7LATS2.assign(Cluster='Cluster 7')
cluster8LATS2=pd.DataFrame({'LATS2':gene_expression_data.loc[geneLATS2,cluster_members[7]]})
cluster8LATS2=cluster8LATS2.assign(Cluster='Cluster 8')
totalLATS2expression=pd.concat([cluster1LATS2,cluster2LATS2,cluster3LATS2,cluster4LATS2,cluster5LATS2,cluster6LATS2,cluster7LATS2,cluster8LATS2],axis=0)

geneNFIB= "NFIB_9_14314045_14315045"
cluster1NFIB=pd.DataFrame({'NFIB':methylationdata.loc[geneNFIB, cluster_members[0]]})
cluster1NFIB=cluster1NFIB.assign(Cluster='Cluster 1')
cluster2NFIB=pd.DataFrame({'NFIB':methylationdata.loc[geneNFIB, cluster_members[1]]})
cluster2NFIB=cluster2NFIB.assign(Cluster='Cluster 2')
cluster3NFIB=pd.DataFrame({'NFIB':methylationdata.loc[geneNFIB, cluster_members[2]]})
cluster3NFIB=cluster3NFIB.assign(Cluster='Cluster 3')
cluster4NFIB=pd.DataFrame({'NFIB':methylationdata.loc[geneNFIB, cluster_members[3]]})
cluster4NFIB=cluster4NFIB.assign(Cluster='Cluster 4')
cluster5NFIB=pd.DataFrame({'NFIB':methylationdata.loc[geneNFIB, cluster_members[4]]})
cluster5NFIB=cluster5NFIB.assign(Cluster='Cluster 5')
cluster6NFIB=pd.DataFrame({'NFIB':methylationdata.loc[geneNFIB, cluster_members[5]]})
cluster6NFIB=cluster6NFIB.assign(Cluster='Cluster 6')
cluster7NFIB=pd.DataFrame({'NFIB':methylationdata.loc[geneNFIB, cluster_members[6]]})
cluster7NFIB=cluster7NFIB.assign(Cluster='Cluster 7')
cluster8NFIB=pd.DataFrame({'NFIB':methylationdata.loc[geneNFIB, cluster_members[7]]})
cluster8NFIB=cluster8NFIB.assign(Cluster='Cluster 8')
totalNFIBexpression=pd.concat([cluster1NFIB, cluster2NFIB, cluster3NFIB, cluster4NFIB, cluster5NFIB, cluster6NFIB, cluster7NFIB, cluster8NFIB], axis=0)

geneLRRC49= "LRRC49_15_71183781_71184781"
cluster1LRRC49=pd.DataFrame({'LRRC49':methylationdata.loc[geneNFIB, cluster_members[0]]})
cluster1LRRC49=cluster1LRRC49.assign(Cluster='Cluster 1')
cluster2LRRC49=pd.DataFrame({'LRRC49':methylationdata.loc[geneNFIB, cluster_members[1]]})
cluster2LRRC49=cluster2LRRC49.assign(Cluster='Cluster 2')
cluster3LRRC49=pd.DataFrame({'LRRC49':methylationdata.loc[geneNFIB, cluster_members[2]]})
cluster3LRRC49=cluster3LRRC49.assign(Cluster='Cluster 3')
cluster4LRRC49=pd.DataFrame({'LRRC49':methylationdata.loc[geneNFIB, cluster_members[3]]})
cluster4LRRC49=cluster4LRRC49.assign(Cluster='Cluster 4')
cluster5LRRC49=pd.DataFrame({'LRRC49':methylationdata.loc[geneNFIB, cluster_members[4]]})
cluster5LRRC49=cluster5LRRC49.assign(Cluster='Cluster 5')
cluster6LRRC49=pd.DataFrame({'LRRC49':methylationdata.loc[geneNFIB, cluster_members[5]]})
cluster6LRRC49=cluster6LRRC49.assign(Cluster='Cluster 6')
cluster7LRRC49=pd.DataFrame({'LRRC49':methylationdata.loc[geneNFIB, cluster_members[6]]})
cluster7LRRC49=cluster7LRRC49.assign(Cluster='Cluster 7')
cluster8LRRC49=pd.DataFrame({'LRRC49':methylationdata.loc[geneNFIB, cluster_members[7]]})
cluster8LRRC49=cluster8LRRC49.assign(Cluster='Cluster 8')
totalRRC49expression=pd.concat([cluster1LRRC49, cluster2LRRC49, cluster3LRRC49, cluster4LRRC49, cluster5LRRC49, cluster6LRRC49, cluster7LRRC49, cluster8LRRC49], axis=0)

geneARHGAP29= "ARHGAP29_1_94713316_94714316"
cluster1ARHGAP29=pd.DataFrame({'ARHGAP29':methylationdata.loc[geneNFIB, cluster_members[0]]})
cluster1ARHGAP29=cluster1ARHGAP29.assign(Cluster='Cluster 1')
cluster2ARHGAP29=pd.DataFrame({'ARHGAP29':methylationdata.loc[geneNFIB, cluster_members[1]]})
cluster2ARHGAP29=cluster2ARHGAP29.assign(Cluster='Cluster 2')
cluster3ARHGAP29=pd.DataFrame({'ARHGAP29':methylationdata.loc[geneNFIB, cluster_members[2]]})
cluster3ARHGAP29=cluster3ARHGAP29.assign(Cluster='Cluster 3')
cluster4ARHGAP29=pd.DataFrame({'ARHGAP29':methylationdata.loc[geneNFIB, cluster_members[3]]})
cluster4ARHGAP29=cluster4ARHGAP29.assign(Cluster='Cluster 4')
cluster5ARHGAP29=pd.DataFrame({'ARHGAP29':methylationdata.loc[geneNFIB, cluster_members[4]]})
cluster5ARHGAP29=cluster5ARHGAP29.assign(Cluster='Cluster 5')
cluster6ARHGAP29=pd.DataFrame({'ARHGAP29':methylationdata.loc[geneNFIB, cluster_members[5]]})
cluster6ARHGAP29=cluster6ARHGAP29.assign(Cluster='Cluster 6')
cluster7ARHGAP29=pd.DataFrame({'ARHGAP29':methylationdata.loc[geneNFIB, cluster_members[6]]})
cluster7ARHGAP29=cluster7ARHGAP29.assign(Cluster='Cluster 7')
cluster8ARHGAP29=pd.DataFrame({'ARHGAP29':methylationdata.loc[geneNFIB, cluster_members[7]]})
cluster8ARHGAP29=cluster8ARHGAP29.assign(Cluster='Cluster 8')
totalARHGAP29expression=pd.concat([cluster1ARHGAP29, cluster2ARHGAP29, cluster3ARHGAP29, cluster4ARHGAP29, cluster5ARHGAP29, cluster6ARHGAP29, cluster7ARHGAP29, cluster8ARHGAP29], axis=0)


LRRC49= "ENSG00000137821.7"
cluster1LRRC49=pd.DataFrame({'LRRC49':gene_expression_data.loc[LRRC49, cluster_members[0]]})
cluster1LRRC49=cluster1LRRC49.assign(Cluster='Cluster 1')
cluster2LRRC49=pd.DataFrame({'LRRC49':gene_expression_data.loc[LRRC49, cluster_members[1]]})
cluster2LRRC49=cluster2LRRC49.assign(Cluster='Cluster 2')
cluster3LRRC49=pd.DataFrame({'LRRC49':gene_expression_data.loc[LRRC49, cluster_members[2]]})
cluster3LRRC49=cluster3LRRC49.assign(Cluster='Cluster 3')
cluster4LRRC49=pd.DataFrame({'LRRC49':gene_expression_data.loc[LRRC49, cluster_members[3]]})
cluster4LRRC49=cluster4LRRC49.assign(Cluster='Cluster 4')
cluster5LRRC49=pd.DataFrame({'LRRC49':gene_expression_data.loc[LRRC49, cluster_members[4]]})
cluster5LRRC49=cluster5LRRC49.assign(Cluster='Cluster 5')
cluster6LRRC49=pd.DataFrame({'LRRC49':gene_expression_data.loc[LRRC49, cluster_members[5]]})
cluster6LRRC49=cluster6LRRC49.assign(Cluster='Cluster 6')
cluster7LRRC49=pd.DataFrame({'LRRC49':gene_expression_data.loc[LRRC49, cluster_members[6]]})
cluster7LRRC49=cluster7LRRC49.assign(Cluster='Cluster 7')
cluster8LRRC49=pd.DataFrame({'LRRC49':gene_expression_data.loc[LRRC49, cluster_members[7]]})
cluster8LRRC49=cluster8LRRC49.assign(Cluster='Cluster 8')
totalLRRC49expression=pd.concat([cluster1LRRC49, cluster2LRRC49, cluster3LRRC49, cluster4LRRC49, cluster5LRRC49, cluster6LRRC49, cluster7LRRC49, cluster8LRRC49], axis=0)

plt.figure()
ax = sns.boxplot(y='DNA',x='value',hue='Cluster',data=meltedFormatCpG,color='grey', width=.8)
plt.setp(ax.get_xticklabels(), rotation=45)
plt.legend(bbox_to_anchor=(0.5, -0.1), loc=0, borderaxespad=0.)
plt.show()
"""
