import pandas as pd
from CreateGenomicInteractionNetworks import CreateGenomicInteractionNetwork
from FilterGeneGeneInteractions import FilterGeneGeneInteractions
from mapTranscriptionFactors import mapTranscriptionFactors

pathToFeatures="C:/Users/tim.kuijpers/Desktop/2019 CCLE analysis/Manuscript/Features cluster 7/GenesToConstructNetwork.txt"
pathToCpGfeatures="C:/Users/tim.kuijpers/Desktop/2019 CCLE analysis/Manuscript/Features cluster 7/Cluster7_FeatureCpGNames.csv"
featuresCluster7=pd.read_csv(pathToFeatures,sep="\t")
cpgFeaturesCluster7=pd.read_csv(pathToCpGfeatures,sep=";")
featuresCluster7GeneSymbols=list(featuresCluster7.Gene)
cpgFeaturesCluster7Symbol=list(cpgFeaturesCluster7.GeneSymbol.unique())

# Create the object for the genomic interaction network
genomicInteractionNetworkCluster=CreateGenomicInteractionNetwork()
genomicInteractionNetworkCluster.set_genes_as_nodes(featuresCluster7GeneSymbols)
genomicInteractionNetworkCluster.get_cpg_gene_interactions(option='FromGeneList',genelist=cpgFeaturesCluster7Symbol)
genomicInteractionNetworkCluster.get_gene_gene_interactions()
genomicInteractionNetworkCluster.get_transcription_factor_regulation(pathToTFlib='E:/GIT repos/NetworkSearchDatabase/transcription_catalogue.txt',sepa='\t',encoding='unicode_escape',genes=featuresCluster7GeneSymbols)
genomicInteractionNetworkCluster.format_gene_gene_networks()
genomicInteractionNetworkCluster.built_the_network()
