metaDataCluster1=cellAnnotation[cellAnnotation['CCLE_ID'].isin(cluster_members[0])]
metaDataCluster2=cellAnnotation[cellAnnotation['CCLE_ID'].isin(cluster_members[1])]
metaDataCluster3=cellAnnotation[cellAnnotation['CCLE_ID'].isin(cluster_members[2])]
metaDataCluster4=cellAnnotation[cellAnnotation['CCLE_ID'].isin(cluster_members[3])]
metaDataCluster5=cellAnnotation[cellAnnotation['CCLE_ID'].isin(cluster_members[4])]
metaDataCluster6=cellAnnotation[cellAnnotation['CCLE_ID'].isin(cluster_members[5])]
metaDataCluster7=cellAnnotation[cellAnnotation['CCLE_ID'].isin(cluster_members[6])]
metaDataCluster8=cellAnnotation[cellAnnotation['CCLE_ID'].isin(cluster_members[7])]

columnnames=['CCLE_ID','Site_Primary','Histology','Hist_Subtype1','Disease','Original.Source.of.Cell.Line','Gender','type','Race']
metaDataCluster1_ForPaper=metaDataCluster1[columnnames]
metaDataCluster2_ForPaper=metaDataCluster2[columnnames]
metaDataCluster3_ForPaper=metaDataCluster3[columnnames]
metaDataCluster4_ForPaper=metaDataCluster4[columnnames]
metaDataCluster5_ForPaper=metaDataCluster5[columnnames]
metaDataCluster6_ForPaper=metaDataCluster6[columnnames]
metaDataCluster7_ForPaper=metaDataCluster7[columnnames]
metaDataCluster8_ForPaper=metaDataCluster8[columnnames]


mutationsCluster1=microsatellite[microsatellite['CCLE_ID'].isin(cluster_members[0])]
mutationsCluster2=microsatellite[microsatellite['CCLE_ID'].isin(cluster_members[1])]
mutationsCluster3=microsatellite[microsatellite['CCLE_ID'].isin(cluster_members[2])]
mutationsCluster4=microsatellite[microsatellite['CCLE_ID'].isin(cluster_members[3])]
mutationsCluster5=microsatellite[microsatellite['CCLE_ID'].isin(cluster_members[4])]
mutationsCluster6=microsatellite[microsatellite['CCLE_ID'].isin(cluster_members[5])]
mutationsCluster7=microsatellite[microsatellite['CCLE_ID'].isin(cluster_members[6])]
mutationsCluster8=microsatellite[microsatellite['CCLE_ID'].isin(cluster_members[7])]



pcbCombinedDataFrame=createdataframeclusters(input=[popExposureCluster1,popExposureCluster2,popExposureCluster3],columnNames=['PCB.170'])
combinedPcbDataFrame=meltdataframetolongformat(input=pcbCombinedDataFrame,idVars=['Cluster'],valueVars=['PCB.170'],varName='PCBs')
creategroupprofileplot(inputData=combinedPcbDataFrame,xVar='Cluster',yVar='value',hueVar='PCBs',title='PCB exposure per group',xlabel='PCBs',ylabel='Exposure')