# -*- coding: utf-8 -*-
"""
Author: David Angeles
Date: 
"""

from __future__ import division, print_function, absolute_import
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.special
import scipy.misc
from scipy import stats
from scipy import optimize
from scipy import misc
import seaborn as sns


input_dir= "../input/"
dictionary= pd.read_csv(input_dir+'smalldictionary.txt', sep= ',')

##==============================================================================
## If you don't have the complete motifAssociations file, uncomment and run this
##==============================================================================
#fileName1= 'INTERGENICDHS_DREME_E0.05_FIMO_0.05_OUT.annotated.uniq.txt'
#fileName2= 'INTRONDHS_DREME_E0.05_FIMO_0.05_OUT.annotated.uniq.txt'
#fileName3= 'NONCODINGDHS_DREME_E0.05_FIMO_0.05_OUT.annotated.uniq.txt'
#fileName4= 'PROMOTERDHS_DREME_E0.05_FIMO_0.05_OUT.annotated.uniq.txt'
#
#df1= pd.read_csv(input_dir+fileName1, sep= '\t')
#df2= pd.read_csv(input_dir+fileName2, sep= '\t')
#df3= pd.read_csv(input_dir+fileName3, sep= '\t')
#df4= pd.read_csv(input_dir+fileName4, sep= '\t')
#
#df1.drop(['chromosome', 'startDHS', 'stopDHS', 'FIMO', 'chromosomeDHS', 'startDHS.1', 'stopDHS.1', 'genomicPosition', 'distance', 'name', 'human'], 1, inplace= True)
#df2.drop(['chromosome', 'startDHS', 'stopDHS', 'FIMO', 'chromosomeDHS', 'startDHS.1', 'stopDHS.1', 'genomicPosition', 'distance', 'name', 'human'], 1, inplace= True)
#df3.drop(['chromosome', 'startDHS', 'stopDHS', 'FIMO', 'chromosomeDHS', 'startDHS.1', 'stopDHS.1', 'genomicPosition', 'distance', 'name', 'human'], 1, inplace= True)
#df4.drop(['chromosome', 'startDHS', 'stopDHS', 'FIMO', 'chromosomeDHS', 'startDHS.1', 'stopDHS.1', 'genomicPosition', 'distance', 'name', 'human'], 1, inplace= True)
#
#df1.to_csv(input_dir+'motifAssociationList.txt', sep= ',', index= False)
#with open(input_dir+'motifAssociationList.txt', 'a') as myfile:
#    for i in range(len(df2)):
#        myfile.write(df2.loc[i][0]+','+df2.loc[i][1]+'\n')
#    for i in range(len(df3)):
#        myfile.write(df3.loc[i][0]+','+df3.loc[i][1]+'\n')
#    for i in range(len(df4)):
#        myfile.write(df4.loc[i][0]+','+df4.loc[i][1]+'\n')
#    myfile.close()


#df= pd.read_csv(input_dir+'motifAssociationList.txt', sep= ',')
#print(len(df))
#print(len(df.groupby('MOTIF')))
#
##Uncomment this if you don't have the file motifAssociationListWithAnatomyTerms
##==============================================================================
## Remove all genes not in the dictionary.
## Also, drop any isoforms. 
##==============================================================================
#clean= lambda x:  dictionary.loc[dictionary.wbid == x].wbid.any() == False
#notInDic= df.wbid.map(clean)
#df= df[~notInDic]
#df.drop_duplicates(inplace= True)
#
dictionary.set_index('wbid', drop= True, inplace= True)
            
#with open(input_dir+'motifAssociationListWithAnatomyTerms.txt', 'w') as myfile:
#    header= ",".join(str(x) for x in df.columns)+','
#    header2=','.join(str(x) for x in dictionary.columns) 
#    myfile.write(header+header2+'\n')
#    for i in range(len(df)):
#        #terms associated with the anatomy
#        terms= ",".join(str(x) for x in dictionary.loc[df.iloc[i][1]])
#        string= df.iloc[i][0] +','+df.iloc[i][1]+','
#        myfile.write(string+terms+'\n')
#    myfile.close()

df= pd.read_csv(input_dir+'motifAssociationListWithAnatomyTerms.txt', sep= ',')

#Remove all motifs with less than 5 associated genes in the dictionary
df= df.groupby('MOTIF').filter(lambda x: len(x) > 5)
df.reset_index()
#Group by into motifs
groupDF= df.groupby('MOTIF')
#print(len(groupDF))
#
#print(groupDF.size())





##==============================================================================
## Plot some simple covariance matrices
##==============================================================================
#def axisLabelerMotifs(x, pos):
#    return str(groupDF.size().index[x])
#def axisLabelerTissues(x, pos):
#    return str(df.columns[x+2])
##Plot the dictionary covariance matrix:
#hmatrix=dictionary.cov()
#
#fig, ax = plt.subplots(figsize=(20, 20))
#ax.set_title('Covariance Between Labels in Dictionary')
#ax = plt.gca()
#plt.xticks(np.arange(0,len(df.columns)-2, 1.0), rotation= 'vertical')
#ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(axisLabelerTissues))
#plt.yticks(np.arange(0,len(df.columns)-2, 1.0))
#ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(axisLabelerTissues))
#plt.draw()
#plt.imshow(hmatrix.values)
#plt.ylabel('Tissue')
#plt.xlabel('Tissue')
#ax.set_aspect('equal')
#plt.colorbar(orientation= 'horizontal')
#fig.savefig('../output/Raw Covariance of Dictionary.png')
#plt.close()
#
##Plot the raw covariance matrix
#hmatrix=df.cov()
#
#fig, ax = plt.subplots(figsize=(20, 20))
#ax.set_title('Covariance Between Labels in Motif Database')
#ax = plt.gca()
#plt.xticks(np.arange(0,len(df.columns)-2, 1.0), rotation= 'vertical')
#ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(axisLabelerTissues))
#plt.yticks(np.arange(0,len(df.columns)-2, 1.0))
#ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(axisLabelerTissues))
#plt.draw()
#plt.imshow(hmatrix.values)
#plt.ylabel('Tissue')
#plt.xlabel('Tissue')
#ax.set_aspect('equal')
#plt.colorbar(orientation= 'horizontal')
#fig.savefig('../output/Raw Covariance Matrix of Motif Database.png')
#plt.close()
#
##Plot the difference between raw covariance matrix of dictionary and df
#hmatrix=np.abs(df.cov()-dictionary.cov())/dictionary.cov()
#
#fig, ax = plt.subplots(figsize=(20, 20))
#ax.set_title('Difference in Covariance Values Between Database and Dictionary')
#ax = plt.gca()
#plt.xticks(np.arange(0,len(df.columns)-2, 1.0), rotation= 'vertical')
#ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(axisLabelerTissues))
#plt.yticks(np.arange(0,len(df.columns)-2, 1.0))
#ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(axisLabelerTissues))
#plt.draw()
#plt.imshow(hmatrix.values)
#plt.ylabel('Tissue')
#plt.xlabel('Tissue')
#ax.set_aspect('equal')
#plt.colorbar(orientation= 'horizontal')
#fig.savefig('../output/Difference of Covariances Between Motifs and Dictionary.png')
#plt.close()
#
##Plot the raw excess matrix
#hmatrix= (groupDF.sum().transpose())/groupDF.size()
#
#fig, ax = plt.subplots(figsize=(20, 20))
#ax.set_title('Sum of Tissue Labels/Number of Genes Associated With Motif for each motif. Motifs with <5 associated genes not plotted')
#ax = plt.gca()
#plt.xticks(np.arange(0,len(groupDF), 1.0), rotation= 'vertical')
#ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(axisLabelerMotifs))
#plt.yticks(np.arange(0,len(df.columns)-2, 1.0))
#ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(axisLabelerTissues))
#plt.draw()
#plt.imshow(hmatrix.values)
#plt.ylabel('Tissue Label')
#plt.xlabel('Motif')
#ax.set_aspect('equal')
#plt.colorbar(orientation= 'horizontal')
#fig.savefig('../output/Tissue Counts Pero Motif Normalized By Motif Size of Filtered Motifs.png')
#plt.close()
#
##Plot the uncorrected excess matrix, correcting for the excess of each and every single tissue
#hmatrix= (groupDF.sum().transpose())/len(df)
#
#fig, ax = plt.subplots(figsize=(20, 20))
#ax.set_title('Sum of Tissue Labels/Total Number of Genes in Database for each motif. Motifs with <5 associated genes not plotted or included in database')
#ax = plt.gca()
#plt.xticks(np.arange(0,len(groupDF), 1.0), rotation= 'vertical')
#ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(axisLabelerMotifs))
#plt.yticks(np.arange(0,len(df.columns)-2, 1.0))
#ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(axisLabelerTissues))
#plt.draw()
#plt.imshow(hmatrix.values)
#plt.ylabel('Tissue Label')
#plt.xlabel('Motif')
#ax.set_aspect('equal')
#plt.colorbar(orientation= 'horizontal')
#fig.savefig('../output/Tissue Counts Pero Motif Normalized By Total Database Size of Filtered Motifs.png')
#plt.close()





#==============================================================================
# Get ready for the randomizer. Margie has a few motifs she'd like to test, so let's take a look at them. 
#==============================================================================

#vector= [ 'AAAATTCMAAA', 'ACTACAAACTAC', 'CTTGTACGGAA', 'CTYCAGCTCC', 'TATTTYAAAAA', 'CGCGACGCR',  'GCRGCCGACA', 'ACCGCRMCGC', 'GCTGCTGCY', 'ACAGAACCGTGG', 'CAACGATGCTC', 'CGCGCAAATGA', 'CGYCAAGGCAC']
vector= [ 'AAAATTCMAAA', 'ACTACAAACTAC', 'CTYCAGCTCC', 'CGCGACGCR',  'GCRGCCGACA', 'ACCGCRMCGC', 'GCTGCTGCY', 'ACAGAACCGTGG', 'CAACGATGCTC']

handle= []
H= {}
for element in vector:
    handle.append(groupDF.get_group(element).sum()[2:])
    H[element]= groupDF.get_group(element).sum()[2:]

maxIndex= []
for element in handle:
    maxIndex.append(sorted(range(len(element[:])), key=lambda i: element[i])[-3:])

j= 0
for element in maxIndex:
    for i in element:
        print(vector[j], df.columns[i+2], handle[j][i])
    j+=1;


#==============================================================================
# Random Simulator
#==============================================================================
#Random simulator:
def axisLabelerTissues(x, pos):
    return str(df.columns[x+1])

    

nInitial= 10**5
index= np.array(df.index)
data= np.zeros(shape= (len(vector), nInitial, len(handle[0])))

entry= {}
i= 0;
for j in vector:
    entry[j]= i
    i+=1;

exact= np.array([])
for j in vector:
    n= nInitial
    sizeOfCurrentGroup= groupDF.get_group(j).shape[0]  
    n= n*sizeOfCurrentGroup;
        
    choose= np.random.choice(index, n)
    sumoflabels= np.array([])
    l= 0;
    for k in choose:
        if sumoflabels.size== 0:
            sumoflabels= np.array(df.loc[k][2:])
        else:
            if l%sizeOfCurrentGroup != 0:
                sumoflabels+= np.array(df.loc[k][2:])
            else:        
                data[entry[j], l/sizeOfCurrentGroup-1, :]= sumoflabels
                sumoflabels= np.array(df.loc[k][2:])
                print('l', l/sizeOfCurrentGroup)
                
                
    with open('../output/randomizationResultsMargieTissueEnrichment'+j+'.txt', 'w') as myfile:
        myfile.write('Motif,Tissue,FractionOfHitsWithSimilarOrGreaterEnrichmentOrDepletion,Enrichment\n')
        for j in vector:
            for i in range(len(handle[0])-1):
                dummy= data[entry[j], :, i]
                perc= float(len(dummy[dummy> H[j][i]])/len(dummy))
                if perc < 0.05:
                    myfile.write(j+ ","+  df.columns[i+2]+ ","+str(perc)+ ',Enriched\n')
        l+=1


#for j in vector: 
#    sns.set_style("whitegrid")
#    fig, ax = plt.subplots(figsize=(20, 20))
#    ax.set_title('Randomization distributions (violinplots) and actual values of each  anatomic tissue for motif {0}'.format(j))
#
#    ax = plt.gca()
#
#    plot1= sns.violinplot(data[entry[j], :, :])
#    plot2= plt.plot(np.arange(1, len(df.columns)-1), handle[vector.index(j)][:-1], 'k--')
#
#
#    plt.xticks(np.arange(1,len(df.columns)-1, 1.0), rotation= 'vertical')
#    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(axisLabelerTissues))
#
#
#    plt.draw()
#    plt.ylabel('Frequency')
#    plt.xlabel('Tissue')
#    ax.set_aspect('equal')
#    fig.savefig('Randomization distributions (violinplots) and actual values of each anatomic tissue for motif {0}.png'.format(j))
#    plt.close()
    
#count the percentage of hits that are below the hits in each motif

            
##Take a slice of data equal to 1,000 points and plot it as a jitter plot:
#for j in vector:
#    new= np.copy(data[entry[j], :100, :])
#    #add a little bit of jitter to these points
#    cov= np.diag(0.005*np.ones(len(handle[0])-1))
#    xpos= np.arange(1, len(handle[0]))
#    for i in range(len(new)):
#        new[i, :]= np.random.multivariate_normal(new[i, :], cov)
#    
#    for i in range(len(xpos)):
#        x = np.random.normal(xpos[i], 0.07, size=len(new))
#        plt.plot(x, new[:, i], 'ko', alpha= 0.3)
#    plt.plot(xpos, handle[entry[j]][:-1], 'bo', ms=3)
#    plt.plot(xpos, handle[entry[j]][:-1], 'ro', ms=15, mfc='none', mec='red')
#    plt.xticks(np.arange(1, len(xpos)+1, 1.0), rotation= 'vertical')
#    ax = plt.gca()        
#    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(axisLabelerTissues))
#    plt.ylabel('Counts')
#    plt.xlim(0, 26)
#    plt.savefig('../output/'+j+'JitterPlot.pdf')
#    plt.close()
#
#
#
#
#
#
#
#
#    