# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 14:07:04 2016
A script to analyze some of the genes identified by Engelmann et al 2011,
A Comprehensive Analysis of Gene Expression Changes Provoked by Bacterial
and Fungal Infection in C elegans.

@author: dangeles
"""
import hypergeometricTests as hgt #the library to be used
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os

#tissue dictionary to use
tissue_df= pd.read_csv('../input/WS252AnatomyDictionary/annot50.thresh0.9.methodavg.csv')
if 'C. elegans Cell and Anatomy WBbt:0000100' in tissue_df.columns:
    tissue_df.drop('C. elegans Cell and Anatomy WBbt:0000100', axis= 1, inplace= True)

#rna-seq datasets to reanalyze
#Bacterial
dfLum= pd.read_csv('../input/Engelmann/luminescens_Engelmann_2011.csv')
dfMarc= pd.read_csv('../input/Engelmann/marcescens_Engelmann_2011.csv')
dfFaec= pd.read_csv('../input/Engelmann/faecalis_Engelmann_2011.csv')

#Fungal
dfHarp= pd.read_csv('../input/Engelmann/harposporium_Engelmann_2011.csv')
dfCon= pd.read_csv('../input/Engelmann/coniospora_Engelmann_2011.csv')


#dfSin= pd.read_csv('../input/Engelmann/sinha_2012.csv')

#dictionary of gene names
names= pd.read_csv('../input/Engelmann/c_elegans.PRJNA13758.WS241.livegeneIDs.unmaprm.txt', 
                   sep= '\t',comment= '#')

#make the directories to place the analyses
dirEngelmann= '../output/Engelmann/'
dirAnalysis= '../output/Engelmann/Analysis'
dirGraphs= '../output/Engelmann/Graphs'
DIRS= [dirEngelmann, dirAnalysis, dirGraphs]

#Make all the necessary dirs if they don't already exist
for d in DIRS:
    if not os.path.exists(d):
        os.makedirs(d)


#Selector functions to draw WBIDs from 
f= lambda x, y: (names.HumanReadable.isin(x[y]))
g= lambda x, y: (names.GeneName.isin(x[y]))

#==============================================================================
#==============================================================================
# # Engelmann Analysis
#==============================================================================
#==============================================================================
Lnames= ['Otorhabdus luminescens', 'Serratia marcescens', 'Enterococcus faecalis',
         'Harsposporium sp', 'Drechmeria coniospora']
Ldf= [dfLum, dfMarc, dfFaec, dfHarp, dfCon]
Ldirection= ['Infection_upregulated', 'Infection_downregulated']

#direction specific
n_genes= []
for i, df in enumerate(Ldf):
    fname= Lnames[i]
    for direction in Ldirection:
        ind= g(df[df[direction] == 1.0], 'SequenceNameGene')
        x= names[ind].WBID
        
        print('---------')
        print(fname + ' ' + direction)
        print('Number of genes submitted for analysis ', len(x))
        y= tissue_df[tissue_df.wbid.isin(x)].wbid.unique().shape[0]
        print('Number of genes used for analysis ', y)
        print('\n')
        
        df_res, unused= hgt.enrichment_analysis(x, tissue_df, show= False)
        hgt.plot_enrichment_results(df_res, title= '{0}, {1}'.format(fname, direction), 
                                    dirGraphs= dirGraphs)
        n_genes.append(['{0}, {1}'.format(fname, direction), y, len(unused)])

print('---------')
print('---------')
print('---------')
print('---------')

#direction agnostic
for i, df in enumerate(Ldf):
    fname= Lnames[i]
    ind1= g(df[ df[Ldirection[0]] == 1.0],  'SequenceNameGene')
    ind3= g(df[ df[Ldirection[1]] == 1.0],'SequenceNameGene' )
    x= names[(ind1) | (ind3)].WBID
    
    print(fname)
    print('Number of genes submitted for analysis ', len(x))
    y= tissue_df[tissue_df.wbid.isin(x)].wbid.unique().shape[0]
    print('Number of genes used for analysis ', y)
    print('\n')
    
    df_res, unused= hgt.enrichment_analysis(x, tissue_df, show= False)
    hgt.plot_enrichment_results(df_res, title= '{0}'.format(fname), dirGraphs= dirGraphs)
##==============================================================================
##==============================================================================
## # Sinha, 2012 Analysis
##==============================================================================
##==============================================================================
#
#sin_ps= ['pval_Bthu', 'pval_Saur', 'pval_Smar', 'pval_Xnem']
#sin_fs= ['log2FC_Bthu','log2FC_Saur','log2FC_Smar','log2FC_Xnem']
#
#print('---------')
#print('---------')
#print('---------')
#print('---------')
#print('Sinha Begins Here')
#print('---------')
#print('---------')
#print('---------')
#print('---------')
#
##direction agnostic
#for pval in sin_ps:
#    
#    fname= pval[5:]    
#    ind= g(dfSin[dfSin[pval] < 0.05], 'Gene ID')
#    ind2= f(dfSin[dfSin[pval] < 0.05], 'Gene ID')
#    x= names[ind | ind2].WBID
#    
#    print(fname)
#    print('Number of genes submitted for analysis ', len(x))
#    y= tissue_df[tissue_df.wbid.isin(x)].wbid.unique().shape[0]
#    print('Number of genes used for analysis ', y)
#    print('\n')
#    df_res, unused= hgt.enrichment_analysis(x, tissue_df, show= False)
#    hgt.plot_enrichment_results(df_res, 
#        title= 'Sinha 2102, {0}'.format(fname), dirGraphs= dirGraphs)    
#
#print('---------')
#print('---------')
#print('---------')
#print('---------')
#print('---------')
#
##direction specific
#for i, pval in enumerate(sin_ps):
#    for i in np.arange(-1, 2 ,2):
#        if i < 0:
#            fname= pval[5:] + ' DownRegulated'
#        else:
#            fname= pval[5:] + ' UpRegulated'
#            
#            
#        sig = (dfSin[pval] < 0.05)
#        fval= (i*dfSin[sin_fs[i]] > 0)
#        ind= g(dfSin[sig  & fval], 'Gene ID')
#        ind2= f(dfSin[dfSin[pval] < 0.05], 'Gene ID')
#        x= names[ind | ind2].WBID
#        
#        print(fname)
#        print('Number of genes submitted for analysis ', len(x))
#        y= tissue_df[tissue_df.wbid.isin(x)].wbid.unique().shape[0]
#        print('Number of genes used for analysis ', y)
#        print('\n')
#        
#        df_res, unused= hgt.enrichment_analysis(x, tissue_df, show= False)
#        hgt.plot_enrichment_results(df_res, 
#            title= 'Sinha, 2012 {0}'.format(fname), dirGraphs= dirGraphs)    

 



























