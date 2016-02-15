# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 14:07:04 2016
A script to analyze the genes identified by Engelmann et al 2011,
A Comprehensive Analysis of Gene Expression Changes Provoked by Bacterial
and Fungal Infection in C elegans.

Note: I've used only the RNA-seq data they generated
@author: dangeles
"""

import hypergeometricTests as hgt #the library to be used
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


tissue_df= pd.read_csv('../input/WS252AnatomyDictionary/25.csv')

#Bacterial
dfLum= pd.read_csv('../input/luminescens_Engelmann_2011.csv')
dfMarc= pd.read_csv('../input/marcescens_Engelmann_2011.csv')
dfFaec= pd.read_csv('../input/faecalis_Engelmann_2011.csv')

#Fungal
dfHarp= pd.read_csv('../input/harposporium_Engelmann_2011.csv')
dfCon= pd.read_csv('../input/coniospora_Engelmann_2011.csv')

#dictionary of gene names
names= pd.read_csv('../input/c_elegans.PRJNA13758.WS241.livegeneIDs.unmaprm.txt', sep= '\t',comment= '#')

#make the directories to place the analyses
dirEngelmann= '../output/Engelmann/'
dirAnalysis= '../output/Engelmann/Analysis'
dirGraphs= '../output/Engelmann/Graphs'

DIRS= [dirEngelmann, dirAnalysis, dirGraphs]
#Make all the necessary dirs if they don't already exist
for d in DIRS:
    if not os.path.exists(d):
        os.makedirs(d)


f= lambda x: (names.HumanReadable.isin(x.SequenceNameGene))
g= lambda x: (names.GeneName.isin(x.SequenceNameGene))

Lnames= ['Otorhabdus luminescens', 'Serratia marcescens', 'Enterococcus faecalis',
         'Harsposporium sp', 'Drechmeria coniospora']
Ldf= [dfLum, dfMarc, dfFaec, dfHarp, dfCon]
Ldirection= ['Infection_upregulated', 'Infection_downregulated']

#direction specific
for i, df in enumerate(Ldf):
    fname= Lnames[i]
    for direction in Ldirection:
        ind= g(df[df[direction] == 1.0])
        ind2= f(df[df[direction]== 1.0])
        x= names[ind | ind2].WBID
        print(fname + direction)
        print(len(x))
        hgt.implement_hypergmt_enrichment_tool(fname + direction,x, tissue_df)
        

#direction agnostic
for i, df in enumerate(Ldf):
    fname= Lnames[i]
    ind1= g(df[ df[Ldirection[0]] == 1.0] )
    ind2= f(df[ df[Ldirection[0]]== 1.0] )
    ind3= g(df[ df[Ldirection[1]] == 1.0] )
    ind4= f(df[ df[Ldirection[1]]== 1.0] )
    x= names[(ind1 | ind2) | (ind3 | ind4)].WBID
    print(fname)
    print(len(x))
    df_res= hgt.implement_hypergmt_enrichment_tool(fname,x, tissue_df)
    hgt.plotting_and_formatting(df_res, ytitle= '{0}'.format(fname), dirGraphs= dirGraphs)
    
    
    
    
    
    
    
    
    
    
    
    
    
    