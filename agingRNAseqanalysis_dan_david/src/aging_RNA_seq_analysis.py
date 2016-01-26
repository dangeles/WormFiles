# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 08:56:36 2016

@author: davidangeles
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

Investigation of the new dictionary
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import hypergeometricTests as hgt
import os


mag= 2 #value of beta from regression
qval= .1 #qvalue from regression
qvalEn= 0.05 #q value for enrichment analysis (tissues)


os.chdir('./')
#gene_lists from sleuth
dfBetaA= pd.read_csv("../input/table_agebeta_genes.csv")
dfBetaG= pd.read_csv("../input/table_genobeta_genes.csv")
dfBetaAG= pd.read_csv("../input/table_genocrossagebeta_genes.csv")

#tissue dictionary-- please cite David Angeles et al TEA publication (forthcoming)
#if using the enrichment tool 
tissue_df= pd.read_csv("../input/dictionary.csv")

#slice all the relevant gene names out
#remove all isoforms!
namesBetaA= \
        dfBetaA[\
        (dfBetaA.qval < qval) & ((dfBetaA.b > mag))]\
        .ens_gene.unique()
namesBetaG= \
        dfBetaG[\
        (dfBetaG.qval < qval) & ((dfBetaG.b > mag))]\
        .ens_gene.unique()
namesBetaAG= \
        dfBetaAG[\
        (dfBetaAG.qval < qval) & ((dfBetaAG.b > mag))]\
        .ens_gene.unique()
        
namesBetaAneg= \
        dfBetaA[\
        (dfBetaA.qval < qval) & ((dfBetaA.b < -mag))]\
        .ens_gene.unique()
namesBetaGneg= \
        dfBetaG[\
        (dfBetaG.qval < qval) & ((dfBetaG.b < -mag))]\
        .ens_gene.unique()
namesBetaAGneg= \
        dfBetaAG[\
        (dfBetaAG.qval < qval) & ((dfBetaAG.b < -mag))]\
        .ens_gene.unique()
        

namesA0= \
        dfBetaA[\
        (dfBetaA.qval < qval) & (dfBetaAG.b > 0)]\
        .ens_gene.unique()
namesG0= \
        dfBetaG[\
        (dfBetaG.qval < qval)& (dfBetaAG.b > 0)]\
        .ens_gene.unique()
namesAG0= \
        dfBetaAG[\
        (dfBetaAG.qval < qval)& (dfBetaAG.b > 0)]\
        .ens_gene.unique()
        
        
namesA0neg= \
        dfBetaA[\
        (dfBetaA.qval < qval) & (dfBetaAG.b < 0)]\
        .ens_gene.unique()
namesG0neg= \
        dfBetaG[\
        (dfBetaG.qval < qval)& (dfBetaAG.b < 0)]\
        .ens_gene.unique()
namesAG0neg= \
        dfBetaAG[\
        (dfBetaAG.qval < qval)& (dfBetaAG.b < 0)]\
        .ens_gene.unique()

#place all the relevant gene lists in this array
array_of_arrays= [namesBetaA,
                  namesBetaG,
                  namesBetaAG,
                  namesBetaAneg,
                  namesBetaGneg,
                  namesBetaAGneg,
                  namesA0,
                  namesG0,
                  namesAG0,
                  namesA0neg,
                  namesG0neg,
                  namesAG0neg
                  ]

#Run the whole thing:


aname0= 'Age Beta> {0}'.format(mag)
aname1= 'Genotype Beta> {0}'.format(mag)
aname2= 'Age::Genotype Beta> {0}'.format(mag)

aname3= 'Age Beta< -{0}'.format(mag)
aname4= 'Genotype Beta< -{0}'.format(mag)
aname5= 'Age::Genotype Beta< -{0}'.format(mag)

aname6= 'Age Beta> 0'
aname7= 'Genotype Beta> 0'
aname8= 'Age::Genotype Beta> 0'

aname9= 'Age Beta< 0'
aname10= 'Genotype Beta< 0'
aname11= 'Age::Genotype Beta< 0'

array_of_anames= [aname0,aname1,aname2,
                  aname3,aname4,aname5,
                  aname6,aname7,aname8,
                  aname9,aname10,aname11]

fname0= 'EnrichmentAnalysisAge_BetaGreaterThan{0}_qval_{1}.csv'.format(mag, qvalEn)
fname1= 'EnrichmentAnalysisGenotype_BetaGreaterThan{0}_qval_{1}.csv'.format(mag, qvalEn)
fname2= 'EnrichmentAnalysisAgeCrossGenotype_BetaGreaterThan{0}_qval_{1}.csv'.format(mag, qvalEn)

fname3= 'EnrichmentAnalysisAgePositive_BetaLessThanNegative{0}_qval_{1}.csv'.format(mag, qvalEn)
fname4= 'EnrichmentAnalysisGenotype_BetaLessThanNegative{0}_qval_{1}.csv'.format(mag, qvalEn)
fname5= 'EnrichmentAnalysisAgeCrossGenotype_BetaLessThanNegative{0}_qval_{1}.csv'.format(mag, qvalEn)

fname6= 'EnrichmentAnalysisAge_AllPosBeta_qval_{0}.csv'.format(qvalEn)
fname7= 'EnrichmentAnalysisGenotype_AllPosBeta_qval_{0}.csv'.format(qvalEn)
fname8= 'EnrichmentAnalysisAgeCrossGenotype_AllPosBeta_qval_{0}.csv'.format(qvalEn)

fname9= 'EnrichmentAnalysisAge_AllNegBeta_qval_{0}.csv'.format(qvalEn)
fname10= 'EnrichmentAnalysisGenotype_AllNegBeta_qval_{0}.csv'.format(qvalEn)
fname11= 'EnrichmentAnalysisAgeCrossGenotype_AllNegBeta_qval_{0}.csv'.format(qvalEn)

array_of_fnames= [fname0, fname1, fname2,
                  fname3, fname4, fname5, 
                  fname6, fname7, fname8, 
                  fname9, fname10, fname11]
                  
array_of_strings= ['namesBetaA', 'namesBetaG', 'namesBetaAG',
                'namesBetaAneg', 'namesBetaGneg', 'namesBetaAGneg',
                'namesA0', 'namesG0', 'namesAG0',
                'namesA0neg', 'namesG0neg', 'namesAG0neg'
                ]
                
dirLists= '../output/Gene_lists_for_analysis'
if not os.path.exists(dirLists):
    os.makedirs(dirLists)
    
dirGraphs= '../output/Graphs'
if not os.path.exists(dirLists):
    os.makedirs(dirGraphs)
    

n_bars= 15 #number of results to plot bars for
#run the whole thing for each dataset
for i, list_of_genes in enumerate(array_of_arrays):        
    
    #run the tissue enrichment tool 
    aname= array_of_anames[i] #name of the analysis
    fname= array_of_fnames[i] #filename to store as
    
    print(aname)
    print(len(list_of_genes))    
    
    df_analysis= \
    hgt.implement_hypergmt_enrichment_tool(aname, list_of_genes,\
                                    tissue_df, qvalEn, f_unused= fname)
    
    with open('../output/EnrichmentAnalysisResults/'+fname, 'w') as f:
        f.write(aname)
        df_analysis.to_csv('../output/EnrichmentAnalysisResults/'+fname)
    
    hgt.plotting_and_formatting(df_analysis, ytitle= aname)
    #close
#    plt.close()
    
#    #save genes with betas significantly different from zero in a list format
#    k= array_of_strings[i][5:]
#            
#    filename= dirLists+'/'+ 'gene_list_{0}_beta{1}_q{2}.csv'.format(k, mag, qvalEn)
#    #print('filename', filename)
#    with open(filename, 'w') as f:
#        for gene in list_of_genes:
#            f.write(gene)
#            f.write('\n')
#    f.close()
    






