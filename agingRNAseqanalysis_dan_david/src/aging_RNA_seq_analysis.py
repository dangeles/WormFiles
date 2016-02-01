# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 08:56:36 2016

@author: davidangeles
"""

# -*- coding: utf-8 -*-


import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import hypergeometricTests as hgt
import os
import mpl_toolkits.mplot3d

import pyrnaseq_graphics as rsq

# Package to perform PCA
import sklearn.datasets
import sklearn.decomposition

sns.set_context("notebook")
mag= 2 #value of beta from regression
qval= .1 #qvalue from regression
qvalEn= 0.05 #q value for enrichment analysis (tissues)
                
dirLists= '../output/Gene_lists_for_analysis'
if not os.path.exists(dirLists):
    os.makedirs(dirLists)
    
dirGraphs= '../output/Graphs'
if not os.path.exists(dirLists):
    os.makedirs(dirGraphs)
    

os.chdir('./')
#gene_lists from sleuth
dfBetaA= pd.read_csv("../input/table_agebeta_genes.csv")
dfBetaG= pd.read_csv("../input/table_genobeta_genes.csv")
dfBetaAG= pd.read_csv("../input/table_genocrossagebeta_genes.csv")

#gold standard datasets
dfDaf12= pd.read_csv('../input/daf12genes.csv')
dfDaf16= pd.read_csv('../input/daf16genes.csv')
dfLund= pd.read_csv('../input/lund_data.csv', header= None, names=['gene'])
dfEckley= pd.read_csv('../input/eckley_data.csv', header= None, names=['gene'])
dfMurphyUp= pd.read_csv('../input/murphy_data_lifespan_extension.csv')
dfMurphyDown= pd.read_csv('../input/murphy_data_lifespan_decrease.csv')
dfHalaschek= pd.read_csv('../input/Halaschek-Wiener_data.csv')

#gpcrs
dfGPCR= pd.read_csv('../input/putative_sensory_gpcrs.csv', header= None, names=['gene'])

#gpcr is going to go into a gold standard fxn so add an 'origin' colmn
dfGPCR['origin']= 'putative sensory gpcrs'

#place all the gold standards in a single dataframe:
dfDaf12['origin']= 'daf-12'
dfDaf16['origin']= 'daf-16'
dfEckley['origin']= 'Eckley'
dfLund['origin']= 'Lund'
dfMurphyUp['origin']= 'MurphyExt'
dfMurphyDown['origin']= 'MurphyDec'
dfHalaschek['origin']= 'Halaschek'
frames= [dfDaf12, dfDaf16, dfEckley, dfLund, dfMurphyDown, dfMurphyUp, dfHalaschek]
dfGoldStandard = pd.concat(frames)

#from wormbase
dfLifespanGenes= pd.read_csv('../input/lifespan gene list complete.csv')

#dfPAN= pd.read_csv()

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
                  namesAG0neg]    
#==============================================================================
# Tissue Enrichment Batch Analysis
#==============================================================================
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
#run the enrichment analysis  for each dataset
for i, list_of_genes in enumerate(array_of_arrays):
    #run the tissue enrichment tool 
    aname= array_of_anames[i] #name of the analysis
    fname= array_of_fnames[i] #filename to store as
    
    #print analysis name so you know what's gong on
    print(aname)
    print(len(list_of_genes))    
    df_analysis= \
    hgt.implement_hypergmt_enrichment_tool(aname, list_of_genes,\
                                    tissue_df, qvalEn, f_unused= fname)
    #save results to csv file
    df_analysis.to_csv('../output/EnrichmentAnalysisResults/'+fname, index= False)
    #reopen the file and add a comment with relevant info for the file
    line= '#' + aname+'\n'
    rsq.line_prepender('../output/EnrichmentAnalysisResults/'+fname, line)
    
    #plot top fold change tissues
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
#==============================================================================
#==============================================================================
# Tissue spread analyses
#==============================================================================
#color vector:
colors= ['#696969','#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00']
tissues= ['labial', 'extracellular', 'tail', 'dopaminergic' ]
#tissues= ['mu_int', 'sex organ', 'excretory', 'gonadal' ]
#tissues= ['gonad', 'sensillum', 'intestine', 'sex organ' ]
#tissues= ['head', 'tail', 'embryonic' ]
#tissues= ['muscle', 'coel', 'hyp' ]
#tissues= ['sperm', 'extracellular', 'tail', 'dopaminergic' ]
df_exp= rsq.organize(tissues, tissue_df)

rsq.explode(qval, dfBetaA, df_exp, colors= colors, title= 'Aging', 
        xlab= r'$\beta_{\mathrm{Aging}}$')
rsq.explode(qval, dfBetaG, df_exp, colors, title= 'Genotype', \
        xlab= r'$\beta_{\mathrm{Genotype}}$')
rsq.explode(qval, dfBetaAG, df_exp, colors, title= 'Aging::Genotype', \
        xlab= r'$\beta_{\mathrm{Aging::Genotype}}$')

colors2= ['#ffff33','#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00']
rsq.kegg(qval, dfBetaA, df_exp, colors2, title= 'Aging', \
        savename= '../output/Graphs/aging_tissue_specific_kde.png',\
        xlab= r'$\beta_{\mathrm{Aging}}$')
rsq.kegg(qval, dfBetaG, df_exp, colors2, title= 'Genotype', \
        savename= '../output/Graphs/genotype_tissue_specific_kde.png',\
        xlab= r'$\beta_{\mathrm{Genotype}}$')
rsq.kegg(qval, dfBetaAG, df_exp, colors2, title= 'Aging::Genotype', \
        savename= '../output/Graphs/agingxgenotype_tissue_specific_kde.png',\
        xlab= r'$\beta_{\mathrm{Aging::Genotype}}$')
#==============================================================================
#Where does gold standard data fall?   
#==============================================================================
#colors
colors= ['#999999','#ffff33','#f781bf','#ff7f00','#984ea3','#4daf4a', '#377eb8',
         '#e41a1c', '#a65628']

rsq.explode_goldstandards(qval, dfBetaA, dfGoldStandard, colors= colors, title= 'Aging', \
        savename= '../output/Graphs/aging_with_goldstandards_volcplot.png',\
        xlab= r'$\beta_{\mathrm{Aging}}$', xlim= [-12,12], ylim= [10**-4, 10**2])
rsq.explode_goldstandards(qval, dfBetaG, dfGoldStandard, colors= colors, title= 'Genotype', \
        savename= '../output/Graphs/genotype_with_goldstandards_volcplot.png',\
        xlab= r'$\beta_{\mathrm{Genotype}}$', xlim= [-12,12], ylim= [10**-4, 10**2])
rsq.explode_goldstandards(qval, dfBetaAG, dfGoldStandard, colors= colors, title= 'Aging::Genotype', \
        savename= '../output/Graphs/agingcrossgenotype_with_goldstandards_volcplot.png',\
        xlab= r'$\beta_{\mathrm{Aging::Genotype}}$', xlim= [-12,12], ylim= [10**-4, 10**2])
        
colors= ['#999999','#ffff33',
         '#e41a1c']
rsq.explode_goldstandards(qval, dfBetaA, dfGPCR, colors= colors, title= 'Aging', \
        savename= '../output/Graphs/aging_with_gpcrs_volcplot.png',\
        xlab= r'$\beta_{\mathrm{Aging}}$', xlim= [-12,12], ylim= [10**-4, 10**2])
rsq.explode_goldstandards(qval, dfBetaG, dfGPCR, colors= colors, title= 'Genotype', \
        savename= '../output/Graphs/genotype_with_gpcrs_volcplot.png',\
        xlab= r'$\beta_{\mathrm{Genotype}}$', xlim= [-12,12], ylim= [10**-4, 10**2])
rsq.explode_goldstandards(qval, dfBetaAG, dfGPCR, colors= colors, title= 'Aging::Genotype', \
        savename= '../output/Graphs/agingcrossgenotype_with_gpcrs_volcplot.png',\
        xlab= r'$\beta_{\mathrm{Aging::Genotype}}$', xlim= [-12,12], ylim= [10**-4, 10**2])



#figure out how many genes in dfLIfespan show up in this analysis
f= lambda x: (dfBetaA.ens_gene.isin(x)) & (dfBetaA.qval < .1)    
ind= f(dfLifespanGenes.gene.values)
m= dfBetaA[ind].ens_gene.values

f= lambda x: (dfBetaG.ens_gene.isin(x)) & (dfBetaG.qval < .1)    
ind= f(dfLifespanGenes.gene.values)
m= np.append(m, dfBetaG[ind].ens_gene.values)

f= lambda x: (dfBetaAG.ens_gene.isin(x)) & (dfBetaAG.qval < .1)    
ind= f(dfLifespanGenes.gene.values)
m= np.append(m, dfBetaAG[ind].ens_gene.values)

m= list(set(m))    
with open('../output/lifespan_genes_that_show_up.csv', 'w') as f:
    f.write('WBID\n')
    for gene in m:
        f.write(gene)
        f.write('\n')
    f.close()


Ldf= [dfBetaA, dfBetaG, dfBetaAG] #list of dataframes
dfnames= ['Age', 'Genotype', 'Age::Genotype']
colors= ['#377eb8','#e41a1c','#4daf4a']
fnames= ['../output/Graphs/positive_aging.png','../output/Graphs/negative_aging.png',
         '../output/Graphs/variable_aging.png','../output/Graphs/unannotated_aging.png']

rsq.kegg_compareall_byval(qval, Ldf, dfLifespanGenes, colors, savenames= fnames,
                   dfnames= dfnames, xscale= 'symlog', save= True)



