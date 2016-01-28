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


mag= 2 #value of beta from regression
qval= .1 #qvalue from regression
qvalEn= 0.05 #q value for enrichment analysis (tissues)


os.chdir('./')
#gene_lists from sleuth
dfBetaA= pd.read_csv("../input/table_agebeta_genes.csv")
dfBetaG= pd.read_csv("../input/table_genobeta_genes.csv")
dfBetaAG= pd.read_csv("../input/table_genocrossagebeta_genes.csv")

dfDaf12= pd.read_csv('../input/daf12genes.csv')
dfDaf16= pd.read_csv('../input/daf16genes.csv')

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
    


#Daf-12 associated genes
ndaf12= dfDaf12.shape[0]
ndaf16= dfDaf16.shape[0]
print('no. of daf-12 genes in list: {0}'.format(ndaf12))
print('no. of daf-16 genes in list: {0}'.format(ndaf16))
print('\n')

nA= dfBetaA[dfBetaA.qval < .1].shape[0]
nG= dfBetaG[dfBetaG.qval < .1].shape[0]
nAG= dfBetaAG[dfBetaAG.qval < .1].shape[0]

print('no. of sig. genes in dfA: {0}'.format(nA))
print('no. of sig. genes in dfG: {0}'.format(nG))
print('no. of sig. genes in dfAG: {0}'.format(nAG))

indA12= (dfBetaA.ens_gene.isin(dfDaf12.gene)) & (dfBetaA.qval < .1)
indA16= (dfBetaA.ens_gene.isin(dfDaf16.gene)) & (dfBetaA.qval < .1)
ndaf12A= dfBetaA.b[indA12].shape[0]
ndaf16A= dfBetaA.b[indA16].shape[0]

print('no. of daf-12 genes with qval < {0:.2} in aging'.format(qval))
print(ndaf12A)
print('frac of genes {0:.2}'.format(ndaf12A/ndaf12))
print('no. of daf-16 genes with qval < {0:.2} in aging'.format(qval))
print(ndaf16A)
print('frac of genes {0:.2}'.format(ndaf16A/ndaf16))

dfBetaA.b[indA12].plot('kde', lw= 4, color= 'b')
dfBetaA.b[indA16].plot('kde', lw= 4)
plt.xlim(-5, 5)


indG12= (dfBetaG.ens_gene.isin(dfDaf12.gene)) & (dfBetaG.qval < .1)
indG16= (dfBetaG.ens_gene.isin(dfDaf16.gene)) & (dfBetaG.qval < .1)
ndaf12G= dfBetaG.b[indG12].shape[0]
ndaf16G= dfBetaG.b[indG16].shape[0]

print('no. of daf-12 genes with qval < {0:.2} in genotype'.format(qval))
print(ndaf12G)
print('frac of genes in dbase {0:.2}'.format(ndaf12G/nG))
print('no. of daf-16 genes with qval < {0:.2} in genotype'.format(qval))
print(ndaf16G)
print('frac of genes in dbase {0:.2}'.format(ndaf16G/nG))

dfBetaG.b[indG12].plot('kde', lw= 4, color= 'b')
dfBetaG.b[indG16].plot('kde', lw= 4)
plt.xlim(-5,5)


indAG12= (dfBetaAG.ens_gene.isin(dfDaf12.gene)) & (dfBetaAG.qval < .1)
indAG16= (dfBetaAG.ens_gene.isin(dfDaf16.gene)) & (dfBetaAG.qval < .1)


ndaf12AG= dfBetaAG.b[indAG12].shape[0]
ndaf16AG= dfBetaAG.b[indAG16].shape[0]

print('no. of daf-12 genes with qval < {0:.2} in aging::genotype'.format(qval))
print(ndaf12AG)
print('frac of genes in dbase {0:.2}'.format(ndaf12AG/nAG))
print('no. of daf-16 genes with qval < {0:.2} in aging::genotype'.format(qval))
print(ndaf16AG)
print('frac of genes in dbase {0:.2}'.format(ndaf16AG/nAG))

dfBetaAG.b[indAG12].plot('kde')
dfBetaAG.b[indAG16].plot('kde')
plt.xlim(-5, 5)


def wbid_extractor(tissue_df, main_tissue):
    """
    Given a string 'main_tissue', find all columns that
    have a substring equal to it in tissue_df. Then,
    extract all the wbids that are expressed in any of these
    columns and return a non-redundant list of these wbids
    """
    if type(main_tissue) != str:
        raise ValueError('please input a string in main tissue')
        
    matching = [s for s in tissue_df.columns if main_tissue in s]
    names= []
    for i in matching:
        if len(names) == 0:
            names= tissue_df.wbid[tissue_df[i]==1].values
        else:
            names= np.append(names, tissue_df.wbid[tissue_df[i]==1].values)
    names= list(set(names))
    
    return names
    
neurons= wbid_extractor(tissue_df, 'neuron')
intestine= wbid_extractor(tissue_df, 'intestine')
germline= wbid_extractor(tissue_df, 'germ')
muscle= wbid_extractor(tissue_df, 'muscle')
hyp= wbid_extractor(tissue_df, 'hyp')
sperm= wbid_extractor(tissue_df, 'sperm')
male= wbid_extractor(tissue_df, 'male')

def organize(names, tissue_df):
    
    index= [None]*len(names)
    for i, value in names:
        
    df_selected_tissues= df.DataFrame(index= len(tissue_df), columns= columns)
    df_selected_tissues.wbid= tissue_df.wbid
    df_selected_tissues[columns[i]][df_selected_tissues.wbid == index[i]]= 1
    df.fillna(0, inplace= True)

all_tiss= list(set(np.append(neurons, intestine)))
all_tiss= list(set(np.append(all_tiss, germline)))
all_tiss= list(set(np.append(all_tiss, muscle)))
all_tiss= list(set(np.append(all_tiss, hyp)))
all_tiss= list(set(np.append(all_tiss, sperm)))
all_tiss= list(set(np.append(all_tiss, male)))


xnotsig= dfBetaA[(dfBetaA.qval > .1)].b
ynotsig= dfBetaA[(dfBetaA.qval > .1)].qval
xnotisssig= dfBetaA[(~dfBetaA.ens_gene.isin(all_tiss)) & (dfBetaA.qval < .1)].b
ynotisssig= dfBetaA[(~dfBetaA.ens_gene.isin(all_tiss)) & (dfBetaA.qval < .1)].qval

def volcano_plot_tissue(tissue, q, df, ax, label, col= 'b'):
    """
    """
    f= lambda x: (df.ens_gene.isin(x)) & (df.qval < q)
    ind= f(tissue) 
    x= df[ind].b
    y= df[ind].qval
    plt.gca().plot(x, -np.log(y), 'o', color= col, ms= 4.5, alpha= .6, label= label)


#color vector:
colors= ['#999999','#a65628','#f781bf','#ffff33','#ff7f00','#984ea3','#4daf4a','#377eb8', '#e41a1c']

fig, ax= plt.subplots()
plt.plot(xnotsig, -np.log(ynotsig), 'o', color=colors[0], ms= 4.5, alpha= .6, label= 'not sig')
plt.plot(xnotisssig, -np.log(ynotisssig), 'o', color=colors[1], ms=4.5, alpha= .6, label= 'sig, no tissue')

#plot all the points not associated with a tissue
volcano_plot_tissue(neurons, .1, dfBetaA, label= 'neurons', col= colors[2], ax= ax)
volcano_plot_tissue(intestine, .1, dfBetaA, label= 'intestine', col= colors[3], ax= ax)
volcano_plot_tissue(germline, .1, dfBetaA, label= 'germline', col= colors[4], ax= ax)
volcano_plot_tissue(muscle, .1, dfBetaA, label= 'muscle', col= colors[5], ax= ax)
volcano_plot_tissue(hyp, .1, dfBetaA, label= 'hyp', col= colors[6], ax= ax)
volcano_plot_tissue(sperm, .1, dfBetaA, label= 'sperm', col= colors[7], ax= ax)
volcano_plot_tissue(male, .1, dfBetaA, label= 'male', col= colors[8], ax= ax)
plt.legend()













    
    