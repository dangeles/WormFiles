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

#sort by target_id
dfBetaA.sort_values('target_id', inplace= True)
dfBetaG.sort_values('target_id', inplace= True)
dfBetaAG.sort_values('target_id', inplace= True)

#gold standard datasets
dfDaf12= pd.read_csv('../input/daf12genes.csv')
dfDaf16= pd.read_csv('../input/daf16genes.csv')
dfLund= pd.read_csv('../input/lund_data.csv', header= None, names=['gene'])
dfEckley= pd.read_csv('../input/eckley_data.csv', header= None, names=['gene'])
dfMurphyUp= pd.read_csv('../input/murphy_data_lifespan_extension.csv')
dfMurphyDown= pd.read_csv('../input/murphy_data_lifespan_decrease.csv')
dfHalaschek= pd.read_csv('../input/Halaschek-Wiener_data.csv')

#gpcrs
dfGPCR= pd.read_csv('../input/all_gpcrs.csv')
dfICh= pd.read_csv('../input/select_ion_transport_genes.csv')
dfAxon= pd.read_csv('../input/axonogenesis_genes.csv')
dfNP= pd.read_csv('../input/neuropeptides.csv')

#gpcr is going to go into a gold standard fxn so add an 'origin' colmn
dfGPCR['origin']= 'gpcrs'
dfICh['origin']= 'select ion transport genes'
dfAxon['origin']= 'axonogenesis genes'
dfNP['origin']= 'neuropeptide genes'
frames= [dfGPCR, dfICh, dfAxon, dfNP]
dfTargets= pd.concat(frames)

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
# KDE for select tissues
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


colors2= ['#ffff33','#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00']
rsq.kegg(qval, dfBetaA, df_exp, colors2, title= 'Aging', \
        savename= '../output/Graphs/aging_tissue_specific_kde.png',\
        xlab= r'$\beta_{\mathrm{Aging}}$', ylab= 'Density')
rsq.kegg(qval, dfBetaG, df_exp, colors2, title= 'Genotype', \
        savename= '../output/Graphs/genotype_tissue_specific_kde.png',\
        xlab= r'$\beta_{\mathrm{Genotype}}$', ylab= 'Density')
rsq.kegg(qval, dfBetaAG, df_exp, colors2, title= 'Aging::Genotype', \
        savename= '../output/Graphs/agingxgenotype_tissue_specific_kde.png',\
        xlab= r'$\beta_{\mathrm{Aging::Genotype}}$', ylab= 'Density')
#==============================================================================
#Volcano plots for gold standards
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

#==============================================================================
# Volcano plots for target genes
#==============================================================================
colors= ['#999999','#ffff33','#f781bf','#ff7f00','#984ea3','#4daf4a', '#377eb8',
         '#e41a1c', '#a65628']
rsq.explode_goldstandards(qval, dfBetaA, dfTargets, colors= colors, title= 'Aging', \
        savename= '../output/Graphs/aging_with_targets_volcplot.png',\
        xlab= r'$\beta_{\mathrm{Aging}}$', xlim= [-12,12], ylim= [10**-4, 10**2])
rsq.explode_goldstandards(qval, dfBetaG, dfTargets, colors= colors, title= 'Genotype', \
        savename= '../output/Graphs/genotype_with_targets_volcplot.png',\
        xlab= r'$\beta_{\mathrm{Genotype}}$', xlim= [-12,12], ylim= [10**-4, 10**2])
rsq.explode_goldstandards(qval, dfBetaAG, dfTargets, colors= colors, title= 'Aging::Genotype', \
        savename= '../output/Graphs/agingcrossgenotype_with_targets_volcplot.png',\
        xlab= r'$\beta_{\mathrm{Aging::Genotype}}$', xlim= [-12,12], ylim= [10**-4, 10**2])


dfLifespanGenes.columns= ['gene', 'origin']
rsq.explode_goldstandards(qval, dfBetaA, dfLifespanGenes, colors= colors, title= 'Aging', \
        savename= '../output/Graphs/aging_lifespan_volcplot.png',\
        xlab= r'$\beta_{\mathrm{Aging}}$', xlim= [-12,12], ylim= [10**-4, 10**2])

rsq.explode_goldstandards(qval, dfBetaG, dfLifespanGenes, colors= colors, title= 'Genotype', \
        savename= '../output/Graphs/genotype_lifespan_volcplot.png',\
        xlab= r'$\beta_{\mathrm{Genotype}}$', xlim= [-12,12], ylim= [10**-4, 10**2])

rsq.explode_goldstandards(qval, dfBetaAG, dfLifespanGenes, colors= colors, title= 'Aging::Genotype', \
        savename= '../output/Graphs/agingcrossgenotype_lifespan_volcplot.png',\
        xlab= r'$\beta_{\mathrm{Aging::Genotype}}$', xlim= [-12,12], ylim= [10**-4, 10**2])

#==============================================================================
# Count how many lifespan genes show up
#==============================================================================
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
#==============================================================================
# KDE for Lifespan Genes from Wormbase
#==============================================================================
rsq.kegg_compareall_byval(qval, Ldf, dfLifespanGenes, colors, savenames= fnames,
                   dfnames= dfnames, xscale= 'symlog', ylab= 'Density', save= True)
                   
#==============================================================================
# KDE for Gene Targets
#==============================================================================
fnames= ['../output/Graphs/gpcrs.png','../output/Graphs/ion_transporters.png',
         '../output/Graphs/axonogenesis_genes.png', '../output/Graphs/neuropeptide_genes.png']
colors= ['#984ea3','#4daf4a', '#377eb8',
         '#e41a1c', '#a65628']
dfTargets.columns= ['gene', 'effect']
rsq.kegg_compareall_byval(qval, Ldf, dfTargets, colors, savenames= fnames,
                   dfnames= dfnames, xscale= 'symlog', ylab= 'Density', save= True)

f= lambda x, y: (x.qval < qval) & (x.ens_gene.isin(dfTargets[dfTargets.effect== y].gene))           

gpcrsA= dfBetaA[f(dfBetaA, 'gpcrs')]
gpcrsG= dfBetaG[f(dfBetaG, 'gpcrs')]
gpcrsAG= dfBetaAG[f(dfBetaAG, 'gpcrs')]

ionA= dfBetaA[f(dfBetaA, 'select ion transport genes')]
ionG= dfBetaG[f(dfBetaG, 'select ion transport genes')]
ionAG= dfBetaAG[f(dfBetaAG, 'select ion transport genes')] #ten of these genes are in axonG as well

axonA= dfBetaA[f(dfBetaA, 'axonogenesis genes')]
axonG= dfBetaG[f(dfBetaG, 'axonogenesis genes')]
axonAG= dfBetaAG[f(dfBetaAG, 'axonogenesis genes')]

gpcrsA[['ens_gene', 'b']].to_csv('../output/gpcrsAging.csv', header= ['GPCRS in Aging', 'reg. value'])
gpcrsG[['ens_gene', 'b']].to_csv('../output/gpcrsGenotype.csv', header= ['GPCRS in Genotype', 'reg. value'])
gpcrsAG[['ens_gene', 'b']].to_csv('../output/gpcrsAgingCrossGenotype.csv', 
                        header= ['GPCRS in Aging::Genotype', 'reg. value'])

ionA[['ens_gene', 'b']].to_csv('../output/ion_transport_genesAging.csv', 
                     header= ['ion transport genes in Aging', 'reg. value'])
ionG[['ens_gene', 'b']].to_csv('../output/ion_transport_genesGenotype.csv', 
                     header= ['ion transport genes in Genotype', 'reg. value'])
ionAG[['ens_gene', 'b']].to_csv('../output/ion_transport_genesAgingCrossGenotype.csv', 
                      header= ['ion transport genes in Aging::Genotype', 'reg. value'])

axonA[['ens_gene', 'b']].to_csv('../output/axonogenesisAging.csv', 
                      header= ['axonogenesis genes in Aging', 'reg. value'])
axonG[['ens_gene', 'b']].to_csv('../output/axonogenesisGenotype.csv', 
                      header= ['axonogenesis genes in Genotype', 'reg. value'])
axonAG[['ens_gene', 'b']].to_csv('../output/axonogenesisAgingCrossGenotype.csv', 
                       header= ['axonogenesis genes in Aging::Genotype', 'reg. value'])



#==============================================================================
# Identify candidates for RNAi
# Good targets: small qval, large positive b vals, not previously described
#==============================================================================
#aging 
dfBetaA.sort_values('b', inplace= True)
ind= (~dfBetaA[dfBetaA.qval < qval].ens_gene.isin(dfGoldStandard.gene)) & \
     (~dfBetaA[dfBetaA.qval < qval].ens_gene.isin(dfLifespanGenes.gene))
     
  
dfBetaA[dfBetaA.qval < qval][ind].tail(55).to_csv('../output/CandidatesAging.csv')


#genotype 
dfBetaG.sort_values('b', inplace= True)
ind= (~dfBetaG[dfBetaG.qval < qval].ens_gene.isin(dfGoldStandard.gene)) & \
     (~dfBetaG[dfBetaG.qval < qval].ens_gene.isin(dfLifespanGenes.gene))
     
  
dfBetaG[dfBetaG.qval < qval][ind].tail(25).to_csv('../output/CandidatesGenotype.csv')
     
#interactions    
dfBetaAG.sort_values('b', inplace= True)
ind= (~dfBetaAG[dfBetaAG.qval < qval].ens_gene.isin(dfGoldStandard.gene)) & \
     (~dfBetaAG[dfBetaAG.qval < qval].ens_gene.isin(dfLifespanGenes.gene)) &\
     (dfBetaA.ix[dfBetaAG[dfBetaAG.qval < qval].index].b > dfBetaA.b.quantile(.75))
     
  
dfBetaAG[dfBetaAG.qval < qval][ind].tail(55).to_csv('../output/CandidatesAgingXGenotype.csv')

#==============================================================================
# Analysis of discordant betas -- i.e., which genes have interaction effects
# That are suggestive of lifespan importance
#==============================================================================
#remove from each dataframe all the genes that aren't in all others
dfBetaA= dfBetaA[dfBetaA.target_id.isin(dfBetaG.target_id)]
dfBetaA= dfBetaA[dfBetaA.target_id.isin(dfBetaAG.target_id)]
dfBetaG= dfBetaG[dfBetaG.target_id.isin(dfBetaA.target_id)]
dfBetaG= dfBetaG[dfBetaG.target_id.isin(dfBetaAG.target_id)]
dfBetaAG= dfBetaAG[dfBetaAG.target_id.isin(dfBetaG.target_id)]
dfBetaAG= dfBetaAG[dfBetaAG.target_id.isin(dfBetaA.target_id)]

dfBetaA.set_index('target_id', inplace= True)
dfBetaG.set_index('target_id', inplace= True)
dfBetaAG.set_index('target_id', inplace= True)

ind= (dfBetaA.qval < qval) & (dfBetaG.qval < qval) & (dfBetaAG.qval < qval)

plt.plot(dfBetaA.b, dfBetaG.b, 'o', alpha= .3)
plt.plot(np.absolute(dfBetaA[ind].b+dfBetaG[ind].b), np.absolute(dfBetaAG[ind].b), 'o', alpha= .3)
#plt.plot(dfBetaA[ind].b, dfBetaAG[ind].b, 'o', alpha= .3)
plt.xlim(-12, 12)
plt.ylim(-12, 12)
plt.plot(dfBetaA[ind].b, dfBetaA[ind].b*dfBetaG[ind].b, 'o', alpha= .3)


ind1= (dfBetaA[ind].b < 0) & (dfBetaG[ind].b <0)
ind2= (dfBetaA[ind].b > 0) & (dfBetaG[ind].b >0)
ind3= dfBetaAG[ind].b > 0

ind4= (ind1 & ind3)
dfBetaA[ind][ind4].ens_gene.to_csv('../output/discordant_betas_age::genotype_less_0.csv', index=False)


listform= dfBetaA[ind][ind4].ens_gene.unique()
df_analysis= \
    hgt.implement_hypergmt_enrichment_tool('aging::genotype less than 0, discordant', listform,\
                                    tissue_df, qvalEn)
#save results to csv file
df_analysis.to_csv('../output/EnrichmentAnalysisResults/aging::genotype less than 0, discordant_Enrich', index= False)
#reopen the file and add a comment with relevant info for the file
#aname= agingANDgenotype_sig_Enrich
#line= '#' + aname+'\n'
#rsq.line_prepender('../output/EnrichmentAnalysisResults/'+fname, line)
#
#plot top fold change tissues
hgt.plotting_and_formatting(df_analysis, ytitle= 'aging::genotype less than 0, discordant Enrich')





tissues= ['neuron', 'muscle']
df_exp= rsq.organize(tissues, tissue_df)

colors2= ['#ffff33','#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00']
rsq.kegg(qval, dfBetaA[ind], df_exp, colors2, title= 'Aging',\
        xlab= r'$\beta_{\mathrm{Aging}}$', ylab= 'Density')
rsq.kegg(qval, dfBetaG[ind], df_exp, colors2, title= 'Genotype',\
        xlab= r'$\beta_{\mathrm{Genotype}}$', ylab= 'Density')
rsq.kegg(qval, dfBetaAG[ind], df_exp, colors2, title= 'Aging::Genotype',\
        xlab= r'$\beta_{\mathrm{Aging::Genotype}}$', ylab= 'Density')
             




#Make figures to show what regressions look like.....
def example(B, x, y):
    """ 
    just a function to show a contour plot example for a linear regression with
    interactions
    """
    b1, b2, b3= B #unpack
    xx, yy= np.meshgrid(x, y)    #mesh
    out= b1*xx + b2*yy + b3*xx*yy #model
    #output a normalized model that always has max absolute value 1 and
    #min of 0
    return (out+np.abs(np.min(out)))/np.max(np.abs(out))
    
x= np.linspace(0, 10)
y= np.linspace(0, 10)
xx, yy= np.meshgrid(x, y)
#set up the figure
fig, (ax1, ax2, ax3)= plt.subplots(ncols= 3, sharey= True, figsize= (16, 4))
cmap = sns.light_palette("navy", reverse=True, as_cmap= True)

#positive interaction
regress= example([1, 1, .5], x, y)
cs= ax1.contourf(xx, yy, regress, cmap= cmap, alpha= 0.7)
ax1.set_title('Linear Model w +ive interactions')
ax1.set_ylabel('y', fontsize= 15)

#no interaction
regress= example([1, 1, 0], x, y)
ax2.contourf(xx, yy, regress, cmap= cmap, alpha= 0.7)
ax2.set_title('Linear Model w/o interactions')
ax2.set_xlabel('x', fontsize= 15)

#negative interaction
regress= example([1, 1, -.5], x, y)
ax3.contourf(xx, yy, regress, cmap= cmap, alpha= 0.7)
ax3.set_title('Linear Model w -ive interactions')
fig.colorbar(cs, ax=ax3, shrink=0.9)
#save
plt.savefig('../output/ExamplesRegression.png')
    
    
    