# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import tissue_enrichment_analysis as tea
import os
import numpy as np


sns.set_context("paper")


mag= 2 #value of beta from regression
qval= .1 #qvalue from regression
qvalEn= 0.05 #q value for enrichment analysis (tissues)
                
#dirLists= '../output/Gene_lists_for_analysis'
#if not os.path.exists(dirLists):
#    os.makedirs(dirLists)
    
#dirGraphs= '../output/Graphs'
#if not os.path.exists(dirLists):
#    os.makedirs(dirGraphs)

os.chdir('./')
#pos beta means high old adults
dfBetaA= pd.read_csv("../input/agebeta_wt.csv")
dfBetaA.dropna(inplace= True)
#pos beta means high in fog2
dfBetaG= pd.read_csv("../input/genotypebeta_wt.csv")
dfBetaG.dropna(inplace= True)
#pos beta means high in fog2-aged
dfBetaS= pd.read_csv("../input/spermbeta_wt.csv")
dfBetaS.dropna(inplace= True)
#likelihood ratio test results
dfLRT= pd.read_csv("../input/lrt.csv")
dfLRT.dropna(inplace= True)


#load results from aging analysis to find any genes that overlap


#sort by target_id
dfBetaA.sort_values('target_id', inplace= True)
dfBetaG.sort_values('target_id', inplace= True)
dfBetaS.sort_values('target_id', inplace= True)
dfLRT.sort_values('target_id', inplace= True)

#fetch tissue dictionary
tissue_df= tea.fetch_dictionary()

#make directories
dirEn= '../output/graphs/enrichment'
#dirVolc= '../output/graphs/volcano'
dirRNAi= '../output/rnai/'

dirs= [dirEn, dirRNAi]
for d in dirs:
    if not os.path.exists(d):
        os.makedirs(d)
        
#==============================================================================
# Identify candidates for RNAi
# Good targets: small qval, large positive b vals, not previously described
#==============================================================================
def exclude(df, excluded_genes, col):
    ind= (~df[col].isin(excluded_genes))
    return df[ind]

def find_molecular_targets(df, to_be_removed, cols='ens_gene', x= 'b', q= 0.1):
    """
    Given a dataframe df, return a new dataframe that:
    Doesn't have WBIDs present in the exclude series
    Has only genes that have q value < qval
    """
    if cols in [str, int, float]:
        cols= [cols]
    
    df.sort_values(x, inplace= True)
    sig= (df.qval < q) #take only sig genes
    
    temp= df[sig].copy() #remove all non-sig genes and make a temp copy
    if isinstance(to_be_removed, list):
        for i, excluded_gene_list in enumerate(to_be_removed):
            temp= exclude(temp, excluded_gene_list, cols[i])
    return temp

def direction_specific_tissue_analysis(anames, fnames, df, Lind, genes='ens_gene'):
    """
    Given a single dataframe, perform a tissue analysis for genes
    using the selection indices in Lind
    """
    if not isinstance(anames, list):
        raise ValueError('anames must be a list!')
    if not isinstance(fnames, list):
        raise ValueError('fnames must be a list!')
    if len(anames) != len(fnames):
        raise ValueError('fnames and anames must match length')
        
    for i, aname in enumerate(anames):
        ind= inds[i]
        df_results, unused= tea.enrichment_analysis(df[ind][genes], tissue_df, qvalEn, show= False)
        tea.plot_enrichment_results(df_results, title= aname, dirGraphs= dirEn)
        plt.close()


#set the path to sve the rnai candidates
#and make sure to exclude known genes from analysis
path= '../output/RNAi Candidates/'                
cols= ['ens_gene', 'target_id']
#==============================================================================
# 
#==============================================================================
#also exclude genes that have significant betas in any other list
x= dfBetaG[dfBetaG.qval<qval].target_id
y= dfBetaS[dfBetaS.qval<qval].target_id
excluded2= pd.Series(x.append(y).unique())

#enrichment on genes associated only with aging
aging_set= find_molecular_targets(dfBetaA, [excluded2], ['target_id'], q=qval)
aname1='aging_up'; fname1= 'aging_up.csv'
aname2='aging_down'; fname2= 'aging_down.csv'
anames= [aname1, aname2]
fnames= [fname1, fname2]
inds= [(aging_set.b>0), (aging_set.b<0)]
direction_specific_tissue_analysis(anames, fnames, aging_set, inds)
#==============================================================================
# 
#==============================================================================
#now for genotype, same thing
x= dfBetaA[dfBetaA.qval<qval].target_id
y= dfBetaS[dfBetaS.qval<qval].target_id
excluded2= pd.Series(x.append(y).unique())
excluded= [excluded2]

genotype_set= find_molecular_targets(dfBetaG, [excluded2], ['target_id'], q=qval)
aname1='genotype up'; fname1= 'genotype up.csv'
aname2='genotype down'; fname2= 'genotype down.csv'
anames= [aname1, aname2]
fnames= [fname1, fname2]
inds= [(genotype_set.b>0), (genotype_set.b<0)]
direction_specific_tissue_analysis(anames, fnames, genotype_set, inds)
#==============================================================================
#   
#==============================================================================
#interactions    
#more than linear increase in age, up in fog2 more than wt
#only take genes that go way up during age
dfBetaS.sort_values('b', inplace= True)

x= dfBetaA[dfBetaA.qval<qval].target_id
y= dfBetaG[dfBetaG.qval<qval].target_id
ind= (dfBetaS.target_id.isin(dfLRT[dfLRT.qval > qval].target_id)) & (dfBetaS.qval < qval)
z= dfBetaS[ind].target_id
excluded= pd.Series(np.unique(np.concatenate([x, y, z])))


interaction_set= find_molecular_targets(dfBetaS, [excluded], ['target_id'], q=qval)
aname1='sperm hollow'; fname1= 'genotype up.csv'
aname2='sperm full'; fname2= 'genotype down.csv'
anames= [aname1, aname2]
fnames= [fname1, fname2]
inds= [(interaction_set.b>0), (interaction_set.b<0)]
direction_specific_tissue_analysis(anames, fnames, interaction_set, inds)
  
df1= pd.read_csv('../../agingRNAseq/output/RNAi Candidates/CandidatesAge_HighInOld.csv')
df2= pd.read_csv('../../agingRNAseq/output/RNAi Candidates/CandidatesAge_LowInOld.csv')
df3= pd.read_csv('../../agingRNAseq/output/RNAi Candidates/CandidatesAgingXGenotypeLessThanExp.csv')
df4= pd.read_csv('../../agingRNAseq/output/RNAi Candidates/CandidatesAgingXGenotypeMoreThanExp.csv')
df5= pd.read_csv('../../agingRNAseq/output/RNAi Candidates/CandidatesGenotype_HighInOld.csv')
df6= pd.read_csv('../../agingRNAseq/output/RNAi Candidates/CandidatesGenotype_LowInOld.csv')

ldf= [df1, df2, df3, df4, df5, df6]
genes= np.array([])
for df in ldf:
    x= df['ens_gene'].values
    genes= np.concatenate([genes, x])

ind= interaction_set.ens_gene.isin(genes)

interaction_set[~ind].tail(55).to_csv(
'../output/RNAi/Sperm_Hollow_New.csv')
interaction_set[~ind].head(55).to_csv(
'../output/RNAi/Sperm_Full_new.csv')
interaction_set[ind]

interaction_set.tail(55).to_csv(
'../output/RNAi/Sperm_Hollow_Complete.csv')
interaction_set.head(55).to_csv(
'../output/RNAi/Sperm_Full_Complete.csv')

        
        