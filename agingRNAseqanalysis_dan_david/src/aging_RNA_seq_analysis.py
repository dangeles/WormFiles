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
import os

#gene_lists from sleuth
dfBetaA= pd.read_csv("../input/table_agebeta_genes.csv")
dfBetaG= pd.read_csv("../input/table_genobeta_genes.csv")
dfBetaAG= pd.read_csv("../input/table_genocrossagebeta_genes.csv")

#tissue dictionary-- please cite David Angeles et al TEA publication (forthcoming)
#if using the enrichment tool 
tissue_df= pd.read_csv("/Users/davidangeles/Downloads/dictionary.csv")

#value of beta and qvalue to slice genes from
mag= 2
qval= .1

#slice all the relevant gene names out
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
# 
#==============================================================================
def pass_list(user_provided, tissue_dictionary):
    """
    A function to check which genes are in the 
    provided list of user provided names and to
    return a dataframe of presence or fail
    
    Takes two vectors, user_provided and tissue_dictionary in that order
    tissue_dictionary is a pandas dataframe that you should know about
    user_provided is a list of gene names
    """
    length= tissue_dictionary.shape[0] #how many genes in the dict
    
    #make an empty dataframe
    present= pd.DataFrame(index= range(length), columns= ['wbid','provided'])
    
    #fill wbid column with all the gene names
    present.wbid= tissue_dictionary.wbid
    
    #go through and pass attendance -- 1 if present, 0 otherwise
    for item in user_provided:
        present.provided[present.wbid==item]= 1
        
    #remove NA's and make them 0
    present.provided= present.provided.fillna(0)
    
    #return df
    return present
    
#==============================================================================
#hgf is short for hypergeometric function
#==============================================================================
def hgf(gene_list, tissue_dictionary, f, dirUnused):
    """
    Given a list of tissues and a gene-tissue dictionary,
    returns a p-dictionary for the enrichment of every tissue
    (a p-dictionary is a vector of length equal to the number
    of tissues in the tissue_dictionary, sorted by value)
    
    The entries of the p-vector are p-values not corrected
    for multiple hypothesis testing. 
    
    gene_list should be a list or list-like
    tissue_dictionary should be a pandas df
    """    
    
    #figure out what genes are in the user provided list
    present= pass_list(gene_list, tissue_dictionary)  
    
    #make a file to let user know what genes were used for the analysis
    present[present.provided == 1].wbid.to_csv(dirUnused+f[:-4]+'_unused_genes.csv', index= False)
    
    
    #slice out only the genes that were present from the user-provided list    
    wanted= present.wbid[present.provided==1]

    #re-index the dictionary s.t. the wbid is the index
    tissue_dictionary= tissue_dictionary.set_index('wbid')
    
    #number of tissues in the dictionary
    sums_of_tissues= tissue_df.sum()[1:] #this object can be identified by 
    #the column names of tissues and excludes gene IDs
    
    #total size of the urn
    total_genes= tissue_dictionary.shape[0] #total genes in the dictionary

    #slice out the rows from tissue_dictionary that came from the user-provided list
    wanted_dictionary= tissue_dictionary.loc[wanted]
    
    #get the total number of labels from each tissue
    wanted_sum= wanted_dictionary.sum()
    
    #get the total number of genes provided by the user that are in the dictionary
    total= wanted.shape[0]    
    
    #make a hash with the p-values for enrichment of each tissue. 
    p_hash= {}
    exp_hash= {} #expected number for each tissue
    for i, name in enumerate(tissue_dictionary.columns.values): 
        #if the total number of genes is zero, return p= 1 for all tissues
        if total == 0:
            p_hash[name]= 1
        else:
            #give entries as total no. of colors
            #total number of balls in urn
            #total number of colors measured
            #total number of balls picked out
            p_hash[name]= \
            stats.hypergeom.sf(\
            wanted_sum[name],\
            total_genes, \
            sums_of_tissues[name],\
            total)
            
            exp_hash[name]= stats.hypergeom.mean(total_genes, sums_of_tissues[name], total)	
            
                    
    #return the p-values, the genes associated with each tissue and the user
    #provided genes associate with each tissue. 
    return p_hash, exp_hash, wanted_dictionary
    
    
#==============================================================================
#     
#==============================================================================
def benjamin_hochberg_stepup(p_vals):
    """
    Given a list of p-values, apply a BH
    FDR correction and return the pertaining q values
    
    alpha is a scalar FDR value
    p_vals should be iterable
    """
    #sort the p_values, but keep the index listed
    index= [i[0] for i in sorted(enumerate(p_vals), key=lambda x:x[1])]
    
    #keep the p_values sorted
    p_vals= sorted(p_vals)
    q_vals= [None]*len(p_vals) #initialize an empty list
    prev_q= 0
    
    #BH Step Up begins here. 
    for i, p in enumerate(p_vals):
        q= len(p_vals)/(i+1)*p #calculate the q_value for the current point
        q= min(q, 1) #if q >1, make it == 1
        q= max(q, prev_q) #preserve monotonicity, since at endpoints the procedure sometimes makes it so that prev_q > new_q
        q_vals[i]= q #store the q_value
        prev_q= q #update the previous q_value
        
    #return q_vals and the index so we can match up each 
    #q_value with its appropriate tissue. 
    return q_vals, index
    
    
#==============================================================================
# 
#==============================================================================
def return_enriched_tissues(p_hash, alpha, analysis_name):
    """
    Given a hash of p-values
    (tissue -> p-values)
    apply an FDR correction and
    return the q_values that pass
    significance along with the 
    tissue type. 
    """
    
    #initialize a list, a hash and a counter
    p_values= []
    index_to_tissue= {}
    k= 0
    #go through the p_hash
    for key, value in p_hash.items():
        #place each value in the list
        p_values.append(value)
        #store the hash key in an array at the same k index as the value is in p_values
        index_to_tissue[k]= key
        #add 1 to the index
        k+=1
    
    #apply FDR, using BH stepup procedure
    q_values, index= benjamin_hochberg_stepup(p_values)
    
    #place everything in a hash
    q_hash= {}
    for i, value in enumerate(q_values):
        j= index_to_tissue[index[i]]
        q_hash[j]= value
    
#    #print the results. This will likely be modified to return a 
#    #file or some such.
#    print(analysis_name+'\n')
#    print\
#    ("q-values less than alpha = {0:.2} are considered statistically significant"\
#    .format(alpha))
#    print("\n\n")
#    print("------------------------")
#    print("tissue,q_value")
#    for key, value in q_hash.items():
#        if value < alpha:
#            print("{0},{1:.3}".format(key, value))
#    print("------------------------\n\n")
    
    return q_hash
#==============================================================================
#     
#==============================================================================    
def implement_hypergmt_enrichment_tool(analysis_name, gene_list, \
    tissue_df= tissue_df, alpha= 0.1, f='EnrichmentAnalysis.csv', \
    dirEnrichment= '../output/EnrichmentAnalysisResults', dirUnused= '../output/UnusedGenes'):
    """
    Calls all the above functions
    
    gene_list: a list of non-redundant gene names
    tissue_df
    alpha: significance threshold, defaults to 0.01
    f: filename for the enrichment analysis
    
    dirEnrichment: directory where the enrichment analysis will be placed
    dirUnusued: directory where the lists of unused genes will be deposited
    
    The directories are useful to specify when users provide multiple analyses 
    in batch
    """
    
    if f[-4:] != '.csv':
        if f[-4:] != '.txt':
            f= f+'.csv'
            
    print('Executing script\n')
    
    #create the directories where the results will go
    if not os.path.exists(dirEnrichment):
        os.makedirs(dirEnrichment)
    if not os.path.exists(dirUnused):
        os.makedirs(dirUnused)
        
    #calculat the enrichment
    p_hash, exp_hash, wanted_dic= hgf(gene_list, tissue_df, f, dirUnused+'/')
    
    #FDR correct
    q_hash= return_enriched_tissues(p_hash, alpha, analysis_name)
    

    #write the results to a file. 
    with open(dirEnrichment+'/'+f, 'w') as file:
        file.write('#Tissue Enrichment Analysis for {0}\n'.format(analysis_name))
        file.write('Tissue,Expected,Observed, Fold Change,Q value\n')
        for tissue, qval in q_hash.items():
            if qval < alpha:
                
                file.write(tissue)
                file.write(',')
                
                expected= exp_hash[tissue]
                file.write('{0:.2}'.format(expected))
                file.write(',')
                
                observed= wanted_dic[tissue].sum()
                file.write(str(observed))
                file.write(',')
                
                file.write('{0:.2}'.format(observed/expected))
                file.write(',')
                
                file.write('{0:.3}'.format(qval))
                file.write('\n')
    
    return q_hash#, p_hash
#==============================================================================
#     
#==============================================================================

def plotting_and_formatting(q, p, want_plots= 'yes', want_files= 'yes', max_n= '15'):
    return
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
# # # # # 
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================




#Run the whole thing:

qvalEn= 0.1

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



#run the whole thing for each dataset
for i, list_of_genes in enumerate(array_of_arrays):
    #run the tissue enrichment tool 
    aname= array_of_anames[i] #name of the analysis
    fname= array_of_fnames[i] #filename to store as
    q= implement_hypergmt_enrichment_tool(aname, list_of_genes, tissue_df, qvalEn, f= fname)
    
    
    #save genes with betas significantly different from zero in a list format
    k= array_of_strings[i][5:]
            
    filename= dirLists+'/'+ 'gene_list_{0}_beta{1}_q{2}.csv'.format(k, mag, qvalEn)
    #print('filename', filename)
    with open(filename, 'w') as f:
        for gene in list_of_genes:
            f.write(gene)
            f.write('\n')
    f.close()
