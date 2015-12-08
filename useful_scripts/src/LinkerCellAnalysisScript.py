# -*- coding: utf-8 -*-
"""
Created on Sun Dec  6 11:56:42 2015

@author: davidangeles

Script for analysis of Mihoko's linker cell data
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, '/Users/davidangeles/github/tissue_enrichment_tool_hypergeometric_test/src/')

import hypergeometricTests.py

df= pd.read_csv('../input/linker_8k_genes_initial_05jun2011.csv', sep= '\t') 

x= np.linspace(0, len(df.columns)-1, len(df.columns)) #column numbers
df.drop(df.columns[map(int, x[5:])], axis= 1, inplace=1)

df.rename(columns= {'nhr-67': 'nhr67'}, inplace= True)
df.L3= df.L3/df.L3.sum()*10**6
df.L4= df.L4/df.L4.sum()*10**6
df.nhr67= df.nhr67/df.nhr67.sum()*10**6
df.N2= df.N2/df.N2.sum()*10**6


fig, ax = plt.subplots()
bins= np.floor(np.sqrt(len(df)))
df.L3.hist(ax=ax, bins=bins, bottom=0.1)
df.L4.hist(ax=ax, bins=bins, bottom=0.1)
df.N2.hist(ax=ax, bins=bins, bottom=0.1)
df.nhr67.hist(ax=ax, bins=bins, bottom=0.1)
ax.set_yscale('log')
ax.set_xscale('log')


#WT linker cell:
df['WT_LC_max']= df[['L3', 'L4']].apply(np.max, 1)
df['ratio']= df['WT_LC_max']/(df['N2']+10**-5)

fig, ax = plt.subplots()
df.ratio[df.ratio > 10**5].hist(ax=ax, bins=bins, bottom=0.1)
ax.set_yscale('log')
ax.set_xscale('log')

names= df.Gene[df.ratio > 10**5].values

f= open('../output/ListOfNamesWithRatio1mill.csv', 'w')
for name in names:
    f.write(name)
    f.write('\n')
f.close()


dfWBID= pd.read_csv('../output/ListOfNamesWithRatio1mill_WBID.csv', header= None, names= ['WBID'])



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

from __future__ import division, print_function, absolute_import
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.special
from scipy import stats
from scipy import optimize
from scipy import misc
import os

#path= 'Users/dangeles/github/WormFiles/tissue_enrichment_hgf/tissue_enrichment/src"
#os.chdir(path)
#gene_file= "../input/controlGE.txt"
#genes1= pd.read_csv(gene_file)
#
#
##since the files we are using include read information
##remove the reads and keep only the gene names
#gene_list1= genes1[genes1.columns[0]].values
#gene_list2= genes2[genes2.columns[0]].values


tissue_df= pd.read_csv("/Users/davidangeles/github/WormFiles/tissue_enrichment_hgf/input/smalldictionary.txt")
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
        present.provided[present.wbid== item]= 1
        
    #remove NA's and make them 0
    present.provided= present.provided.fillna(0)
    
    #return df
    return present
    
#==============================================================================
#hgf is short for hypergeometric function
#==============================================================================
def hgf(gene_list, tissue_dictionary):
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
    for i, name in enumerate(tissue_dictionary.columns.values): 
        if total == 0:
            p_hash[name]= 1
        else:
            p_hash[name]= stats.hypergeom.sf(wanted_sum[name],total_genes, sums_of_tissues[name],total)
        
    #return the p-values. 
    return p_hash
    
    
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
def return_enriched_tissues(p_hash, alpha):
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
    for key, value in p_hash.iteritems():
        #place each value in the list
        p_values.append(value)
        #and the index of the value in the hash with value equal to the tissue
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
    
    #print the results. This will likely be modified to return a 
    #file or some such.
    print("q-values less than alpha = {0:.2} are considered statistically significant".format(alpha))
    print("\n\n")
    print("------------------------")
    print("tissue,q_value")
    for key, value in q_hash.iteritems():
        if value < alpha:
            print("{0},{1:.3}".format(key, value))
    print("------------------------\n\n")

        
    return q_hash
#==============================================================================
#     
#==============================================================================    

def implement_hypergmt_enrichment_tool(gene_list, tissue_df= tissue_df, alpha= 0.01):
    """
    Calls all the above functions
    
    gene_list: a list of non-redundant gene names
    tissue_df
    alpha: significance threshold, defaults to 0.01
    """
    
    
    print('Executing script\n')
    p_hash= hgf(gene_list, tissue_df)
    
    q_hash= return_enriched_tissues(p_hash, alpha)
    
    return q_hash, p_hash
#==============================================================================
#     
#==============================================================================
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
q, p= implement_hypergmt_enrichment_tool(names, tissue_df)
#q2= implement_hypergmt_enrichment_tool(gene_list2, tissue_df)
