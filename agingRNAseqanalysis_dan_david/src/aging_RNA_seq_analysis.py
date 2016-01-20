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


#==============================================================================
#Formatting
#==============================================================================
def fixer_upper(df1, df2, df3, gene_col_name= 'gene', reads_col_name= 'reads'):
    """
    A function that takes 3 dataframes from different files, appends their columns,
    sets the wbids as the gene name
    and calculates a mean and a std column
    returns a df with 5 columns (3 measurements, 1mean, 1 std)
    
    WARNINGS: DROPS NAS
    """
    df1.set_index(gene_col_name, inplace= True)
    df2.set_index(gene_col_name, inplace= True)
    df3.set_index(gene_col_name, inplace= True)
    
    df1.sort_index(inplace= True)
    df2.sort_index(inplace= True)
    df3.sort_index(inplace= True)    
    
    df1[reads_col_name+'2']= df2[reads_col_name]
    df1[reads_col_name+'3']= df3[reads_col_name]
    
    df1['reads_mean']= df1.mean(axis= 1)
    df1['reads_std']= df1.std(axis= 1)
    
    df1.dropna(inplace= True)
    return df1
    
    
    
    
    
    
    
    
    
    
#==============================================================================
# Load files
#==============================================================================
    

dfJKAged= pd.read_csv("../input/JK_aged_adult_1.count",\
 header= None, names= ['gene', 'reads'], sep= '\t')
dfJKAged2= pd.read_csv("../input/JK_aged_adult_2.count",\
 header= None, names= ['gene', 'reads'], sep= '\t')
dfJKAged3= pd.read_csv("../input/JK_aged_adult_3.count",\
 header= None, names= ['gene', 'reads'], sep= '\t')

dfJKYoung= pd.read_csv("../input/JK_young_adult_1.count",\
 header= None, names= ['gene', 'reads'], sep= '\t')
dfJKYoung2= pd.read_csv("../input/JK_young_adult_2.count",\
 header= None, names= ['gene', 'reads'], sep= '\t')
dfJKYoung3= pd.read_csv("../input/JK_young_adult_3.count",\
 header= None, names= ['gene', 'reads'], sep= '\t')

dfN2Aged= pd.read_csv("../input/N2_aged_adult_1.count",\
 header= None, names= ['gene', 'reads'], sep= '\t')
dfN2Aged2= pd.read_csv("../input/N2_aged_adult_1.count",\
 header= None, names= ['gene', 'reads'], sep= '\t')
dfN2Aged3= pd.read_csv("../input/N2_aged_adult_1.count", header= None, names= ['gene', 'reads'], sep= '\t')

dfN2Young= pd.read_csv("../input/N2_young_adult_1.count",\
 header= None, names= ['gene', 'reads'], sep= '\t')
dfN2Young2= pd.read_csv("../input/N2_young_adult_1.count",\
 header= None, names= ['gene', 'reads'], sep= '\t')
dfN2Young3= pd.read_csv("../input/N2_young_adult_1.count",\
 header= None, names= ['gene', 'reads'], sep= '\t')

dfJKAged= fixer_upper(dfJKAged, dfJKAged2, dfJKAged3)
dfJKYoung= fixer_upper(dfJKYoung, dfJKYoung2, dfJKYoung3)
dfN2Aged= fixer_upper(dfN2Aged, dfN2Aged2, dfN2Aged3)
dfN2Young= fixer_upper(dfN2Young, dfN2Young2, dfN2Young3)


#define column names 
genotype_cols= ['WTY','WTA','JKY','JKAged']
#place all files into a single dataframe
df= pd.concat([dfN2Young, dfN2Aged, dfJKAged, dfN2Aged], \
    keys= genotype_cols, levels= ['genotype', 'reads'], axis= 1)


#gene lengths:
df_gl= pd.read_csv('../input/c_elegans_gene_lengths_PRJNA13758.txt')
df_gl.set_index('WBID',inplace = True)
df_gl.sort_index(inplace= True)


def tpm_calc(df= df, dfgl= df_gl):
    """
    WARNING: DF1 AND DFGL MUST BE PREINDEXED AND PRESORTED BY GENEID      
    A function to calculate transcripts per million for rna-seq data
    """
    
    if df.shape[0] == dfgl[0]:
        t= '_tpm'
        for i in np.arange(5):
            if i > 0 and i < 3:
                r= 'reads'+str(int(i))
            elif i == 0:
                r= 'reads'
            elif i == 3:
                r= 'reads_mean'
            else:
                r= 'reads_std'
            df[r+t]= df[r]/dfgl['length'] #normalize for read length
            df[r+t]= df[r+t]/df[r+t].sum()*10**6 #normalize
    else:
        print('Dataframes are not the same size')
        if df.shape[0] > dfgl[0]:
            print('Some genes are missing from the length df')
        else:
            print('dfgl is too long, please check')
        
        if df.shape[0]== dfgl[0] +5:
            print('this df might still have 5 lines of code not associated with \
            a WBID. Check to make sure they are not preprocessing info')
    


def plot_tha_hist(df, n_obs, plot_name, xname, yname):
    fig, ax = plt.subplots()
    bins2= np.logspace(-10, 10, n_obs)
    ax.hist(df, bins= bins2)
    ax.set_xscale('log')
    ax.set_title(plot_name)
    ax.set_xlabel(xname)
    ax.set_ylabel(yname)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')













dfJKAged['JKA_N2A']=(.01+dfJKAged.avg)/(.01+dfN2Aged.avg)
dfJKAged['JKA_N2Y']=(.01+dfJKAged.avg)/(.01+dfN2Young.avg)
dfJKAged['JKA_JKY']=(.01+dfJKAged.avg)/(.01+dfJKYoung.avg)
dfJKAged['JKY_N2A']=(.01+dfJKYoung.avg)/(.01+dfN2Aged.avg)
dfJKAged['JKY_N2Y']=(.01+dfJKYoung.avg)/(.01+dfN2Young.avg)
dfJKAged['N2A_N2Y']=(.01+dfN2Aged.avg)/(.01+dfN2Young.avg)


n= 2
names= dfJKAged.index[dfJKAged.JKA_N2A>10**n].values
names2= dfJKAged.index[dfJKAged.JKA_N2Y>10**n].values
names3= dfJKAged.index[dfJKAged.JKA_JKY>10**n].values
names4= dfJKAged.index[dfJKAged.JKY_N2A>10**n].values
names5= dfJKAged.index[dfJKAged.JKY_N2Y>10**n].values
names6= dfJKAged.index[dfJKAged.N2A_N2Y>10**n].values
#fix name formatting to get into WBID format
#all WBIDs are 14 chars long, but you could use lambda functions to remove
#formatting as well. but i'm too lazy to do this in proper program form
#g= lambda x: x.find('|')
#end= map(g, names)


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


tissue_df= pd.read_csv("/Users/davidangeles/Downloads/dictionary.csv")
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
    print(analysis_name+'\n')
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

def implement_hypergmt_enrichment_tool(analysis_name, gene_list, tissue_df= tissue_df, alpha= 0.01):
    """
    Calls all the above functions
    
    gene_list: a list of non-redundant gene names
    tissue_df
    alpha: significance threshold, defaults to 0.01
    """
    
    
    
    print('Executing script\n')
    p_hash= hgf(gene_list, tissue_df)
    
    q_hash= return_enriched_tissues(p_hash, alpha, analysis_name)
    
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
q, p= implement_hypergmt_enrichment_tool('JKAged/N2A 10**{0}'.format(n), names, tissue_df, 0.1)
q2, p= implement_hypergmt_enrichment_tool('JKAged/N2Y 10**{0}'.format(n), names2, tissue_df, 0.1)
q3, p= implement_hypergmt_enrichment_tool('JKAged/JKY 10**{0}'.format(n), names3, tissue_df, 0.1)
q4, p= implement_hypergmt_enrichment_tool('JKY/N2A 10**{0}'.format(n), names4, tissue_df, 0.1)
q5, p= implement_hypergmt_enrichment_tool('JKY/NY2 10**{0}'.format(n), names5, tissue_df, 0.1)
q6, p= implement_hypergmt_enrichment_tool('N2A/N2Y 10**{0}'.format(n), names6, tissue_df, 0.1)



#==============================================================================
# plotting functions
#==============================================================================
def plot_tha_seqs(dfx, dfy, name, namex, namey, plot_line= True, xmax= None, ymax= None):
    fig, ax = plt.subplots()
    ax.scatter(dfx, dfy, s= 1, c= 'b', alpha= 0.3 )
    if xmax == None:
        xmax= dfx.max()
    if ymax == None:
        ymax= dfy.max()
    ax.set_xlim(10**-2, xmax) 
    ax.set_ylim(10**-1, ymax)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_title(name)
    ax.set_xlabel(namex)
    ax.set_ylabel(namey)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    if plot_line:
        plt.gca().plot([10**-2, xmax], [10**-2, xmax], '-', lw= 2)
        

plot_tha_seqs(df['reads_mean'], df['reads_mean']/dfJKAged['reads_mean'], '', 'fem-1 aged', 'n2aged')

