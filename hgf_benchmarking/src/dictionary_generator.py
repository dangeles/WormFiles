# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from urllib.request import urlopen
import simplejson
import json
import numpy as np
import pandas as pd
#main solr url
solr_url = 'http://wobr.caltech.edu:8082/solr/anatomy/';

#queries as lambda functions
#query for terms. Finds terms that have x or more annotating genes
query_terms= lambda x: 'select?qt=standard&indent=on&wt=json&version=2.2&fl=id&start=0&rows=0&q=document_category:bioentity&facet=true&facet.field=regulates_closure&facet.limit=-1&facet.mincount={0}&facet.sort=count&fq=source:%22WB%22&fq=-qualifier:%22not%22'.format(x)

#query for relationships. given a wbbt ID, find the nodes connected to it. 
query_relation= lambda x: "select?qt=standard&fl=topology_graph_json&version=2.2&wt=json&indent=on&rows=1&q=id:%22{0}%22&fq=document_category:%22ontology_class%22".format(x)

#query for number of genes. given a wbbt ID, find genes assoc. with it. 
query_genes= lambda x: "select?qt=standard&indent=on&wt=json&version=2.2&fl=id&start=0&rows=10000&q=document_category:bioentity&fq=source:%22WB%22&fq=-qualifier:%22not%22&fq=regulates_closure:%22{0}%22".format(x)

query_readable= "select?qt=standard&fl=id,annotation_class_label&version=2.2&wt=json&indent=on&rows=100000&q=id:*&fq=document_category:ontology_class&fq=-is_obsolete:true"

#==============================================================================
# timeit using query_terms: 134ms
# using query_relation:      77ms
#==============================================================================
def solr_query(query, main= solr_url):
    """
    Given a query, append it to the main url. Open URL and
    use simplejson to load the results
    """
    conn= urlopen(solr_url+ query)
    return simplejson.load(conn)
#==============================================================================
# find_all_nodes timeit: 114ms
#==============================================================================
def find_all_nodes_with(x):
    """
    Given a minimum annotation filter, find all the
    wbbt ids that have more than or equal to that number
    of gene annotations
    """
    if type(x) not in [int, float]:
        raise ValueError('x must be an int or float')
    if type(x) == float:
        x= int(x)
        
    #get all tissue terms with at least 1 entry 
    rsp_terms = solr_query(query_terms(x))
    #extract all wbbts from query
    i=0
    wbbt= []
    for k in enumerate(rsp_terms['facet_counts']['facet_fields']['regulates_closure']):
        if i%2 == 0:
            wbbt.append(k[1]) 
        i+=1
    
    return wbbt
#==============================================================================
# find_family timeit: 81.5ms
#==============================================================================
def find_family(wbbt):
    """
    For a given node, find it's complete family
    and return two arrays, 'fathers_array' and 'daughters_array'
    """
    
    if type(wbbt) != str:
        raise ValueError('wbbt must be a string')
    if wbbt[0:5] != 'WBbt:':
        raise ValueError('wbbt must start with WBbt:')    
    if len(wbbt) != 12:
        raise ValueError('wbbt must be 12 letters long')
    
    #get the json object
    rsp_rlshp= solr_query(query_relation(wbbt))  
    fathers_array= [] #initialize father array
    daughters_array= [] #initialize daughters array
    
    #extract the array with all the right information
    array_of_rlshps= rsp_rlshp['response']['docs'][0]
    
    #go through the array, turning each line into a dictionary
    #these mini-dictionaries contain the edges between nodes
    for j in json.loads(array_of_rlshps['topology_graph_json'])['edges']:        
        #if the object isn't the same as the wbbt, object is parent to wbbt        
        #if object is same as wbbt, wbbt is parent to subject
        if wbbt != j['obj']:
            fathers_array.append(j['obj'])
        else:
            daughters_array.append(j['sub'])

    return fathers_array, daughters_array
#==============================================================================
# find_all_families timeit (limit= 20): 1.83s --> for all nodes with 1 gene == 1.3hrs
#==============================================================================
def find_all_families(array, flux_control= False, limit= 20):
    """
    Given an array, find all the sets of daughters and parents
    and return these sets in two hashes, parents and daughters.
    
    parents -- key = daughter, value= parents (in list form)
    daughters -- key= parent, value= daughters (in list form)
    """
    parents= {} #key is daughter node, value is a list of parents
    daughters= {} #key is parent node, value is a list of daughters

    i=0
    for wbbt in array:        
        fathers_array, daughters_array= find_family(wbbt)
        
        if len(fathers_array) > 0:
            parents[wbbt]= fathers_array
        if len(daughters_array) >0:
            daughters[wbbt]= daughters_array
        
        if flux_control:
            if i > limit:
                break
        i+=1
        next

    return parents, daughters    
#==============================================================================
# 
#==============================================================================
def find_complete_daughters(wbbt, daughters_all, daughters_filter):
    """"
    wbbt -- the array containing the wbbts where filtered daughters are
    daughters_all -- the complete dictionary of daughters
    daughters_filter -- the dictionary of daughters post-filtering
    
    given daughters all, daughters_filter, find the set of all dads
    that is complete 
    """
    kill_dad= []    
    
    #check all the dads in the complete daughter set
    for dad in iter(daughters_all):
        #if that dad still exists after filtering
        if dad in daughters_filter:
            #check if his daughter set is the same as in the complete set
            #if so, kill him
            daughters_found= []
            for daughter in daughters_filter[dad]:
                if daughter in wbbt:
                    daughters_found.append(daughter)
                
            if len(daughters_all[dad]) == len(daughters_found):
                kill_dad.append(dad)
        
    return kill_dad
#==============================================================================
#     
#==============================================================================
def trim(array, remove):
    """
    Given two arrays, find all the entries in the 2nd array present in the 1st
    and remove them from array 
    """
    if type(array) != list:
        raise ValueError('array must be a list')
    if remove == str:
        remove= [remove]
    if type(remove) != list:
        raise ValueError('remove must be a list')
        
    #make sure you get the union of the list
    remove= set(remove)
    array= set(array)
    trimmed= list(array - remove)
    return trimmed
#==============================================================================
# Filter all nodes where all daughters are present
# timeit x=25, min_annot= 20: 3.8s
#==============================================================================
def filter_by_daughters(x, min_annot= 20, flux_control= True, limit= 20):
    """
    Given a number 'x', find all the nodes with at least 'x' annotations.
    Further trim that set of nodes by identifying succesful parents who
    have a complete progeny set. I.e., parents whose daughters all have
    at least 'x' annotations get killed and dropped from the set
    """    
    if type(x) != int:
        raise ValueError('x must be an int')
    if x < min_annot:
        raise ValueError('x must be larger than min_annot')

    #first find all nodes with 1
    wbbt_all= find_all_nodes_with(min_annot)
    wbbt_filter= find_all_nodes_with(x)
    
    #find parents, daughters
    parents_all, daughters_all= find_all_families(wbbt_all, flux_control, limit= limit)
    parents_filter, daughters_filter= find_all_families(wbbt_filter, flux_control, limit= limit)
    
    #check all the dads in the complete daughter set
    kill_dad= find_complete_daughters(wbbt_filter, daughters_all, daughters_filter)
    
    #remove doomed parents    
    wbbt_filter_final= trim(wbbt_filter, kill_dad)
    
    return wbbt_filter_final, kill_dad
#==============================================================================
# 
#==============================================================================
def find_genes(wbbt):
    """
    For a given wbbt, find the genes associated with it
    """
    rsp_genes= solr_query(query_genes(wbbt))
    #extract the array with all the right information
    array_of_genes= rsp_genes['response']['docs']
    #go through the array, turning each line into a dictionary
    genes= []
    for entry in array_of_genes:
        genes.append(entry['id'][3:]) #remove WB: from the string
    
    return genes
#==============================================================================
# 
#==============================================================================
def find_all_genes(array, flux_control= False, limit= 20):
    """
    For a given array, find all genes associated with each wbbt in the array
    and return an annotations hash where key= wbbt and value= list_of_genes
    """
    annotations= {}    
    i=0 
    
    if array in [int, float, str]:
        array= [array]
    for wbbt in array:
        annotations[wbbt]= find_genes(wbbt)
        if flux_control:
            if i > limit:
                break
            i+=1
    return annotations
#==============================================================================
# 
#==============================================================================
def calculate_similarity(sisters, annotations, thresh= 1.0, method= 'avg'):
    """
    Given an array of sisters and an annotations hash, figure out the
    similarity of the sisters to each other. If they pass the threshold
    'thresh' of similarity, they are too similar and should be removed
    in which case a value of 1 (remove) will be returned
    sisters -- array of wbbts
    annotations hash -- from find_all_genes
    thresh -- scalar between 0, 1
    method -- one of 'avg' or 'any'
    
    outputs: 
    1 - remove
    0 - don't remove
    """
    #tests    
    if type(thresh) not in [float, int]:
        raise ValueError('thresh variable must be float between 0 and 1.0')  
    if thresh > 1.0:
        thresh= 1.0            
    if method not in ['avg', 'any']:
        raise ValueError('method must be one of \'avg\', \'any\'')
    #if one sister (or empty) don't flag for removal
    if len(sisters) <= 1:
        return 0
        
    complete_set= []
    similarity= np.zeros(len(sisters))    
    #go through each sister and extract the genes assoc. with her
    for i, sister in enumerate(sisters):
        similarity[i]= len(annotations[sister])
        complete_set= complete_set+annotations[sister]  

    #calculate simialrity
    similarity= similarity/len(list(set(complete_set)))
    #apply the correct method and return the answer. 
    if method == 'avg':
        if similarity.mean() > thresh:
            return 1
        else:
            return 0
    else:
        if len(np.where(similarity > thresh)[0]) > 0:
            return 1
        else:
            return 0
#==============================================================================
# 
#==============================================================================
def kill_by_similarity(sisters, thresh, method= 'avg'):
    """
    Given an array, find all the sets of daughters and parents
    and return these sets in two hashes, parents and daughters.
    
    parents -- key = daughter, value= parents (in list form)
    daughters -- key= parent, value= daughters (in list form)
    """
    kill= 0 #daughters to remove from list
    genes= find_all_genes(sisters)
    kill= calculate_similarity(sisters, genes, thresh, method= method)
    return kill   
#==============================================================================
# 
#==============================================================================
def filter_by_similarity(array, method= 'avg', flux_control= False, limit= 20, 
                         filter_by_sim= False, thresh= 1.0):
    """
    Given an array, find all the sets of daughters and parents
    and return these sets in two hashes, parents and daughters.
    
    parents -- key = daughter, value= parents (in list form)
    daughters -- key= parent, value= daughters (in list form)
    """
    #warnings and bells
    if filter_by_sim not in [True, False]:    
        raise ValueError('filter_by_sim must be a boolean parameter True or False')
    if type(thresh) not in [int, float]:
        raise ValueError('thresh must be a float between 0 and 1.0')
    if thresh > 1.0:
        thresh= 1.0
    if flux_control not in [True, False]:
        raise ValueError('flux_control must be one of True or False')
    if len(array) == 0:
        raise ValueError('array is empty')
    
    parents= {} #key is daughter node, value is a list of parents
    daughters= {} #key is parent node, value is a list of daughters
    
    i=0
    remove= [] #where daughters who fail are kept
    for i, wbbt in enumerate(array):
        #solr query
        fathers_array, daughters_array= find_family(wbbt)
    
        if len(fathers_array) > 0:
            parents[wbbt]= fathers_array
        if len(daughters_array) >0:
            daughters[wbbt]= daughters_array
            
        if filter_by_sim:
            #solr query
            #if they are meant to be killed, add them to remove array
            if kill_by_similarity(daughters_array, thresh, method):
                remove= remove+daughters_array

        if flux_control:
            if i > limit:
                break
            i+=1    

    #remove the offending nodes    
    trimmed= trim(array, remove)
    
    #test
#    print('original {0}, trimmed {1}'.format(len(array), len(trimmed)))
            
    return trimmed, parents, daughters
#==============================================================================
# 
#==============================================================================
def filter_by_daughters_and_sim(x, thresh= 1.0, filter_by_sim= False, min_annot= 20, 
                                flux_control= True, limit= 20, method= 'avg'):
    """
    Given a number 'x', find all the nodes with at least 'x' annotations.
    Further trim that set of nodes by identifying succesful parents who
    have a complete progeny set. I.e., parents whose daughters all have
    at least 'x' annotations get killed and dropped from the set
    """
    if type(x) != int:
        raise ValueError('x must be an int')
    
    if x < min_annot:
        raise ValueError('x must be larger than min_annot')
        
    #first find all nodes with 1
    wbbt_all= find_all_nodes_with(min_annot)
    wbbt_x= find_all_nodes_with(x)
    
    #all parents
    parents_all, daughters_all= find_all_families(wbbt_all, flux_control, limit= limit)
    
    #filtered by similarity
    wbbt_filter, parents_filter, daughters_filter= \
    filter_by_similarity(wbbt_x, thresh= thresh, filter_by_sim= filter_by_sim, 
                         method= method,
                         flux_control= flux_control,limit= limit)
        
    
    #check all the dads in the complete daughter set
    kill_dad= find_complete_daughters(wbbt_filter, daughters_all, daughters_filter)
    #remove the doomed fathers    
    wbbt_filter_final= trim(wbbt_filter, kill_dad)
    
    #return complete list
    return wbbt_filter_final, kill_dad
#==============================================================================
#     
#==============================================================================
#==============================================================================
# Filter all nodes with less than 25 annotations
#==============================================================================
#build dictionary
def build_dictionary(tissues):
    #given a list of tissues, find the genes associated with each tissue and 
    #place them in a vector..... 
    
#    print(tissues)
    gene_dict= find_all_genes(tissues) #now all your tissues are in a hash
    genes= []
    
    for tissue in iter(gene_dict):
        genes= genes + gene_dict[tissue]
        
    genes= list(set(genes))
    
    mat= np.zeros(shape= (len(genes), len(tissues)))
    for i, gene in enumerate(genes):
        for j, tissue in enumerate(tissues):
            if gene in gene_dict[tissue]:
                mat[i, j]= 1
    
    cols= tissues
    df= pd.DataFrame(mat, columns= cols)
    df.insert(0, 'wbid', genes)
    return df
#==============================================================================
# 
#==============================================================================
def fetch_human_readable(query):
    
    rsp_names=solr_query(query)   
    
    #process names to return a hash in the form 'id': 'name'
    array_of_names= rsp_names['response']['docs']
    #go through the array, turning each line into a dictionary
    names= {}
    for entry in array_of_names:
        names[entry['id']]= entry['annotation_class_label'] 
    
    return names
#==============================================================================
# 
#==============================================================================
def implement(x, fname, filter_by_sim= True, thresh= .9, method= 'any', min_annot= 20, 
              flux_control= False,  **kwargs):
    #don't let threshold go below one
    if thresh*x < 1:
        thresh= 1- 1.5/float(x)
        print('warning, thresh was set too high, so the value was changed to {0:.2}'.format(thresh))
        
    if flux_control:
        print('warning, you are in debug mode. please make sure you want to be in this mode')
        
    tissues, killed= filter_by_daughters_and_sim(x, thresh= thresh, 
                                           filter_by_sim= filter_by_sim, 
                                           min_annot= min_annot, method= method, 
                                           flux_control= flux_control)
        
    #get names
    names= fetch_human_readable(query_readable)       
    df= build_dictionary(tissues)
    
    #build the new columns:
    col= {}
    for tid in tissues:
        col[tid]= names[tid]+' '+tid
    
    y= df.columns
    ordered= ['wbid']
    for column in y:
        if column != 'wbid':
            ordered.append(col[column])
    
    df.columns= ordered
    
    #drop the root term 0000100 bc it contains all the terms -- so it will
    #always test significant. 
    df.drop('C. elegans Cell and Anatomy WBbt:0000100', axis= 1, inplace= True)

    df.to_csv(fname, index= False)
    return df
    
x_array= [100, 50, 25]
thresh_array= [1, .96, .9]
methods= ['avg', 'any']
path= '../input/WS252AnatomyDictionary/'
    
for x in x_array:
    for thresh in thresh_array:
        for method in methods:
            fname= path + 'annot{}.thresh{}.method{}.csv'.format(x, thresh, method)
            df= implement(x, thresh= thresh, method= 'any', fname= fname)
            










#reference values:
ref, killed1= filter_by_daughters(25, min_annot=20, flux_control=False)
x10, killed= filter_by_daughters_and_sim(25, thresh= .1, filter_by_sim= True, 
                                        min_annot= 20, method= 'any', flux_control=False)
x90, killed= filter_by_daughters_and_sim(25, thresh= .9, filter_by_sim= True, 
                                        min_annot= 20, method= 'any', flux_control=False)
x100, killed= filter_by_daughters_and_sim(25, thresh= 1.0, filter_by_sim= True, 
                                        min_annot= 20, method= 'any', flux_control=False)


        
df= build_dictionary(x90)         
df.to_csv('../input/WS252AnatomyDictionary/annot25sim90.csv', index= False)






























            
flp= 'WBbt:0006828'
flpl= 'WBbt:0004798'        
flpr= 'WBbt:0004797'         
asil= 'WBbt:0003888'
asir= 'WBbt:0003887'
asi= 'WBbt:0005666'
f= [flpr, flpl]
a= [asir, asil]

g= find_all_genes(a)
calculate_similarity(a, g, thresh= .1, method= 'any')
kill_by_similarity(a, thresh= .1, method= 'any')  


#test find_family
#find the parent
parent= 'WBbt:0006783'      
daughter= 'WBbt:0008596'
parents , daughters= find_family(daughter)
if parent not in parents:
    raise ValueError('Parent WBbt:0006783 is not in list of parents for \
    node WBbt:0008596. Test failed')

#find the daughter
parent='WBbt:0006783'        
daughter= 'WBbt:0008596'
parents , daughters= find_family(parent)
if daughter not in daughters:
    raise ValueError('Daughter WBbt:0008596 is not in list of parents for \
    node WBbt:0006783. Test failed')

#test find_all_families
parent= 'WBbt:0006783'      
daughter= 'WBbt:0008596'
parents , daughters= find_all_families([daughter])
if parent not in parents[daughter]:
    raise ValueError('Parent WBbt:0006783 is not in list of parents for \
    node WBbt:0008596. Test failed')

parent= 'WBbt:0006783'      
daughter= 'WBbt:0008596'
parents , daughters= find_all_families([parent])
if daughter not in daughters[parent]:
    raise ValueError('Parent WBbt:0006783 is not in list of daughters for \
    node WBbt:0008596. Test failed')


#test 1
#filter_by_daughters_and_sim should give the same answer as filter_by_daughters
#reference values:
ref, killed1= filter_by_daughters(25, min_annot=20, flux_control=True)
x, killed= filter_by_daughters_and_sim(25, min_annot=20 )
            
            
if len(x) != len(ref):
    raise ValueError('filter_by_daughters_and_sim is not giving the\
    same answer as filter_by_daughters when not filtering for similarity')
         
#test 2
#if filter_by_daughters_and_sim has more nodes than filter_by_daughters, something is wrong
x, killed= filter_by_daughters_and_sim(25, thresh= .01, filter_by_sim= True, min_annot= 20)
if len(x) > len(ref):
    raise ValueError('Error: filter_by_daughters_and_sim is returning  more nodes\
    than filter_by_daughters when  filtering for similarity')

x2, killed= filter_by_daughters_and_sim(25, thresh= .01, filter_by_sim= True, 
                                        min_annot= 20, method= 'any')
if len(x2) > len(x):
    raise ValueError('Error: filter_by_daughters_and_sim is returning  more nodes\
    when using \'any\' than when using \'avg\'')
            
            
x, killed= filter_by_daughters_and_sim(25, thresh= .9, filter_by_sim= True, 
                                        min_annot= 20, method= 'any')
if len(x2) > len(x):
    raise ValueError('Error: filter_by_daughters_and_sim is returning  more nodes\
    when using \'any\' than when using \'avg\'')
            
            
#most stringent to not stringent filters. 
#check that nodes disappear in a predictable fashion.
#reference values:
ref, killed1= filter_by_daughters(25, min_annot=20, flux_control=False)
x, killed= filter_by_daughters(25, min_annot=20, flux_control=False)
x10, killed= filter_by_daughters_and_sim(25, thresh= .1, filter_by_sim= True, 
                                        min_annot= 20, method= 'any', flux_control=False)
x90, killed= filter_by_daughters_and_sim(25, thresh= .9, filter_by_sim= True, 
                                        min_annot= 20, method= 'any', flux_control=False)
x100, killed= filter_by_daughters_and_sim(25, thresh= 1.0, filter_by_sim= True, 
                                        min_annot= 20, method= 'any', flux_control=False)
            
if len(x100) < len(ref):
    raise ValueError('filter_by_daughters_and_sim doesn\'t reduce to filter by daughters when\
    threshold is set to 1')            
