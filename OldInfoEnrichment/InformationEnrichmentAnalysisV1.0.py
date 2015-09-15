# -*- coding: utf-8 -*-
"""
A script to perform exploratory analysis of data via information theory tissue enrichment. 
A script that encodes a super basic enrichment tool for tissues in elegans. 
Author: David Angeles
Date: 6 Nov 2014
"""

from __future__ import division, print_function, absolute_import
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.special
import EnrichmentTools as et
from scipy import stats
from scipy import optimize
from scipy import misc
import seaborn as sns
#
#
##Load the dictionary
#dictionary= pd.read_csv("../input/gene_tissue_expression.txt", sep= "\t", comment= "#")
##Load the control database
#control= pd.read_csv("../input/controlGE.txt", sep= ",", comment= "#")
##Experimental
#experimental= pd.read_csv("../input/experimentalGE.txt", sep= ",", comment= "#")
#
##rename the columns
#control.columns= ['u', 'counts', 'id']
#experimental.columns= ['u', 'counts', 'id']
#dictionary.columns= ['wbid', 'GeneName', 'tissues']
#
##drop the first column of control and experimental -- it's useless!
#control.drop(control.columns[0], axis= 1, inplace=True)
#experimental.drop(experimental.columns[0], axis= 1, inplace=True)
#
##drop the GeneName column of dictionary, it also is useless. We work only with WBID
#dictionary= dictionary.loc[:, ['wbid', 'tissues']]
#
##==============================================================================
## Term extraction -- figure out how many distinct terms there are
##==============================================================================
##The dictionary should have 3 columns, GeneID, GeneName, and tissues. 
##Let's figure out what tissues are inside the dictionary
##put all the terms in a redundant list
#terms= []; 
#for element in dictionary.tissues:
#    #elements may have more than one tissue listed separated by commas
#    listElements= element.split(',')
#    for item in listElements:
#        terms.append(item);
#
##==============================================================================
## Over-representation within the dictionary
##==============================================================================
#i= 0;#initialize a counter
#representation= [[]];
#for element in set(terms):
#    representation[i].append(element);
#    count= terms.count(element);
#    representation[i].append(count);
#    representation.append([])
#    i+=1
#
##==============================================================================
## Dictionary of terms, non-redundant:
##==============================================================================
#terms= list(set(terms)) 
#
#
##==============================================================================
## Clean the databases; remove any entries not in the dictionary
##==============================================================================
##Control
#for name in control.id:
#    existence= dictionary.loc[dictionary.wbid == name].wbid.any()
#    if existence == False:
#        control= control[control.id != name];
##Experimental
#for name in experimental.id:
#    existence= dictionary.loc[dictionary.wbid == name].wbid.any()
#    if existence == False:
#        experimental= experimental[experimental.id != name]
#        
##Re-set index   
#control.reset_index(inplace= True)
#experimental.reset_index(inplace= True)
#
##==============================================================================
## Shannon Information Calculations
##==============================================================================
##Make a matrix of 0's and 1's, where each row is a tissue type
##And the columns are specific genes. 
##The specific genes MUST BE SORTED IN DESCENDING ORDER
##Control
#binaryControl= np.zeros(shape= (len(terms), len(control)))        
#controlTissues=[]
#controlRankTissues= {}
#for name in control.id:
#    currentTerms= dictionary.loc[dictionary.wbid == name].tissues.str.split(',') #make an array
#    currentTerms= currentTerms.iloc[0];
#    for label in terms:
#        if label in currentTerms:
#            index1= terms.index(label)
#            index2= control.index[control.id == name]
#            binaryControl[index1, index2]= 1;
#            controlTissues.append(label)
#            
#            if label in controlRankTissues.keys():
#                rank= control[control.id == name].index.values
#                controlRankTissues[label]= np.append(controlRankTissues[label], rank[0])
#            else:
#                rank= control[control.id == name].index.values
#                controlRankTissues[label]= np.array(rank[0])
#            
#            currentTerms.remove(label)
#            #currentTerms.remove(label)
#            #if the array is empty, you're done:
#            if not currentTerms:
#                break;
##Experiment
#binaryExperimental= np.zeros(shape= (len(terms), len(experimental)))
#expTissues= []      
#expRankTissues= {}
##make a new dataframecolumn that corresponds to a term. 
#for name in experimental.id:
#    currentTerms= dictionary.loc[dictionary.wbid == name].tissues.str.split(',') #make an array
#    currentTerms= currentTerms.iloc[0];
#    for label in terms:
#        if label in currentTerms:
#            index1= terms.index(label)
#            index2= experimental.index[experimental.id == name]
#            binaryExperimental[index1, index2]= 1;
#            expTissues.append(label)
#            
#            if label in expRankTissues.keys():
#                rank= experimental[experimental.id == name].index.values
#                expRankTissues[label]= np.append(expRankTissues[label], rank[0])
#            else:
#                rank= experimental[experimental.id == name].index.values
#                expRankTissues[label]= np.array(rank)
#                
#            currentTerms.remove(label)
#            if not currentTerms:
#                break
#
##==============================================================================
## Information of the Systems
##==============================================================================
#maxInfo= np.log(len(terms))
#infoControl= et.information(binaryControl, len(terms))
#infoExp= et.information(binaryExperimental, len(terms))
#
##Calculate the normalization number, N:
#NCont= binaryControl.shape[0]*binaryControl.shape[1]
#NExp= binaryExperimental.shape[0]*binaryExperimental.shape[1]
#
##==============================================================================
## Entropy Calculation
##==============================================================================
##Calculate the entropy of each sample with a variable size window. 
#XCont, H_Control= et.slider(binaryControl, len(terms))
#XExp, H_Experimental= et.slider(binaryExperimental, len(terms))
#
##==============================================================================
## Calculate some basic stats of the libraries in question 
##==============================================================================
#controlTissueCounts, controlTissues= et.countLabels(controlTissues)
#controlTissueCounts= et.completeAllTissues(terms, controlTissues, controlTissueCounts)
#
#expTissueCounts, expTissues= et.countLabels(expTissues)
#expTissueCounts= et.completeAllTissues(terms, expTissues, expTissueCounts)
#
#==============================================================================
# Sliding Window Calculation.
# Counts the number of occurences of strings of length n with value 1:
# So: looks for strings of the form 111....1 with exactly n ones.  
#==============================================================================
#XC, D_Control= deltaSlider(binaryControl, len(terms))
#XE, D_Exp= deltaSlider(binaryExperimental, len(terms))

#
#==============================================================================
# Result summary:
#==============================================================================
found, diff, fdiff, delta, error= et.summary(H_Control, H_Experimental, XCont, XExp, terms, NCont, infoControl, NExp, infoExp, controlTissueCounts, expTissueCounts, controlRankTissues, expRankTissues, len(control), len(experimental), fold1= False, search= 'none', interval3= '3')


print("\nmaxI ", maxInfo)
print("infoC ", infoControl)
print("infoExp ", infoExp)
print("infoC fraction {:f} ".format(infoControl/maxInfo))
print("infoExp fraction {:f}".format(infoExp/maxInfo))

#integral1= scipy.integrate.simps(H_Control[:, 521], dx= dx1)/NCont/infoControl
#integral2= scipy.integrate.simps(H_Experimental[:, 521], dx= dx2)/NExp/infoExp
#0 787  interval

#==============================================================================
# Plot the results
#==============================================================================
if error == 0:     
    import EnrichmentTools as et
    et.plotter(found, controlRankTissues, expRankTissues, len(control), len(experimental), fileName= 'violinplot.png')

#    et.shannonplotter(H_Control, H_Experimental, XCont, XExp, NCont, infoControl, NExp, infoExp,terms, found)



#
#totalEntropyC=np.array([H_Control[0,:].sum()])
#for i in range(1, len(H_Control)):
#    totalEntropyC= np.append(totalEntropyC, H_Control[i,:].sum())
#totalEntropyE=np.array([H_Experimental[0,:].sum()])
#for i in range(1, len(H_Experimental)):
#    totalEntropyE= np.append(totalEntropyE, H_Experimental[i,:].sum())
#plt.plot(XCont[:-3], totalEntropyC[:-3]/infoControl/NCont, '.')
#plt.plot(XExp[:-3], totalEntropyE[:-3]/infoExp/NExp, '-')
#
#
#
