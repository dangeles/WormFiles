# -*- coding: utf-8 -*-
"""
Created on Mon Nov  3 13:28:07 2014
A script to re-format a couple of files that were given me by Margie into the adequate format.
These scripts are SWAG1 and SWAG2 gene expression profiles. 
@author: davidaalbores
"""
from __future__ import division, print_function, absolute_import
import pandas as pd

#Load the gene name dictionary
names= pd.read_csv("../input/c_elegans.PRJNA13758.WS241.livegeneIDs.unmaprm.txt", sep= "\t", comment= "#")

#Load the control database
control= pd.read_csv("../input/AVAcontrol.txt", sep= "\t", comment= "#")
control.columns= ['gene', 'reads']
##Experimental
#experimental= pd.read_csv("../input/glp-4dissectedgut(SWAG1).txt", sep= "\t", comment= "#")

#Discard the first line of the dataframes; it's useless. 
#control= control[1:]
#experimental= experimental[1:]

#Make a new dataframe that only contains the counts and the gene name, dropping any NAs
control= control.loc[:, ['gene', 'reads']].dropna()
#experimental= experimental.loc[:, ['SWAG1', 'gene']].dropna()

#change the names into WB IDs
for i in range(0, len(control)):
    #get the gene associated with locus i
    g= control.iloc[i].gene
    if names[names.GeneName == g].empty == False:
       #make a variable x that extracts the worm ID from the name df. this is a
    #dataframe itself!
        print('here')
        x= names[names.GeneName == g].WBID 
        #Replace the gene name in the control array        
        control= control.replace(g, x.iloc[0])
    else:
        if names[names.WBID == g].empty == True:
            #There's going to be a problem.
            print(names[names.WBID == g])
            control= control.replace(g, float('nan'))
        else:
            print(names[names.WBID == g])

##change the names into WB IDs
#for i in range(0, len(experimental)):
#    #get the gene associated with locus i
#    g= experimental.iloc[i].gene
#    if names[names.GeneName == g].empty == False:
#       #make a variable x that extracts the worm ID from the name df. this is a
#    #dataframe itself!
#        x= names[names.GeneName == g].WBID 
#        #Replace the gene name in the experimental array        
#        experimental= experimental.replace(g, x.iloc[0])
#    else:
#        #There's going to be a problem.
#        experimental= experimental.replace(g, float('nan'))

    
control= control.dropna()
#experimental= experimental.dropna()


control.to_csv("../input/controlAVA131014.txt", sep=',')
#experimental.to_csv("experimentalGE.txt", sep=',')




