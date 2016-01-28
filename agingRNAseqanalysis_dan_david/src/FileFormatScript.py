# -*- coding: utf-8 -*-
"""
Created on Mon Nov  3 13:28:07 2014
A script to re-format a file
@author: davidaalbores
"""
from __future__ import division, print_function, absolute_import
import pandas as pd

#Filenames!
fileInput= '../input/glp-4dissectedgut(SWAG1).txt'; sepInput= '\t'; comInput= '#';
newFile= '../input/experimentalGE.txt'
errorFile= '../output/ErrorFile.txt'
summaryFile= '../output/Summary.txt'

#Filename and location of the dataframe you want to translate, 
#kind of separator used and comment marker
filename= '../input/c_elegans.PRJNA13758.WS241.livegeneIDs.unmaprm.txt'
separator= "\t"
comments= "#"

#Load the gene name dictionary
names= pd.read_csv(filename, sep= separator, comment= comments)

#Load the control database
df= pd.read_csv(fileInput, sep= sepInput, comment= comInput)
df= df.loc[:, ['gene', 'reads'] ]
#df.columns= ['gene', 'reads']

nInitial= len(df)
nulls= []

#change the names into WB IDs
for i in range(0, len(df)):
    #get the gene associated with locus i
    g= df.iloc[i].gene
    if names[names.GeneName == g].empty == False:
        #make a variable x that extracts the worm ID from the name df. this is a
        #dataframe itself!
        x= names[names.GeneName == g].WBID 
        #Replace the gene name in the control array        
        df= df.replace(g, x.iloc[0])
    else:
        if names[names.WBID == g].empty == True:
            if names[names.HumanReadable == g].empty == False and g.lower() != 'nan':
                x= names[names.HumanReadable == g].WBID
                df= df.replace(g, x)
            else:
                #There's going to be a problem.
                nulls.append([g, df.iloc[0].reads])
                df= df.replace(g, float('nan'))

df= df.dropna()

#Write to main file
df.to_csv(newFile, sep=',', index= False)

#Write error file
f= open(errorFile, 'w')
f.write("#This file stores any genes that could not be translated into a WBID\n")
f.write("gene_Name\,reads\n")
for element in nulls:
    string= str(element[0])+","+str(element[1])+"\n"
    f.write(string)
f.close()

#Write summary
f1= open(summaryFile, 'w')
f1.write("#This file gives a summary of the results of running FileFormatScript.py\n")
n= len(df)
f1.write("Original list contained {0} entries. List with WBIDs contains {1} entries".format(nInitial, n))
f1.close()

