# -*- coding: utf-8 -*-
"""
Spyder Editor


@david angeles
dangeles@caltech.edu
"""

import hypergeometricTests as hgt #the library to be used
import pandas as pd
import os
import importlib as imp
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
imp.reload(hgt)


pd.set_option('display.float_format', lambda x:'%f'%x) 


#this script generates a few directories. 
dirOutput= '../output/'
dirSummaries= '../output/SummaryInformation/'
dirAnalyses= '../output/Analyses'
dirHGT25= '../output/HGT25Results/'
dirHGT25_any= '../output/HGT25_any_Results/'
dirHGT25_avg= '../output/HGT25_avg_Results/'
dirHGT50= '../output/HGT50Results/'

#open the relevant file
path_sets= '../input/genesets_golden/'
path_dicts= '../input/WS252AnatomyDictionary/'


#Make the dataframe to know how many tissues are being tested
a= [10, 25, 50, 100]
b= [None]*4
tissue_numbers_df= pd.DataFrame(index= a, columns= ['No. of Tissues'])


#Make all the necessary dirs if they don't already exist
if not os.path.exists(dirSummaries):
    os.makedirs(dirSummaries)

if not os.path.exists(dirHGT25):
    os.makedirs(dirHGT25)
    
if not os.path.exists(dirHGT25_any):
    os.makedirs(dirHGT25_any)
    
if not os.path.exists(dirHGT25_avg):
    os.makedirs(dirHGT25_avg)
        
if not os.path.exists(dirHGT50):
    os.makedirs(dirHGT50)

if not os.path.exists(dirAnalyses):
    os.makedirs(dirAnalyses)

#Make the file that will hold the summaries and make the columns. 
with open(dirSummaries+'ExecutiveSummary.csv', 'w') as fSum:
    fSum.write('#Summary of results from all benchmarks\n')
    fSum.write('DictUsed,EnrichmentSetUsed,TissuesTested,GenesTested,TissuesReturned,AvgFold,AvgQ\n')

#==============================================================================
#==============================================================================
# # Perform the bulk of the analysis, run every single dictionary on every set
#==============================================================================
#==============================================================================
i= 0
#look in the dictionaries
for folder in os.walk(path_dicts):
    #open each one
    for f_dict in folder[2]:
        
        tissue_df= pd.read_csv(path_dicts+f_dict)
        dname= int(f_dict[:-4])
        ntiss= len(tissue_df.columns)
        tissue_numbers_df.loc[dname, 'No. of Tissues']= ntiss
        
        #open each enrichment set
        for fodder in os.walk(path_sets):
            for f_set in fodder[2]:
                df= pd.read_csv(path_sets + f_set)
                test= df.gene.values
                ntest= len(test)
                short_name= f_set[16:len(f_set)-16]
                
                df_analysis= \
                hgt.implement_hypergmt_enrichment_tool(short_name, test, tissue_df, alpha= 0.05)
                
                nana= len(df_analysis)                
                avf= df_analysis['Fold Change'].mean()
                avq= df_analysis['Q value'].mean()
                s= '{0},{1},{2},{3},{4},{5},{6}'.format(dname, f_set, ntiss, ntest, nana, avf,avq)
                with open(dirSummaries+'ExecutiveSummary.csv', 'a+') as fSum:
                    fSum.write(s)
                    fSum.write('\n')
    
#Print summary to csv
df_summary= pd.read_csv(dirSummaries+'ExecutiveSummary.csv', comment= '#')

#some entries contain nulls. before I remove them, I can inspect them
df_summary.isnull().any()
indexFold = df_summary['AvgFold'].index[df_summary['AvgFold'].apply(np.isnan)]
indexQ = df_summary['AvgQ'].index[df_summary['AvgQ'].apply(np.isnan)]
df_summary.ix[indexFold[0]]
df_summary.ix[indexQ[5]]


#kill all nulls!
df_summary.dropna(inplace= True)
#calculate fraction of tissues that tested significant in each run
df_summary['fracTissues']= df_summary['TissuesReturned']/df_summary['TissuesTested']


#==============================================================================
#==============================================================================
# # Plot summary graphs
#==============================================================================
#==============================================================================

#KDE of the fraction of all tissues that tested significant
df_summary[df_summary.DictUsed==10]['fracTissues'].plot('kde', label= '10', color= '#1b9e77', lw= 5)
df_summary[df_summary.DictUsed==25]['fracTissues'].plot('kde', label= '25', color= '#d95f02', lw= 5)
df_summary[df_summary.DictUsed==50]['fracTissues'].plot('kde', label= '50', color= '#7570b3', lw= 5)
df_summary[df_summary.DictUsed==100]['fracTissues'].plot('kde', label= '100', color= '#e7298a', lw= 5)
df_summary[df_summary.DictUsed==200]['fracTissues'].plot('kde', label= '200', color= '#66a61e', lw= 5)
plt.xlabel('Fraction of all tissues that tested significant')
plt.xlim(0, 1)
plt.title('KDE Curves for all dictionaries, benchmarked on all gold standards')
plt.legend(['10', '25', '50', '100', '200'])
plt.savefig(dirSummaries+'fractissuesKDE.png')
plt.close()

#KDE of the fraction of avgq
df_summary[df_summary.DictUsed==10]['AvgQ'].plot('kde', label= '10', color= '#1b9e77', lw= 5)
df_summary[df_summary.DictUsed==25]['AvgQ'].plot('kde', label= '25', color= '#d95f02', lw= 5)
df_summary[df_summary.DictUsed==50]['AvgQ'].plot('kde', label= '50', color= '#7570b3', lw= 5)
df_summary[df_summary.DictUsed==100]['AvgQ'].plot('kde', label= '100', color= '#e7298a', lw= 5)
df_summary[df_summary.DictUsed==200]['AvgQ'].plot('kde', label= '200', color= '#66a61e', lw= 5)
plt.xlabel('AvgQ value')
plt.xlim(0, 0.05)
plt.title('KDE Curves for all dictionaries, benchmarked on all gold standards')
plt.legend(['10', '25', '50', '100', '200'])
plt.savefig(dirSummaries+'avgQKDE.png')
plt.close()

#KDE of the fraction of avgFold
df_summary[df_summary.DictUsed==10]['AvgFold'].plot('kde', label= '10', color= '#1b9e77', lw= 5)
df_summary[df_summary.DictUsed==25]['AvgFold'].plot('kde', label= '25', color= '#d95f02', lw= 5)
df_summary[df_summary.DictUsed==50]['AvgFold'].plot('kde', label= '50', color= '#7570b3', lw= 5)
df_summary[df_summary.DictUsed==100]['AvgFold'].plot('kde', label= '100', color= '#e7298a', lw= 5)
df_summary[df_summary.DictUsed==200]['AvgFold'].plot('kde', label= '200', color= '#66a61e', lw= 5)
plt.xlabel('Avg Fold Change value')
plt.xlim(0, 15)
plt.title('KDE Curves for all dictionaries, benchmarked on all gold standards')
plt.legend(['10', '25', '50', '100', '200'])
plt.savefig(dirSummaries+'avgFoldChangeKDE.png')
plt.close()



def line_prepender(filename, line):
    """
    Given a filename, opens it and prepends the line 'line' 
    at the beginning o fthe file
    """
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)
        
        
#==============================================================================
#==============================================================================
# # Detailed analysis of 25 and 50 genes per node dictionaries
#==============================================================================
#==============================================================================
tissue_df= pd.read_csv('../input/WS252AnatomyDictionary/25.csv')
for fodder in os.walk(path_sets):
   for f_set in fodder[2]:
       df= pd.read_csv(path_sets + f_set)
       short_name= f_set[16:len(f_set)-16]
       test= df.gene.values
       df_analysis= hgt.implement_hypergmt_enrichment_tool(test, tissue_df, aname= short_name)
       df_analysis.to_csv(dirHGT25+short_name+'.csv')
       line= '#' + short_name+'\n'
       line_prepender(dirHGT25+short_name+'.csv', line)
       hgt.plotting_and_formatting(df_analysis, ytitle= short_name, dirGraphs= dirHGT25)
               
tissue_df= pd.read_csv('../input/WS252AnatomyDictionary/50.csv')
for fodder in os.walk(path_sets):
   for f_set in fodder[2]:
       df= pd.read_csv(path_sets + f_set)
       short_name= f_set[16:len(f_set)-16]
       test= df.gene.values
       df_analysis= hgt.implement_hypergmt_enrichment_tool(short_name, test, tissue_df)
       df_analysis.to_csv(dirHGT50+short_name+'.csv')
       line= '#' + short_name+'\n'
       line_prepender(dirHGT50+short_name+'.csv', line)
       hgt.plotting_and_formatting(df_analysis, ytitle= short_name, dirGraphs= dirHGT50)

tissue_df= pd.read_csv('../input/WS252AnatomyDictionary/25_any.csv')
for fodder in os.walk(path_sets):
   for f_set in fodder[2]:
       df= pd.read_csv(path_sets + f_set)
       short_name= f_set[16:len(f_set)-16]
       test= df.gene.values
       df_analysis= hgt.implement_hypergmt_enrichment_tool(short_name, test, tissue_df)
       df_analysis.to_csv(dirHGT25_any+short_name+'.csv')
       line= '#' + short_name+'\n'
       line_prepender(dirHGT50+short_name+'.csv', line)
       hgt.plotting_and_formatting(df_analysis, ytitle= short_name, dirGraphs= dirHGT25_any)

tissue_df= pd.read_csv('../input/WS252AnatomyDictionary/25_avg.csv')
for fodder in os.walk(path_sets):
   for f_set in fodder[2]:
       df= pd.read_csv(path_sets + f_set)
       short_name= f_set[16:len(f_set)-16]
       test= df.gene.values
       df_analysis= hgt.implement_hypergmt_enrichment_tool(short_name, test, tissue_df)
       df_analysis.to_csv(dirHGT25_avg+short_name+'.csv')
       line= '#' + short_name+'\n'
       line_prepender(dirHGT50+short_name+'.csv', line)
       hgt.plotting_and_formatting(df_analysis, ytitle= short_name, dirGraphs= dirHGT25_avg)

