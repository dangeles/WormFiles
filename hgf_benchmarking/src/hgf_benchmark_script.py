# -*- coding: utf-8 -*-
"""
Spyder Editor


@david angeles
dangeles@caltech.edu
"""

import hypergeometricTests as hgt
import pandas as pd
import os
import importlib as imp
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
imp.reload(hgt)


pd.set_option('display.float_format', lambda x:'%f'%x) 

dirOutput= '../output/'
dirSummaries= '../output/SummaryInformation/'
dirAnalyses= '../output/Analyses'

#open the relevant file
path_sets= '../input/genesets_golden/'
path_dicts= '../input/WS252AnatomyDictionary/'


a= [10, 25, 50, 100]
b= [None]*4
tissue_numbers_df= pd.DataFrame(index= a, columns= ['No. of Tissues'])


with open(dirSummaries+'ExecutiveSummary.csv', 'w') as fSum:
    fSum.write('#Summary of results from all benchmarks\n')
    fSum.write('DictUsed,EnrichmentSetUsed,TissuesTested,GenesTested,TissuesReturned,AvgFold,AvgQ\n')
    


if not os.path.exists(dirSummaries):
    os.makedirs(dirSummaries)
    
if not os.path.exists(dirAnalyses):
    os.makedirs(dirAnalyses)

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
    
#df_analysis.to_csv('../output/'+f+f_dict[:-4])
df_summary= pd.read_csv(dirSummaries+'ExecutiveSummary.csv', comment= '#')

#some entries contain nulls. before I remove them, I can inspect them
df_summary.isnull().any()
indexFold = df_summary['AvgFold'].index[df_summary['AvgFold'].apply(np.isnan)]
indexQ = df_summary['AvgQ'].index[df_summary['AvgQ'].apply(np.isnan)]
df_summary.ix[indexFold[0]]
df_summary.ix[indexQ[5]]
#kill all nulls!
df_summary.dropna(inplace= True)

df_summary['fracTissues']= df_summary['TissuesReturned']/df_summary['TissuesTested']
bins= np.floor(np.sqrt(len(df_summary[df_summary.DictUsed==10])))

groupbydict= df_summary.groupby('DictUsed')

groupbydict.plot(x= 'fracTissues', kind= 'kde')

df_summary[df_summary.DictUsed==10].hist('TissuesReturned', bins= bins)
df_summary[df_summary.DictUsed==25].hist('TissuesReturned', bins= bins)
df_summary[df_summary.DictUsed==50].hist('TissuesReturned', bins= bins)
df_summary[df_summary.DictUsed==100].hist('TissuesReturned', bins= bins)
df_summary[df_summary.DictUsed==200].hist('TissuesReturned', bins= bins)

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


df_summary[df_summary.DictUsed==10]['TissuesTested'].unique()
df_summary[df_summary.DictUsed==25]['TissuesTested'].unique()
df_summary[df_summary.DictUsed==50]['TissuesTested'].unique()
df_summary[df_summary.DictUsed==100]['TissuesTested'].unique()
df_summary[df_summary.DictUsed==200]['TissuesTested'].unique()

df_summary[df_summary.DictUsed==10].hist('AvgQ', bins= bins)
df_summary[df_summary.DictUsed==25].hist('AvgQ', bins= bins)
df_summary[df_summary.DictUsed==50].hist('AvgQ', bins= bins)
df_summary[df_summary.DictUsed==100].hist('AvgQ', bins= bins)
df_summary[df_summary.DictUsed==200].hist('AvgQ', bins= bins)
        
singleton= ['WBGene00001187']
tissue_df= pd.read_csv('../input/dictionary.csv')

#x='syncytium(WBbt:0008074)'
#n= 0
#print(stats.hypergeom.sf(n, tg, s[x], 1))
#print(scipy.misc.comb(s[x], n) * scipy.misc.comb(tg - s[x], 1-n) / scipy.misc.comb(tg, 1))

