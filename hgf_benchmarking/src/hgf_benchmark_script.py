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
import re
imp.reload(hgt)


pd.set_option('display.float_format', lambda x:'%f'%x) 


#this script generates a few directories. 
dirOutput= '../output/'
dirSummaries= '../output/SummaryInformation/'
dirAnalyses= '../output/Analyses'
dirHGT25_any= '../output/HGT25_any_Results/'
dirHGT25_avg= '../output/HGT25_avg_Results/'
dirHGT50_any= '../output/HGT50_any_Results/'
dirHGT50_avg= '../output/HGT50_avg_Results/'
dirHGT100_any= '../output/HGT100_any_Results/'

DIRS= [dirOutput, dirSummaries, dirHGT25_any, dirHGT25_avg,
       dirHGT50_any, dirHGT50_avg, dirHGT100_any]
#open the relevant file
path_sets= '../input/genesets_golden/'
path_dicts= '../input/WS252AnatomyDictionary/'


#Make all the necessary dirs if they don't already exist
for d in DIRS:
    if not os.path.exists(d):
        os.makedirs(d)

#Make the file that will hold the summaries and make the columns. 
with open(dirSummaries+'ExecutiveSummary.csv', 'w') as fSum:
    fSum.write('#Summary of results from all benchmarks\n')
    fSum.write('NoAnnotations,Threshold,Method,EnrichmentSetUsed,TissuesTested,GenesSubmitted,TissuesReturned,GenesUsed,AvgFold,AvgQ\n')

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
        if f_dict == '.DS_Store':
            continue
        tissue_df= pd.read_csv(path_dicts+f_dict)
        
        #tobedropped when tissue dictionary is corrected
        tissue_df.drop('C. elegans Cell and Anatomy WBbt:0000100', axis= 1, inplace= True)
        annot, thresh= re.findall(r"[-+]?\d*\.\d+|\d+", f_dict)    
        annot= int(annot); thresh= float(thresh) #typecasting
        method= f_dict[-7:-4]
        
        ntiss= len(tissue_df.columns)
        
        #open each enrichment set
        for fodder in os.walk(path_sets):
            for f_set in fodder[2]:
                df= pd.read_csv(path_sets + f_set)
                test= df.gene.values
                ntest= len(test)
                short_name= f_set[16:len(f_set)-16]
                
                df_analysis, unused= \
                hgt.enrichment_analysis(test, tissue_df, alpha= 0.05, show= False)
                
                nana= len(df_analysis) #len of results
                nun= len(unused) #number of genes dropped
                avf= df_analysis['Fold Change'].mean()
                avq= df_analysis['Q value'].mean()
                s= '{0},{1},{2},{3},{4},{5},{6},{7},{8},{9}'.format(
                annot, thresh, method, f_set, ntiss, ntest, nana, ntest-nun, avf,avq)
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

df_summary.sort_values(['NoAnnotations', 'Threshold', 'Method'], inplace= True)
#==============================================================================
#==============================================================================
# # Plot summary graphs
#==============================================================================
#==============================================================================
sel= lambda x, y, z: (df_summary.NoAnnotations == x) & (df_summary.Threshold == y) & (df_summary.Method == z)


#KDE of the fraction of all tissues that tested significant
cols= ['#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e'] #used with varying colors
ls= ['-', '--', ':'] # used with varying thresh
thresh= df_summary.Threshold.unique()
NoAnnotations= df_summary.NoAnnotations.unique()

def resplot(column, method= 'any'):
    for j, annots in enumerate(NoAnnotations):
        for i, threshold in enumerate(thresh):
            df_summary[sel(annots, threshold, method)][column].plot('kde', 
                            color= cols[j], ls= ls[i], lw= 4, 
                            label= 'Annotation Cut-off: {0}, Threshold: {1}'.format(annots, threshold))



resplot('fracTissues')
plt.xlabel('Fraction of all tissues that tested significant')
plt.xlim(0, 1)
plt.title('KDE Curves for all dictionaries, benchmarked on all gold standards')
plt.legend()
plt.savefig(dirSummaries+'fractissuesKDE_method=any.png')
plt.close()

resplot('AvgQ', method= 'avg')
plt.xlabel('Fraction of all tissues that tested significant')
plt.xlim(0, 0.05)
plt.title('KDE Curves for all dictionaries, benchmarked on all gold standards')
plt.legend()
plt.savefig(dirSummaries+'avgQKDE_method=avg.png')
plt.close()


resplot('AvgQ')
plt.xlabel('AvgQ value')
plt.xlim(0,.05)
plt.title('KDE Curves for all dictionaries, benchmarked on all gold standards')
plt.legend()
plt.savefig(dirSummaries+'avgQKDE_method=any.png')
plt.close()

#KDE of the fraction of avgFold
resplot('AvgFold')
plt.xlabel('Avg Fold Change value')
plt.xlim(0, 15)
plt.title('KDE Curves for all dictionaries, benchmarked on all gold standards')
plt.legend()
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
def walker(tissue_df, directory):
   tissue_df.drop('C. elegans Cell and Anatomy WBbt:0000100', axis= 1, inplace= True)
   for fodder in os.walk(path_sets):
       for f_set in fodder[2]:
           df= pd.read_csv(path_sets + f_set)
           short_name= f_set[16:len(f_set)-16]
           test= df.gene.values
           df_analysis, unused= hgt.enrichment_analysis(test, tissue_df, show= False)
           if df.empty == False:
               df_analysis.to_csv(directory+short_name+'.csv')
               line= '#' + short_name+'\n'
               line_prepender(directory+short_name+'.csv', line)
               hgt.plot_enrichment_results(df_analysis, title= short_name, dirGraphs= directory)
           
        
tissue_df= pd.read_csv('../input/WS252AnatomyDictionary/annot25.thresh0.9.methodany.csv')
walker(tissue_df, dirHGT25_any)
               
tissue_df= pd.read_csv('../input/WS252AnatomyDictionary/annot50.thresh0.9.methodany.csv')
walker(tissue_df, dirHGT50_any)


tissue_df= pd.read_csv('../input/WS252AnatomyDictionary/annot100.thresh0.9.methodany.csv')
walker(tissue_df, dirHGT100_any)


tissue_df= pd.read_csv('../input/WS252AnatomyDictionary/annot25.thresh0.9.methodavg.csv')
walker(tissue_df, dirHGT25_avg)


with open('../output/SummaryInformation/TissueNumbers.csv', 'w') as f:
    f.write('No. Of Annotations,Threshold,Method,No. Of Tissues in Dictionary\n')
    for annots in NoAnnotations:
        for threshold in thresh:
            for method in ['any', 'avg']:
                x= df_summary[sel(annots, threshold, method)].TissuesTested.unique()[0]
                f.write('{0},{1},{2},{3}\n'.format(annots, threshold, method, x))





