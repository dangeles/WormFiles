# -*- coding: utf-8 -*-
"""
Created on Sat Nov  7 18:07:32 2015

A script to analyze motion data from single worm track data

Accepts: WormTracker data. But be warned! You must manually re-write all NaNs
because the script that generates the files inserts random spaces in there.

Outputs: Bunch of graphs and summary of the statistics. 

@author: davidangeles
"""

import pandas as pd
import numpy as np
import scipy.stats as scs
import matplotlib.pyplot as plt
import time
## dd/mm/yyyy format

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 30}
plt.rc('font', **font)

worm_name= 'C. elegans hermaphrodites'

new_names= ['worm','rpm','DurF',	'DurB','DurFB','DurS','DistF','DistB','DistFB','SpdF','SpdB','SpdFB','Revs']
units= {'rpm': 'reversals/minute', 
        'DurB': 's',
        'DurF': 's',
        'DurFB': 's',
        'DurS': 's',
        'DistF': 'mm',
        'DistB': 'mm',
        'DistFB': 'mm',
        'SpdF': 'mm/s',
        'SpdB': 'mm/s',
        'SpdFB': 'mm/s',
        'Revs': 'total reversals'}
        
df= pd.read_csv('../input/n2_motion.csv', na_values= 'NaN')

old_names= df.columns.values
df.rename(columns=dict(zip(old_names, new_names)), inplace=True)


#No worm comparison
for name in new_names[1:]:
    
    data= df[name].dropna()
    
    if len(data) > 7:
        x= scs.gaussian_kde(df[name].values)
        
        xaxis = np.linspace(df[name].min()-df[name].std(), 
                            df[name].max()+df[name].std(),
                            len(df[name])*2)
                            
        plt.close()
        fig,ax= plt.subplots(figsize= (16,16))
        plt.plot(xaxis, x.evaluate(xaxis), linewidth= 2)
        ax.set_title('Gaussian KD Estimate for ' + worm_name + ' ' + name, y=1.05)
        ax.set_xlabel(name + ', '+ units[name])
        ax.set_ylabel('Normalized Density, '+ units[name]+'^-1')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.set_xlim(xaxis[0], xaxis[len(xaxis)-1])
        plt.savefig('../output/'+worm_name+'_'+name+time.strftime("%d_%m_%Y")+'.png')
        



df.describe().to_csv('../output/' + worm_name + '_summary_'+time.strftime("%d_%m_%Y")+'.csv')







        