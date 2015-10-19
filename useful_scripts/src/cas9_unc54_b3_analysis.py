# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 17:53:47 2015

@author: davidangeles
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 21:48:17 2015
A script to analyze the survival and brood size data of female P. red worms
injected with cas9.

@author: davidangeles
"""


from __future__ import division, print_function, absolute_import
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import scipy.special
from scipy import stats
from scipy import optimize
from scipy import misc
from datetime import date
import os
import sys
import itertools as itt
import re
from numpy.lib.stride_tricks import as_strided as ast
import pylab as pyl
import os 
import lifelines as lf

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 30}
plt.rc('font', **font)



maxI= 3#last day for which brood sizes are known


strain_made= 'cas-9 unc-54 b3'
k= 0
k2= 0
k3= 0
df_lifespan= pd.read_csv(
'../input/cas9_unc54_injection_survival_rate_06_october_2015.csv', sep= ',')

df_brood_size= pd.read_csv(
'../input/cas9-unc54-b3-broodsize-10-06-2015.csv', sep= ',')
df_bags= pd.read_csv(
'../input/cas9-unc54-b3-bags-10-06-2015.csv', sep= ',')
df_mutants= pd.read_csv(
'../input/cas9-unc54-b3-mutants-10-06-2015.csv', sep= ',')

#==============================================================================
# 
#==============================================================================
#manipulate the data set into the format we need, i.e., total days lived
#and whether the worm was censored or not. 
#-1 = censored, 0= dead, 1= alive
def when_did_it_die(P0, days):
    """
    Given a P0 and an array of days in which it was alive, censored or dead,
    figure out its lifespan and the number of days it was alive'
    
    Note, the lifespan is considered to last up to the last observation.
    I.e., a worm that was observed to live 3 days, and died on day 4
    would have a lifespan of 4 days.     
    
    If a worm wasn't observed to die, it is considered censored for all practical
    purposes
    
    P0= a unique identifier, typically a number, always a scalar or a string
    days= an array or array-like of floats or ints consisting of -1, 0, or 1
    
    returns:
    Worm Identifier, Lifespan, Death Observed(1) or Censored (0)
    """
    counter=0
    for day in days:
        if day == 0: #worm died, add last day and return
#            counter+=1
            return [P0, counter, 1]
        elif day == -1: #worm was censored, add last day and return
#            counter+=1
            return [P0, counter, 0]
        elif day == 1:
            counter+=1
    return [P0, counter, 0]
#==============================================================================
# 
#==============================================================================
lifespan_data= [[]]*df_lifespan.shape[0]
counter= 0
for p0 in df_lifespan.P0_Number:
    life= df_lifespan[df_lifespan.P0_Number == p0].values[0][1:]
    curr_row= when_did_it_die(p0, life)
    lifespan_data[counter]= curr_row
    counter+=1

df_lf= pd.DataFrame(lifespan_data, columns= ['p0', 'Lifespan', 'Event'])


#Kaplan Meier Estimate of Survival
from lifelines import KaplanMeierFitter
kmf = KaplanMeierFitter(alpha= .95)
kmf.fit(df_lf.Lifespan, df_lf.Event) # more succiently, kmf.fit(T,E)


kmf.survival_function_
kmf.median_
ytcs= [0, .5, 1, 1.1]
fig,ax1= plt.subplots(figsize= (16,16))
ax1= kmf.plot(ax= ax1)
ax1.set_ylim(0, 1.1)
ax1.set_title('Survival Curve')
ax1.set_xlabel('Days since Injection')
ax1.set_ylabel('Fraction Surviving')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')
plt.savefig('../output/Kaplan-Meier Curve for {0} injected females.png'.format(strain_made))

#Hazard estimates. First a cumulative hazard, then a corrected, normalized hazard
from lifelines import NelsonAalenFitter
naf = NelsonAalenFitter()

naf.fit(df_lf.Lifespan, df_lf.Event, label= \
'Nelson Aalen Cumulative Hazard Estimate')

#ax1= naf.plot(ax= ax1)

#Smooth to calculate the hazard parameter, Lambda(t)
plt.clf()
fig,ax1= plt.subplots(figsize= (16, 16))

naf.fit(df_lf.Lifespan, df_lf.Event, label= 'Injected P0')
b= 2 #bandwidth
ax1= naf.plot_hazard(ax=ax1, bandwidth= b)
ax1.set_title('Hazard Curve Curve | bandwidth = {0}'.format(b))
ax1.set_xlabel('Days Post Injection')
ax1.set_ylabel('Hazard')
ax1.set_ylim(0, .25)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')
plt.savefig('../output/Nelson Aalen Cumulative Hazard Estimate for {0} injected females.png'.format(strain_made))


plt.clf()
fig, ax1= plt.subplots(figsize= (16, 16))
plt.xlim(0, df_lf.Lifespan.max())
plt.ylim(1,df_lf.p0.max())
plt.xlabel('Days Post-Injection')
plt.ylabel('P0 identifier')
plt.title('Lifespans of Injected Worms')
ax1= lf.plotting.plot_lifetimes(df_lf.Lifespan, df_lf.Event)
fig.savefig('../output/Lifetime Plot for {0} injected females.png'.format(strain_made))


#==============================================================================
#==============================================================================
# # Brood size analysis
#==============================================================================
#==============================================================================
plt.clf()
ax1= plt.subplots(figsize= (16, 16))
ax1= df_brood_size.Day1[df_lifespan.Day1==1].plot(style= 'bo', label= 'Day 1', alpha= 0.8)
df_brood_size.Day2[df_lifespan.Day2==1].plot(style= 'ro', label= 'Day 2', alpha= 0.8)
df_brood_size.Day3[df_lifespan.Day2==1].plot(style= 'go', label= 'Day 3', alpha= 0.8)
ax1.set_ylim(-5, 70)
ax1.set_xlim(-2, 45)
ax1.set_xlabel('P0 identifier')
ax1.set_ylabel('Brood Size')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')
ax1.legend(loc= 'upper right', fancybox=True, framealpha=0.5, numpoints=1)


#==============================================================================
# 
#==============================================================================
fig, ax1= plt.subplots(figsize= (16, 16))
for i in range(1, maxI +1):
    s= 'Day'+str(int(i))
    df_brood_size[s][df_lifespan[s]==1].hist(label= s, alpha= 0.8)
ax1.set_xlabel('Brood Size')
ax1.set_ylabel('Frequency')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')
ax1.legend(loc= 'upper right', fancybox=True, framealpha=0.5, numpoints=1)

#==============================================================================
# 
#==============================================================================
plt.clf()
t= np.empty(maxI)
perctles= np.empty((maxI, 3))
mode= np.empty((maxI))
mean= np.empty((maxI))
for i in range(1, maxI+1):
    s= 'Day'+str(int(i))
    t[i-1]= i
    brood= df_brood_size[s][df_lifespan[s]==1].values
    perctles[i-1, :]= np.percentile(brood, [25, 50, 75])
    mode[i-1]=df_brood_size[s][df_lifespan[s]==1].mode() 
    mean[i-1]=df_brood_size[s][df_lifespan[s]==1].mean() 

#Plot
plt.clf()
fig, ax1= plt.subplots(figsize= (16, 16))
plt.plot(t, mode, 'bv', markersize= 24, label= 'mode')
plt.plot(t, mean, 'rx', markersize= 30, mew= 3, label= 'mean')
plt.errorbar(t, perctles[:,1], yerr= [perctles[:,1]-perctles[:,0], 
                         perctles[:,2]-perctles[:,1]], label= 'median',
                         fmt= 'go', markersize= 24)
ax1.set_ylim(-5, 1.25*np.max(perctles))
ax1.set_xlim(0, maxI+1)
ax1.set_xlabel('Days Post-Injection')
ax1.set_ylabel('Brood Size per P0')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')
ax1.legend(loc= 'upper right', fancybox=True, framealpha=0.5, numpoints=1)
ax1.set_title('Brood Size per P0 vs Time')
plt.savefig('../output/Mean,Mode,Median of Brood Size per {0} P0 through time'.format(strain_made))
#==============================================================================
# 
#==============================================================================

if k== 0:
    temp= list(df_brood_size)
    temp.remove('P0_Number')
#    temp.remove('Day2')
#    temp.remove('Day3')
#    temp.remove('Day4')
#
    df_brood_size['total_brood']= df_brood_size[temp].sum(axis=1)
    df_brood_size['Lifespan']= df_lf['Lifespan']
    df_brood_size['Event']= df_lf['Event']
    k=1

plt.clf()
fig, ax1= plt.subplots(figsize= (16, 16))

for i in range(0, 2):
    if i == 0:
        s= 'Alive'
        col= 'g'
    else:
        s= 'Dead'
        col= 'b'
#    ax1= plt.gca()
    #Add a little gaussian noise to lifespan:
    x= df_brood_size.Lifespan[df_brood_size.Event==i].values
    x_noise= np.random.normal(0, np.max(x)/100, len(x))
    ax1.scatter(x+x_noise,
                df_brood_size.total_brood[df_brood_size.Event==i], s= 250, c= col, alpha= 0.6,
                label= s)
ax1.set_ylim(-5, 1.25*np.max(df_brood_size.total_brood))
ax1.set_xlim(0, int(maxI)+1)
ax1.set_xlabel('P0 Lifespan (Days Post-Injection)')
ax1.set_ylabel('Brood Size per P0')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')
ax1.legend(loc= 'upper right', fancybox=True, framealpha=0.5, numpoints=1, scatterpoints= 1)
ax1.set_title('Brood Size per P0 vs Lifespan of P0')
plt.savefig('../output/BroodSize vs Lifespan {0} injection.png'.format(strain_made))





p0s_in_bags= df_bags.P0_Number.values
total_size= df_brood_size.total_brood\
[df_brood_size.P0_Number.isin(p0s_in_bags)].values


if k2== 0:
    temp= list(df_bags)
    temp.remove('P0_Number')
    df_bags['total_bags']= df_bags[temp].sum(axis=1).values
    df_bags['total_brood']= total_size
    df_bags['Lifespan']= df_lf['Lifespan'][df_lf.p0.isin(p0s_in_bags)].values
    df_bags['Event']= df_lf['Event'][df_lf.p0.isin(p0s_in_bags)].values
    k2=1
    
df_bags['inv_total_brood']= 1/df_bags.total_brood.values  
df_bags['ratio']= df_bags['total_bags']/df_bags['total_brood']




p0s_in_mutants= df_mutants.P0_Number.values
if k3== 0:
    temp= list(df_mutants)
    temp.remove('P0_Number')
    df_mutants['total_mutants']= df_mutants[temp].sum(axis=1).values
    df_mutants['total_brood']= total_size
    df_mutants['Lifespan']= df_lf['Lifespan'][df_lf.p0.isin(p0s_in_mutants)].values
    df_mutants['Event']= df_lf['Event'][df_lf.p0.isin(p0s_in_mutants)].values
    k3=1
    
df_mutants['inv_total_brood']= 1/df_mutants.total_brood.values  
df_mutants['ratio']= df_mutants['total_mutants']/df_mutants['total_brood']

df_mutants= df_mutants.replace([np.inf, -np.inf], np.nan)


xlim_max= 10**3
fig, ax= plt.subplots(figsize = (16, 16))

plt.plot(df_bags['total_brood']+1, df_bags['ratio'], 'o', markersize= 15,\
label= 'bags') #plus one so log plot doesn't freak
plt.plot(df_mutants['total_brood']+1, df_mutants['ratio'], 'o', markersize= 13, label= 'other_mutants') #plus one so log plot doesn't freak


plt.hlines(np.median(df_bags['total_bags']/df_bags['total_brood']),0, xlim_max, 
           linestyles= '-', label= 'median no. bags per p0')

#plt.hlines(np.std(df_bags['total_bags']/df_bags['total_brood'])+
#           np.mean(df_bags['total_bags']/df_bags['total_brood']),0, xlim_max, 
#           linestyles= '--', label= 'average bags per p0')



ax.set_xlim(1, xlim_max)
ax.set_xscale('log')
ax.set_xlabel('Total Brood Size')
ax.set_ylabel('Bags/Total Brood Size')
ax.set_ylim(-.01, .3)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.legend(loc= 'upper right', fancybox=True, framealpha=0.5, numpoints=1, scatterpoints= 1)
ax.set_title('Jackpot Plot',  y=1.08)
plt.savefig('../output/Jackpot plot, {0}.png'.format(strain_made))

pd.set_option('display.float_format', lambda x: '%.2f' % x)
print(df_bags[['total_brood', 'total_bags', 'ratio']].dropna().sort('ratio'))
pd.set_option('display.float_format', lambda x: '%.7f' % x)






#total 'edited' visible larvae: 118
x= df_bags.P0_Number[df_bags.total_bags > 0].values
y= df_mutants.P0_Number[df_mutants.total_mutants>0].values
z= np.append(x, y)
z= np.unique(z)

print(len(z))
print(len(df_bags))
print(len(z)/len(df_bags))



xlim_max= 10**3
fig, ax= plt.subplots(figsize = (16, 16))

plt.plot(df_bags['total_brood']+1, np.sqrt(df_bags['ratio']+df_mutants['ratio'])*df_bags['total_brood'], 'o', markersize= 15,\
label= 'max_edit_freq') #plus one so log plot doesn't freak
plt.plot(df_bags['total_brood']+1, (df_bags['ratio']+df_mutants['ratio'])*df_bags['total_brood'], 'go', markersize= 13,\
label= 'min_edit_freq') #plus one so log plot doesn't freak


ax.set_xlim(1, xlim_max)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Total Brood Size')
ax.set_ylabel('Projected Edited Larvae')
ax.set_ylim(1, xlim_max)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.legend(loc= 'upper right', fancybox=True, framealpha=0.5, numpoints=1, scatterpoints= 1)
ax.set_title('Upper and Lower Bounds of Edit Frequency',  y=1.08)
lims = [
    np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
    np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
]

# now plot both limits against eachother
ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
plt.savefig('../output/Editing Bounds Estimate {0}.png'.format(strain_made))