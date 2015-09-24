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
        'size'   : 22}
plt.rc('font', **font)


k= 0
df_lifespan= pd.read_csv(
'../input/cas9_unc54_injection_survival_rate_13_september_2015.csv', sep= ';')

df_brood_size= pd.read_csv(
'../input/brood_size_cas9_unc54_13_sept_2015.csv', sep= ';')
df_bags= pd.read_csv(
'../input/unc54_cas9_injection_13sept2015_bags_of_worms_found.csv', sep= ';')

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
fig,ax1= plt.subplots(figsize= (16, 16))
ax1= kmf.plot(ax= ax1)
ax1.set_ylim(0, 1.1)
ax1.set_title('Survival Curve')
ax1.set_xlabel('Days since Injection')
ax1.set_ylabel('Fraction Surviving')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')
plt.savefig('../output/Kaplan-Meier Curve for pred-unc-54 cas9 injected females.png')

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
plt.savefig('../output/Nelson Aalen Cumulative Hazard Estimate for pred-unc-54 cas9 injected females.png')


plt.clf()
fig, ax1= plt.subplots(figsize= (16, 16))
plt.xlim(0, df_lf.Lifespan.max())
plt.ylim(1,df_lf.p0.max())
plt.xlabel('Days Post-Injection')
plt.ylabel('P0 identifier')
plt.title('Lifespans of Injected Worms')
ax1= lf.plotting.plot_lifetimes(df_lf.Lifespan, df_lf.Event)
fig.savefig('../output/Lifetime Plot for pred-unc-54 cas9 injected females (generated 21 Sep 2015).png')


#==============================================================================
#==============================================================================
# # Brood size analysis
#==============================================================================
#==============================================================================
maxI= 4 #last day for which brood sizes are known


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
plt.savefig('../output/Mean,Mode,Median of Brood Size per P0 through time')
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
        s= 'Dead'
        col= 'g'
    else:
        s= 'Alive'
        col= 'b'
#    ax1= plt.gca()
    #Add a little gaussian noise to lifespan:
    x= df_brood_size.Lifespan[df_brood_size.Event==i].values
    x_noise= np.random.normal(0, np.max(x)/100, len(x))
    ax1.scatter(x+x_noise,
                df_brood_size.total_brood[df_brood_size.Event==i], s= 100, c= col, alpha= 0.6,\
                label= s)
ax1.set_ylim(-5, 1.25*np.max(df_brood_size.total_brood))
ax1.set_xlim(0, 7)
ax1.set_xlabel('P0 Lifespan (Days Post-Injection)')
ax1.set_ylabel('Brood Size per P0')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')
ax1.legend(loc= 'upper right', fancybox=True, framealpha=0.5, numpoints=1, scatterpoints= 1)
ax1.set_title('Brood Size per P0 vs Lifespan of P0')
plt.savefig('../output/BroodSize vs Lifespan, cas9 pred-unc54 injection.png')




brooddays= df_brood_size.Day1[df_brood_size.P0_Number.isin(df_bags.P0_Number.values)].values\
 + df_brood_size.Day2[df_brood_size.P0_Number.isin(df_bags.P0_Number.values)].values
total_size= df_brood_size.total_brood\
[df_brood_size.P0_Number.isin(df_bags.P0_Number.values)].values+1 #plus one so log plot doesn't freak
bags= df_bags.Day1+df_bags.Day2
fig, ax= plt.subplots(figsize = (8, 8))
plt.plot(total_size, bags/total_size, 'o', markersize= 15)
ax.set_xscale('log')
ax.set_ylim(-.1, np.max(bags/total_size)+.1)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
