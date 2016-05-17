# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 10:47:12 2016

A script to find the overlap between differentially expressed genes in fly, Zebrafish
and chick. Returns the overlapping genes in flybase id format.

@author: dangeles
"""
import pandas as pd
import numpy as np
import os

print(os.getcwd())
df_fish1= pd.read_csv('../input/48fish_to_fly_refseqID_ensembl.txt')
df_fish2= pd.read_csv('../input/36fish_to_fly_refseqID_ensembl.txt')
df_fish3= pd.read_csv('../input/48fish_to_fly_predictedID_ensembl.txt')
df_fish4= pd.read_csv('../input/36fish_to_fly_predictedID_ensembl.txt')

df_fly= pd.read_csv('../input/border_cell_data.csv')
df_chick= pd.read_csv('../input/chick_to_fly_ensembl.txt')

print(df_fly.head())
#drop all nans from everything
df_fish1.dropna(inplace= True)
df_fish2.dropna(inplace= True)
df_fish3.dropna(inplace= True)
df_fish4.dropna(inplace= True)

df_fly.dropna(inplace= True)
df_chick.dropna(inplace= True)


#fly and chick may have repeated fish IDS. drop them
df_fly.drop_duplicates('Acc#', inplace= True)
df_chick.drop_duplicates('Fruitfly Ensembl Gene ID', inplace= True)
df_fish1.drop_duplicates('Fruitfly Ensembl Gene ID', inplace= True)
df_fish2.drop_duplicates('Fruitfly Ensembl Gene ID', inplace= True)
df_fish3.drop_duplicates('Fruitfly Ensembl Gene ID', inplace= True)
df_fish4.drop_duplicates('Fruitfly Ensembl Gene ID', inplace= True)

#gather all the ensemble transcripts into a single list
Ldff= [df_fish1, df_fish2, df_fish3, df_fish4]

overlap_all= np.array([])
for df in Ldff:
    genes= df['Fruitfly Ensembl Gene ID'].values
    if len(overlap_all) == 0:
        overlap_all= genes
    else:
        overlap_all= np.concatenate([overlap_all, genes])

#remove any overlaps within fish -- these are spurious results
fish_genes= overlap_all.copy()
overlap_all= np.unique(overlap_all)
overlap_fly_fish= overlap_all
overlap_chick_fish= np.array([])

#append chick ids
genes_chick= df_chick['Fruitfly Ensembl Gene ID'].values
overlap_all= np.concatenate([overlap_all, genes_chick])

#append fly IDS
genes_fly= df_fly['Fruitfly Ensembl Gene ID'].values
overlap_all= np.concatenate([overlap_all, genes_fly])

unq, unq_idx, unq_cnt = np.unique(overlap_all, return_inverse=True, return_counts=True)
unq_cnt = np.bincount(unq_idx)
cnt_mask = unq_cnt > 2
unique_overlapped_all = unq[cnt_mask]


unq, unq_idx, unq_cnt = np.unique(fc_overlap, return_inverse=True, return_counts=True)
unq_cnt = np.bincount(unq_idx)
cnt_mask = unq_cnt > 1
unique_fly_chick_overlap = unq[cnt_mask]

unq, unq_idx, unq_cnt = np.unique(overlap_fly_fish, return_inverse=True, return_counts=True)
unq_cnt = np.bincount(unq_idx)
cnt_mask = unq_cnt > 1
unique_fly_fish_overlap = unq[cnt_mask]


unq, unq_idx, unq_cnt = np.unique(overlap_chick_fish, return_inverse=True, return_counts=True)
unq_cnt = np.bincount(unq_idx)
cnt_mask = unq_cnt > 1
unique_overlap_chick_fish = unq[cnt_mask]

ovr= [unique_overlapped_all, unique_fly_chick_overlap, unique_fly_fish_overlap, overlap_chick_fish]

#write to files:
with open('../output/overlap_all.txt', 'w') as f:
    f.write('This file contains the genes that overlapped in all comparisons\n')
    for value in unique_overlapped_all:
        f.write(value+'\n')

with open('../output/overlap_fish_fly.txt', 'w') as f:
    f.write('This file contains the genes that overlapped between fish and fly\n')
    for value in unique_fly_fish_overlap:
        f.write(value+'\n')

with open('../output/overlap_fish_chick.txt', 'w') as f:
    f.write('This file contains the genes that overlapped between fish and chick\n')
    for value in unique_overlap_chick_fish:
        f.write(value+'\n')

with open('../output/overlap_chick_fly.txt', 'w') as f:
    f.write('This file contains the genes that overlapped between fly and chick\n')
    for value in unique_fly_chick_overlap:
        f.write(value+'\n')



print('hi')
print('goodbye')
