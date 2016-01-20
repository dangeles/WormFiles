# -*- coding: utf-8 -*-
"""
Spyder Editor

A script to obtain all gene lengths in c. elegans


@david angeles
dangeles@caltech.edu
"""

import pandas as pd


fname= '../input/Caenorhabditis_elegans.WBcel235.rel79.cdna.all.fa'

wbids=[]
gene_lengths= []

with open(fname, 'r') as f:
    i= 0
    gene= ''
    for line in f:
        if line[0] == '>':
            start= line.find('gene:') + len('gene:')
            end= start + 15
            wbid= line[start:end].rstrip()
            wbids.append(wbid)
            if i != 0:
                gene= gene.rstrip()
                gene_lengths.append(len(gene))
            else:
                i+=1
            gene= ''
        else:
            gene= gene+line.rstrip()
    gene_lengths.append(len(gene))

cols= ['WBID', 'length']
data= list(zip(wbids, gene_lengths))
df= pd.DataFrame(data, columns= cols)

df.to_csv('../output/c_elegans_gene_lengths_PRJNA13758.txt', index = False)