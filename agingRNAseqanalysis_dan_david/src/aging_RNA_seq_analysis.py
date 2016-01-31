# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 08:56:36 2016

@author: davidangeles
"""

# -*- coding: utf-8 -*-


import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import hypergeometricTests as hgt
import os


mag= 2 #value of beta from regression
qval= .1 #qvalue from regression
qvalEn= 0.05 #q value for enrichment analysis (tissues)
                
dirLists= '../output/Gene_lists_for_analysis'
if not os.path.exists(dirLists):
    os.makedirs(dirLists)
    
dirGraphs= '../output/Graphs'
if not os.path.exists(dirLists):
    os.makedirs(dirGraphs)
    

os.chdir('./')
#gene_lists from sleuth
dfBetaA= pd.read_csv("../input/table_agebeta_genes.csv")
dfBetaG= pd.read_csv("../input/table_genobeta_genes.csv")
dfBetaAG= pd.read_csv("../input/table_genocrossagebeta_genes.csv")

dfDaf12= pd.read_csv('../input/daf12genes.csv')
dfDaf16= pd.read_csv('../input/daf16genes.csv')

dfLund= pd.read_csv('../input/lund_data.csv', header= None, names=['gene'])
dfEckley= pd.read_csv('../input/eckley_data.csv', header= None, names=['gene'])
dfMurphyUp= pd.read_csv('../input/murphy_data_lifespan_extension.csv')
dfMurphyDown= pd.read_csv('../input/murphy_data_lifespan_decrease.csv')


dfDaf12['origin']= 'daf-12'
dfDaf16['origin']= 'daf-16'
dfEckley['origin']= 'Eckley'
dfLund['origin']= 'Lund'
dfMurphyUp['origin']= 'MurphyExt'
dfMurphyDown['origin']= 'MurphyDec'
frames= [dfDaf12, dfDaf16, dfEckley, dfLund, dfMurphyDown, dfMurphyUp]
#frames= [dfEckley, dfLund]
dfGoldStandard = pd.concat(frames)

#from wormbase
dfLifespanGenes= pd.read_csv('../input/lifespan gene list complete.csv')

#dfPAN= pd.read_csv()

#tissue dictionary-- please cite David Angeles et al TEA publication (forthcoming)
#if using the enrichment tool 
tissue_df= pd.read_csv("../input/dictionary.csv")

#slice all the relevant gene names out
#remove all isoforms!
namesBetaA= \
        dfBetaA[\
        (dfBetaA.qval < qval) & ((dfBetaA.b > mag))]\
        .ens_gene.unique()
namesBetaG= \
        dfBetaG[\
        (dfBetaG.qval < qval) & ((dfBetaG.b > mag))]\
        .ens_gene.unique()
namesBetaAG= \
        dfBetaAG[\
        (dfBetaAG.qval < qval) & ((dfBetaAG.b > mag))]\
        .ens_gene.unique()
        
namesBetaAneg= \
        dfBetaA[\
        (dfBetaA.qval < qval) & ((dfBetaA.b < -mag))]\
        .ens_gene.unique()
namesBetaGneg= \
        dfBetaG[\
        (dfBetaG.qval < qval) & ((dfBetaG.b < -mag))]\
        .ens_gene.unique()
namesBetaAGneg= \
        dfBetaAG[\
        (dfBetaAG.qval < qval) & ((dfBetaAG.b < -mag))]\
        .ens_gene.unique()
        

namesA0= \
        dfBetaA[\
        (dfBetaA.qval < qval) & (dfBetaAG.b > 0)]\
        .ens_gene.unique()
namesG0= \
        dfBetaG[\
        (dfBetaG.qval < qval)& (dfBetaAG.b > 0)]\
        .ens_gene.unique()
namesAG0= \
        dfBetaAG[\
        (dfBetaAG.qval < qval)& (dfBetaAG.b > 0)]\
        .ens_gene.unique()
        
        
namesA0neg= \
        dfBetaA[\
        (dfBetaA.qval < qval) & (dfBetaAG.b < 0)]\
        .ens_gene.unique()
namesG0neg= \
        dfBetaG[\
        (dfBetaG.qval < qval)& (dfBetaAG.b < 0)]\
        .ens_gene.unique()
namesAG0neg= \
        dfBetaAG[\
        (dfBetaAG.qval < qval)& (dfBetaAG.b < 0)]\
        .ens_gene.unique()

#place all the relevant gene lists in this array
array_of_arrays= [namesBetaA,
                  namesBetaG,
                  namesBetaAG,
                  namesBetaAneg,
                  namesBetaGneg,
                  namesBetaAGneg,
                  namesA0,
                  namesG0,
                  namesAG0,
                  namesA0neg,
                  namesG0neg,
                  namesAG0neg
                  ]




def wbid_extractor(tissue_df, main_tissue):
    """
    Given a string 'main_tissue', find all columns that
    have a substring equal to it in tissue_df. Then,
    extract all the wbids that are expressed in any of these
    columns and return a non-redundant list of these wbids
    """
    if type(main_tissue) != str:
        raise ValueError('please input a string in main tissue')
        
    matching = [s for s in tissue_df.columns if main_tissue in s]
    names= []
    for i in matching:
        if len(names) == 0:
            names= tissue_df.wbid[tissue_df[i]==1].values
        else:
            names= np.append(names, tissue_df.wbid[tissue_df[i]==1].values)
    names= list(set(names))
    
    return names


def organize(names, tissue_df):
    """
    Pick your favourite tissues, place them in a list
    provide the tissue dictionary and they will be assembled
    into a dataframe of all tissues that include that word 
    """
    #guarantee its iterable
    if type(names) in [str]:
        names= [names]
    
    names= ['wbid']+names
        
    df= pd.DataFrame(index= np.arange(0, len(tissue_df)), columns= names)    
    df.wbid= tissue_df.wbid
    for i, value in enumerate(names):
        genes_in_tissue= wbid_extractor(tissue_df, value)
        df[value][df.wbid.isin(genes_in_tissue)]= 1

    df.fillna(0, inplace= True)
    df1= pd.melt(df, id_vars='wbid', value_vars=tissues)
    return df1

def fix_axes(**kwargs):
    """
    Makes modifications to axes to ensure proper plotting.
    """ 
    title= kwargs.pop('title',  '')
    savename= kwargs.pop('savename', '')
    xlab= kwargs.pop('xlab', r'$\beta$')
    ylab= kwargs.pop('ylab', r'log$_{10}Q$')
    yscale= kwargs.pop('yscale', 'log')
    xscale= kwargs.pop('xscale', 'symlog')
    xlim= kwargs.pop('xlim', [])
    ylim= kwargs.pop('xlim', [])
    loc= kwargs.pop('loc', 1)
    
    plt.legend(fontsize= 12, loc= loc)
    plt.title(title)
    plt.xlabel(xlab, fontsize= 15)
    plt.ylabel(ylab, fontsize= 15)
    plt.xscale(xscale)
    plt.yscale(yscale)
    
    if len(xlim) != 0:
        plt.xlim(xlim)
    if len(ylim) != 0:
        plt.ylim(ylim)
    
    if savename:
        plt.savefig(savename)

def volcano_plot_tissue(tissue, q, dfplot, dfindex, ax, label, col= 'b', a= .8):
    """
    Plots all the tissue specific genes,i.e. all genes that appear in one and only
    one 'tissue'
    """
    g= lambda x:((dfindex.value == 1) & (dfindex.variable == x))\
        & (~dfindex[dfindex.value == 1].duplicated('wbid')) 
    f= lambda x: (dfplot.ens_gene.isin(x)) & (dfplot.qval < q)
    
    gene_selection= g(tissue)
    genes_to_plot= dfindex[gene_selection].wbid
    
    ind= f(genes_to_plot)
    x= dfplot[ind].b
    y= dfplot[ind].qval
    plt.gca().plot(x, -np.log10(y), 'o', color= col, ms= 6, alpha= a, label= label)    


def explode(q, dfvals, dftiss, colors, **kwargs):
    """
    A function that generates all the relevant volcano plots
    """
    a= kwargs.pop('a', .7)
    
    ind1= (dftiss.value==0) | (dftiss[dftiss.value==1].duplicated('wbid'))
    ind2= (dfvals.ens_gene.isin(dftiss[ind1].wbid)) & (dfvals.qval < q)
    
    xnotisssig= dfvals[ind2].b
    ynotisssig= dfvals[ind2].qval
    
    fig, ax= plt.subplots()
    plt.plot(xnotisssig, -np.log10(ynotisssig), 'o', \
    color=colors[0], ms=6, alpha= a, label= 'all others')
        
    #plot all the points not associated with a tissue
    for i, value in enumerate(tissues):
        volcano_plot_tissue(value, q, dfvals, dftiss, label= value,\
        col= colors[i+2], ax= ax, a=a)
        
def volcano_plot_goldstandards(q, dfplot, dfgenes, ax, colors, a= 1, **kwargs):
    """
    Plots all the tissue specific genes,i.e. all genes that appear in one and only
    one 'tissue'
    """
    f= lambda x: (dfplot.ens_gene.isin(x))# & (dfplot.qval < q)    
    
    nvals= len(dfgenes.origin.unique())
    ncolors= len(colors)
    if  nvals > ncolors:
        raise ValueError('Please provide as many colors as there are datasets. {0} {1}'
        .format(ncolors, nvals))
    
    for i, origin in enumerate(dfgenes.origin.unique()):
        
        ind= f(dfgenes[dfgenes.origin == origin].gene.values)
        x= dfplot[ind].b
        y= dfplot[ind].qval
        
        ngsig= len(dfplot[ind & (dfplot.qval < .1)].ens_gene.unique()) #no. of genes showing up in assay
        tg= len(dfgenes[dfgenes.origin == origin].gene.unique()) #no. of genes in dataset
        label= '{0} {1}= {2}, tot= {3}'.format(origin, r'$n_{sig}$', ngsig, tg)

        plt.gca().plot(x, -np.log10(y), 'o', color= colors[i], ms= 6, alpha= a, label= label)  
    
def explode_goldstandards(q, dfplot, dfgenes, colors, **kwargs):
    """
    A function that generates all the relevant volcano plots
    """
    a= kwargs.pop('a', .7)
    loc= kwargs.pop('loc', 'lower right')
    
    ind1= (~dfplot.ens_gene.isin(dfgenes.gene)) & (dfplot.qval < q) #sig genes not in any given df
    ind2= (~dfplot.ens_gene.isin(dfgenes.gene)) & (dfplot.qval > q) #nonsig genes
    
    xnotsig= dfplot[ind2].b
    ynotsig= dfplot[ind2].qval
    
    xsig= dfplot[ind1].b
    ysig= dfplot[ind1].qval
    
    nnotsig= len(dfplot[ind2].ens_gene.unique())
    fig, ax= plt.subplots()
    plt.plot(xnotsig, -np.log10(ynotsig), 'o', \
    color=colors[0], ms=6, alpha= .1, label= 'not sig n= {0}'.format(nnotsig))
    
    nsig= len(dfplot[ind2].ens_gene.unique())
    plt.plot(xsig, -np.log10(ysig), 'o', \
    color=colors[1], ms=6, alpha= .2, label= 'sig n= {0}'.format(nsig))
    
    #plot all the points not associated with a tissue
    volcano_plot_goldstandards(q, dfplot, dfgenes, colors= colors[2:], ax= ax, a=a)
    
    fix_axes(loc= loc, **kwargs)
    leg= ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    leg.get_frame().set_facecolor('#00FFCC')

def kde_tissue(tissue, q, dfplot, dfindex, ax, label, col= 'b'):
    """
    Plots all the tissue specific genes,i.e. all genes that appear in one and only
    one 'tissue'
    """
    g= lambda x:((dfindex.value == 1) & (dfindex.variable == x))\
       # & (~dfindex[dfindex.value == 1].duplicated('wbid')) 
    f= lambda x: (dfplot.ens_gene.isin(x)) & (dfplot.qval < q)
    
    gene_selection= g(tissue)
    
    
    genes_to_plot= dfindex[gene_selection].wbid
    
    ind= f(genes_to_plot)
    x= dfplot[ind].b
    if len(x) > 15:
        sns.kdeplot(x, color= col,label= label+' n= {0}'.format(len(x)), ax= ax, 
                    lw= 5, cut=0.5)        
        if len(x) <= 20:
            sns.rugplot(x, color= col, ax= ax, height= .07, lw= 2)

def kegg(q, dfvals, dftiss, colors, **kwargs):
    """
    A function that generates all the relevant volcano plots
    """
    
    #set scale parameters
    yscale= kwargs.pop('yscale', 'linear')
    xscale= kwargs.pop('xscale', 'linear')
    xlim= kwargs.pop('xlim', [-8,8])
    ylim= kwargs.pop('ylim', [0, .5])
    
    ind1= (dftiss.value==0)# | (dftiss[dftiss.value==1].duplicated('wbid'))
    ind2= (dfvals.ens_gene.isin(dftiss[ind1].wbid)) & (dfvals.qval < q)
    
    xnotisssig= dfvals[ind2].b
    
    fig, ax= plt.subplots()
    sns.kdeplot(xnotisssig, \
        color=colors[0], label= 'all others n= {0}'.format(len(xnotisssig)), ax= ax, \
        lw= 5, cut= 0.5)
    plt.axvline(0, ls= '--', color= 'black', lw= 3)
    #plot all the points not associated with a tissue
    for i, value in enumerate(tissues):
        kde_tissue(value, .1, dfvals, dftiss, label= value,\
        col= colors[i+2], ax= ax)
        
    fix_axes(xscale= xscale, yscale= yscale, xlim= xlim, ylim= ylim, **kwargs)


def kde_value(value, q, dfplot, dfindex, ax, label, col= 'b', min_length= 10, rug_length= 20):
    """
    Plots all the value specific genes,i.e. all genes that appear in one and only
    one 'tissue'
    """
    g= (dfindex.effect == value)
    f= lambda x: (dfplot.ens_gene.isin(x)) & (dfplot.qval < q)
        
    
    genes_to_plot= dfindex[g].gene
    
    ind= f(genes_to_plot)
    x= dfplot[ind].b
    
    if len(x) > min_length:
        sns.kdeplot(x, color= col,label= label+' n= {0}'.format(len(x)), ax= ax, 
                    lw= 5, cut=0.5)    

        if len(x) < rug_length:
            sns.rugplot(x, color= col, ax= ax, height= .1, lw= 2)
    else:
        print('too few values to plot {0}'.format(label+' n= {0}'.format(len(x))))

def kegg_compare_byval(value, q, Ldf, dfindex, colors, **kwargs):
    """
    Given a list of dataframes, Ldf, and a list of target genes dfindex, compare the 
    distributions of genes within dfindex throughout every list for genes with
    trait 'value'. 
    Ldf= a list of dataframes. must have df.ens_gene, df.b, df.qval exactly as written
    dfindex= a list of genes to select, must have df.gene, df.effect
    value= a trait associated with genes in dfindex, (an entry in df.effect)
    colors= an array of len(Ldf) of colors
    """
    dfnames= kwargs.pop('dfnames', ['']*len(Ldf))
    xlim= kwargs.pop('xlim', [-10,10])
    ylim= kwargs.pop('ylim', [0, 1])
    zeroline= kwargs.pop('zeroline', True)
    xscale= kwargs.pop('xscale', 'linear')
    yscale= kwargs.pop('yscale', 'linear')
    
    if len(Ldf) < len(colors):
        raise ValueError('Please provide as many colors as dataframes')
    
    if len(dfindex[dfindex.effect == value]) == 0:
        raise ValueError('Value \'{0}\' is not containted within dfindex'.format(value))
    
    if dfnames:
        if len(Ldf) != len(dfnames):
            raise ValueError('dfnames must be the same length as Ldf')
    
    fig, ax= plt.subplots()
    for i, df in enumerate(Ldf):
        kde_value(value, q, df, dfindex, ax, dfnames[i], colors[i])
    
    if zeroline:
        plt.axvline(0, ls= '--', color= 'black', lw= 2.5)
        
    fix_axes(xlim= xlim, ylim= ylim, xscale= xscale, yscale= yscale, **kwargs)


def kegg_compareall_byval(q, Ldf, dfindex, colors, **kwargs):
    """
    Given a list of dataframes Ldf, and a list of selection genes with criteria
    make all the plots of interest
    """
    vals= dfindex.effect.unique()
    savenames= kwargs.pop('savenames', ['']*len(vals))
    titles= kwargs.pop('titles', ['']*len(vals))
    
    if len(titles) <= len(vals):
        raise ValueError('There are not enough titles for plots'
        .format(len(titles), len(vals)))
        
    if len(savenames) != len(vals):
        raise ValueError('list savenames ({0}) must be the same length as Ldf ({1})'
        .format(len(savenames), len(Ldf)))    
        
    for i, value in enumerate(vals):
        savename= savenames[i]
        title= titles[i]
        kegg_compare_byval(value, q, Ldf, dfindex, colors, savename= savename, 
                           title= title, **kwargs)
    
#==============================================================================
# Tissue Enrichment Batch Analysis
#==============================================================================
aname0= 'Age Beta> {0}'.format(mag)
aname1= 'Genotype Beta> {0}'.format(mag)
aname2= 'Age::Genotype Beta> {0}'.format(mag)

aname3= 'Age Beta< -{0}'.format(mag)
aname4= 'Genotype Beta< -{0}'.format(mag)
aname5= 'Age::Genotype Beta< -{0}'.format(mag)

aname6= 'Age Beta> 0'
aname7= 'Genotype Beta> 0'
aname8= 'Age::Genotype Beta> 0'

aname9= 'Age Beta< 0'
aname10= 'Genotype Beta< 0'
aname11= 'Age::Genotype Beta< 0'
array_of_anames= [aname0,aname1,aname2,
                  aname3,aname4,aname5,
                  aname6,aname7,aname8,
                  aname9,aname10,aname11]

fname0= 'EnrichmentAnalysisAge_BetaGreaterThan{0}_qval_{1}.csv'.format(mag, qvalEn)
fname1= 'EnrichmentAnalysisGenotype_BetaGreaterThan{0}_qval_{1}.csv'.format(mag, qvalEn)
fname2= 'EnrichmentAnalysisAgeCrossGenotype_BetaGreaterThan{0}_qval_{1}.csv'.format(mag, qvalEn)

fname3= 'EnrichmentAnalysisAgePositive_BetaLessThanNegative{0}_qval_{1}.csv'.format(mag, qvalEn)
fname4= 'EnrichmentAnalysisGenotype_BetaLessThanNegative{0}_qval_{1}.csv'.format(mag, qvalEn)
fname5= 'EnrichmentAnalysisAgeCrossGenotype_BetaLessThanNegative{0}_qval_{1}.csv'.format(mag, qvalEn)

fname6= 'EnrichmentAnalysisAge_AllPosBeta_qval_{0}.csv'.format(qvalEn)
fname7= 'EnrichmentAnalysisGenotype_AllPosBeta_qval_{0}.csv'.format(qvalEn)
fname8= 'EnrichmentAnalysisAgeCrossGenotype_AllPosBeta_qval_{0}.csv'.format(qvalEn)

fname9= 'EnrichmentAnalysisAge_AllNegBeta_qval_{0}.csv'.format(qvalEn)
fname10= 'EnrichmentAnalysisGenotype_AllNegBeta_qval_{0}.csv'.format(qvalEn)
fname11= 'EnrichmentAnalysisAgeCrossGenotype_AllNegBeta_qval_{0}.csv'.format(qvalEn)
array_of_fnames= [fname0, fname1, fname2,
                  fname3, fname4, fname5, 
                  fname6, fname7, fname8, 
                  fname9, fname10, fname11]               
array_of_strings= ['namesBetaA', 'namesBetaG', 'namesBetaAG',
                'namesBetaAneg', 'namesBetaGneg', 'namesBetaAGneg',
                'namesA0', 'namesG0', 'namesAG0',
                'namesA0neg', 'namesG0neg', 'namesAG0neg'
                ]
#run the enrichment analysis  for each dataset
for i, list_of_genes in enumerate(array_of_arrays):        
    
    #run the tissue enrichment tool 
    aname= array_of_anames[i] #name of the analysis
    fname= array_of_fnames[i] #filename to store as
    
    print(aname)
    print(len(list_of_genes))    
    df_analysis= \
    hgt.implement_hypergmt_enrichment_tool(aname, list_of_genes,\
                                    tissue_df, qvalEn, f_unused= fname)
    with open('../output/EnrichmentAnalysisResults/'+fname, 'w') as f:
        f.write(aname)
        df_analysis.to_csv('../output/EnrichmentAnalysisResults/'+fname)
    
    hgt.plotting_and_formatting(df_analysis, ytitle= aname)
    #close
#    plt.close()
    
#    #save genes with betas significantly different from zero in a list format
#    k= array_of_strings[i][5:]
#            
#    filename= dirLists+'/'+ 'gene_list_{0}_beta{1}_q{2}.csv'.format(k, mag, qvalEn)
#    #print('filename', filename)
#    with open(filename, 'w') as f:
#        for gene in list_of_genes:
#            f.write(gene)
#            f.write('\n')
#    f.close()
#==============================================================================
# Analysis of some gold standard sets within our data
#==============================================================================
#Daf-12 associated genes
ndaf12= dfDaf12.shape[0]
ndaf16= dfDaf16.shape[0]
nA= dfBetaA[dfBetaA.qval < .1].shape[0]
nG= dfBetaG[dfBetaG.qval < .1].shape[0]
nAG= dfBetaAG[dfBetaAG.qval < .1].shape[0]
#da12 and daf16 genes with sig age betas
indA12= (dfBetaA.ens_gene.isin(dfDaf12.gene)) & (dfBetaA.qval < .1)
indA16= (dfBetaA.ens_gene.isin(dfDaf16.gene)) & (dfBetaA.qval < .1)
ndaf12A= dfBetaA.b[indA12].shape[0]
ndaf16A= dfBetaA.b[indA16].shape[0]
#da12 and daf16 genes with sig genotype betas
indG12= (dfBetaG.ens_gene.isin(dfDaf12.gene)) & (dfBetaG.qval < .1)
indG16= (dfBetaG.ens_gene.isin(dfDaf16.gene)) & (dfBetaG.qval < .1)
ndaf12G= dfBetaG.b[indG12].shape[0]
ndaf16G= dfBetaG.b[indG16].shape[0]
#da12 and daf16 genes with sig age:: genotype betas
indAG12= (dfBetaAG.ens_gene.isin(dfDaf12.gene)) & (dfBetaAG.qval < .1)
indAG16= (dfBetaAG.ens_gene.isin(dfDaf16.gene)) & (dfBetaAG.qval < .1)
ndaf12AG= dfBetaAG.b[indAG12].shape[0]
ndaf16AG= dfBetaAG.b[indAG16].shape[0]
#==============================================================================
# 
#==============================================================================
#color vector:
colors= ['#696969','#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00']
tissues= ['labial', 'extracellular', 'tail', 'dopaminergic' ]
#tissues= ['mu_int', 'sex organ', 'excretory', 'gonadal' ]
#tissues= ['gonad', 'sensillum', 'intestine', 'sex organ' ]
#tissues= ['head', 'tail', 'embryonic' ]
#tissues= ['muscle', 'coel', 'hyp' ]
#tissues= ['sperm', 'extracellular', 'tail', 'dopaminergic' ]
df_exp= organize(tissues, tissue_df)

explode(qval, dfBetaA, df_exp, colors, title= 'Aging', \
        savename= '../output/Graphs/aging_tissue_specific_volcplot.png',\
        xlab= r'$\beta_{\mathrm{Aging}}$')
explode(qval, dfBetaG, df_exp, colors, title= 'Genotype', \
        savename= '../output/Graphs/genotype_tissue_specific_volcplot.png', \
        xlab= r'$\beta_{\mathrm{Genotype}}$')
explode(qval, dfBetaAG, df_exp, colors, title= 'Aging::Genotype', \
        savename= '../output/Graphs/agingxgenotype_tissue_specific_volcplot.png',\
        xlab= r'$\beta_{\mathrm{Aging::Genotype}}$')

colors2= ['#ffff33','#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00']
kegg(qval, dfBetaA, df_exp, colors2, title= 'Aging', \
        savename= '../output/Graphs/aging_tissue_specific_kde.png',\
        xlab= r'$\beta_{\mathrm{Aging}}$')
kegg(qval, dfBetaG, df_exp, colors2, title= 'Genotype', \
        savename= '../output/Graphs/genotype_tissue_specific_kde.png',\
        xlab= r'$\beta_{\mathrm{Genotype}}$')
kegg(qval, dfBetaAG, df_exp, colors2, title= 'Aging::Genotype', \
        savename= '../output/Graphs/agingxgenotype_tissue_specific_kde.png',\
        xlab= r'$\beta_{\mathrm{Aging::Genotype}}$')
#==============================================================================
#Where does gold standard data fall?   
#==============================================================================
#colors
colors= ['#a65628','#ffff33','#ff7f00','#984ea3','#4daf4a','#377eb8', '#e41a1c', '#f781bf']

explode_goldstandards(qval, dfBetaA, dfGoldStandard, colors= colors)
explode_goldstandards(qval, dfBetaG, dfGoldStandard, colors= colors)
explode_goldstandards(qval, dfBetaAG, dfGoldStandard, colors= colors)



#figure out how many genes in dfLIfespan show up in this analysis
f= lambda x: (dfBetaA.ens_gene.isin(x)) & (dfBetaA.qval < .1)    
ind= f(dfLifespanGenes.gene.values)
m= dfBetaA[ind].ens_gene.values

f= lambda x: (dfBetaG.ens_gene.isin(x)) & (dfBetaG.qval < .1)    
ind= f(dfLifespanGenes.gene.values)

m= np.append(m, dfBetaG[ind].ens_gene.values)

f= lambda x: (dfBetaAG.ens_gene.isin(x)) & (dfBetaAG.qval < .1)    
ind= f(dfLifespanGenes.gene.values)

m= np.append(m, dfBetaAG[ind].ens_gene.values)

m= list(set(m))    
with open('lifespan_genes_that_show_up.csv', 'w') as f:
    f.write('WBID\n')
    for gene in m:
        f.write(gene)
        f.write('\n')
    f.close()


Ldf= [dfBetaA, dfBetaG, dfBetaAG] #list of dataframes
dfnames= ['Age', 'Genotype', 'Age::Genotype']
titles= ['Genes associated positively with lifespan', 
         'Genes associated negatively with lifespan',
         'Genes associated variably with lifespan',
         'Genes associated with lifespan but unannotated',
         ]
colors= ['#377eb8','#e41a1c','#4daf4a']
fnames= ['../output/Graphs/positive_aging.png','../output/Graphs/negative_aging.png',
         '../output/Graphs/variable_aging.png','../output/Graphs/unannotated_aging.png']

kegg_compareall_byval(qval, Ldf, dfLifespanGenes, colors, savenames= fnames,
                   dfnames= dfnames, xscale= 'symlog', titles=  titles)



