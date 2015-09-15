# -*- coding: utf-8 -*-
"""
Created on Fri Nov  7 17:43:10 2014

@author: davidaalbores
"""

from __future__ import division, print_function, absolute_import
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.special
from scipy import stats
from scipy import optimize
from scipy import misc
from scipy.integrate import simps
import seaborn as sns

pp= 5*10**-3


#==============================================================================
# 
#==============================================================================
#==============================================================================
# Calculate the information of each term.
#==============================================================================
def information(matrix, totalLabels):
    """
    """
    #given a matrix, sum over each of it's rows and calculate the info.
    I= 0; N= matrix.shape[0]*matrix.shape[1]
    #sum up the info
    for i in range(totalLabels):
        n= matrix[i].sum()
        if n != 0:
            I+= -n/N*np.log(n/N)/np.log(2) 
    
    return I
#==============================================================================
# 
#==============================================================================
#==============================================================================
# 
#==============================================================================
def entropy(matrix, totalLabels, separatorPos, vector= False):
    """
    """
    #Matrix is an m x n matrix where the first dimension given by shape is the number of terms in use currently and the second dimension is the number of genes that was measured. 
    if vector == False:
        A= matrix.shape[0]*separatorPos;
        B= matrix.shape[0]*(matrix.shape[1]-separatorPos)
        H= np.zeros(totalLabels)
    
        for i in range(totalLabels):
            a= matrix[i][:separatorPos].sum()
            b= matrix[i][separatorPos:].sum()
            if a != 0:
                H[i]+= -a*np.log(a/A)*np.log(2)
                
            if b != 0:
                H[i]+= - b*np.log(b/B)*np.log(2)

    else: #you were given a vector, treat it as such
        A= separatorPos;
        B= len(matrix) - separatorPos
        H= 0
        
        a= matrix[:separatorPos].sum()
        b= matrix[separatorPos:].sum()
        print(a, " a")
        print(b, " b" )
#        print(A, "sepA")
#        print(a, "a")
#        print(b, "b")
        if a != 0:
            H+= -a*np.log(a/A)*np.log(2)
            
        if b != 0:
            H+= - b*np.log(b/B)*np.log(2)
        
    return H
#==============================================================================
# 
#==============================================================================
#==============================================================================
# 
#==============================================================================
def slider(matrix, totalLabels, maxSep= 0, vector= False):
    """
    """
    if vector== False:
        maxSep= int(np.floor((matrix.shape[1])-2))
    
    H= np.zeros(shape= (maxSep-1, totalLabels))
    for separator in range(1, maxSep):
        H[separator-1]= entropy(matrix, totalLabels, separator, vector);
        
    X= np.array(range(0, len(H)))/maxSep;

    return X, H
#==============================================================================
#     
#==============================================================================
#==============================================================================
#     
#==============================================================================
def deltaWindow(matrix, totalLabels, sepSize):
    """
    """
    D= np.zeros(totalLabels)
    for i in range(totalLabels):
        x= np.floor(pd.rolling_sum(matrix[i, :], sepSize)[sepSize-1:]/sepSize)
        x.astype(int)
        D[i]= x.sum()
    return D
#==============================================================================
#     
#==============================================================================
def deltaSlider(matrix, totalLabels):
    """
    """
    maxSlider= 3;
    maxSlider= int(maxSlider)
    D= np.zeros(shape=(maxSlider, totalLabels))    
    for separator in range(1, maxSlider+1):
        D[separator-1]= deltaWindow(matrix, totalLabels, separator)
    
    X= np.array(range(1, maxSlider+1))    
    return X, D
#==============================================================================
#     
#==============================================================================
#==============================================================================
# 
#==============================================================================
def threshold(H1, H2, X1, X2, terms, delta, N1, info1, N2, info2, counts1, counts2, fold= True, plotter = False, interval= 'all'):    
    """
    """

    percentageCounts= pp
    
    found= []
    diff= []
    foldDiff= []
    
    a1= 0; b1= len(H1[:, 1])-1
    a2= 0; b2= len(H2[:, 1])-1

    if interval== '1':
        #only itegrate over first quartile:
        a1= 0;
        b1= len(H1[:, 1])/4
        a2= 0;
        b2= len(H2[:, 1])/4
    elif interval == '2':
        a1= len(H1[:, 1])/4
        b1= len(H1[:, 1])/2
        a2= len(H2[:, 1])/4
        b2= len(H2[:, 1])/2
    elif interval == '3':
        a1= len(H1[:, 1])/2
        b1= len(H1[:, 1])/4*3
        a2= len(H2[:, 1])/2
        b2= len(H2[:, 1])/4*3
    elif interval == '4':
        a1= len(H1[:, 1])/4*3
        b1= len(H1[:, 1])-1
        a2= len(H2[:, 1])/4*3
        b2= len(H2[:, 1])-1
    
    a1= int(np.floor(a1)); b1= int(np.floor(b1))
    a2= int(np.floor(a2)); b2= int(np.floor(b2))


    dx1= np.diff(X1)[0]
    dx2= np.diff(X2)[0]
    
    for i in range(len(terms)):
        integral1= scipy.integrate.simps(H1[a1:b1, i], dx= dx1)/N1/info1
        integral2= scipy.integrate.simps(H2[a2:b2, i], dx= dx2)/N2/info2
        if integral1 != 0 and integral2 != 0: #only consider stuff that is present in both 
            if fold == False:
                difference= (integral2 - integral1)
            else:
                difference= (integral2 - integral1) / integral1  
                        
            if np.absolute(difference) > delta:
                if counts1[i] > percentageCounts*counts1.sum() and counts2[i]>percentageCounts*counts2.sum():
                    found.append(terms[i])
                    diff.append(integral2-integral1)
                    foldDiff.append((integral2 - integral1)/integral1)
#                    print(found)

#                    if terms[i] == 'head':
#                        print(a1, b1,  " interval")
#                        print(len(H1[:, i]))
#                        print(integral1)
#                        print(a2, b2,  " interval2")
#                        print(len(H2[:, i]))
#                        print(integral1)
                    
                    if plotter:                
                        plt.figure()
                        plt.plot(X1, H1[:, i]/N1/info1, 'r', label= r'Control')
                        plt.plot(X2, H2[:, i]/N2/info2, 'b', label= r'Experimental')
                        plt.title(terms[i])
                        plt.xlabel(r'Separator Location')
                        plt.ylabel(r'H/I')
                        plt.legend(numpoints=1,loc=5)
    
                    
    return found, diff, foldDiff
#==============================================================================
# 
#==============================================================================
def reporter(found, diff, fDiff, counts1, counts2, rank1, rank2, terms, totalRanks1, totalRanks2, sort= True):
    
    
    if sort == True:
        sortedIndices= [i[0] for i in sorted(enumerate(fDiff), key=lambda x:x[1])]
    else:
        sortedIndices= [i[0] for i in sorted(enumerate(diff), key=lambda x:x[1])]
        
    indicator= -1
    
    for elements in sortedIndices:
        
        if diff[elements] > 0 and indicator != 1:
            indicator= 0;
        
        
        if indicator == 0:
            print("\n\nUPREGULATED:\n\n")
            print("Tissue\tDifference\tFoldChange\tControlCounts\tExperimentalCounts\tMin. Rank. Control\tMedian Rank Control\t Max. Rank Control\Min. Rank. Exp.\tMedian Rank Control\tMax. Rank Control\tMann Whitney P-Value\tSignificance")
            indicator= 1;
        
        print(indicator)
            
        if indicator == -1:
            print("\n\nDOWNREGULATED\n\n")
            print("Tissue\tDifference\tFoldChange\tControlCounts\tExperimentalCounts\tMin. Rank. Control\tMedian Rank Control\t Max. Rank Control\Min. Rank. Exp.\tMedian Rank Control\tMax. Rank Control\tMann Whitney P-Value\tSignificance")
            indicator= 0
        
        #find the index of the found term:
        index= terms.index(found[elements])
        #get the count
        c1= counts1[index]; c2= counts2[index]
        
        if found[elements] in rank1.keys() and found[elements] in rank2.keys():
            medianRank1= np.median(rank1[found[elements]])/totalRanks1;
            medianRank2= np.median(rank2[found[elements]])/totalRanks2; 
            meanRank1= np.mean(rank1[found[elements]])/totalRanks1;
            meanRank2= np.mean(rank2[found[elements]])/totalRanks2
            minRank1= np.min(rank1[found[elements]])/totalRanks1
            minRank2= np.min(rank2[found[elements]])/totalRanks2
            maxRank1= np.max(rank1[found[elements]])/totalRanks1
            maxRank2= np.max(rank2[found[elements]])/totalRanks2
        
            length1= 0; length2= 0;
            q= True
            
            while q:
                try:
                    length1= len(rank1[found[elements]])
                    q= False;
                except TypeError:
                    length1= 1;
                    q= False
                    
                        
            q= True;
            while q:
                try:
                    length2= len(rank2[found[elements]])
                    q= False;
                except TypeError:
                    length2= 1;
                    q= False
                    
            if length1 > 20 and length2> 20:
                mannWhitney= scipy.stats.mannwhitneyu(rank1[found[elements]], rank2[found[elements]])

                #use bonferroni correct:
                if mannWhitney[1] < 0.05/len(found):
                    significance= 'YES'
                else:
                    significance= 'NO'
        
                print("{0}\t{1:.3f}\t{2:.1f}\t{3:.4f}\t{4:.4f}\t{5:.2}\t{6:.2}\t{7:.2}\t{8:.2}\t{9:.2}\t{10:.2}\t{11:.2}\t{12:.2}\t{13:.2}\t{14}".format(found[elements], diff[elements], fDiff[elements], c1, c2, minRank1, medianRank1, meanRank1, maxRank1, minRank2, medianRank2, meanRank2, maxRank2, mannWhitney[1], significance))
            else:
                print("{0}\t{1:.3f}\t{2:.1f}\t{3:.4f}\t{4:.4f}\t{5:.2}\t{6:.2}\t{7:.2}\t{8:.2}\t{9:.2}\t{10:.2}\t{11:.2}\t{12:.2}\t{13:.2}\t{14}".format(found[elements], diff[elements], fDiff[elements], c1, c2, minRank1, medianRank1, meanRank1, maxRank1, minRank2, medianRank2, meanRank2, maxRank2, 'NA', 'NA'))

        else:
            if found[elements] in rank1.keys():
                print("{0}\t{1:.3f}\t{2:.1f}\t{3:.4f}\t{4:.4f}\t{5:.2}\t{6:.2}\t{7:.2}\t{8:.2}\t{9:.2}\t{10:.2}\t{11:.2}\t{12:.2}\t{13:.2}\t{14}".format(found[elements], diff[elements], fDiff[elements], c1, c2, minRank1, medianRank1, meanRank1, maxRank1, 'NA', 'NA', 'NA', 'NA', 'NA', 'NA'))
            else: 
                print("{0}\t{1:.3f}\t{2:.1f}\t{3:.4f}\t{4:.4f}\t{5:.2}\t{6:.2}\t{7:.2}\t{8:.2}\t{9:.2}\t{10:.2}\t{11:.2}\t{12:.2}\t{13:.2}\t{14}".format(found[elements], diff[elements], fDiff[elements], c1, c2, 'NA', 'NA', 'NA', 'NA', minRank2, medianRank2, meanRank2, maxRank2, 'NA', 'NA'))


#==============================================================================
# 
#==============================================================================
#==============================================================================
# 
#==============================================================================
def countLabels(tissues, threshold= 0, junk= False):
    """
    A function that takes as input an array of tissues and a scalar threshold
    value in order to generate a number of outputs.
    It generates an array of the labels; an array with the number of times any
    one label appears; it makes a total count calculcation and it figures out
    what the most common tissue is and how many times it's called. 
    """
    countsArray= np.array([]); countsLabels=[];
    totalCounts= 0;
    maxCount=0; maxTissue= 0;
    total= 0;

    for element in set(tissues):
        count= tissues.count(element)
        if count > maxCount:
            maxCount= count;
            maxCount= element
        if count >= threshold:
            countsLabels.append(element);
            countsArray= np.append(countsArray, count)
            totalCounts+=count;
    if junk:
        return countsArray, countsLabels, maxCount, maxTissue, total
    else:
        return countsArray, countsLabels
#==============================================================================
# 
#==============================================================================
#==============================================================================
# 
#==============================================================================
def completeAllTissues(complete, labels, countsIncomplete):
    """
    A function that given a complete list of terms, a list of terms present in an  experiment (labels), the counts in the experiment for each label and the fraction of the total that label represents, returns the counts and fractional amount associated with a label in an array with length equal to the length of the complete list of terms. 
    """
    counts= np.array([])
    for element in complete:
        #get the index of the element in the control or experimental set
        index= 0
        if element in labels:
            index= labels.index(element)
        else:
                #element is only in the exp set.
                #Hence, the expected number of counts for that tissue type is zero.
                #Since this can't happen though, let it be some somewhat small number
                #but not too small or otherwise it's going to bias the whole thing.
                #this is arguably the most arbitrary part of it all. 
            index= complete.index(element)
            counts= np.append(counts, 0)
            continue;
        counts= np.append(counts, countsIncomplete[index])
    
    return counts
#==============================================================================
# 
#==============================================================================
#==============================================================================
# 
#==============================================================================
def thresholdFinder(H1, H2, X1, X2, terms, N1, info1, N2, info2, counts1, counts2, keyword, fold= True, interval2= 'all'):    
    """
    """    
    delta = 0
    
#    
#    
#    
    def differencer(H1, H2, X1, X2, terms, keyword, delta, N1, info1, N2, info2, counts1, counts2, fold1= True,  interval= 'all'):
        
        percent= pp
        a1= 0; b1= len(H1[:, 1])-1
        a2= 0; b2= len(H2[:, 1])-1
        
        if interval== '1':
            #only itegrate over first quartile:
            a1= 0;
            b1= len(H1[1,:])/4
            a2= 0;
            b2= len(H2[:,1])/4
        elif interval == '2':
            a1= len(H1[:,1])/4
            b1= len(H1[:,1])/2
            a2= len(H2[:,1])/4
            b2= len(H2[:,1])/2
        elif interval == '3':
            a1= len(H1[:,1])/2
            b1= len(H1[:,1])/4*3
            a2= len(H2[:,1])/2
            b2= len(H2[:,1])/4*3
        elif interval == '4':
            a1= len(H1[:,1])/4*3
            b1= len(H1[:,1])-1
            a2= len(H2[:,1])/4*3
            b2= len(H2[:,1])-1
        else:
            a1= a1; b1= b1;
    
        a1= int(np.floor(a1)); b1= int(np.floor(b1))
        a2= int(np.floor(a2)); b2= int(np.floor(b2))
            
            
        dx1= np.diff(X1)[0]
        dx2= np.diff(X2)[0]

    
        index= terms.index(keyword)
        integral1= scipy.integrate.simps(H1[a1:b1, index], dx= dx1)/N1/info1
        integral2= scipy.integrate.simps(H2[a2:b2, index], dx= dx2)/N2/info2

        if integral1 != 0 and integral2 != 0: #only consider stuff that is present in biasoth
            if fold == False:
                difference= (integral2 - integral1)
            else:
                difference= (integral2 - integral1) / integral1                
            if np.absolute(difference) > delta :#and H1[:,i].max>0.1*N1*info1 and   H2[:, i].max >0.1*N2*info2
                if counts1[index] > percent*counts1.sum() and counts2[index]>percent*counts2.sum():
                    if keyword == terms[index]:
                        return 1;
        return 0;
#        
#        
#        
        
    inspector = 1; 
    iterations= 0
    maxIterations= 1000;
    error= 0;
    while inspector != 0:
        x= differencer(H1, H2, X1, X2, terms, keyword, delta, N1, info1, N2, info2, counts1, counts2, fold1= fold, interval= interval2)
        if iterations > maxIterations:
            if x != 1:
                delta-= 10**-3;
                if iterations > 10*3*maxIterations:
                    print("Your tissue cannot be found; please pick a different tissue and try again. The error may be due to low copy reads")
                    error= 1
                    delta= 0;
                    break;
            else:
                inspector= 0;
                break;
        else:
            if x == 1: 
                delta+= 10**-2 #your delta is too low
            else: 
                delta-= 10**-3
        iterations+=1
    
    
        
    return delta, error
    
#==============================================================================
#     
#==============================================================================
def summary(H_Control, H_Experimental, XCont, XExp, terms, NCont, infoControl, NExp, infoExp, controlTissueCounts, expTissueCounts, controlRankTissues, expRankTissues, cLength, eLength, fold1= True, search= 'none', plotter1= False, interval3= 'all'):
    """    
    """
    if search.lower() == 'none':
        if fold1 == True:
            delta= 0.1;
            error= 0;
        else:
            delta= 0.0001; 
            error= 0;
    else:
        delta, error= thresholdFinder(H_Control, H_Experimental, XCont, XExp, terms, NCont, infoControl, NExp, infoExp, controlTissueCounts, expTissueCounts, search, fold= fold1, interval2= interval3)

    if error == 1:
        print('Could not find delta. Terminating program')
        return 0, 0, 0, delta, 1

    found, diff, fdiff= threshold(H_Control, H_Experimental, XCont, XExp, terms, delta, NCont, infoControl, NExp, infoExp, controlTissueCounts, expTissueCounts, fold= fold1, plotter= plotter1, interval= interval3)
    print (delta) 
    
    if not found:
        print("Didn't find anything that seemed to change. Terminating program")
        return 0, 0, 0, delta, 1;

    reporter(found, diff, fdiff, controlTissueCounts, expTissueCounts, controlRankTissues, expRankTissues, terms, cLength, eLength)
    
    return found, diff, fdiff, delta, 0
#==============================================================================
# 
#==============================================================================
#==============================================================================
# 
#==============================================================================
def plotter(found, rank1, rank2, c, e, fileName= 'None'):
    """
    found -- an array of names
    rank1, rank2 -- hashes with keys names (in found) and values that are arrays of ranks
    c, e -- normalization factors.
    fileName -- a filename if you want to save this. 
    """
    #Close everything
    plt.clf()
    plt.close('all')
    #Declare some variables
    tickNames= []
    ytot= np.array([])
    NormControl= c
    NormExp= e
    
    #Make an array to parse into violinplots, want y1 next to y2 for every y1andy2
    for name in found:
        #get the value of y
        y = rank1[name]/NormControl
        y2= rank2[name]/NormExp
        #append to a numpy array
        ytot= np.append(ytot, [y, y2])
        #make the labels
        tickNames.append(name)
        tickNames.append(name+' Exp')
    
    #plot
    #Make the axis and figure
    fig, ax1 = plt.subplots(figsize=(10,6))
    sns.set_context(rc={'lines.markeredgewidth': 0.1})
    sns.set(style="whitegrid")
    
    #Plot, label and make pretty.
    plot= sns.violinplot(ytot, inner= 'points', color= 'Paired')
    sns.despine(offset= 10, trim= True)
    plot.set_title('Violinplot of tissue rank data')
    plot.set_ylabel('Normalized Rank')
    plot.set_xticklabels(tickNames, rotation=45, fontsize= 12)
    
    if fileName != 'None':
        fig.savefig('../output/' +fileName, bbox_inches='tight')
#==============================================================================
# 
#==============================================================================
#==============================================================================
# 
#==============================================================================
def shannonplotter(H1, H2, X1, X2, N1, info1, N2, info2, terms, found, filename= 'none'):
    
    for element in found:
        index= terms.index(element)
        fig, ax= plt.subplots(figsize= (10, 6))
        plt.title('Shannon Entropy Curves for ' + element)
        plt.plot(X1, H1[:, index]/info1/N1, 'r', label= r"Control")
        plt.plot(X2, H2[:, index]/info2/N2, 'b', label= r"Experimental")
        plt.xlabel(r'Separator Location')
        plt.ylabel(r'H/I')
        plt.legend(numpoints=1,loc=5)

        if filename != "none":
            fig.savefig('../output/'+filename)
    
    
    
    
