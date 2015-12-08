# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 09:00:12 2015

A script to automatically select cas9 guideRNAs in a given exon,
as well as make the gfp homology arms (30bp) and the sgRNA primers,
using the same material as in Schwartz and Sternberg, 2014

Note: all coordinates are given in terms of the original 5' strand provided

@author: davidangeles
"""
import re as re
import pandas as pd
import math
import numpy as np

#given an exon, and a PAM motif, make sgRNAs, their primers and primers for 
#30bp homology to gfp
exon_name= 'sqt-1'
pam= 'ggngg' #desired motif to look for
guide_begins= 2 #position in pam where the sgRNA should start. 0 if pam= 'ngg', 2 if pam='ggngg'
gfpF= 'GGATCAGGTGGAGGAGGTGGAG'
gfpR= 'CTATTTGTATAGTTCATCCATGCCATGTG'

preference= 'append' #choose between 'append', 'replace' and 'none
#the gibson sequences to be appended to the oligo
gibson5prime= 'GATCCCCCGGGCTGCAGGAATTCATTTAGGTGACACTATA'
gibson3prime= 'GACTAGCCTTATTTTAACTTGCTATTTCTAGCTCTAAAAC'


#the name of the file that will contain the primers
out_file= 'cas9_primers_unc54_exon1.txt'
#choose whether to add a 'g' to the front of the sequence if it doesn't begin
#with a 'g', or whether to replace the first letter of the sequence with a 'g'
#notice: if you desire to append, your preference will only be applied if the
#oligo is less than 60bp long before appending
#none means no g's will be appended or replaced.
preference= 'append' #choose between 'append', 'replace' and 'none
#the gibson sequences to be appended to the oligo
gibson5prime= 'GATCCCCCGGGCTGCAGGAATTCATTTAGGTGACACTATA'
gibson3prime= 'GACTAGCCTTATTTTAACTTGCTATTTCTAGCTCTAAAAC'
exon= 'ATGTCGGTGAACAAGGCGACGATCGGAGCGATTAGCTTCTCGGGCTTGACCCTGCTAGCCTGCTTGTGCGCTATTAGCGTCATCTACTCGGATGTGCAGGCGATCTGGAACGAGCTTGACGCTGAGATGGACGGCTTCAAGTTGTTGGCTGATGACCTTTGGCGCGATATGGTCTCGATGGGTGCTGGTAGCCCGGCTAACCGTGCTCGTCGTCAGGCGTATGGTGGATATGGCGCTTCCGGCTCGAACACTGGCGGTTACGGCGGCTACGGCAGCGGTGGCAACCCGACTGGTCCTGGCGGTCCTTCGGACTCCCCCAACTTCCCTGGCGGCAGCTCTCCCAACGTTCCCGGCGGCAGCAGCCCCATCTTCCCCGGAGCTGGCGGTCCCAACGGCGGCAGCAACAACGGTGGTCCCAACGGTCAGCAGTGCCAGTGCGAAAGCAAGAACACCTGCCCCGAAGCCCCCGCTGGTCCCCCTGGCGAGGCTGGCCCCTCCGGTTTGGACGGTCTCCCCGGTCTGCCTGGTGTCGATGGCAAGGACGCCGAAGACATCCACAACGAGCCTCCTCAGGGCTGCTTCACCTGCCCTCAGGGTCCCGCTGGTCCCCAAGGTCCCTCTGGCAAGCCTGGTATCCGCGGTATGCGCGGCCCCAAGGGCTCCCCTGGCTTCCCCGGTAACGATGGTTTCCCCGGCACCCCCGGAGAACAAGGACCGCCCGGTGCTCCCGGCACTGACGGCAAGGCCGGCACTCGCGGTGACAAGGGTAACGATGCTGAGAAGATCGTCGGTCGCAAGGGCAACCGTGGCCTCCCCGGTCCCCAGGGTGATGAAGGCGAAGAAGGTGATGCTGGCAAGGCTGGTGCCCCTGGTGCTCCCGGTGCCGAAGGTCCCGCTGGCTCTCCCGGCTTCCAAGGCCCCACTGGTATCGATGGTGAAGAAGGTCCCGAAGGTCCCGCTGGCGCTCCCGGCAAGGATGCCGAGgtaagtcatacagttcgctttgatgatcatagctaaacctcttcttttcagTACTGCAAGTGCCCTGCCCGCAACGGCGCTGGTGGAGCCGGTGGCCGCAACGGTGGTGGAGCCGGCGGCAACGGCGGTGCTGGCCGTCACGGTGGTCGCGGCGGAAACGGTGGCGGAGCTGGCGGCGCTGGTGGCAACGATGGCTACAACGGTGGTGCCGGTGATGCTGGCGCTTCGCCCTACCGTCGTCTCCGCATCTAA' #in frame!

exon= exon.lower()
size_pam= len(pam)
#==============================================================================
# 
#==============================================================================
def pam_to_regex(pam):
    """
    Given a string, pam, re-write it as a very simple regex.
    """    
    genetic_code= {'a': 'a', 't': 't', 'c':'c', 'g':'g', 'n': '[atgc]'}
    regex= ''
    for letter in list(pam):
        regex+= genetic_code[letter]
    return regex
#==============================================================================
# 
#==============================================================================
#==============================================================================
# 
#==============================================================================
def conjugate_strand(strand):
    """
    Given a dna strand, re-write it as it's reverse complement
    """
    basepair= {'a':'t', 't':'a', 'g':'c', 'c':'g'}
    rev_strand= ''.join(basepair[s] for s in strand[::-1] if s in basepair.keys())
    return rev_strand
#==============================================================================
# 
#==============================================================================
#==============================================================================
# the following algorithm is from John G Doench et al, 
# 'Rational design of highly active sgRNAs for CRISPR-Cas9-mediated gene 
# inactivation
#==============================================================================
def calc_score(s):
    """
    Given a 20-mer, figure out the score for a given sgRNA
    
    Input:
    s, a sgRNA 20bp long
    
    Author: John G Doench et al
    """
    s= s.upper()
    if len(s) > 20:
        print('Warning, sgRNA is more than 20bp long. Exiting')
        return -np.inf()
#    nuc_hash = {'A':0, 'T':1, 'C':2, 'G':3}    
    s_list= list(s)
    #base score given by the broad    
    score = 0.597636154
    gc = s.count('G')+s.count('C')
    
    #gc energies?
    gc_low = -0.202625894
    gc_high = -0.166587752
    
    #modify score based on gc content
    if gc < 10:
        gc_val = abs(gc-10)
        score = score+(gc_val*gc_low)
    elif gc > 10:
        gc_val = gc-10
        score = score+(gc_val*gc_high)
#    print (score)
    #single nucleotide energy matrix
    #rows[1-30]cols['ATCG']
    sing_nuc_hash = {'G2':-0.275377128,'A3':-0.323887456,'C3':0.172128871,\
                    'C4':-0.100666209,'C5':-0.20180294, \
                    'G5':0.245956633,'A6':0.036440041,'C6':0.098376835,'C7':-0.741181291,\
                    'G7':-0.393264397,'A12':-0.466099015,'A15':0.085376945,'C15':-0.013813972,\
                    'A16':0.272620512,'C16':-0.119022648,'T16':-0.285944222,'A17':0.097454592,\
                    'G17':-0.17554617,'C18':-0.345795451,'G18':-0.678096426,'A19':0.22508903,\
                    'C19':-0.507794051,'G20':-0.417373597,'T20':-0.054306959,'G21':0.379899366,\
                    'T21':-0.090712644,'C22':0.057823319,'T22':-0.530567296,'T23':-0.877007428,\
                    'C24':-0.876235846,'G24':0.278916259,'T24':-0.403102218,'A25':-0.077300704,\
                    'C25':0.287935617,'T25':-0.221637217,'G28':-0.689016682,'T28':0.117877577,\
                    'C29':-0.160445304,'G30':0.386342585}
#    score_mat = np.matrix('0 0 0 0;0 0 0 -0.275377128;-0.323887456 0 0.172128871 0;0 0 -0.100666209 0;0 0 -0.20180294 0.245956633;0.036440041 0 0.098376835 0;0 0 -0.741181291 -0.393264397;0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0;-0.466099015 0 0 0;0 0 0 0;0 0 0 0;0.085376945 0 -0.013813972 0;0.272620512 -0.285944222 -0.119022648 0;0.097454592 0 0 -0.17554617;0 0 -0.345795451 -0.678096426;0.22508903 0 -0.507794051 0;0 -0.054306959 0 -0.417373597;0 -0.090712644 0 0.379899366;0 -0.530567296 0.057823319 0;0 -0.877007428 0 0;0 -0.403102218 -0.876235846 0.278916259;-0.077300704 -0.221637217 0.287935617 0;0 0 0 0;0 0 0 0;0 0.117877577 0 -0.689016682;0 0 -0.160445304 0;0 0 0 0.386342585')               
    #double nucleotide energy matrix
    dinuc_hash = {'GT2':-0.625778696,'GC5':0.300043317,'AA6':-0.834836245,'TA6':0.760627772,\
    'GG7':-0.490816749,'GG12':-1.516907439,'TA12':0.7092612,'TC12':0.496298609,\
    'TT12':-0.586873894,'GG13':-0.334563735,'GA14':0.76384993,'GC14':-0.53702517,\
    'TG17':-0.798146133,'GG19':-0.66680873,'TC19':0.353183252,'CC20':0.748072092, \
    'TG20':-0.367266772,'AC21':0.568209132,'CG21':0.329072074,'GA21':-0.836456755,\
    'GG21':-0.782207584,'TC22':-1.029692957,'CG23':0.856197823,'CT23':-0.463207679,\
    'AA24':-0.579492389,'AG24':0.649075537,'AG25':-0.077300704,'CG25':0.287935617,\
    'TG25':-0.221637217,'GT27':0.117877577,'GG29':-0.697740024}
    
    #go through each position 
    for i,nuc in enumerate(s_list):
        key = nuc+str(i+1)
        #assign an energy per a nucleotide according to the energy matrix
        if sing_nuc_hash.has_key(key):
            nuc_score = sing_nuc_hash[key] 
        else:
            nuc_score = 0
#        nuc_score = score_mat[i,nuc_hash[nuc]]
        #add to the score
        score+= nuc_score
        #i <30 in the published algorithm
        if i<19:
            #assign the dinucleotide energies
            dinuc = nuc+s[i+1]+str(i+1)
            if dinuc in dinuc_hash.keys():
                score+= dinuc_hash[dinuc]
    
    partial_score = math.e**-score
    final_score = 1/(1+partial_score)
    return final_score
#==============================================================================
#==============================================================================
# #     
#==============================================================================
#==============================================================================
    
    
def PAM_finder(pam, exon):
    """
    """
    pam= pam_to_regex(pam)
    #in the exon, find all ocurrences of PAM motif
    pam_5prime= [motif.start() for motif in re.finditer(pam, exon)]
    #now search 3prime
    rev_exon= conjugate_strand(exon)
    #return ocurrences of pam motif in coordinates for the template strand provided
    pam_3prime= [len(exon)-size_pam-motif.start() \
                for motif in re.finditer(pam, rev_exon)]
    
    return pam_5prime, pam_3prime

#==============================================================================
# 
#==============================================================================
#once all sgRNAs have been found, select the guide, annotating their position

def make_sgRNA_list(pam, exon):
    """
    """
    
    pam_5prime, pam_3prime= PAM_finder(pam, exon)
    
    #array contains:
    #strand direction, motif coordinate, sequence of sgRNA,
    #a binary variable indicating whether the guide begins with a g or not and
    #a score for the sgRNA sequence 
    guide_list= [[]*5]*(len(pam_5prime)+len(pam_3prime))
    counter= 0
    
    #add all the 5prime guides
    for motif in pam_5prime:
        #add strand sgRNA belongs to, motif coordinate and sequence to array in that order
        #extract the guideRNA from the sequence
        sgRNA= exon[motif+guide_begins-20: motif+guide_begins]
        
        #if the guide isn't 20bp, it's because you've reached the end of the string,
        #so don't make those guides. 
        
        if len(sgRNA) == 20:
            score= calc_score(sgRNA)
            begins_with_g= [1 if sgRNA[0] == 'g' else 0]
            guide_list[counter]= [1, motif, sgRNA, begins_with_g[0], score]
        counter+=1

    for motif in pam_3prime:
        #extract the guideRNA from the sequence
        rev_sgRNA= exon[motif+size_pam-guide_begins: motif+size_pam-guide_begins+20]
        
        #reverse complement it
        sgRNA= conjugate_strand(rev_sgRNA)
        #if the guide isn't 20bp, it's because you've reached the end of the string,
        #so don't make those guides. 
        if len(sgRNA) == 20:
            score= calc_score(sgRNA)
            #check if it ends with a c, since this guide is against the opposite strand
            begins_with_g= [1 if sgRNA[0] == 'g' else 0]
            #add strand sgRNA belongs to, motif coordinate and sequence to array in that order
            guide_list[counter]= [-1, motif, sgRNA, begins_with_g[0], score]
        counter+=1
    #place in a dataframe for ease of use
    df= pd.DataFrame(guide_list,
                     columns= ['strand', 'coordinate', 'sequence',
                     'begins_with_g', 'score']) 
    df.dropna(inplace= True)
        
    return df
#==============================================================================
# 
#==============================================================================
def trim_and_score(df, n= 0, q= .75):
    """
    """
    
    if n== 0 and q > 0:
        #if there are enough guide rnas, keep only the ones that start with a g. 
        if len(df[df.begins_with_g == 1]) > 1:
            
            df= df[df.begins_with_g == 1]
            return df[df.score > df.score.quantile(q)] 
            
        else:
            
            return df[df.score > df.score.quantile(q)] 
            
            
    elif n > 0:
        
        if n >= len(df):
            return df
        
        if len(df[df.begins_with_g == 1]) > 1:
            df= df[df.begins_with_g == 1]
            l= float(len(df))
            nq= 1- n/l            
#            print(nq, 'first if')
            if n>= len(df):
                return df
            else:
                attempt= df[df.score > df.score.quantile(nq)]
                
                while len(attempt) != n:
                    nq= nq-1/l
                    attempt= df[df.score > df.score.quantile(nq)] 
                return df[df.score > df.score.quantile(nq)]
            
        else:
            l= float(len(df))
            nq= 1-n/l
#            print(nq, 'seond')
            attempt= df[df.score > df.score.quantile(nq)]
            
            while len(attempt) != n:
                nq= nq-1/l
                attempt= df[df.score > df.score.quantile(nq)] 
            return df[df.score > df.score.quantile(nq)]

#==============================================================================
#     
#==============================================================================
def make_homology_arms(max_oligo_len, primerF, primerR, df, size_pam, exon):
    """
    """
    #design the gfp homology arms: 
    arm_lengthF= max_oligo_len-len(primerF)
    arm_lengthR= max_oligo_len-len(primerR)
        
    microhomology= [[]*3]*len(df[['strand', 'coordinate']].values)
    counter= 0
    for strand, coordinate in df[['strand', 'coordinate']].values:
        
        coordinate= int(coordinate)
        strand= int(strand)
        if strand == 1:
            #where cas9 cuts
            cut_site= coordinate-3-coordinate%3
            #select the dna beyond the pam site -- we want it removed so the guide
            #doesn't cut!        
            end_pam= coordinate+size_pam + 3-(coordinate+size_pam)%3 
            
            #make sure the arm is in frame the whole time. 
            in_frame= 3-(cut_site-arm_lengthF)%3
            
            beginF= cut_site-arm_lengthF+in_frame
            endR= end_pam + arm_lengthR - (end_pam + arm_lengthR)%3
            
            armF= exon[beginF: cut_site]
            armR= exon[end_pam: endR]
        else:
            cut_site= coordinate+size_pam+3-(coordinate+size_pam)%3
            end_pam= coordinate-size_pam +3 - (coordinate-size_pam)%3
            in_frame= (end_pam-arm_lengthF)%3
            
            beginF= end_pam - arm_lengthF + in_frame
            endR= cut_site+arm_lengthR-arm_lengthR%3
            
            armF= exon[beginF: end_pam]
            armR= exon[cut_site: endR]
            
        armF= armF+primerF.upper()
        armR= conjugate_strand(armR)+primerR.upper()
        
        if len(armF) < 60 and len(armR) < 60:
            microhomology[counter]= [coordinate, armF, armR]
            counter+=1
                
    df_primer= pd.DataFrame(microhomology, 
                            columns= ['coordinate', 'armF', 'armR'])
    df_primer.dropna(inplace= True)
        
    return df_primer
   
#==============================================================================
#Append a 'g' for sp6 RNA pol if it is necessary
#==============================================================================
    #for each sequence, check that it starts with a 'GG' or 'GA' and
    #that the user actually wants to have a 'g' substituted in or out. 
    #. Else, change whatever
    #you need to change to get it to start that way. 
def sp6(df, preference):
#==============================================================================
#     
#==============================================================================
    def sp6_seq(seq, preference):
        """
        Given an sgRNA sequence, make sure the sequence starts with
        a 'g' either by replacing or appending a nucleotide
        seq= a string
        preference= a string, either 'none', 'replace' or 'append'
        """
        if seq != 'sequence' and seq[0] != 'g' and preference != 'none':
            #if the oligo is too long, always replace the first letter
            #also replace the letter if ordered by the user
            if preference == 'replace' or len(seq+gibson5prime) >= 60:
            #first letter is not g, so change it:
                seq= seq.replace(seq[0], 'g', 1)
                #if the primer is not too long, and the user wants to append, append
            elif preference == 'append' and len(seq+gibson5prime) < 60:
                seq= 'g'+seq
        
        return seq
#==============================================================================
#                 
#==============================================================================

    for seq in df.sequence:
        seq_g_corrected= sp6_seq(seq, preference)
        df.sequence.replace(seq, seq_g_corrected, inplace= True)
    
    return df

#==============================================================================
# reverse complement, add gibson sequences and add to the final arrays
#==============================================================================
def make_Gibson(df, gibson5prime, gibson3prime):
    """
    seq= sgRNA sequence
    name= name of the sgRNA, a string. 
    rev= reverse array
    forw= forward oligo array
    """
    #reverse complement the sequence
    #reverse complement:
    df['sgRNA_primer_R']= gibson3prime+df['sequence'].apply(conjugate_strand)
    df['len_R']= df['sgRNA_primer_R'].apply(len)
    #append the forward sequence with the gibson sequence to the final array:
    df['sgRNA_primer_F']= gibson5prime+df['sequence']
    df['len_F']= df['sgRNA_primer_F'].apply(len)

    return df

#==============================================================================
# if any primers are too long, warn the user
#==============================================================================
def oligo_check(df, max_oligo_len= 60):
    """
    """
    #if any of the primers are too long, print a warning with the name:
    df= df[df['len_F'] <= max_oligo_len]
    df= df[df['len_R'] <= max_oligo_len]
    return df

#==============================================================================
# output a primer file for use by idt
#==============================================================================
def write_files(df, exon_name, path= '../output/'):
    """
    """
    name_sgrna= path+exon_name+'_sgrna_and_flanking_homology_arms.csv'
    
    columns= ['sequence', 'armF', 'armR', 'strand', 'score']
    header= ['sgRNA', 'homology_arm_F', 'homology_arm_R', 'strand', 'score'] 
    df.to_csv(name_sgrna, columns= columns, header= header)
        
    name_sgrna= path+exon_name+'primer_IDT_order_sheet.csv'
    
    primers= ['sgRNA_primer_F', 'sgRNA_primer_R', 'armF', 'armR']
    with open(path + exon_name+'_IDT_order_sheet.csv', 'w') as f:
        for coor in df.index.values:
            for primer in primers:
                line= exon_name +'_co_'+ str(int(coor)) + primer + ','
                line+= df.ix[coor, primer]+ ',25nm,STD\n'
                f.write(line)
    f.close()
#==============================================================================
# 
#==============================================================================
def master_call(pam, guide_begins, exon_name, exon, primerF, 
                primerR, preference, gibs5, gibs3, q=0, 
                n= 0, max_oligo_len= 60):
    """
    """
    
    
    if q == 0 and n == 0:
        print('Please specify at least one of q or n, thank you!\nTerminating.')
        return
    
    df_sgrna= make_sgRNA_list(pam, exon)
    df_sgrna= trim_and_score(df_sgrna, n, q)
    
    if len(df_sgrna) == 0:
        print('no guideRNAs found, sorry!')
        return 0, 0
        
    df_gfp= make_homology_arms(max_oligo_len, primerF, primerR, 
                               df_sgrna, len(pam), exon)
    df_sgrna= sp6(df_sgrna, preference)
    df_sgrna= make_Gibson(df_sgrna, gibs5, gibs3)
    df_sgrna= oligo_check(df_sgrna, max_oligo_len)
    
    df_gfp['len_F']= df_gfp['armF'].apply(len)
    df_gfp['len_R']= df_gfp['armR'].apply(len)
    
    df_gfp= oligo_check(df_gfp, max_oligo_len)    
    df_sgrna= df_sgrna[df_sgrna.coordinate.isin(df_gfp.coordinate.values)]
    df_sgrna.set_index('coordinate', inplace= True)
    df_gfp.set_index('coordinate', inplace= True)
    df_sgrna= pd.concat([df_sgrna, df_gfp], axis= 1)
    write_files(df_sgrna, exon_name)
    
    return df_sgrna, df_gfp

    
x, y= master_call(pam, guide_begins, exon_name, exon, gfpF, gfpR, 'append',
                  gibson5prime, gibson3prime, n= 2)
                  

