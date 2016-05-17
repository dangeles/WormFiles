"""
A script to generate primers for gibson assembly of cas9 probes.

author: David Angeles, dangeles@caltech.edu
"""
# -*- coding: utf-8 -*-

from __future__ import division, print_function, absolute_import

# location of the file
path = "/Users/davidangeles/Documents/"
# file name
probe_file = "cas9_probes_unc54_pred.csv"
# the name of the file that will contain the primers
out_file = 'cas9_primers_unc54.txt'

# choose whether to add a 'g' to the front of the sequence if it doesn't begin
# with a 'g', or whether to replace the first letter of the sequence with a 'g'
# notice: if you desire to append, your preference will only be applied if the
# oligo is less than 60bp long before appending
# none means no g's will be appended or replaced.
preference = 'append'  # choose between 'append', 'replace' and 'none

# the gibson sequences to be appended to the oligo
gibson5prime = 'GATCCCCCGGGCTGCAGGAATTCATTTAGGTGACACTATA'
gibson3prime = 'GACTAGCCTTATTTTAACTTGCTATTTCTAGCTCTAAAAC'

basepair = {'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}

# ==============================================================================
#
# ==============================================================================
data = []
with open(path+probe_file, 'r') as f:  # read the file in
    for line in f:
        name = line.split(',')[0].strip()
        seq = line.split(',')[1].strip().lower()
        data.append([name, seq])

# ==============================================================================
#
# ==============================================================================
# reverse complement the sequence
forw = []  # this will be the complete forward sequence
rev = []  # this will be the complete reverse sequence
f = 0
k = 0
for line in data:
    name = line[0]
    seq = line[1]

    # ==============================================================================
    # Append a 'g' for sp6 RNA pol if it is necessary
    # ==============================================================================
    # for each sequence, check that it starts with a 'GG' or 'GA' and
    # that the user actually wants to have a 'g' substituted in or out
    # Else, change whatever
    # you need to change to get it to start that way.
    if seq != 'sequence' and seq[0] != 'g' and preference != 'none':
        # if the oligo is too long, always replace the first letter
        # also replace the letter if ordered by the user
        if preference == 'replace' or len(seq+gibson5prime) >= 60:
            # first letter is not g, so change it:
            seq = seq.replace(seq[0], 'g', 1)
        # if the primer is not too long, and the user wants to append, append
        elif preference == 'append' and len(seq+gibson5prime) < 60:
            seq = 'g'+seq

# ==============================================================================
# reverse complement, add gibson sequences and add to the final arrays
# ==============================================================================
    # reverse complement the sequence
    if f != 0:
        # reverse complement:
        rev_seq = ''.join(basepair[s] for s in seq[::-1]
                          if s in basepair.keys())
        # add gibson sequence
        rev_seq = rev_seq + gibson3prime
        # append to the reverse array
        rev.append([name+'_R', rev_seq])
        # append the fwd sequence with the gibson sequence to the final array:
        seq = seq+gibson5prime
        forw.append([name+'_F', seq])
    else:
        rev.append([name + '_R', seq])
        forw.append([name + '_F', seq])
        f = 1

# ==============================================================================
# if any primers are too long, warn the user
# ==============================================================================
# if any of the primers are too long, print a warning with the name:
for index, value in enumerate(forw):
    name, seq = value[0], value[1]
    if len(seq) > 60:
        print("Primer {0} is too long, check your input sequence".format(name))
for index, value in enumerate(rev):
    name, seq = value[0], value[1]
    if len(seq) > 60:
        print("Primer {0} is too long, check your input sequence".format(name))

# ==============================================================================
# output a primer file for use by idt
# ==============================================================================
with open(path + out_file, 'w') as f:
    for index, value in enumerate(forw):
        if index != 0:
            line = ",".join(value) + ',25nm,STD\n'
            f.write(line)
    for index, value in enumerate(rev):
        if index != 0:
            line = ",".join(value) + ',25nm,STD\n'
            f.write(line)
