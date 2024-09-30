#!/usr/bin/env python
"""
Module of MeStanG v0.1
@author: Daniel Ramos Lopez
Metagenomic Standards Generator (MeStanG) for HTS Nanopore datasets 
"""

import sys
import random
import numpy as np

def ratio_h(input_r, input_pr):
    host_reads = int(input_r)
    prob = ['p', 'h']
    for i in range(len(input_pr)):
        if input_pr[i] >= 0.5:
            sys.stderr.write("\nPlease provide valid pathogen ratios\n")
            sys.exit(1)
        input_pr[i] = np.random.choice(prob, input_r, p=[input_pr[i], 1 - input_pr[i]]).tolist().count('p')
        host_reads -= input_pr[i]
        if host_reads <= 0:
            sys.stderr.write("\nError: There are more pathogens than host itself!\n"
                             "Ratios might be too high for the population\n")
            sys.exit(1)
    input_pr.append(host_reads)
    return input_pr
    

def ratio_e(input_r, input_pr):
    prob = ['r', 'd']
    while True:
        reads = int(input_r)
        member_pr = input_pr[:]
        member_pr[random.randint(0, len(member_pr) - 1)] = 0
        for i in range(len(member_pr)):
            if member_pr[i] > 1.0 or (member_pr[i] >= 1.0 and i > 0):
                sys.stderr.write("\nPlease provide valid ratios\n")
                sys.exit(1)
            if member_pr[i] != 0:
                member_pr[i] = np.random.choice(prob, input_r, p=[member_pr[i], 1 - member_pr[i]]).tolist().count('r')
            reads -= member_pr[i]
        if reads > 0:
            break
    member_pr[member_pr.index(0)] = reads       
    return member_pr


def dist_ratio(number, subtaxa, equally):
    nr = []
    if equally:
        number_for_subtaxon = int(number / subtaxa)
        for i in range(subtaxa):
            nr.append(number_for_subtaxon)
            if i == subtaxa - 1:
                for j in range(number % subtaxa):
                    nr[j] += 1
    else:
        number_subtaxon = int(number)
        for i in range(subtaxa):
            if number_subtaxon != 0:
                number_for_subtaxon = np.random.randint(1, number_subtaxon + 1)
                number_subtaxon -= number_for_subtaxon
                nr.append(number_for_subtaxon)
                if i == subtaxa - 1 and number_subtaxon != 0:
                    add = int(number_subtaxon / subtaxa)
                    for j in range(subtaxa):
                        nr[j] += add
                        if j == subtaxa - 1:
                            for k in range(number_subtaxon % subtaxa):
                                nr[k] += 1
        if len(nr) != subtaxa:
            for j in range(subtaxa - len(nr)):
                nr.append(0)
    random.shuffle(nr)
    return nr
    
    
def dist_weight(n_reads, file_l):
    seq_len = 0
    weighted_reads = []
    for i in range(int(len(file_l) / 2)):
        weighted_reads.append(len(file_l[2 * i + 1]))
        seq_len += len(file_l[2 * i + 1])
    for i in range(len(weighted_reads)):
        weighted_reads[i] = int((float(weighted_reads[i]) / float(seq_len)) * n_reads)
        if weighted_reads[i] == 0:
            weighted_reads[i] = 1
    if sum(weighted_reads) < n_reads:
        residual = int(n_reads - sum(weighted_reads))
        for j in range(residual):
            random_entry = random.randint(0, len(weighted_reads) - 1)
            weighted_reads[random_entry] += 1
    return weighted_reads


def dist_ratio_communities(big_ratio, subtaxa):
    nr = []
    number_subtaxon = float(big_ratio)
    for i in range(subtaxa):
        if number_subtaxon != 0:
            number_for_subtaxon = random.uniform(0, number_subtaxon)
            number_subtaxon -= number_for_subtaxon
            nr.append(number_for_subtaxon)
            if i == subtaxa - 1 and number_subtaxon != 0:
                add = float(number_subtaxon / subtaxa)
                for j in range(subtaxa):
                    nr[j] += add
    if len(nr) != subtaxa:
        for j in range(subtaxa - len(nr)):
            nr.append(0)
    random.shuffle(nr)
    return nr
    
