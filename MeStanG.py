#!/usr/bin/env python
"""
MeStanG v1.0
@author: Daniel Ramos Lopez
Metagenomic Standards Generator (MeStanG) 
"""

from __future__ import print_function
from __future__ import with_statement

import multiprocessing as mp
from textwrap import dedent
from time import strftime
import sys
import os
import argparse
import numpy
import random
import string
from Bio import SeqIO
from Bio.Seq import Seq

VERSION = "1.0"
PROGRAM = "MeStanG"
AUTHOR = "Daniel Ramos Lopez"
CONTACT = "aramosl@okstate.edu, adramoslp@outlook.com"

### Distribution and ratio

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
    

### Files to inputs

def files_to_inputs_se(master, support, number, equally):
    files_list = []
    with open(master) as r:
        header_f = r.readline().strip('\n')
        header_f = header_f.split('\t')
        while True:
            s = r.readline().strip('\n')
            if s == '':
                break
            files_list.append(s.split('\t'))
    
    m_ratios = []
    
    if support:
        supp = []
        with open(support) as r:
            header_s = r.readline().strip('\n')
            header_s = header_s.split('\t')
            while True:
                s = r.readline().strip('\n')
                if s == '':
                    break
                supp.append(s.split('\t'))

        if 'ratio' in header_s and 'reads' in header_s:
            sys.stderr.write('\nInsert reads or ratio, not both\n')
            sys.exit(1)
        if 'taxon' not in header_s or 'file' not in header_f or 'taxon' not in header_f:
            sys.stderr.write('\nInsert valid input files, files and taxon must be consistent\n')
            sys.exit(1)
        if 'ratio' in header_s:
            if 'ratio' in header_f:
                sys.stderr.write('\nInsert valid input files, insert ratios only in one file\n')
                sys.exit(1)
        elif 'reads' in header_s:
            if 'reads' in header_f:
                sys.stderr.write('\nInsert valid input files, insert reads only in one file\n')
                sys.exit(1)
        elif 'ratio' not in header_s and 'ratio' not in header_s and 'ratio' not in header_f and 'ratio' not in header_f:
            sys.stderr.write('\nInsert valid input files, insert ratios or reads\n')
            sys.exit(1)

        if 'ratio' in header_s:
            for i in supp:
                m_ratios.append(float(i[header_s.index('ratio')]))
            big_ratios = ratio_e(number, m_ratios)
        else:
            for i in supp:
                m_ratios.append(int(i[header_s.index('reads')]))
            big_ratios = m_ratios[:]

        sys.stdout.write("\n" + strftime("%Y-%m-%d %H:%M:%S") + ": Calculating number of reads of each member of "
                                                                "the community\n")

        ratios = []
        files = []
        header = header_s[:]
        
        for i in range(len(supp)):
            subtaxa_for_taxon = []
            for j in files_list:
                if supp[i][0] in j:
                    subtaxa_for_taxon.append(j)
            sub_ratios = dist_ratio(big_ratios[i], len(subtaxa_for_taxon), equally)
            for j in range(len(subtaxa_for_taxon)):
                if sub_ratios[j] != 0:
                    entry = [subtaxa_for_taxon[j][0]]
                    for k in range(len(supp[i]) - 2):
                        entry.append(supp[i][k + 2])
                    entry.append(supp[i][0])
                    files.append(entry)
                    ratios.append(sub_ratios[j])
        
        header[0] = 'file'
        if 'ratio' in header:
            header.remove('ratio')
        if 'reads' in header:
            header.remove('reads')
        header.append('taxon')
        
    else:
        if 'file' not in header_f:
            sys.stderr.write('\nFile column must be present\n')
            sys.exit(1)
        if 'ratio' in header_f and 'reads' in header_f:
            sys.stderr.write('\nInsert reads or ratio, not both\n')
            sys.exit(1)

        sys.stdout.write("\n" + strftime("%Y-%m-%d %H:%M:%S") + ": Calculating number of reads of each member of "
                                                                "the community\n")
        if 'ratio' in header_f:
            for i in files_list:
                m_ratios.append(float(i[header_f.index('ratio')]))
            ratios = ratio_e(number, m_ratios)
        else:
            for i in files_list:
                m_ratios.append(int(i[header_f.index('reads')]))
            ratios = m_ratios[:]
            
        header = header_f[:]
        
        for i in files_list:
            if 'ratio' in header_f:
                i.pop(header.index('ratio'))
            else:
                i.pop(header.index('reads'))
        files = files_list[:]

        if 'ratio' in header:
            header.remove('ratio')
        if 'reads' in header:
            header.remove('reads')

    return files, header, ratios


def files_to_inputs_sh(master, support, number, equally_p, equally_h):
    files_list = []
    hosts_list = []
    with open(master) as r:
        header_f = r.readline().strip('\n')
        header_f = header_f.split('\t')
        while True:
            s = r.readline().strip('\n')
            if s == '':
                break
            if 'host' in s.split('\t'):
                hosts_list.append(s.split('\t'))
            else:
                files_list.append(s.split('\t'))

    pathogen_ratios = []
    
    if support:
        supp = []
        with open(support) as r:
            header_s = r.readline().strip('\n')
            header_s = header_s.split('\t')
            while True:
                s = r.readline().strip('\n')
                if s == '':
                    break
                supp.append(s.split('\t'))

        if 'ratio' in header_s and 'reads' in header_s:
            sys.stderr.write('\nInsert reads or ratio, not both\n')
            sys.exit(1)
        if 'taxon' not in header_s or 'file' not in header_f or 'taxon' not in header_f:
            sys.stderr.write('\nInsert valid input files, files and taxon must be consistent\n')
            sys.exit(1)
        if 'ratio' in header_s:
            if 'ratio' in header_s and 'ratio' in header_f:
                sys.stderr.write('\nInsert valid input files, insert ratios only in one file\n')
                sys.exit(1)
        elif 'reads' in header_s:
            if 'reads' in header_s and 'reads' in header_f:
                sys.stderr.write('\nInsert valid input files, insert reads only in one file\n')
                sys.exit(1)
        else:
            sys.stderr.write('\nInsert valid input files, insert ratios or reads\n')
            sys.exit(1)

        if 'ratio' in header_s:
            for i in supp:
                pathogen_ratios.append(float(i[header_s.index('ratio')]))
            big_ratios = ratio_h(number, pathogen_ratios)
        else:
            for i in supp:
                pathogen_ratios.append(int(i[header_s.index('reads')]))
            big_ratios = pathogen_ratios[:]
            pathogen_reads = 0
            for i in big_ratios:
                pathogen_reads += i
            big_ratios.append(number - pathogen_reads)

        sys.stdout.write("\n" + strftime("%Y-%m-%d %H:%M:%S") + ": Calculating number of reads of each member of "
                                                                "the community\n")

        ratios = []
        files = []
        header = header_s[:]
        for i in range(len(supp)):
            subtaxa_for_taxon = []
            for j in files_list:
                if supp[i][0] in j:
                    subtaxa_for_taxon.append(j)
            sub_ratios = dist_ratio(big_ratios[i], len(subtaxa_for_taxon), equally_p)
            for j in range(len(subtaxa_for_taxon)):
                if sub_ratios[j] != 0:
                    entry = [subtaxa_for_taxon[j][0]]
                    for k in range(len(supp[i]) - 2):
                        entry.append(supp[i][k + 2])
                    entry.append(supp[i][0])
                    files.append(entry)
                    ratios.append(sub_ratios[j])

        if files == []:
            files = [[0, 0]]
        
        header[0] = 'file'
        if 'ratio' in header:
            header.remove('ratio')
        if 'reads' in header:
            header.remove('reads')
        header.append('taxon')
        hosts_ratios = dist_ratio(big_ratios[len(big_ratios) - 1], len(hosts_list), equally_h)

        for i in hosts_list:
            for j in range(len(files[0]) - 2):
                i.append('-')
            i.append('host')
            i.pop(1)

        files += hosts_list
        ratios += hosts_ratios

    else:
        if 'file' not in header_f:
            sys.stderr.write('\nFile column must be present\n')
            sys.exit(1)
        if 'ratio' in header_f and 'reads' in header_f:
            sys.stderr.write('\nInsert reads or ratio, not both\n')
            sys.exit(1)

        sys.stdout.write("\n" + strftime("%Y-%m-%d %H:%M:%S") + ": Calculating number of reads of each member of "
                                                                "the community\n")
        if 'ratio' in header_f:
            for i in files_list:
                pathogen_ratios.append(float(i[header_f.index('ratio')]))
            ratios = ratio_h(number, pathogen_ratios)
        else:
            for i in files_list:
                pathogen_ratios.append(int(i[header_f.index('reads')]))
            ratios = pathogen_ratios[:]
            pathogen_reads = 0
            for i in ratios:
                pathogen_reads += i
            ratios.append(int(number) - int(pathogen_reads))

        header = header_f[:]
        for i in files_list:
            if 'ratio' in header_f:
                i.pop(header.index('ratio'))
            else:
                i.pop(header.index('reads'))
        files = files_list[:]

        hosts_ratios = dist_ratio(ratios[len(ratios) - 1], len(hosts_list), equally_h)
        ratios.pop()

        for i in hosts_list:
            if 'ratio' in header_f:
                i.pop(header.index('ratio'))
            else:
                i.pop(header.index('reads'))

        files += hosts_list
        ratios += hosts_ratios

        if 'ratio' in header:
            header.remove('ratio')
        if 'reads' in header:
            header.remove('reads')

    return files, header, ratios


def files_to_inputs_re(master, number, equally, master_out):
    taxa_list = []
    with open(master) as r:
        header_t = r.readline().strip('\n')
        header_t = header_t.split('\t')
        while True:
            s = r.readline().strip('\n')
            if s == '':
                break
            taxa_list.append(s.split('\t')[header_t.index('taxon')])

    taxa_list = list(set(taxa_list))

    if 'file' not in header_t or 'taxon' not in header_t:
        sys.stderr.write('\nInsert valid input file\n')
        sys.exit(1)

    r_ratios = dist_ratio_communities(1, len(taxa_list))
    reads_per_member = ratio_e(number, r_ratios)

    with open('temp_' + master_out + '.txt', 'w') as temp:
        temp.write('taxon\treads\n')
        for i in range(len(taxa_list)):
            temp.write(taxa_list[i] + '\t' + str(reads_per_member[i]) + '\n')

    files, header, ratios = files_to_inputs_se(master, 'temp_' + master_out + '.txt', number, equally)
    os.remove('temp_' + master_out + '.txt')
    
    return files, header, ratios


def community_generator(community_list, header_s, community, master_out, number, min_ratio, max_ratio, targets, decoys):
    if community == 't':
        max_reads = int(number / 2)
        while max_reads >= int(number / 2):
            big_ratio = 0
            while big_ratio == 0:
                big_ratio = random.uniform(min_ratio, max_ratio)
            ratios = dist_ratio_communities(big_ratio, len(targets))
            zero_ratios = 0
            max_reads = 0
            for i in range(len(ratios)):
                current_reads = ratio_h(number, [ratios[i]])[0]
                if current_reads == 0:
                    zero_ratios += 1
                max_reads += current_reads
                ratios[i] = current_reads
        if zero_ratios == len(ratios):
            ratios[random.randint(0, len(ratios) - 1)] = 1
        if len(header_s) > 2:
            header_temp = '\t' + '\t'.join(header_s[2:]) + '\n'
        else:
            header_temp = '\n'
        with open('temp_' + master_out + '.txt', 'w') as temp:
            temp.write('taxon\treads' + header_temp)
            for i in range(len(targets)):
                if len(header_s) > 2:
                    params = '\t' + '\t'.join(targets[i][2:]) + '\n'
                else:
                    params = '\n'
                temp.write(targets[i][0] + '\t' + str(ratios[i]) + params)

    if community == 'd':
        max_reads = int(number / 2)
        while max_reads >= int(number / 2):
            big_ratio = 0
            while big_ratio == 0:
                big_ratio = random.uniform(min_ratio, max_ratio)
            ratios = dist_ratio_communities(big_ratio, len(decoys))
            zero_ratios = 0
            max_reads = 0
            for i in range(len(ratios)):
                current_reads = ratio_h(number, [ratios[i]])[0]
                if current_reads == 0:
                    zero_ratios += 1
                max_reads += current_reads
                ratios[i] = current_reads
        if zero_ratios == len(ratios):
            ratios[random.randint(0, len(ratios) - 1)] = 1
        if len(header_s) > 2:
            header_temp = '\t' + '\t'.join(header_s[2:]) + '\n'
        else:
            header_temp = '\n'
        with open('temp_' + master_out + '.txt', 'w') as temp:
            temp.write('taxon\treads' + header_temp)
            for i in range(len(decoys)):
                if len(header_s) > 2:
                    params = '\t' + '\t'.join(decoys[i][2:]) + '\n'
                else:
                    params = '\n'
                temp.write(decoys[i][0] + '\t' + str(ratios[i]) + params)

    if community == 't+d':
        max_reads = int(number / 2)
        while max_reads >= int(number / 2):
            big_ratio = 0
            while big_ratio == 0:
                big_ratio = random.uniform(min_ratio, max_ratio)
            ratios = dist_ratio_communities(big_ratio, len(community_list))
            max_reads = 0
            for i in range(len(ratios)):
                current_reads = ratio_h(number, [ratios[i]])[0]
                max_reads += current_reads
                ratios[i] = current_reads
        zt_count = 0
        for i in range(len(targets)):
            if ratios[community_list.index(targets[i])] == 0:
                zt_count += 1
        if zt_count == len(targets):
            ratios[community_list.index(targets[random.randint(0, len(targets) - 1)])] = 1
        
        zd_count = 0                    
        for i in range(len(decoys)):
            if ratios[community_list.index(decoys[i])] == 0:
                zd_count += 1
        if zd_count == len(decoys):
            ratios[community_list.index(decoys[random.randint(0, len(decoys) - 1)])] = 1
        
        if len(header_s) > 2:
            header_temp = '\t' + '\t'.join(header_s[2:]) + '\n'
        else:
            header_temp = '\n'
        with open('temp_' + master_out + '.txt', 'w') as temp:
            temp.write('taxon\treads' + header_temp)
            for i in range(len(community_list)):
                if len(header_s) > 2:
                    params = '\t' + '\t'.join(community_list[i][2:]) + '\n'
                else:
                    params = '\n'
                temp.write(community_list[i][0] + '\t' + str(ratios[i]) + params)


def files_to_inputs_rh(master, number, equally_p, equally_h, master_out, sp_file, community, min_ratio, max_ratio):
    community_list = []
    targets = []
    decoys = []
    with open(sp_file) as r:
        header_s = r.readline().strip('\n')
        header_s = header_s.split('\t')
        while True:
            s = r.readline().strip('\n')
            if s == '':
                break
            community_list.append(s.split('\t'))
            if community_list[len(community_list) - 1][header_s.index('type')] == 'decoy':
                decoys.append(community_list[len(community_list) - 1])
            if community_list[len(community_list) - 1][header_s.index('type')] == 'target':
                targets.append(community_list[len(community_list) - 1])
    
    if 'taxon' not in header_s or 'type' not in header_s or targets == []:
        sys.stderr.write('\nInsert valid input file\n')
        sys.exit(1)
    
    community_generator(community_list, header_s, community, master_out, number, min_ratio, max_ratio, targets, decoys)
    files, header, ratios = files_to_inputs_sh(master, 'temp_' + master_out + '.txt', number, equally_p, equally_h)
    os.remove('temp_' + master_out + '.txt')

    return files, header, ratios


### File manager

def base_conversion(base):
    base = base.upper()
    if base not in ['A', 'C', 'G', 'T']:
        base_code = {'Y': ['C', 'T'], 'R': ['A', 'G'],
                     'W': ['A', 'T'], 'S': ['G', 'C'],
                     'K': ['T', 'G'], 'M': ['C', 'A'],
                     'D': ['A', 'G', 'T'], 'V': ['A', 'C', 'G'], 'H': ['A', 'C', 'T'], 'B': ['C', 'G', 'T'],
                     'N': ['A', 'T', 'C', 'G']}
        base = random.choice(base_code[base])
    return base


def fasta_to_list(file):
    records = list(SeqIO.parse(file,"fasta"))
    f_list = []
    for i in records:
        f_list.append(">" + str(i.id).replace("_","-"))
        f_list.append(str(i.seq))
    return f_list


def model_loader(file):
    mis = {'A': {'C': 0.25, 'G': 0.5, 'T': 0.25}, 'C': {'A': 0.25, 'G': 0.25, 'T': 0.5}, 'G': {'A': 0.5, 'C': 0.25, 'T': 0.25}, 'T': {'A': 0.25, 'C': 0.5, 'G': 0.25}}
    ins = {'A': {'A': 0.5, 'C': 0.125, 'G': 0.25, 'T': 0.125}, 'C': {'C': 0.5, 'A': 0.125, 'G': 0.125, 'T': 0.25}, 'G': {'G': 0.5, 'A': 0.25, 'C': 0.125, 'T': 0.125}, 'T': {'T': 0.5, 'A': 0.125, 'C': 0.25, 'G': 0.125}}        

    if file:
        model_file = open(file, "r")
        model_file = model_file.readlines()
        for i in range(len(model_file)):
            model_file[i] = model_file[i].strip('\n')
            model_file[i] = model_file[i].split(':')

        mis['A']['C'], mis['C']['A'] = float(model_file[1][1]), float(model_file[1][1])
        mis['A']['G'], mis['G']['A'] = float(model_file[2][1]), float(model_file[2][1])
        mis['A']['T'], mis['T']['A'] = float(model_file[3][1]), float(model_file[3][1])
        mis['C']['G'], mis['G']['C'] = float(model_file[4][1]), float(model_file[4][1])
        mis['C']['T'], mis['T']['C'] = float(model_file[5][1]), float(model_file[5][1])
        mis['G']['T'], mis['T']['G'] = float(model_file[6][1]), float(model_file[6][1])

        ins['A']['A'] = float(model_file[8][1])
        ins['C']['C'] = float(model_file[12][1])
        ins['G']['G'] = float(model_file[15][1])
        ins['T']['T'] = float(model_file[17][1])
        ins['A']['C'], ins['C']['A'] = float(model_file[9][1]), float(model_file[9][1])
        ins['A']['G'], ins['G']['A'] = float(model_file[10][1]), float(model_file[10][1])
        ins['A']['T'], ins['T']['A'] = float(model_file[11][1]), float(model_file[11][1])
        ins['C']['G'], ins['G']['C'] = float(model_file[13][1]), float(model_file[13][1])
        ins['C']['T'], ins['T']['C'] = float(model_file[14][1]), float(model_file[14][1])
        ins['G']['T'], ins['T']['G'] = float(model_file[16][1]), float(model_file[16][1])

    for key in mis:
        if (sum(list(mis[key].values()))) != 1.0:
            print('Mismatch model probabilities do not sum 1, fix the model rates')
            sys.exit(1)
        if (sum(list(ins[key].values()))) != 1.0:
            print('Insertion model probabilities do not sum 1, fix the model rates')
            sys.exit(1)
            
    return mis, ins


def remove_ids(master_out):
    tag = ''.join(random.choices(string.ascii_uppercase + string.ascii_lowercase + string.digits, k=10))
    sys.stdout.write("\n" + strftime("%Y-%m-%d %H:%M:%S") + ": Removing sequence ids\n")
    name_o = "./" + master_out + "/" + master_out
    with open(name_o + "_no_ids.fasta", 'a') as output:
        records = SeqIO.parse(name_o + "_reads.fasta", "fasta")
        seqind = 1
        for record in records:
            record.id = master_out + "_sequence_" + str(seqind) + "_" + tag
            record.description = master_out + "_sequence_" + str(seqind) + "_" + tag
            SeqIO.write(record, output, "fasta")
            seqind += 1
    os.remove(name_o + "_reads.fasta")
    os.rename(name_o + "_no_ids.fasta", name_o + "_reads.fasta")


def dorado_gen(error_prob):
    dist_error = numpy.random.dirichlet(numpy.ones(3),size=1).tolist()[0]
    dist_error = [i * error_prob[1] for i in dist_error]
    error_prob.pop(1)
    return error_prob + dist_error
    

# Based on profiles presented in doi.org/10.3389/fgene.2019.01332, doi.org/10.5281/zenodo.10038673, and doi.org/10.5281/zenodo.10397818
def error_selector(profile, basecaller):
    if profile == 'perfect':
        error_prob = [1, 0, 0, 0]
    else:
        if basecaller == 'guppy':
            if profile == 'virus':
                error_prob = [0.9038, 0.03, 0.0202, 0.046]
            elif profile == 'bacteria_ec':
                error_prob = [0.9095, 0.0302, 0.0197, 0.0406]
            elif profile == 'bacteria_kpn':
                error_prob = [0.8725, 0.0449, 0.0416, 0.0410]
            elif profile == 'bacteria':
                error_prob = [0.891, 0.03755, 0.03065, 0.0408]
            elif profile == 'human':
                error_prob = [0.8990, 0.0286, 0.0246, 0.0478]
            else:
                error_prob = [0.896200, 0.033425, 0.026525, 0.043850]
        if basecaller == 'dorado':
            if profile == 'fast':
                error_prob = [0.91, 0.09]
            elif profile == 'hac':
                error_prob = [0.962, 0.038]
            elif profile == 'sup':
                error_prob = [0.977, 0.023]
            elif profile == 'res':
                error_prob = [0.973, 0.027]
            error_prob = dorado_gen(error_prob)
    return error_prob


def get_custom_profile(custom):
    error_prob=[1, 0, 0, 0]
    error = 0
    with open(custom) as summary:
        while True:
            line = summary.readline().strip('\n')
            if line == '':
                break
            line = line.split('\t')
            if line[0] in ['Mismatches', 'mismatches', 'Mismatch', 'mismatch']:
                error_prob[1] = float(line[1])
            if line[0] in ['Insertions', 'insertions', 'Insertion', 'insertion']:
                error_prob[2] = float(line[1])
            if line[0] in ['Deletions', 'deletions', 'Deletion', 'deletion']:
                error_prob[3] = float(line[1])
            error += float(line[1])
    error_prob[0] = float(error_prob[0] - error)
    if error_prob[0] <= 0:
        print("\nPlease input valid error rates\n")
        sys.exit(1)
    return error_prob
    

def read_summary(file):
    met = []
    with open(file) as summary:
        summary.readline()
        while True:
            line = summary.readline().strip('\n')
            if line == '':
                break
            line = line.split(' ')
            met.append(line[len(line) - 1])
    for i in range(len(met)):
        if i < 4:
            met[i] = int(met[i])
        else:
            met[i] = float(met[i])
    return met


def validate_header(stat, header, files, args, parser, i):
    if 'mean' in header and files[i][header.index('mean')] != '-':
        stat.append(int(files[i][header.index('mean')]))
    else:
        stat.append(args.mean)

    if 'sd_len' in header and files[i][header.index('sd_len')] != '-':
        stat.append(float(files[i][header.index('sd_len')]))
    else:
        stat.append(args.sd_len)

    if 'posrate' in header and files[i][header.index('posrate')] != '-':
        stat.append(float(files[i][header.index('posrate')]))
    else:
        stat.append(args.posrate)

    if 'profile' in header and files[i][header.index('profile')] != '-':
        profile = files[i][header.index('profile')]
    else:
        profile = args.profile

    if 'basecaller' in header and files[i][header.index('basecaller')] != '-':
        basecaller = files[i][header.index('basecaller')]
    else:
        basecaller = args.basecaller

    if 'circular' in header and files[i][header.index('circular')] != '-':
        if files[i][header.index('circular')] == 'true':
            circular = True
        else:
            circular = False
    else:
        circular = args.circular

    if 'custom' in header and files[i][header.index('custom')] != '-':
        custom = files[i][header.index('custom')]
    else:
        custom = args.custom
        
    if 'em_model' in header and files[i][header.index('em_model')] != '-':
        em_model = files[i][header.index('em_model')]
    else:
        em_model = args.em_model

    if 'error_profile' in header and files[i][header.index('error_profile')] != '-':
        if files[i][header.index('error_profile')] == 'true':
            error_profile = True
        else:
            error_profile = False
    else:
        error_profile = args.error_profile

    if 'no_metrics' in header and files[i][header.index('no_metrics')] != '-':
        if files[i][header.index('no_metrics')] == 'true':
            no_metrics = True
        else:
            no_metrics = False
    else:
        no_metrics = args.no_metrics

    if 'unweighted' in header and files[i][header.index('unweighted')] != '-':
        if files[i][header.index('unweighted')] == 'true':
            unweighted = True
        else:
            unweighted = False
    else:
        unweighted = args.unweighted

    if stat[3] and (stat[3] < 0 or stat[3] > 1):
        print("\nPlease input proper posrate value between 0 and 1\n")
        parser.print_help(sys.stderr)
        sys.exit(1)
        
    if basecaller == 'guppy' and profile in ['fast', 'hac', 'sup', 'res']:
        print("\nPlease input proper profile model for Guppy, model for Dorado was given\n")
        parser.print_help(sys.stderr)
        sys.exit(1)
        
    if basecaller == 'dorado' and profile in ['virus', 'bacteria', 'human', 'bacteria_ec', 'bacteria_kpn']:
        print("\nPlease input proper profile model for Dorado, model for Guppy was given\n")
        parser.print_help(sys.stderr)
        sys.exit(1)

    return profile, basecaller, circular, custom, error_profile, no_metrics, unweighted, em_model


def get_abundance_file(name_o, ratios, files, header, total):
    with open(name_o + "_abundance.tsv", 'a') as meta_comp:
        if 'taxon' in header:
            meta_comp.write('Member\tTaxon\tN_Reads\tPercentage\n')
            for j in range(len(files)):
                if ratios[j] != 0:
                    files[j][0] = files[j][0].split("/")[len(files[j][0].split("/")) - 1]
                    member_ratio = (float(ratios[j]) / float(total)) * 100
                    meta_comp.write(
                        '%s\t%s\t%i\t%.3f\n' % (os.path.splitext(files[j][0])[0], files[j][header.index('taxon')],
                                                ratios[j], member_ratio))
        else:
            meta_comp.write('Member\tN_Reads\tPercentage\n')
            for j in range(len(files)):
                if ratios[j] != 0:
                    files[j][0] = files[j][0].split("/")[len(files[j][0].split("/")) - 1]
                    member_ratio = (float(ratios[j]) / float(total)) * 100
                    meta_comp.write(
                        '%s\t%i\t%.3f\n' % (os.path.splitext(files[j][0])[0], ratios[j], member_ratio))


def display_console(file, stat, profile, basecaller, circular, custom, error_profile, no_metrics, unweighted, num_threads, name_o, equally, em_model):
    print("\nRunning error insertion with the following parameters:\n")
    print("input: %s" % file)
    print("number: %i" % stat[0])
    print("mean: %i" % stat[1])
    print("sd_len: %i" % stat[2])
    print("posrate: %.3f" % stat[3])
    print("profile: %s" % profile)
    print("basecaller: %s" % basecaller)
    print("circular: %s" % circular)
    print("custom: %s" % custom)
    print("em_model: %s" % em_model)
    print("error_profile: %s" % error_profile)
    print("no_metrics: %s" % no_metrics)
    print("unweighted: %s" % unweighted)
    print("equally: %s" % equally)
    print("threads: %s" % num_threads)
    with open(name_o + "_parameters.txt", 'a') as parameters:
        parameters.write("input: " + str(file) + "\n" +
                        "number: " + str(stat[0]) + "\n" +
                        "mean: " + str(stat[1]) + "\n" +
                        "sd_len: " + str(stat[2]) + "\n" +
                        "posrate: " + str(stat[3]) + "\n" +
                        "profile: " + str(profile) + "\n" +
                        "basecaller: " + str(basecaller) + "\n" +
                        "circular: " + str(circular) + "\n" +
                        "custom: " + str(custom) + "\n" +
                        "error_profile: " + str(error_profile) + "\n" +
                        "no_metrics: " + str(no_metrics) + "\n" +
                        "unweighted: " + str(unweighted) + "\n" + 
                        "equally: " + str(equally) + "\n\n")


### MeStanG main

def trim_matches(value):
    errors = ['mis', 'ins', 'del']
    if value in errors:
        return True
    else:
        return False


def insert_error(perfect_read, error_prob, out_reads, seq_index, f_name, sense, initial_index,
                   error_profile, out_error, no_metrics, out_metrics, mis, ins):
    numpy.random.seed()
    if not no_metrics:
        error_ls = [0, 0, 0]
    name_error = ['mat', 'mis', 'ins', 'del']
    error = list(filter(trim_matches, numpy.random.choice(name_error, len(perfect_read), p=error_prob).tolist()))
    error_pos = random.sample(range(len(perfect_read)), len(error))
    error_pos.sort()
    tag = ''.join(random.choices(string.ascii_uppercase + string.ascii_lowercase + string.digits, k=5))
    if len(error) > 0 and len(error_pos) > 0:
        profiles = [[error_pos[0]], [error[0]], [1], [], []]
        for i in range(len(error) - 1):
            if profiles[0][len(profiles[0]) - 1] == error_pos[i + 1] - 1 and profiles[1][len(profiles[1]) - 1] == error[i + 1]:
                profiles[2][len(profiles[2]) - 1] += 1
            else:
                profiles[0].append(error_pos[i + 1])
                profiles[1].append(error[i + 1])
                profiles[2].append(1)
        simulated_read = ''
        i_pos = 0
        for i in range(len(profiles[1])):
            ref_seq = ''
            error_seq = ''
            for j in range(profiles[2][i]):
                if profiles[1][i] == 'mis':
                    if perfect_read[profiles[0][i] + j] not in ['A','C','G','T']:
                        base = base_conversion(perfect_read[profiles[0][i] + j])
                    else:
                        base = perfect_read[profiles[0][i] + j]
                    ref_seq += base
                    error_seq += numpy.random.choice(list(mis[base].keys()), p=list(mis[base].values()))
                    if not no_metrics:
                        error_ls[0] = error_ls[0] + 1
                if profiles[1][i] == 'ins':
                    if j == 0:
                        base = perfect_read[profiles[0][i]]
                        ref_seq += base
                        error_seq += base
                    else:
                        base = error_seq[len(error_seq) - 1]
                    if base not in ['A','C','G','T']:
                        base = base_conversion(base)
                    error_seq += numpy.random.choice(list(ins[base].keys()), p=list(ins[base].values()))
                    if not no_metrics:
                        error_ls[1] = error_ls[1] + 1
                if profiles[1][i] == 'del':
                    ref_seq += perfect_read[profiles[0][i]]
                    error_seq += '-'
                    if not no_metrics:
                        error_ls[2] = error_ls[2] + 1
            profiles[3].append(ref_seq)
            profiles[4].append(error_seq)
            simulated_read += perfect_read[i_pos:profiles[0][i]]
            if error_profile:
                out_error.write('%s_%i_%s\t%i\t%i\t%s\t%i\t%s\t%s\n' % (
                        f_name[1:], seq_index, tag, initial_index + profiles[0][i], len(simulated_read), profiles[1][i], profiles[2][i], profiles[3][i], profiles[4][i]))
            if profiles[1][i] != 'del':
                simulated_read += profiles[4][i]
            if profiles[1][i] != 'ins':
                i_pos = profiles[0][i] + profiles[2][i]
            else:
                i_pos = profiles[0][i] + 1
        simulated_read += perfect_read[i_pos:]
    else:
        simulated_read = perfect_read
    
    out_reads.write('%s_%i_%s_%s_%i_%i\n%s\n' % (
                    f_name, seq_index, tag, sense, initial_index, len(simulated_read), simulated_read))
    summary = []

    if not no_metrics and len(error) > 0 and len(error_pos) > 0:
        accuracy = ((float(len(simulated_read)) - error_ls[0]) * 100) / (float(len(simulated_read)) + error_ls[1] + error_ls[2])
        error_rate = (float(error_ls[0] + error_ls[1] + error_ls[2]) * 100) / (float(len(simulated_read)) + error_ls[1] + error_ls[2])
        out_metrics.write('%s_%i_%s\t%s\t%i\t%i\t%i\t%i\t%i\t%.5f\t%.5f\n' %
                          (f_name[1:], seq_index, tag, sense, initial_index, len(simulated_read), error_ls[0], error_ls[1], error_ls[2], accuracy, error_rate))
        summary = [len(simulated_read)]
        summary += error_ls
        summary += [accuracy, error_rate]
    
    if not no_metrics and len(summary) == 0:
        summary = [len(simulated_read)]
        summary += [0, 0, 0]
        summary += [100, 0]
        
    return summary


def get_template_reads(stat, file_l, circular, error_prob, out_reads_name, error_profile, no_metrics, out_error_name, out_metrics_name, num_threads, mis, ins):
    out_reads = open(out_reads_name, "a")
    if error_profile:
        out_error = open(out_error_name, "a")
    else:
        out_error = ''

    if not no_metrics:
        out_metrics = open(out_metrics_name, "a")
        summaries = [0, 0, 0, 0, 0, 0]
    else:
        out_metrics = ''

    for i in range(stat[0]):
        numpy.random.seed()
        if num_threads != 1:
            seq_index = total_reads.value
            total_reads.value += 1
        else:
            seq_index = i
            total_reads.value += 1
        len_read = int(numpy.random.normal(stat[1], float(stat[2])))
        if len_read > len(file_l[1]):
            len_read = len(file_l[1])
        if circular:
            initial_index = numpy.random.randint(0, len(file_l[1]))
            if initial_index > (len(file_l[1]) - len_read):
                perfect_read = (file_l[1] + file_l[1])[initial_index:initial_index + len_read]
            else:
                perfect_read = file_l[1][initial_index:initial_index + len_read]
        else:
            initial_index = numpy.random.randint(0, len(file_l[1]) - len_read + 1)
            perfect_read = file_l[1][initial_index:initial_index + len_read]
        sense = 'F'
        posrate = stat[3]
        if posrate != 1.0:
            if numpy.random.choice(['F', 'R'], p=[posrate, 1 - posrate]) == 'R':
                perfect_read = str(Seq(perfect_read).reverse_complement())
                sense = 'R'
        summary = insert_error(perfect_read, error_prob, out_reads, seq_index + 1, file_l[0], sense, initial_index,
                                 error_profile, out_error, no_metrics, out_metrics, mis, ins)

        if not no_metrics:
            for j in range(len(summaries)):
                summaries[j] = summaries[j] + float(summary[j]) / float(stat[0])

    out_reads.close()
    if error_profile:
        out_error.close()
    if not no_metrics and stat[0] != 0:
        out_summary = open(out_metrics_name + "_summary", "a")
        if stat[0] == 1:
            n_of_reads = " read)"
        else:
            n_of_reads = " reads)"
        out_summary.write("Metrics in terms of averages for " + str(file_l[0][1:]) + " (" + str(stat[0])
                          + n_of_reads + "\n" + "Read length = %i" % int(
                          summaries[0]) + "\n" + "Mismatches = %i" % int(
                          summaries[1]) + "\n" + "Insertions = %i" % int(
                          summaries[2]) + "\n" + "Deletions = %i" % int(
                          summaries[3]) + "\n" + "Accuracy = %.5f" % float(
                          summaries[4]) + "\n" + "Error rate = %.5f" % float(
                          summaries[5]) + "\n\n")
        out_summary.close()
        out_metrics.close()


def denovo(file, stat, profile, basecaller, circular, name_o, error_profile, custom, no_metrics, unweighted, num_threads, em_model):
    global total_reads
    total_reads = mp.Value("i", 0, lock=True)
    out_reads_name = name_o + '_reads.fasta'
    if not custom:
        error_prob = error_selector(profile, basecaller)
    else:
        error_prob = get_custom_profile(custom)
    mis, ins = model_loader(em_model)
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Processing file\n")
    file_l = fasta_to_list(file)
    if len(file_l) == 2:
        unweighted = False
        weighted_reads = [int(stat[0])]
    else:
        if not unweighted:
            weighted_reads = dist_weight(stat[0], file_l)
        else:
            weighted_reads = [int(stat[0])] * int(len(file_l) / 2)
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Error insertion running\n")
    out_error_name = ''
    if error_profile:
        out_error_name = name_o + '_error_profile'
        with open(out_error_name, "a") as out_error:
            if os.stat(out_error_name).st_size == 0:
                out_error.write("Seq_name\tRef_pos\tSeq_pos\terror_type\terror_length\tref_base\tseq_base\n")

    out_metrics_name = ''
    if not no_metrics:
        out_metrics_name = name_o + '_metrics'
        with open(out_metrics_name, "a") as out_metrics:
            if os.stat(out_metrics_name).st_size == 0:
                out_metrics.write("Read\tSense\tS.pos_genome\tbp_ref_genome\tMismatches\tInsertions\tDeletions\tAccuracy\tError_rate\n")

    if num_threads == 1:
        for i in range(len(weighted_reads)):
            stat[0] = weighted_reads[i]
            get_template_reads(stat, file_l[i * 2:i * 2 + 2], circular, error_prob, out_reads_name, error_profile, no_metrics, out_error_name, out_metrics_name, num_threads, mis, ins)
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Finished! reads generated for " + str(file_l[i * 2][1:]) +
                             " >> " + str(int(stat[0])) + "\n")

    else:
        
        if len(weighted_reads) > 1 and not no_metrics:
            m_summaries = []
            for i in range(len(weighted_reads) - 1):
                m_summaries.append(out_metrics_name + "_summary_" + str(i + 1))

        for entry in range(len(weighted_reads)):
            procs = []
            out_reads_subfiles = []
            out_error_subfiles = []
            out_metrics_subfiles = []
            out_summary_subfiles = []
            stat[0] = weighted_reads[entry]
            num_reads = int(stat[0] / num_threads)
            last_reads = int(stat[0] % num_threads)
            stat[0] = int(num_reads)

            for i in range(num_threads):
                out_reads_subfile = out_reads_name + "{}".format(i)
                out_reads_subfiles.append(out_reads_subfile)

                if error_profile:
                    out_error_subfile = out_error_name + "{}".format(i)
                    out_error_subfiles.append(out_error_subfile)
                else:
                    out_error_subfile = ''

                if not no_metrics:
                    out_metrics_subfile = out_metrics_name + "{}".format(i)
                    out_metrics_subfiles.append(out_metrics_subfile)
                    out_summary_subfile = out_metrics_name + "{}".format(i) + "_summary"
                    out_summary_subfiles.append(out_summary_subfile)
                else:
                    out_error_subfile = ''

                if i != num_threads - 1:
                    if no_metrics:
                        out_metrics_subfile = "null"
                    p = mp.Process(target=get_template_reads, args=(stat, file_l[entry * 2:entry * 2 + 2], circular, error_prob, out_reads_subfile, error_profile, no_metrics, out_error_subfile, out_metrics_subfile, num_threads, mis, ins))
                    procs.append(p)
                    p.start()
                else:
                    if no_metrics:
                        out_metrics_subfile = "null"
                    stat[0] = stat[0] + last_reads
                    p = mp.Process(target=get_template_reads, args=(stat, file_l[entry * 2:entry * 2 + 2], circular, error_prob, out_reads_subfile, error_profile, no_metrics, out_error_subfile, out_metrics_subfile, num_threads, mis, ins))
                    procs.append(p)
                    p.start()

            for p in procs:
                p.join()

            with open(out_reads_name, "a") as output_r:
                for fname in out_reads_subfiles:
                    with open(fname) as infile:
                        output_r.write(infile.read())

            for fname in out_reads_subfiles:
                os.remove(fname)

            if error_profile:
                with open(out_error_name, "a") as output_e:
                    for fname in out_error_subfiles:
                        with open(fname) as infile:
                            output_e.write(infile.read())

                for fname in out_error_subfiles:
                    os.remove(fname)

            if not no_metrics:
                with open(out_metrics_name, "a") as output_m:
                    for fname in out_metrics_subfiles:
                        with open(fname) as infile:
                            output_m.write(infile.read())

                for fname in out_metrics_subfiles:
                    os.remove(fname)

                len_s = int(len(out_summary_subfiles))
                not_path = []
                for fname in range(len_s):
                    if not os.path.exists(out_summary_subfiles[fname]):
                        not_path.append(out_summary_subfiles[fname])

                for fname in not_path:
                    out_summary_subfiles.remove(fname)

                if entry == 0:
                    suffix = ''
                else:
                    suffix = "_{}".format(entry)

                with open(out_metrics_name + "_summary" + suffix, "a") as output_s:
                    big_summary = [0, 0, 0, 0, 0, 0]
                    for fname in out_summary_subfiles:
                        met = read_summary(fname)
                        for ent in range(len(met)):
                            big_summary[ent] = big_summary[ent] + met[ent]
                    for ent in range(len(big_summary)):
                        big_summary[ent] = big_summary[ent] / len(out_summary_subfiles)
                    if weighted_reads[entry] == 1:
                        n_of_reads = " read)"
                    else:
                        n_of_reads = " reads)"
                    output_s.write("Metrics in terms of averages for " + file_l[entry * 2][1:] + " (" + str(weighted_reads[entry])
                                   + n_of_reads + "\n" + "Read length = %i" % int(
                                   big_summary[0]) + "\n" + "Mismatches = %i" % int(
                                   big_summary[1]) + "\n" + "Insertions = %i" % int(
                                   big_summary[2]) + "\n" + "Deletions = %i" % int(
                                   big_summary[3]) + "\n" + "Accuracy = %.5f" % float(
                                   big_summary[4]) + "\n" + "Error rate = %.5f" % float(
                                   big_summary[5]) + "\n")

                for fname in out_summary_subfiles:
                    os.remove(fname)

        if len(weighted_reads) > 1 and not no_metrics:
            with open(out_metrics_name + "_summary", "a") as output_s:
                for fname in m_summaries:
                    with open(fname, "r") as infile:
                        output_s.write("\n")
                        output_s.write(infile.read())
                output_s.write("\n\n")

            for fname in m_summaries:
                    os.remove(fname)

        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Finished! reads generated for " + str(file) +
                         " >> " + str(sum(weighted_reads)) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description=dedent('''
        MeStanG
        -----------------------------------------------------------
        Given a set of genomes, assemblies, or contigs
	    generate standard datasets of raw ONT reads
        '''),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-v', '--version', action='version', version='MeStanG ' + VERSION)
    subparsers = parser.add_subparsers(help="You may run on env or host", dest='sample',
                                       description=dedent('''
                For detailed usage of each sample source:
                    MeStanG.py sample -h
                -------------------------------------------------------
                '''))
    
    env = subparsers.add_parser('env', help="Environmental samples")

    subenv = env.add_subparsers(help="You may run env on standard or random mode",
                                dest='mode', description='''
                For detailed usage of each mode:
                    MeStanG.py env mode -h
                ''')

    parser_se = subenv.add_parser('st', help="Generate standard environmental sample")

    parser_se.add_argument('-f', '--file', help='Input file list', required=True)
    parser_se.add_argument('-tx', '--taxon', help='File with organism taxon designation (Optional)'
                          , default=None)
    parser_se.add_argument('-n', '--number', help='Number of reads to be simulated for each sequence (Default = 100)',
                          type=int, default=100)
    parser_se.add_argument('-m', '--mean', help='Average read length (Default = 2000)', type=int, default=2000)
    parser_se.add_argument('-sd', '--sd_len', help='Standard deviation of read length in normal scale (Default = 0)',
                          type=float, default=0)
    parser_se.add_argument('-s', '--posrate',
                          help='Positive sense read rate (Default = 1.0)', type=float,
                          default=1.0)
    parser_se.add_argument('-p', '--profile',
                          help='Basecaller profile error model (Default = fast)',
                          choices=["virus", "bacteria_ec", "bacteria_kpn", "bacteria", "human", "consensus", "perfect", "fast", "hac", "sup", "res"], default='fast')
    parser_se.add_argument('-b', '--basecaller',
                          help='Basecaller error model (Default = dorado)',
                          choices=["dorado", "guppy"], default='dorado')
    parser_se.add_argument('--circular', help='Set genome structure as circular',
                          action='store_true', default=False)
    parser_se.add_argument('-c', '--custom', help='File with custom basecalling error model', default=None)
    parser_se.add_argument('-em', '--em_model', help='File with custom base emission model (Optional)'
                          , default=None)
    parser_se.add_argument('-o', '--output', help='Output directory (Default = MeStanG_output)',
                          default='MeStanG_output')
    parser_se.add_argument('--equally', help='Distribute the ratio of a taxon in its subtaxa equally',
                          action='store_true', default=False)
    parser_se.add_argument('--error_profile', help='Generate error insertion profile output',
                          action='store_true', default=False)
    parser_se.add_argument('--no_metrics', help='Skip generating metrics files',
                          action='store_true', default=False)
    parser_se.add_argument('--no_ids', help='Generate sequences without id',
                          action='store_true', default=False)
    parser_se.add_argument('--unweighted', help='Use the same number of reads for every entry in Multifasta files',
                          action='store_true', default=False)
    parser_se.add_argument('-t', '--num_threads', help='Number of threads'
                            '(Default = 1)', type=int, default=1)

    parser_re = subenv.add_parser('rd', help="Generate randomly distributed environmental sample")

    parser_re.add_argument('-f', '--file', help='Input file list', required=True)
    parser_re.add_argument('-n', '--number', help='Number of reads to be simulated for each sequence (Default = 100)',
                          type=int, default=100)
    parser_re.add_argument('-m', '--mean', help='Average read length (Default = 2000)', type=int, default=2000)
    parser_re.add_argument('-sd', '--sd_len', help='Standard deviation of read length in normal scale (Default = 0)',
                          type=float, default=0)
    parser_re.add_argument('-s', '--posrate',
                          help='Positive sense read rate (Default = 1.0)', type=float,
                          default=1.0)
    parser_re.add_argument('-p', '--profile',
                          help='Basecaller profile error model (Default = fast)',
                          choices=["virus", "bacteria_ec", "bacteria_kpn", "bacteria", "human", "consensus", "perfect", "fast", "hac", "sup", "res"], default='fast')
    parser_re.add_argument('-b', '--basecaller',
                          help='Basecaller error model (Default = dorado)',
                          choices=["dorado", "guppy"], default='dorado')
    parser_re.add_argument('--circular', help='Set genome structure as circular',
                          action='store_true', default=False)
    parser_re.add_argument('-c', '--custom', help='File with custom basecalling error model', default=None)
    parser_re.add_argument('-em', '--em_model', help='File with custom base emission model (Optional)'
                          , default=None)
    parser_re.add_argument('-o', '--output', help='Output directory (Default = MeStanG_output)',
                          default='MeStanG_output')
    parser_re.add_argument('--equally', help='Distribute abundance of a taxon among subtaxa equally',
                          action='store_true', default=False)
    parser_re.add_argument('--error_profile', help='Generate error insertion profile output',
                          action='store_true', default=False)
    parser_re.add_argument('--no_metrics', help='Skip generating metrics files',
                          action='store_true', default=False)
    parser_re.add_argument('--no_ids', help='Generate sequences without id',
                          action='store_true', default=False)
    parser_re.add_argument('--unweighted', help='Use the same number of reads for every entry in Multifasta files',
                          action='store_true', default=False)
    parser_re.add_argument('-t', '--num_threads', help='Number of threads'
                            '(Default = 1)', type=int, default=1)
    
    host = subparsers.add_parser('host', help="Host/pathogen samples")

    subhost = host.add_subparsers(help="You may run host on standard or random mode",
                                  dest='mode', description='''
                    For detailed usage of each mode:
                        MeStanG.py host mode -h
                    ''')
    
    parser_sh = subhost.add_parser('st', help="Generate standard host/pathogen sample")

    parser_sh.add_argument('-f', '--file', help='Input file list', required=True)
    parser_sh.add_argument('-tx', '--taxon', help='File with pathogens taxon designation (Optional)'
                          , default=None)
    parser_sh.add_argument('-n', '--number', help='Number of reads to be simulated for each sequence (Default = 100)',
                          type=int, default=100)
    parser_sh.add_argument('-m', '--mean', help='Average read length (Default = 2000)', type=int, default=2000)
    parser_sh.add_argument('-sd', '--sd_len', help='Standard deviation of read length in normal scale (Default = 0)',
                          type=float, default=0)
    parser_sh.add_argument('-s', '--posrate',
                          help='Positive sense read rate (Default = 1.0)', type=float,
                          default=1.0)
    parser_sh.add_argument('-p', '--profile',
                          help='Basecaller profile error model (Default = fast)',
                          choices=["virus", "bacteria_ec", "bacteria_kpn", "bacteria", "human", "consensus", "perfect", "fast", "hac", "sup", "res"], default='fast')
    parser_sh.add_argument('-b', '--basecaller',
                          help='Basecaller error model (Default = dorado)',
                          choices=["dorado", "guppy"], default='dorado')
    parser_sh.add_argument('--circular', help='Set genome structure as circular',
                          action='store_true', default=False)
    parser_sh.add_argument('-c', '--custom', help='File with custom basecalling error model', default=None)
    parser_sh.add_argument('-em', '--em_model', help='File with custom base emission model (Optional)'
                          , default=None)
    parser_sh.add_argument('-o', '--output', help='Output directory (Default = MeStanG_output)',
                          default='MeStanG_output')
    parser_sh.add_argument('--equally_p', help='Distribute abundance of a pathogen taxon among subtaxa equally',
                          action='store_true', default=False)
    parser_sh.add_argument('--equally_h', help='Distribute abundance of host among subtaxa equally',
                          action='store_true', default=False)
    parser_sh.add_argument('--error_profile', help='Generate error insertion profile output',
                          action='store_true', default=False)
    parser_sh.add_argument('--no_metrics', help='Skip generating metrics files',
                          action='store_true', default=False)
    parser_sh.add_argument('--no_ids', help='Generate sequences without id',
                          action='store_true', default=False)
    parser_sh.add_argument('--unweighted', help='Use the same number of reads for every entry in Multifasta files',
                          action='store_true', default=False)
    parser_sh.add_argument('-t', '--num_threads', help='Number of threads'
                            '(Default = 1)', type=int, default=1)

    parser_rh = subhost.add_parser('rd', help="Generate randomly distributed Host/pathogen sample")

    parser_rh.add_argument('-f', '--file', help='Input file with parameters', required=True)
    parser_rh.add_argument('-tx', '--taxon', help='File with pathogens taxon designation', required=True)
    parser_rh.add_argument('-ct', '--ctype', help='Sample type', choices=["t", "t+d", "d"])
    parser_rh.add_argument('-n', '--number', help='Number of reads to be simulated for each sequence (Default = 100)',
                          type=int, default=100)
    parser_rh.add_argument('-m', '--mean', help='Average read length (Default = 2000)', type=int, default=2000)
    parser_rh.add_argument('-sd', '--sd_len', help='Standard deviation of read length in normal scale (Default = 0)',
                          type=float, default=0)
    parser_rh.add_argument('-maxr', '--max_ratio', help='Maximum non-host ratio (Default = 0.5)',
                            type=float, default=0.5)
    parser_rh.add_argument('-minr', '--min_ratio', help='Minimum non-host ratio (Default = 0)',
                            type=float, default=0)
    parser_rh.add_argument('-s', '--posrate',
                          help='Positive sense read rate (Default = 1.0)', type=float,
                          default=1.0)
    parser_rh.add_argument('-p', '--profile',
                          help='Basecaller profile error model (Default = fast)',
                          choices=["virus", "bacteria_ec", "bacteria_kpn", "bacteria", "human", "consensus", "perfect", "fast", "hac", "sup", "res"], default='fast')
    parser_rh.add_argument('-b', '--basecaller',
                          help='Basecaller error model (Default = dorado)',
                          choices=["dorado", "guppy"], default='dorado')
    parser_rh.add_argument('--circular', help='Set genome structure as circular',
                          action='store_true', default=False)
    parser_rh.add_argument('-c', '--custom', help='File with custom basecalling error model', default=None)
    parser_rh.add_argument('-em', '--em_model', help='File with custom base emission model (Optional)'
                          , default=None)
    parser_rh.add_argument('-o', '--output', help='Output directory (Default = MeStanG_output)',
                          default='MeStanG_output')
    parser_rh.add_argument('--equally_p', help='Distribute abundance of a pathogen taxon among subtaxa equally',
                          action='store_true', default=False)
    parser_rh.add_argument('--equally_h', help='Distribute abundance of host equally',
                          action='store_true', default=False)
    parser_rh.add_argument('--error_profile', help='Generate error insertion profile output',
                          action='store_true', default=False)
    parser_rh.add_argument('--no_metrics', help='Skip generating metrics files',
                          action='store_true', default=False)
    parser_rh.add_argument('--no_ids', help='Generate sequences without id',
                          action='store_true', default=False)
    parser_rh.add_argument('--unweighted', help='Use the same number of reads for every entry in Multifasta files',
                          action='store_true', default=False)
    parser_rh.add_argument('-t', '--num_threads', help='Number of threads'
                            '(Default = 1)', type=int, default=1)


    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    
    if args.sample == "env":
        if args.mode == "st":
            sys.stdout.write("\n" + strftime("%Y-%m-%d %H:%M:%S") + ": Beginning process\n")

            master = args.file
            support = args.taxon
            master_out = args.output
            number = args.number
            equally = args.equally
            no_ids = args.no_ids
            num_threads = int(max(args.num_threads, 1))

            files, header, ratios = files_to_inputs_se(master, support, number, equally)
            total = int(sum(ratios))
            
            os.mkdir(master_out)
            name_o = "./" + master_out + "/" + master_out

            for i in range(len(files)):
                stat = []
                file = files[i][header.index('file')]
                number = ratios[i]
                if number != 0:
                    stat.append(number)
                    profile, basecaller, circular, custom, error_profile, no_metrics, unweighted, em_model = validate_header(stat, header, files, args, parser, i)
                    display_console(file, stat, profile, basecaller, circular, custom, error_profile, no_metrics, unweighted, num_threads, name_o, equally, em_model)
                    denovo(file, stat, profile, basecaller, circular, name_o, error_profile, custom, no_metrics, unweighted, num_threads, em_model)
            
            get_abundance_file(name_o, ratios, files, header, total)

            if no_ids:
                remove_ids(master_out)
            
            sys.stdout.write("\n" + strftime("%Y-%m-%d %H:%M:%S") + ": Finished Process!\n")
            sys.stdout.close()
        
        if args.mode == "rd":
            sys.stdout.write("\n" + strftime("%Y-%m-%d %H:%M:%S") + ": Beginning process\n")

            master = args.file
            support = args.taxon
            master_out = args.output
            number = args.number
            equally = args.equally
            no_ids = args.no_ids
            num_threads = int(max(args.num_threads, 1))

            files, header, ratios = files_to_inputs_re(master, number, equally, master_out)
            total = int(sum(ratios))
            
            os.mkdir(master_out)
            name_o = "./" + master_out + "/" + master_out

            for i in range(len(files)):
                stat = []
                file = files[i][header.index('file')]
                number = ratios[i]                
                if number != 0:
                    stat.append(number)
                    profile, basecaller, circular, custom, error_profile, no_metrics, unweighted, em_model = validate_header(stat, header, files, args, parser, i)
                    display_console(file, stat, profile, basecaller, circular, custom, error_profile, no_metrics, unweighted, num_threads, name_o, equally, em_model)
                    denovo(file, stat, profile, basecaller, circular, name_o, error_profile, custom, no_metrics, unweighted, num_threads, em_model)
            
            get_abundance_file(name_o, ratios, files, header, total)

            if no_ids:
                remove_ids(master_out)
            
            sys.stdout.write("\n" + strftime("%Y-%m-%d %H:%M:%S") + ": Finished Process!\n")
            sys.stdout.close()

    if args.sample == "host":
        if args.mode == "st":
            sys.stdout.write("\n" + strftime("%Y-%m-%d %H:%M:%S") + ": Beginning process\n")

            master = args.file
            support = args.taxon
            master_out = args.output
            number = args.number
            equally_p = args.equally_p
            equally_h = args.equally_h
            no_ids = args.no_ids
            num_threads = int(max(args.num_threads, 1))

            files, header, ratios = files_to_inputs_sh(master, support, number, equally_p, equally_h)
            total = int(sum(ratios))
            
            os.mkdir(master_out)
            name_o = "./" + master_out + "/" + master_out

            for i in range(len(files)):
                stat = []
                file = files[i][header.index('file')]
                number = ratios[i]                
                if number != 0:
                    stat.append(number)
                    profile, basecaller, circular, custom, error_profile, no_metrics, unweighted, em_model = validate_header(stat, header, files, args, parser, i)                    
                    if files[i][header.index('taxon')] == 'host':
                        if equally_h:
                            equally = True
                        else:
                            equally = False
                    else:
                        if equally_p:
                            equally = True
                        else:
                            equally = False
                    display_console(file, stat, profile, basecaller, circular, custom, error_profile, no_metrics, unweighted, num_threads, name_o, equally, em_model)
                    denovo(file, stat, profile, basecaller, circular, name_o, error_profile, custom, no_metrics, unweighted, num_threads, em_model)
            
            get_abundance_file(name_o, ratios, files, header, total)

            if no_ids:
                remove_ids(master_out)
            
            sys.stdout.write("\n" + strftime("%Y-%m-%d %H:%M:%S") + ": Finished Process!\n")
            sys.stdout.close()
        
        if args.mode == "rd":
            sys.stdout.write("\n" + strftime("%Y-%m-%d %H:%M:%S") + ": Beginning process\n")

            master = args.file
            support = args.taxon
            community = args.ctype
            max_ratio = args.max_ratio
            min_ratio = args.min_ratio
            master_out = args.output
            number = args.number
            equally_p = args.equally_p
            equally_h = args.equally_h
            no_ids = args.no_ids
            num_threads = int(max(args.num_threads, 1))

            files, header, ratios = files_to_inputs_rh(master, number, equally_p, equally_h, master_out, support, community, min_ratio, max_ratio)
            total = int(sum(ratios))
            
            os.mkdir(master_out)
            name_o = "./" + master_out + "/" + master_out + "_" + str(community)

            for i in range(len(files)):
                stat = []
                file = files[i][header.index('file')]
                number = ratios[i]
                if number != 0:
                    stat.append(number)
                    profile, basecaller, circular, custom, error_profile, no_metrics, unweighted, em_model = validate_header(stat, header, files, args, parser, i)
                    if files[i][header.index('taxon')] == 'host':
                        if equally_h:
                            equally = True
                        else:
                            equally = False
                    else:
                        if equally_p:
                            equally = True
                        else:
                            equally = False
                    display_console(file, stat, profile, basecaller, circular, custom, error_profile, no_metrics, unweighted, num_threads, name_o, equally, em_model)
                    denovo(file, stat, profile, basecaller, circular, name_o, error_profile, custom, no_metrics, unweighted, num_threads, em_model)
            
            get_abundance_file(name_o, ratios, files, header, total)

            if no_ids:
                remove_ids(master_out)
            
            sys.stdout.write("\n" + strftime("%Y-%m-%d %H:%M:%S") + ": Finished Process!\n")
            sys.stdout.close()

if __name__ == "__main__":
    main()
