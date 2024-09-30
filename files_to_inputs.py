#!/usr/bin/env python
"""
Module of MeStanG v0.1
@author: Daniel Ramos Lopez
Metagenomic Standards Generator (MeStanG) for HTS Nanopore datasets 
"""

import sys
import os
import random
from time import strftime
import distribution_and_ratio as dnr

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
            big_ratios = dnr.ratio_e(number, m_ratios)
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
            sub_ratios = dnr.dist_ratio(big_ratios[i], len(subtaxa_for_taxon), equally)
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
            ratios = dnr.ratio_e(number, m_ratios)
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
            big_ratios = dnr.ratio_h(number, pathogen_ratios)
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
            sub_ratios = dnr.dist_ratio(big_ratios[i], len(subtaxa_for_taxon), equally_p)
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
        hosts_ratios = dnr.dist_ratio(big_ratios[len(big_ratios) - 1], len(hosts_list), equally_h)

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
            ratios = dnr.ratio_h(number, pathogen_ratios)
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

        hosts_ratios = dnr.dist_ratio(ratios[len(ratios) - 1], len(hosts_list), equally_h)
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

    r_ratios = dnr.dist_ratio_communities(1, len(taxa_list))
    reads_per_member = dnr.ratio_e(number, r_ratios)

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
            ratios = dnr.dist_ratio_communities(big_ratio, len(targets))
            zero_ratios = 0
            max_reads = 0
            for i in range(len(ratios)):
                current_reads = dnr.ratio_h(number, [ratios[i]])[0]
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
            ratios = dnr.dist_ratio_communities(big_ratio, len(decoys))
            zero_ratios = 0
            max_reads = 0
            for i in range(len(ratios)):
                current_reads = dnr.ratio_h(number, [ratios[i]])[0]
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
            ratios = dnr.dist_ratio_communities(big_ratio, len(community_list))
            max_reads = 0
            for i in range(len(ratios)):
                current_reads = dnr.ratio_h(number, [ratios[i]])[0]
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
