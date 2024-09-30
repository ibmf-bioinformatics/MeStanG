#!/usr/bin/env python
"""
Module of MeStanG v0.1
@author: Daniel Ramos Lopez
Metagenomic Standards Generator (MeStanG) for HTS Nanopore datasets 
"""

import os
import sys
import random
import numpy
import string
from time import strftime
from Bio import SeqIO

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
        stat.append(int(files[i][header.index('sd_len')]))
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
    print("\nRunning error insertion with following parameters:\n")
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
