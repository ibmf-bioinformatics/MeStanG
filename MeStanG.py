#!/usr/bin/env python
"""
MeStanG v0.1
@author: Daniel Ramos Lopez
Metagenomic Standards Generator (MeStanG) for HTS Nanopore datasets 
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
from Bio.Seq import Seq

import files_to_inputs as fti
import file_manager as fm
import distribution_and_ratio as dnr

VERSION = "0.1"
PROGRAM = "MeStanG"
AUTHOR = "Daniel Ramos Lopez"
CONTACT = "aramosl@okstate.edu, adramoslp@outlook.com"


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
                        base = fm.base_conversion(perfect_read[profiles[0][i] + j])
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
                        base = fm.base_conversion(base)
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
        error_prob = fm.error_selector(profile, basecaller)
    else:
        error_prob = fm.get_custom_profile(custom)
    mis, ins = fm.model_loader(em_model)
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Processing file\n")
    file_l = fm.fasta_to_list(file)
    if len(file_l) == 2:
        unweighted = False
        weighted_reads = [int(stat[0])]
    else:
        if not unweighted:
            weighted_reads = dnr.dist_weight(stat[0], file_l)
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
                        met = fm.read_summary(fname)
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

    parser_se.add_argument('-f', '--file', help='Input file with parameters', required=True)
    parser_se.add_argument('-tx', '--taxon', help='File with organism taxon designation (Optional)'
                          , default=None)
    parser_se.add_argument('-n', '--number', help='Number of reads to be simulated for each sequence (Default = 100)',
                          type=int, default=100)
    parser_se.add_argument('-m', '--mean', help='Average read length (Default = 2000)', type=int, default=2000)
    parser_se.add_argument('-sd', '--sd_len', help='Standard deviation of read length in normal scale (Default = 0)',
                          type=int, default=0)
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
    parser_se.add_argument('--equally', help='Distribute equally the ratio of a taxon in its subtaxa',
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

    parser_re.add_argument('-f', '--file', help='Input file with parameters', required=True)
    parser_re.add_argument('-n', '--number', help='Number of reads to be simulated for each sequence (Default = 100)',
                          type=int, default=100)
    parser_re.add_argument('-m', '--mean', help='Average read length (Default = 2000)', type=int, default=2000)
    parser_re.add_argument('-sd', '--sd_len', help='Standard deviation of read length in normal scale (Default = 0)',
                          type=int, default=0)
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

    parser_sh.add_argument('-f', '--file', help='Input file with parameters', required=True)
    parser_sh.add_argument('-tx', '--taxon', help='File with pathogens taxon designation (Optional)'
                          , default=None)
    parser_sh.add_argument('-n', '--number', help='Number of reads to be simulated for each sequence (Default = 100)',
                          type=int, default=100)
    parser_sh.add_argument('-m', '--mean', help='Average read length (Default = 2000)', type=int, default=2000)
    parser_sh.add_argument('-sd', '--sd_len', help='Standard deviation of read length in normal scale (Default = 0)',
                          type=int, default=0)
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
                          type=int, default=0)
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
    parser_rh.add_argument('--equally_h', help='Distribute abundance of host among subtaxa equally',
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

            files, header, ratios = fti.files_to_inputs_se(master, support, number, equally)
            total = int(sum(ratios))
            
            os.mkdir(master_out)
            name_o = "./" + master_out + "/" + master_out

            for i in range(len(files)):
                stat = []
                file = files[i][header.index('file')]
                number = ratios[i]
                if number != 0:
                    stat.append(number)
                    profile, basecaller, circular, custom, error_profile, no_metrics, unweighted, em_model = fm.validate_header(stat, header, files, args, parser, i)
                    fm.display_console(file, stat, profile, basecaller, circular, custom, error_profile, no_metrics, unweighted, num_threads, name_o, equally, em_model)
                    denovo(file, stat, profile, basecaller, circular, name_o, error_profile, custom, no_metrics, unweighted, num_threads, em_model)
            
            fm.get_abundance_file(name_o, ratios, files, header, total)

            if no_ids:
                fm.remove_ids(master_out)
            
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

            files, header, ratios = fti.files_to_inputs_re(master, number, equally, master_out)
            total = int(sum(ratios))
            
            os.mkdir(master_out)
            name_o = "./" + master_out + "/" + master_out

            for i in range(len(files)):
                stat = []
                file = files[i][header.index('file')]
                number = ratios[i]                
                if number != 0:
                    stat.append(number)
                    profile, basecaller, circular, custom, error_profile, no_metrics, unweighted, em_model = fm.validate_header(stat, header, files, args, parser, i)
                    fm.display_console(file, stat, profile, basecaller, circular, custom, error_profile, no_metrics, unweighted, num_threads, name_o, equally, em_model)
                    denovo(file, stat, profile, basecaller, circular, name_o, error_profile, custom, no_metrics, unweighted, num_threads, em_model)
            
            fm.get_abundance_file(name_o, ratios, files, header, total)

            if no_ids:
                fm.remove_ids(master_out)
            
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

            files, header, ratios = fti.files_to_inputs_sh(master, support, number, equally_p, equally_h)
            total = int(sum(ratios))
            
            os.mkdir(master_out)
            name_o = "./" + master_out + "/" + master_out

            for i in range(len(files)):
                stat = []
                file = files[i][header.index('file')]
                number = ratios[i]                
                if number != 0:
                    stat.append(number)
                    profile, basecaller, circular, custom, error_profile, no_metrics, unweighted, em_model = fm.validate_header(stat, header, files, args, parser, i)                    
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
                    fm.display_console(file, stat, profile, basecaller, circular, custom, error_profile, no_metrics, unweighted, num_threads, name_o, equally, em_model)
                    denovo(file, stat, profile, basecaller, circular, name_o, error_profile, custom, no_metrics, unweighted, num_threads, em_model)
            
            fm.get_abundance_file(name_o, ratios, files, header, total)

            if no_ids:
                fm.remove_ids(master_out)
            
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

            files, header, ratios = fti.files_to_inputs_rh(master, number, equally_p, equally_h, master_out, support, community, min_ratio, max_ratio)
            total = int(sum(ratios))
            
            os.mkdir(master_out)
            name_o = "./" + master_out + "/" + master_out + "_" + str(community)

            for i in range(len(files)):
                stat = []
                file = files[i][header.index('file')]
                number = ratios[i]
                if number != 0:
                    stat.append(number)
                    profile, basecaller, circular, custom, error_profile, no_metrics, unweighted, em_model = fm.validate_header(stat, header, files, args, parser, i)
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
                    fm.display_console(file, stat, profile, basecaller, circular, custom, error_profile, no_metrics, unweighted, num_threads, name_o, equally, em_model)
                    denovo(file, stat, profile, basecaller, circular, name_o, error_profile, custom, no_metrics, unweighted, num_threads, em_model)
            
            fm.get_abundance_file(name_o, ratios, files, header, total)

            if no_ids:
                fm.remove_ids(master_out)
            
            sys.stdout.write("\n" + strftime("%Y-%m-%d %H:%M:%S") + ": Finished Process!\n")
            sys.stdout.close()

if __name__ == "__main__":
    main()
