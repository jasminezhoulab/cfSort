import sys, gzip, re
import numpy as np
from scipy.sparse import lil_matrix
import pandas as pd
from datetime import datetime
from itertools import combinations
import collections
from collections import Counter
from operator import itemgetter

def convert_int2int_dict_to_str(d, sep=','):
    keys_list = sorted(d.keys())
    int_keys_str = sep.join(list(map(str, keys_list)))
    int_values_str = sep.join(['%d'%d[k] for k in keys_list])
    return((int_keys_str, int_values_str))

def convert_str2int_dict_to_str(d, sep=','):
    keys_list = sorted(d.keys())
    keys_str = sep.join(keys_list)
    int_values_str = sep.join(['%d'%d[k] for k in keys_list])
    return((keys_str, int_values_str))

def convert_methylation_str_list_to_distribution(meth_str_list):
    distribution_meth_counts = {}
    distribution_alpha_values = {}
    max_len = 0
    for m in meth_str_list:
        meth_count = m.count('1')
        alpha_str = '%.3g'%(meth_count / float(len(m)))
        if meth_count not in distribution_meth_counts:
            distribution_meth_counts[meth_count] = 1
        else:
            distribution_meth_counts[meth_count] += 1
        if alpha_str not in distribution_alpha_values:
            distribution_alpha_values[alpha_str] = 1
        else:
            distribution_alpha_values[alpha_str] += 1
        if max_len<len(m):
            max_len = len(m)
    return((distribution_meth_counts, distribution_alpha_values, max_len))

def summarize_mary_file_binary_meth_values_for_distribution_file(input_methy_reads_binning_file, output_distribution_file):
    with gzip.open(input_methy_reads_binning_file,'rt') as fin, gzip.open(output_distribution_file, 'wt') as fout:
        next(fin) # skip the first header line
        # distribution_of_marker = {'marker_index': None, 'num_cpg': 0,
        #                           'distribution': {'strand+':{'num_read': 0, 'unique_meth_counts_to_freq':{}},
        #                                            'strand-':{'num_read': 0, 'unique_meth_counts_to_freq':{}}
        #                                            }
        #                           }
        fout.write('marker_index\tmax_num_cpg\tnum_read\tunique_alpha_values\tread_freq_of_unique_alpha_values\tunique_meth_counts\tread_freq_of_unique_meth_counts\n')
        distribution_of_marker = {'marker_index': -1,
                                  'max_num_cpg': 0,
                                  'methy_str_list': [],
                                  'unique_meth_counts_to_freq': {},
                                  'unique_alpha_values_to_freq': {},
                                  }
        for line in fin:
            if 'NA' in line:
                continue
            marker_index, _, meth_string, _, _, _ = line.rstrip().split('\t')
            if marker_index != distribution_of_marker['marker_index']:
                # A new marker begins, we need to print the old marker
                if distribution_of_marker['marker_index']!=-1:
                    num_reads = len(distribution_of_marker['methy_str_list'])
                    distribution_of_marker['unique_meth_counts_to_freq'], distribution_of_marker['unique_alpha_values_to_freq'], distribution_of_marker['max_num_cpg'] = convert_methylation_str_list_to_distribution(distribution_of_marker['methy_str_list'])
                    # Print 'distribution_of_marker' to output
                    unique_meth_count_str, read_freq_of_unique_meth_counts_str = convert_int2int_dict_to_str(distribution_of_marker['unique_meth_counts_to_freq'])
                    unique_alpha_values_str, read_freq_of_unique_alpha_values_str = convert_str2int_dict_to_str(
                        distribution_of_marker['unique_alpha_values_to_freq'])
                    fout.write('%s\t%d\t%d\t%s\t%s\t%s\t%s\n'%(distribution_of_marker['marker_index'],
                                                               distribution_of_marker['max_num_cpg'],
                                                               num_reads,
                                                               unique_alpha_values_str,
                                                               read_freq_of_unique_alpha_values_str,
                                                               unique_meth_count_str,
                                                               read_freq_of_unique_meth_counts_str,
                                                           ))
                # Clear the old marker info, and initialize read counting for a new marker
                distribution_of_marker = {'marker_index': marker_index,
                                          'max_num_cpg': 0,
                                          'methy_str_list': [],
                                          'unique_meth_counts_to_freq': {},
                                          'unique_alpha_values_to_freq': {},
                                          }
            distribution_of_marker['methy_str_list'].append(meth_string)

        # At the end of the input file
        # write the last marker of the file
        if distribution_of_marker['marker_index'] != -1:
            num_reads = len(distribution_of_marker['methy_str_list'])
            distribution_of_marker['unique_meth_counts_to_freq'], distribution_of_marker['unique_alpha_values_to_freq'], \
            distribution_of_marker['max_num_cpg'] = convert_methylation_str_list_to_distribution(
                distribution_of_marker['methy_str_list'])
            # Print 'distribution_of_marker' to output
            unique_meth_count_str, read_freq_of_unique_meth_counts_str = convert_int2int_dict_to_str(
                distribution_of_marker['unique_meth_counts_to_freq'])
            unique_alpha_values_str, read_freq_of_unique_alpha_values_str = convert_str2int_dict_to_str(
                distribution_of_marker['unique_alpha_values_to_freq'])
            fout.write('%s\t%d\t%d\t%s\t%s\t%s\t%s\n' % (distribution_of_marker['marker_index'],
                                                         distribution_of_marker['max_num_cpg'],
                                                         num_reads,
                                                         unique_alpha_values_str,
                                                         read_freq_of_unique_alpha_values_str,
                                                         unique_meth_count_str,
                                                         read_freq_of_unique_meth_counts_str,
                                                         ))


debug = False
# debug = True
if debug:
    input_reads_binning_file = './input.data/binary_meth_values_files/debug.binary_meth_values.txt.gz'
    output_alpha_value_distribution_file = './debug.alpha_distr.txt.gz'
else:
    input_reads_binning_file = sys.argv[1]
    output_alpha_value_distribution_file = sys.argv[2]

print('input: %s'%input_reads_binning_file)
print('output: %s'%output_alpha_value_distribution_file)
summarize_mary_file_binary_meth_values_for_distribution_file(input_reads_binning_file, output_alpha_value_distribution_file)
print('py_done')
