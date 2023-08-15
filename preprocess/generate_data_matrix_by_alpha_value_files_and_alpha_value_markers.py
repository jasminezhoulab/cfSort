#
# Python 3.x
#
# Deployed to the working directory: /u/project/xjzhou/wenyuan/data/reads_alpha_distribution_pipeline
#
import sys, gzip, re, string
import numpy as np
from scipy import sparse
from scipy.sparse import lil_matrix, csc_matrix
import pandas as pd
from datetime import datetime
from itertools import combinations
import collections
from collections import Counter
from operator import itemgetter
import seaborn as sns

def read_lines(file):
    lines = []
    if file == 'stdin':
        f = sys.stdin
    elif file.endswith('gz'):
        f = gzip.open(file, 'rt')
    else:
        f = open(file, 'rt')
    for line in f:
        lines.append(line.rstrip())
    if file != 'stdin':
        f.close()
    return(lines)

def write_lines(fout, lines, header_line=None):
    if header_line is not None:
        fout.write('%s\n'%header_line)
    for line in lines:
        fout.write('%s\n'%line)

def count_nonzero_nonnan_max_of_rows(mat):
    _, ncol = mat.shape
    mat_type = str(type(mat))
    if ('lil_matrix' in mat_type) or ('csc_matrix' in mat_type) or ('csr_matrix' in mat_type) or (
            'coo_matrix' in mat_type) or ('bsr_matrix' in mat_type):
        mat_np_array = mat.toarray()
        nonzero_counts_of_rows = np.count_nonzero(mat_np_array > 0, axis=1)
        nonnan_counts_of_rows = np.count_nonzero(~np.isnan(mat_np_array), axis=1)
        max_of_rows = np.nanmax(mat_np_array, axis=1)
    elif 'numpy.ndarray' in mat_type:
        nonzero_counts_of_rows = np.count_nonzero(mat>0, axis=1)
        nonnan_counts_of_rows = np.count_nonzero(~np.isnan(mat), axis=1)
        max_of_rows = np.nanmax(mat, axis=1)
    nonzero_frac_of_rows = nonzero_counts_of_rows / float(ncol)
    nonnan_frac_of_rows = nonnan_counts_of_rows / float(ncol)
    return((nonzero_counts_of_rows, nonzero_frac_of_rows, nonnan_counts_of_rows, nonnan_frac_of_rows, max_of_rows))

def write_int_matrix_with_row_labels(mat, row_names_list, fid):
    nrow, ncol = mat.shape
    for i in range(nrow):
        profile_str = '\t'.join(['%.3g' % v for v in mat[i, :]])
        profile_str = profile_str.replace('nan', 'NA')
        fid.write('%s\t%s\n'%(row_names_list[i], profile_str))

def normalize_matrix_by_rows(mat, rows_factor):
    mat_type = str(type(mat))
    if ('lil_matrix' in mat_type) or ('csc_matrix' in mat_type) or ('csr_matrix' in mat_type) or (
            'coo_matrix' in mat_type) or ('bsr_matrix' in mat_type):
        new_mat = mat.toarray() * rows_factor
    else:
        new_mat = mat * rows_factor
    return(new_mat)

def parse_markers_with_dynamic_alpha_thresholds(file):
    markers_list = {'marker_index':[], 'alpha_threshold':[], 'header_line':None, 'lines':[]}
    if file.endswith('gz'):
        f = gzip.open(file, 'rt')
    elif file == 'stdin':
        f = sys.stdin
    else:
        f = open(file)
    markers_list['header_line'] = next(f).rstrip() # skip the header line
    for line in f:
        items = line.rstrip().split('\t')
        marker_index = int(items[0])
        markers_list['marker_index'].append( marker_index )
        markers_list['alpha_threshold'].append( float(items[1]) )
        markers_list['lines'].append( line.rstrip() )
    if file != 'stdin':
        f.close()
    return(markers_list)


def get_number_of_reads_meeting_criterion(unique_alpha_values, read_freq_of_unique_alpha_values, alpha_threshold, direction='<='):
    read_num = 0
    if direction=='<=':
        read_num = np.sum(read_freq_of_unique_alpha_values[unique_alpha_values <= alpha_threshold])
    elif direction == '>=':
        read_num = np.sum(read_freq_of_unique_alpha_values[unique_alpha_values >= alpha_threshold])
    elif direction == '<':
        read_num = np.sum(read_freq_of_unique_alpha_values[unique_alpha_values < alpha_threshold])
    elif direction == '>':
        read_num = np.sum(read_freq_of_unique_alpha_values[unique_alpha_values > alpha_threshold])
    elif direction == '==':
        read_num = np.sum(read_freq_of_unique_alpha_values[unique_alpha_values == alpha_threshold])
    return(read_num)


def generate_plasma_sample_profile_with_given_markers_and_dynamic_alpha_threshold_by_parsing_alpha_value_distribution_file_quick_version(file, markers_list, marker_type):
    n_markers = len(markers_list['marker_index'])
    profile = np.empty(n_markers, dtype=float)
    profile[:] = np.nan
    total_reads = 0
    current_marker_pointer_in_markers_list = 0
    with gzip.open(file, 'rt') as f:
        next(f) # skip the header line
        line = f.readline()
        while line != "":
            marker_index, max_num_cpg, num_read, unique_alpha_values_str, read_freq_of_unique_alpha_values_str, _, _ = line.rstrip().split('\t')
            marker_index = int(marker_index)
            if current_marker_pointer_in_markers_list >= n_markers:
                break
            if marker_index > markers_list['marker_index'][current_marker_pointer_in_markers_list]:
                while True:
                    current_marker_pointer_in_markers_list += 1
                    if current_marker_pointer_in_markers_list >= n_markers:
                        break
                    if marker_index <= markers_list['marker_index'][current_marker_pointer_in_markers_list]:
                        break
                continue
            if current_marker_pointer_in_markers_list >= n_markers:
                break
            if marker_index < markers_list['marker_index'][current_marker_pointer_in_markers_list]:
                line = f.readline()
                continue
            # now 'marker_index == markers_list_int[current_marker_pointer_in_markers_list]'
            # idx = markers_list['marker_index'].index(marker_index)
            idx = current_marker_pointer_in_markers_list
            Ta = markers_list['alpha_threshold'][idx]
            current_marker_pointer_in_markers_list += 1 # for the next round comparison
            unique_alpha_values = np.array(list(map(float, unique_alpha_values_str.split(','))))
            read_freq_of_unique_alpha_values = np.array(
                list(map(int, read_freq_of_unique_alpha_values_str.split(','))))
            # identify tumor reads in tumor
            if 'hyper' in marker_type:
                # tumor reads in tumor & plasma are defined as 'alpha>Ta';
                read_cov_of_criterion = get_number_of_reads_meeting_criterion(unique_alpha_values,
                                                                              read_freq_of_unique_alpha_values,
                                                                              Ta,
                                                                              '>')
            elif 'hypo' in marker_type:
                # tumor reads in tumor & plasma are defined as 'alpha<Ta';
                read_cov_of_criterion = get_number_of_reads_meeting_criterion(unique_alpha_values,
                                                                              read_freq_of_unique_alpha_values,
                                                                              Ta,
                                                                              '<')
            profile[idx] = read_cov_of_criterion
            total_reads += int(num_read)
    return( (profile, total_reads) )


def generate_matrix_for_many_plasma_samples_with_given_markers_and_dynamic_alpha_threshold_by_parsing_alpha_value_distribution_files(input_alpha_value_distribution_files_list,
                                                                                                                                     markers_list,
                                                                                                                                     marker_type):
    n_sample = len(input_alpha_value_distribution_files_list)
    total_reads_list = []
    for i in range(n_sample):
        input_sample_alpha_distr_file = input_alpha_value_distribution_files_list[i]
        print('  (%d/%d) %s'%(i+1, n_sample, input_sample_alpha_distr_file), end="\t", flush=True)
        profile, total_reads = generate_plasma_sample_profile_with_given_markers_and_dynamic_alpha_threshold_by_parsing_alpha_value_distribution_file_quick_version(input_sample_alpha_distr_file,
                                                                                                                                                       markers_list,
                                                                                                                                                       marker_type)
        print(datetime.now(), flush=True)
        if i==0:
            # data = csr_matrix((n_sample, len(profile)), dtype=np.single).toarray()
            data = lil_matrix((n_sample, len(profile)), dtype=np.single)
        data[i,:] = profile
        total_reads_list.append(total_reads)
    return( (data, np.array(total_reads_list)) )


marker_type = sys.argv[1]
input_markers_file_with_dynamic_alpha_threshold = sys.argv[2]
input_plasma_samples_list_file = sys.argv[3]
output_raw_matrix_file = sys.argv[4]
output_normalized_matrix_file = ''
if len(sys.argv) == 6:
    output_normalized_matrix_file = sys.argv[5]

normalize_type = 'None'  # default is not to use normalization.
if len(output_normalized_matrix_file) > 0:
    all_normalized_types_list = 'cpm.bpm'.split('.')
    for t in all_normalized_types_list:
        if t in output_normalized_matrix_file:
            normalize_type = t
            break
    if normalize_type == 'None':
        sys.stderr.write('Error: argument output_normalized_matrix_file (see below) does not contain [%s].\nExit.\n'%(', '.join(all_normalized_types_list)))


print('REQUIRE: marker_index in both markers_file and alpha_values_distribution_files are sorted in increasing order of integeter!')

print('Parse samples list file\n  %s'%input_plasma_samples_list_file, flush=True)
samples_info = {'samples':[], 'alpha_values_files':[]}
samples_info['samples'] = read_lines(input_plasma_samples_list_file)
for s in samples_info['samples']:
    samples_info['alpha_values_files'].append( '%s'%(s) )
print('  #samples: %d'%len(samples_info['samples']), flush=True)

print('Parse markers list file\n  %s'%input_markers_file_with_dynamic_alpha_threshold, flush=True)
markers_list = parse_markers_with_dynamic_alpha_thresholds(input_markers_file_with_dynamic_alpha_threshold)
print('  #markers: %d'%len(markers_list['lines']))

print('Parse alpha value distribution file of each sample', flush=True)
data, total_reads_of_samples = generate_matrix_for_many_plasma_samples_with_given_markers_and_dynamic_alpha_threshold_by_parsing_alpha_value_distribution_files(samples_info['alpha_values_files'],
                                                                                                                                                                markers_list,
                                                                                                                                                                marker_type)

print('  #original_markers: %d'%len(markers_list['lines']))

print('Output raw matrix file\n  matrix: %s' % (output_raw_matrix_file), flush=True)
with gzip.open(output_raw_matrix_file, 'wt') as fout:
    write_int_matrix_with_row_labels(data.toarray(), samples_info['samples'], fout)

if len(output_normalized_matrix_file) > 0:
    print('\n==========\nNormalize raw matrix and output ...\n')
    if normalize_type == 'cpb':
        rows_factor = 10**9 / total_reads_of_samples
        print('Normalize row factor: %s (10^9 / total_read_count)' % normalize_type, flush=True)
    elif normalize_type == 'cpm':
        rows_factor = 10**6 / total_reads_of_samples
        print('Normalize row factor: %s (10^6 / total_read_count)' % normalize_type, flush=True)
    else:
        sys.stderr.write('Error: normalize_type should be cpb (count per billion reads) or cpm (count per million reads)\nExit.\n')
        sys.exit(-1)
    rows_factor[rows_factor==np.inf] = 1 # For those rows whose samples' total reads are zero, row_factor is set to one.
    rows_factor = rows_factor.reshape((len(rows_factor), 1)) # Convert to array of size nrow X 1
    data_normalized = normalize_matrix_by_rows(data, rows_factor)

    print('Output normalized matrix file\n  matrix: %s' % (output_normalized_matrix_file), flush=True)
    with gzip.open(output_normalized_matrix_file, 'wt') as fout:
        write_int_matrix_with_row_labels(data_normalized, samples_info['samples'], fout)

print('Done.')
