#ifndef UTILS_H
#define UTILS_H

#include <sstream>
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include "data_types.h"
#include "matrix.h"

using namespace std;

void multi_vec_by_number(vector<double>& vec, double v);
bool assign_vector_zeros(vector<double>& vec);
bool str_ends_with(string& str, string& suffix);
void make_complete_chromosomes(vector<string> & all_chrs);

void read_list_strings_from_file(string input_file, vector<string>& strs);
void read_two_columns_of_list_strings_from_file(string input_file, vector<string>& strs1, vector<string>& strs2);
void read_two_columns_of_list_strings_from_file(string input_file, vector<string>& strs1, vector<string>& strs2,
	map<string, vector<string> > & map_str1tostrs2, string delimit);

void read_tissue_markers_txt_file(string tissue_markers_file, int value_column_start_index, int num_tissue_types, Bins2Values & bins2values, vector<string> & value_names);
void read_tissue_markers_gz_file(string tissue_markers_file, int value_column_start_index, int num_tissue_types, Bins2Values & bins2values, vector<string> & value_names);

unsigned long calc_read_probability_by_marker2beta(string reads_binning_file, Bins2Values & marker2beta, Matrix_Double& reads_likelihoods, Bins2UnsignedIntegers& marker2rowindexes, Bins2Value& marker2ambiguousreadcounts, vector<int>&Rm, vector<int>& Rl, double likelihood_ratio_cutoff);

void get_reads_methy_data_from_reads_binning_file(string reads_binning_file,
	vector<int> & num_cpg_sites, vector<int> & num_methy_cpg_sites);

void print_vec_of_uint(ostream& of, vector<unsigned int> & v);
void print_vec_of_ulong(ostream& of, vector<unsigned long> & v);
void print_map_of_strings(ostream& of, map<string, vector<string> > & map_str1tostrs2);
void print_str_vectors(ofstream& o, vector<string> & s);
void strings2floats(string str, vector<float> & vec, string delimit=",");
float min(float a, float b);
float max(float a, float b);
float meta_p(vector<float> & probs);

#endif
