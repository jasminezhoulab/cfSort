#include <cstdlib>
#include <sstream>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <boost/algorithm/string.hpp>
#include <boost/assign.hpp>
#include "data_types.h"

using namespace std;
using namespace boost;

// 
// bins (or features) annotation format
//
//1. file of bins (features) annotation: "biomarkers.all_bins"
//Each line is a bin (or feature). All columns are delimited by TAB. There is one header line.
//Column 1: marker_index, 1-based index. Those bins that are not markers, are indexed as 0.
//Column 2: chr
//Column 3: start coordinate of read (1-base)
//Column 4: end coordinate of read (1-base). The range of the bin is [start, end)
//Column 5: marker_type. "I" is Marker_type_I, "II" is Marker_type_II, "-" is the complementary bin, only facilitating the searching.
//
// The following is an example file
//
//marker_index    chr     start   end     marker_type
//0       chr1    1       855266  -
//1       chr1    855266  855766  II
//0       chr1    855766  969796  -
//2       chr1    969796  970296  II
//0       chr1    970296  1099044 -
//3       chr1    1099044 1099544 II
//0       chr1    1099544 1109315 -
//4       chr1    1109315 1109815 II
//
void read_bins_annot_file(string input_bins_annot_file, Bins_end_coord & bins_end_coord,
	Bins_index & bins_index, Bins_info & bins_info, bool has_header_line=true)
{
	ifstream fin;
	fin.open(input_bins_annot_file.c_str());
	if (fin.fail()){
        cerr << "Error: Unable to open " << input_bins_annot_file << " in read_bins_annot_file()" << endl;
        exit(EXIT_FAILURE);
    }
	string line;
	if (has_header_line)
		// skip the first header line
		getline(fin, line);
	unsigned long i=0;
	string old_chr;
	while (!fin.eof()) {
		getline(fin, line);
		if (line.empty()) {
			// this is the last line of the file
			break;
		}
		//cout << line << endl;

		vector<string> strs1;
		split(strs1, line, is_any_of("\t"));
		int bin_index = atoi(strs1[0].c_str());
		string chr = strs1[1];
		//int start_coord = atoi(strs1[2].c_str()); // start coordinate
		int end_coord = atoi(strs1[3].c_str()); // end coordinate
		if (i==0) {
			// this is the first bin of all the genome, so we initialize old_chr
			old_chr = chr;
		}
		if (chr.compare(old_chr)!=0) {
			// This is the 1st bin of the new chromosome
			old_chr = chr;
			bins_end_coord.insert(make_pair(chr, vector<unsigned int>()));
			bins_index.insert(make_pair(chr, vector<int >()));
			bins_info.insert(make_pair(chr, vector<string >()));
		}
		vector<unsigned int> & coords = bins_end_coord[chr];
		coords.push_back( end_coord );
		vector<int> & indexes = bins_index[chr];
		indexes.push_back( bin_index );
		vector<string> & infos = bins_info[chr];
		infos.push_back( line );
		i++;
	}
	cerr << "#bins=" << i << endl;
}

int get_num_of_non_void_bins(Bins_index & bins_index, vector<int> & returned_markers_index)
{
	int num_of_non_void_bins = 0;
	Bins_index::iterator it;
	for (it=bins_index.begin(); it!=bins_index.end(); ++it) {
		vector<int> & bins_of_chr = it->second;
		for (int i=0; i<bins_of_chr.size(); i++)
			if (bins_of_chr[i] > 0) {
				num_of_non_void_bins++;
				returned_markers_index.push_back( bins_of_chr[i] );
			}
	}
	return num_of_non_void_bins;
}

// We find out the index of the bin with the range [bin_start_coord, bin_end_coord), where both "bin_start_coord" and "bin_end_coord" are 0-base.
// Returned bin_internal_index is 0-base. If not found, return -1
int find_exact_bin(Bins_end_coord & bins_end_coord, string chr, unsigned int bin_start_coord, unsigned int bin_end_coord) {
	int bin_internal_index=-1;
	vector<unsigned int> & coords_bins_of_chr = bins_end_coord[chr]; // a vector of end coordinates (1-base) of this chr.
	vector<unsigned int>::iterator bin_it = find( coords_bins_of_chr.begin(), coords_bins_of_chr.end(), bin_end_coord);
	if (bin_it!=coords_bins_of_chr.end()) {
		// found
		bin_internal_index = bin_it-coords_bins_of_chr.begin();
	}
	return bin_internal_index;
}

// We find out the index of the bin with the range [bin_start_coord, bin_end_coord), where the input paramter "position" (1-base) falls into this bin.
// Returned bin_internal_index is 0-base. If not found, return -1
int find_bin_of_position(Bins_end_coord & bins_end_coord, string chr, unsigned int position) {
	int bin_internal_index=-1;
	if ( bins_end_coord.find(chr) == bins_end_coord.end() ) {
		// chr name is not found in binning system
		bin_internal_index = -1;
	} else {
		vector<unsigned int> & coords_bins_of_chr = bins_end_coord[chr]; // a vector of end coordinates (1-base) of this chr.
		vector<unsigned int>::iterator bin_it = lower_bound( coords_bins_of_chr.begin(), coords_bins_of_chr.end(), position);
		if (position==*bin_it) bin_it++;
		bin_internal_index = bin_it - coords_bins_of_chr.begin(); // bin_internal_index is 0-base
		if (bin_internal_index==coords_bins_of_chr.size()) {
			//cerr << "position(1-base): " << position << " doesn't exist in " << chr << endl;
			bin_internal_index=-1;
		}
	}
	return bin_internal_index;
}

// Each bin is in the range [ bins_end_coord[i-1], bins_end_coord[i] )
// Given a query region, we want to know which bin has overlap with this query region. If there is overlap,
// return (1) bin index, and (2) the overlap length
// ongoing devevloping
int find_overlap_bin(Bins_end_coord & bins_end_coord, string query_region_chr, unsigned int query_region_start_coord,
	unsigned int query_region_end_coord, int & overlap_length)
{
	int bin_internal_index = -1;
	overlap_length = -1;
	if ( bins_end_coord.find(query_region_chr) != bins_end_coord.end() ) {
		// chr name is found in binning system
		bin_internal_index = find_bin_of_position(bins_end_coord, query_region_chr, query_region_start_coord);
		if (bin_internal_index != -1) {
			unsigned int bin_end_coord = bins_end_coord[query_region_chr][bin_internal_index];
			if ( query_region_end_coord > bin_end_coord )
				overlap_length = bin_end_coord - query_region_start_coord + 1;
			else
				overlap_length = query_region_end_coord - query_region_start_coord + 1;
		}
	}
	return bin_internal_index;
}

void print_uint_vec( ostream& os, vector<unsigned int>& v, int len )
{
	int i;
	if (len>v.size() || len==0) len=v.size();
	if (len==0 || v.size()==0) {
		os << "[" << "]";
	} else {
		os << "[";
		for (i=0; i<len-1; i++)
			os << v[i] << ",";
		os << v[i] << "]";
	}
}

void print_bins( ostream& os, Bins_end_coord & bins_end_coord, Bins_index & bins_index, Bins_info & bins_info) {
	Bins_end_coord::iterator it;
	for (it=bins_end_coord.begin(); it!=bins_end_coord.end(); ++it) {
		string chr = it->first;
		vector<unsigned int> coords=it->second;
		vector<int> indexes=bins_index[chr];
		vector<string> infos=bins_info[chr];
		int n_bins = coords.size();
		for (int i=0; i<n_bins; i++) {
			os << indexes[i] << "\t" << chr << "\t" << coords[i] << "\t'" << infos[i] << "'"  << endl;
		}
	}
}

/*
void print_bins_fullinfo( Bins_FullInfo & bins_fullinfo ) {
	Bins_FullInfo::iterator: it;
	for (it=bins.begin(); it!=bins.end(); ++it) {
		cout << it.first << "\t";
		print_int_vec(cout, it.second, it.second.size());
		cout << endl;
	}
}
*/

// Bins2Values: a map of bin_index -> a vector of values. bin_index is always 1-base
void create_Bins2Values(int num_bins, int num_of_values, double init_value, Bins2Values & bins2values)
{
	for (int bin_index=1; bin_index<=num_bins; bin_index++) {
		vector<double> values;
		for (int i=0; i<num_of_values; i++)
			values.push_back( init_value );
		bins2values[bin_index] = values;
	}
}

// Bins2Values: a map of marker_index -> a vector of values.
void create_Bins2Values(vector<int> markers_index, int num_of_values, double init_value, Bins2Values & bins2values)
{
	int num_bins = markers_index.size();
	for (int ibin=0; ibin<num_bins; ibin++) {
		vector<double> values;
		for (int i=0; i<num_of_values; i++)
			values.push_back( init_value );
		bins2values[markers_index[ibin]] = values;
	}
}

void print_Bins2Values(Bins2Values & bins2values)
{
	cout.precision(15);
	Bins2Values::iterator it;
	for (it=bins2values.begin(); it!=bins2values.end(); ++it) {
		vector<double> & values = it->second;
		cout << it->first;
		for (int i=0; i<values.size(); i++)
			cout << "\t" << values[i];
		cout << endl;
	}
}

void print_Bins2UnsignedIntegers(Bins2UnsignedIntegers& bins2values)
{
	Bins2UnsignedIntegers::iterator it;
	for (it=bins2values.begin(); it!=bins2values.end(); ++it) {
		vector<unsigned int> & values = it->second;
		cout << it->first;
		for (int i=0; i<(int)values.size(); i++)
			cout << "\t" << values[i];
		cout << endl;
	}
}

// when optional_write==TRUE, we assume there are two values associated with each bin
void write_Bins2Values(Bins2Values & bins2values, vector<string> & columns_names,
	string output_file, bool optional_write)
{
	ofstream out;
	out.open(output_file.c_str());
	if (out.fail()){
		cerr << "Error: Unable to write " << output_file << " in write_Bins2Value()" << endl;
		exit(EXIT_FAILURE);
	}
	int i;
	for (i=0; i<columns_names.size()-1; i++)
		out << columns_names[i] << "\t";
	out << columns_names[i] << endl;

	out.precision(15);
	Bins2Values::iterator it;
	for (it=bins2values.begin(); it!=bins2values.end(); ++it) {
		vector<double> & values = it->second;
		out << it->first;
		if (optional_write) {
			// assume there are at least two values associated with each bin
			// for example, when we have three associated values, they can be 
			// (1) methylation_count
			// (2) unmethylation_count
			// (3) number of reads
			double n = values[0] + values[1];
			double v;
			if (n==0) v=0;
			else v=values[0]/n;
			out << "\t" << v << "\t" << n;
		}
		for (i=0; i<values.size(); i++)
			out << "\t" << values[i];
		out << endl;
	}
	out.close();
}

// Bins2Value: a map of bin_index -> a value. bin_index is always 1-base
void create_Bins2Value(int num_bins, double init_value, Bins2Value & bins2value)
{
	for (int bin_index=1; bin_index<=num_bins; bin_index++)
		bins2value[bin_index] = init_value;
}

ostream& operator<<(ostream& out, Bins2Value& bins2value) {
	out << "bin_index" << "\t" << "value" << endl;
	out.precision(15);
	Bins2Value::iterator it;
	for (it=bins2value.begin(); it!=bins2value.end(); ++it)
		out << it->first << "\t" << it->second << endl;
	return(out);
}

void write_Bins2Value(Bins2Value & bins2value, string output_file)
{
	ofstream out;
	out.open(output_file.c_str());
	if (out.fail()){
		cerr << "Error: Unable to write " << output_file << " in write_Bins2Value()" << endl;
		exit(EXIT_FAILURE);
	}
	out << "bin_index" << "\t" << "value" << endl;
	out.precision(15);
	Bins2Value::iterator it;
	for (it=bins2value.begin(); it!=bins2value.end(); ++it)
		out << it->first << "\t" << it->second << endl;
	out.close();
}

