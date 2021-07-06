#ifndef DATA_TYPES_H
#define DATA_TYPES_H

using namespace std;

//
// definitions of data types
//
typedef map<string, vector<unsigned int> > Bins_end_coord; // a map: chr -> a vector of end coordinate of each bin. coordinate is 1-base
typedef map<string, vector<int> > Bins_index; // a map: chr -> a vector of bin index. Note some bins are void/complementary, so their indexes are 0
typedef map<string, vector<string> > Bins_info; // a map: chr -> a vector of bins' info. Each bin's info is a string

typedef map<string, map<int,map<int,double> > > Wig_bins_data; // chr -> (bin_internal_index -> (position -> value)). This is a three-layer map

typedef unsigned int GENOME_POSITION;

// Bins2Values: a map for bin_index -> a vector of value. bin_index is always 1-base
typedef map<int, vector<double> > Bins2Values;
typedef map<int, vector<unsigned int> > Bins2UnsignedIntegers;

// Bins2PairedValues: a map for bin_index -> two vectors of values. bin_index is always 1-base
typedef map<int, pair<vector<double>,vector<double> > > Bins2PairedValues;

// Bins2PositionValuePairs: a map for bin_index -> a pair of (position vector, value vector). bin_index is always 1-base
typedef map<unsigned int, pair<vector<GENOME_POSITION>,vector<double> > > Bins2PositionValuePairs;

// Bins2Value: a map of bin_index -> a value. bin_index is always 1-base
typedef map<int, double> Bins2Value;

//
// functions for these data types
//
void read_bins_annot_file(string input_bins_annot_file, Bins_end_coord & bins_end_coord,
	Bins_index & bins_index, Bins_info & bins_info, bool has_header_line);
int get_num_of_non_void_bins(Bins_index & bins_index, vector<int> & returned_markers_index);
int find_exact_bin(Bins_end_coord & bins_end_coord, string chr, unsigned int bin_start_coord, unsigned int bin_end_coord);
int find_bin_of_position(Bins_end_coord & bins_end_coord, string chr, unsigned int position);
int find_overlap_bin(Bins_end_coord & bins_end_coord, string query_region_chr, unsigned int query_region_start_coord,
	unsigned int query_region_end_coord, int & overlap_length);
void print_bins( ostream& os, Bins_end_coord & bins_end_coord, Bins_index & bins_index, Bins_info & bins_info);

void create_Bins2Values(int num_bins, int num_of_values, double init_value, Bins2Values & bins2values);
void create_Bins2Values(vector<int> markers_index, int num_of_values, double init_value, Bins2Values & bins2values);
void write_Bins2Values(Bins2Values & bins2values, vector<string> & columns_names, string output_file, bool optional_write=false);
void print_Bins2Values(Bins2Values & bins2values);
void print_Bins2UnsignedIntegers(Bins2UnsignedIntegers& bins2values);

void create_Bins2Value(int num_bins, double init_value, Bins2Value & bins2value);
ostream& operator<<(ostream& out, Bins2Value& bins2value);
void write_Bins2Value(Bins2Value & bins2value, string output_file);

#endif

