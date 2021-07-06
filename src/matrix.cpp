#include <cstdlib>
#include <cmath>
#include <sstream>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <cassert>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/assign.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <bits/stdc++.h>  // for functions: min_element() and max_element()
#include "matrix.h"
#include "utils.h"

using namespace std;
using namespace boost;

#include "template_utils.cpp"

Matrix_Double::Matrix_Double(const unsigned int nrow_, const unsigned int ncol_, const double init_value)
{
	nrow = nrow_;
	ncol = ncol_;
	empty = false;
	mat.reserve(nrow_);
	for (int i=0; i<nrow; i++) {
		vector<double> temp;
		temp.reserve(ncol);
		for (int j=0; j<ncol; j++)
			temp.push_back(init_value);
		mat.push_back( temp );
	}
	empty = mat.empty();
}

// a line represents a row of the matrix. Each element in this row is delimited by a TAB.
void Matrix_Double::process_one_line_of_mat_file(string & line, int num_header_column)
{
	vector<string> items;
	split(items, line, is_any_of("\t"));
	vector<double> v;
	for (int i=0; i<items.size(); i++) {
		if (i==0) row_labels.push_back( atoi(items[0].c_str()) ); // the first column is the label of this row.
		if (i<num_header_column) continue; // skip the first num_header_column header columns
		v.push_back( atof(items[i].c_str()) );
	}
	mat.push_back( v );
}

Matrix_Double::Matrix_Double(const string & mat_file,
	bool header_line=false, int num_header_column=1)
{
	unsigned long i_line=0;
	// input is a plain text file
	istream * in=&cin; // default is stdin
	ifstream fin;
	if ( !(mat_file.compare("stdin")==0)) {
		fin.open(mat_file.c_str());
		if (fin.fail()){
			cerr << "Error: Unable to open " << mat_file << " in Matrix_Double()" << endl;
			exit(EXIT_FAILURE);
		}
		in = &fin;
	}
	string line;
	while (!(*in).eof()) {
		getline((*in), line);
		//cerr << "Line: " << line << endl;
		if (header_line) {
			if (i_line==0) {
				// skip the header line
				i_line++;
				continue;
			}
		}
		if (line.empty()) {
			// this is the last line of the file
			break;
		}
		process_one_line_of_mat_file(line, num_header_column);
		i_line++;
	}
	if ( !(mat_file.compare("stdin")==0)) {
		fin.close();
	}
	empty = mat.empty();
	nrow = mat.size();
	if (nrow>=1) ncol = mat[0].size();
	else ncol = 0;
}

// return the row index of this vector
long Matrix_Double::append_row_vector(vector<double>& v_, int row_label_) {
	if (v_.size()!=ncol) {
		cerr << "Error of Matrix_Double::append_row_vector: the appended vector size (" << v_.size() << ") does not match the number of column in the matrix (" << ncol << ")!\nExit." << endl;
		exit(EXIT_FAILURE);
	}
	vector<double> v;
	for (int i=0; i<v_.size(); i++) {
		v.push_back( v_[i] );
	}
	mat.push_back( v );
	nrow = mat.size();
	row_labels.push_back(row_label_);
	empty=false;
	return(nrow-1);
}

//
// if max(v_)/min(v_) >= min_threshold_maxv_minv, then we append this vector, otherwise ignore it.
//
// return row index of this appended vector, if successfully appending; return -1 otherwise.
//
long Matrix_Double::append_row_vector_with_filter(vector<double>& v_, int row_label_, double min_threshold_maxv_minv) {
	if (v_.size()!=ncol) {
		cerr << "Error of Matrix_Double::append_row_vector_with_filter: the appended vector size (" << v_.size() << ") does not match the number of column in the matrix (" << ncol << ")!\nExit." << endl;
		exit(EXIT_FAILURE);
	}
	double min_ = *min_element(v_.begin(), v_.end());
	double max_ = *max_element(v_.begin(), v_.end());
	//cout << "max: " << max_ << ", min: " << min_ << ", ratio: " << max_/min_ << endl;
	if (min_==0) {
		if (max_==0) return(-1);
	} else {
		double ratio=max_/min_;
		if (ratio<min_threshold_maxv_minv) return(-1);
	}
	long row_index = append_row_vector(v_, row_label_);
	return(row_index);
}

bool Matrix_Double::get_element(const int i, const int j, double & v)
{
	if (!empty) {
		if (i<nrow && j<ncol) {
			v = mat[i][j];
			return true;
		} else {
			cerr << "Warning: i=" << i << " and j=" << j << " exceed matrix size!" << endl;
			return false;
		}
	} else {
		v = 0;
		return false;
	}
}

bool Matrix_Double::set_element(const int i, const int j, double v)
{
	if (!empty) {
		if (i<nrow && j<ncol) {
			mat[i][j] = v;
			return true;
		} else {
			cerr << "Warning: i=" << i << " and j=" << j << " exceed matrix size!" << endl;
			return false;
		}
	} else {
		return false;
	}
}

bool Matrix_Double::get_column_sum(const int j, double & v)
{
	if (!empty) {
		if (j<ncol) {
			v = 0;
			for (int k=0; k<nrow; k++) {
				v += mat[k][j];
			}
			return true;
		} else {
			cerr << "Warning: j=" << j << " exceeds matrix column size (" << ncol << ")!" << endl;
			return false;
		}
	} else {
		v = 0;
		return false;
	}
}

bool Matrix_Double::get_row_sum(const int i, double & v)
{
	if (!empty) {
		if (i<nrow) {
			v = 0;
			for (int k=0; k<ncol; k++) {
				v += mat[i][k];
			}
			return true;
		} else {
			cerr << "Warning: i=" << i << " exceeds matrix row size (" << nrow << ")!" << endl;
			return false;
		}
	} else {
		v = 0;
		return false;
	}
}

void Matrix_Double::get_unique_row_labels(vector<int> & uniq_labels)
{
	map<int, int> counts;
	for (int i=0; i<row_labels.size(); i++)
		if (counts.find(row_labels[i]) != counts.end())
			counts[row_labels[i]] += 1; // found this label, then increment count
		else {
			counts.insert(make_pair(row_labels[i],1)); // not found this label, then count it once
			//counts[row_labels[i]] = 1;
		}
	map<int,int>::iterator it;
	for (it=counts.begin(); it!=counts.end(); ++it)
		uniq_labels.push_back(it->first);
}

void Matrix_Double::print_with_additional_column_of_Bins2Value(ostream & os, Bins2Value& additional_column_data)
{
	// obtain unique row labels
	vector<int> unique_row_labels;
	int row_label_;
	for (int i=0; i<nrow; i++) {
		row_label_ = row_labels[i];
		unique_row_labels.push_back( row_label_ );
	}
	Bins2Value::iterator iter;
	for (iter=additional_column_data.begin(); iter!=additional_column_data.end(); iter++) {
		row_label_ = iter->first;
		if (!exist_row_label(row_label_)) {
			unique_row_labels.push_back( row_label_ );
		}
	}
	sort(unique_row_labels.begin(), unique_row_labels.end());
	// print both 'mat' and 'additional_column_data' as the last column
	for (int i=0; i<unique_row_labels.size(); i++) {
		row_label_ = unique_row_labels[i];
		os << row_label_;
		int row_index = get_row_index(row_label_);
		if (row_index!=-1) {
			for (int j=0; j<ncol; j++) {
				os << "\t" << mat[row_index][j];
			}
		} else {
			for (int j=0; j<ncol; j++) {
				os << "\t0";
			}
		}
		if (additional_column_data.find(row_label_)!=additional_column_data.end()) {
			os << "\t" << additional_column_data[row_label_];
		} else {
			os << "\t0";
		}
		os << endl;
	}
}

ostream & operator<<(ostream & os, Matrix_Double & mat)
{
	double v;
	int row_label;
	if (!mat.isempty()) {
		for (int i=0; i<mat.get_row_num(); i++) {
			int j;
			mat.get_row_label(i, row_label);
			os << row_label << "\t";
			for (j=0; j<mat.get_column_num()-1; j++) {
				mat.get_element(i,j,v);
				os << v << "\t";
			}
			mat.get_element(i,j,v);
			os << v << endl;
		}
	}
	return os;
}

double objective_em_supervise(Matrix_Double & p, vector<double> & theta)
{
	unsigned int ncol = p.get_column_num();
	unsigned int nrow = p.get_row_num();
	double v;
	double obj = 0;
	for (int i=0; i<nrow; i++) {
		double sum = 0;
		for (int j=0; j<ncol; j++) {
			p.get_element(i,j,v);
			sum += theta[j]*v;
		}
		obj += log(sum);
	}
	return obj;
}
//
// There are T known tissues and 1 unknown tissue (described by a double vector "m")
//
// input:
//   p is a matrix of N X T, where N is number of reads and T is number of known tissue
//
// output:
//   theta (model parameters), a vector with "T" elements.
//   q is a matrix of N X T, the tissue-specific posterior probabilty of each read
//
void em_supervise(Matrix_Double & p, int max_iter, vector<double> & theta, Matrix_Double& q)
{
	cout.precision(15);
	cerr.precision(15);
	unsigned int ncol = p.get_column_num();
	unsigned int nrow = p.get_row_num(); // nrow is number of tissues.
	theta.resize(ncol, 0); // assign (num_tissues) space and initialize to 0s.
	// initialize model parameters as uniform distribution
	// alternatively, model parameters can be random numbers by satisfying the crition
	//    (1) \sum_{i=1}^{ncol}{theta_i}=1
	for (int j=0; j<ncol; j++) {
		theta[j] = 1/(double)ncol;
	}
	// create and initialize q with the same size of p and with all elements initialized as 0
	//Matrix_Double q(nrow, ncol, 0);
	//cerr << "iter 0\t" << objective_em_supervise(p, theta) << endl;
	double v1, v2;
	for (int iter=0; iter<max_iter; iter++) {
		cerr << iter+1 << "," ;
		// E-step: estimate q
		for (int i=0; i<nrow; i++) {
			double sum = 0;
			for (int j=0; j<ncol; j++) {
				p.get_element(i,j,v1);
				v2 = theta[j]*v1;
				q.set_element( i, j, v2 );
				sum += v2;
			}
			for (int j=0; j<ncol; j++) {
				q.get_element(i,j,v2);
				v2 /= sum;
				q.set_element( i, j, v2 );
			}
		}
		// M-step: estimate theta
		for (int j=0; j<ncol; j++) {
			double sum=0;
			for (int i=0; i<nrow; i++) {
				q.get_element(i,j,v2);
				sum += v2;
			}
			theta[j] = sum / nrow;
		}
		// for debug
		//cerr << "iter " << (iter+1) << "\t" << objective_em_supervise(p, theta) << endl;
	}
	// last E-step: estimate q
	for (int i=0; i<nrow; i++) {
		double sum = 0;
		for (int j=0; j<ncol; j++) {
			p.get_element(i,j,v1);
			v2 = theta[j]*v1;
			q.set_element( i, j, v2 );
			sum += v2;
		}
		for (int j=0; j<ncol; j++) {
			q.get_element(i,j,v2);
			v2 /= sum;
			q.set_element( i, j, v2 );
		}
	}
	cerr << endl;
}

//
// There are T known tissues and 1 unknown tissue (described by a double vector "m")
//
// input:
//   p is a matrix of N X T, where N is number of reads and T is number of known tissue
//     p.row_label is an int vector of N X 1, each element is the marker Id (1-base) of a read.
//   Rm is a vector of N X 1, each element is number of valid methylated CpG sites in a read
//   Rl is a vector of N X 1, each element is number of all valid CpG sites in a read
//
// output:
//   theta is a vector of (T+1) X 1, where T is number of known tissue. This vector is already allocated with space of (ncol+1) units.
//   m is a vector of M X 1, where M is number of markers.
//   q is a matrix of N X T, the tissue-specific posterior probabilty of each read
//   q_unknown is a vector of N X 1, the posterior probabilty of each read for the unknown class. It does not need to be allocated space before this function calling.
//
void em_semisupervise(Matrix_Double & p, vector<int> & Rm, vector<int> & Rl,
	int max_iter, vector<double> & theta, Matrix_Double& q, vector<double>& q_unknown,
	vector<double> & m)
{
	cout.precision(15);
	cerr.precision(15);
	unsigned int ncol = p.get_column_num(); // T, number of known tissues
	unsigned int nrow = p.get_row_num(); // N, number of reads

	theta.resize(ncol+1, 0); // assign (num_tissues+1) space and initialize to 0s.

	vector<int> uniq_marker_Ids;
	p.get_unique_row_labels(uniq_marker_Ids);
	int nMarker = uniq_marker_Ids.size(); // M, number of markers that all N reads cover. Some markers may not be covered by these N reads.
	m.resize(nMarker, 0); // assign nMarker space and initialize to 0s.
	map<int,int> markerId2Index; // marker index is 0-based.
	for (int index=0; index<nMarker; index++) {
		markerId2Index.insert(make_pair(uniq_marker_Ids[index],index));
	}

	// initialize model parameters theta as uniform distribution, and m as 0.5
	// alternatively, theta and m can be random numbers, by satisfying:
	//    (1) \sum_{i=1}^{ncol+1}{theta_i}=1
	//    (2) 0 <= m_k <=1 for all k=1,2,...,#marker_covered_by_all_reads
	for (int j=0; j<ncol+1; j++) {
		theta[j] = 1/(double)(ncol+1);
	}
	for (int k=0; k<nMarker; k++) {
		m[k] = 0.5;
	}
	for (int i=0; i<nrow; i++)
		q_unknown.push_back(0);

	// EM algorithm
	for (int iter=0; iter<max_iter; iter++) {
		cerr << iter+1 << "," ;
		// E-step: estimate q (for T known tissues) and q_unknown (for one unknown tissue)
		double v1, v2, likelihood_unknown_class;
		int markerId, markerIndex;
		for (int i=0; i<nrow; i++) {
			// process the first T known tissues
			double sum = 0;
			for (int j=0; j<ncol; j++) {
				p.get_element(i,j,v1);
				v2 = theta[j]*v1;
				q.set_element( i, j, v2 );
				sum += v2;
			}
			// process the last unknown tissue
			p.get_row_label(i, markerId); // get markerId (1-base) of read i
			markerIndex = markerId2Index[markerId];
			likelihood_unknown_class = pow(m[markerIndex],Rm[i]) * pow(1-m[markerIndex],Rl[i]-Rm[i]);
			q_unknown[i] = theta[ncol]*likelihood_unknown_class;
			sum += q_unknown[i];

			// update q and q_unknown
			for (int j=0; j<ncol; j++) {
				q.get_element(i,j,v2);
				v2 /= sum;
				q.set_element( i, j, v2 );
			}
			q_unknown[i] /= sum;
		}

		// M-step 1: estimate unknown class's methylation level (m) of each marker
		double sum1=0, sum2=0;
		int curr_marker_index=-1, prev_marker_index=-1; // all marker index is 0-base
		for (int i=0; i<nrow; i++) {
			// for each row (or read)
			p.get_row_label(i, markerId); // get the marker_index (0-base) of read i
			curr_marker_index = markerId2Index[markerId];
			if (curr_marker_index != prev_marker_index) {
				// this is the 1st read of a new marker
				// we need to summarize m value of previous marker
				if (prev_marker_index!=-1) { // current marker is not the 1st marker
					m[prev_marker_index] = (sum2!=0 ? sum1/sum2 : 0);
					sum1 = 0;
					sum2 = 0;
				}
				prev_marker_index = curr_marker_index;
			}
			sum1 += q_unknown[i]*Rm[i];
			sum2 += q_unknown[i]*Rl[i];
		}
		m[curr_marker_index] = (sum2!=0 ? sum1/sum2 : 0); // estimate m of the last marker

		//cout << "iter=" << iter << ", " << "m=" << endl; // for debug
		//print_vec(cout, m, "\n"); // for debug
		//cout << endl; // for debug

		// M-step 2: estimate theta, which has ncol+1 values
		double sum;
		for (int j=0; j<ncol; j++) {
			sum=0;
			for (int i=0; i<nrow; i++) {
				q.get_element(i,j,v2);
				sum += v2;
			}
			theta[j] = sum / nrow;
		}
		sum=0;
		for (int i=0; i<nrow; i++)
			sum += q_unknown[i];
		theta[ncol] = sum / nrow;

		//cout << "round " << iter+1 << ", theta="; // for debug
		//print_vec(cout, theta, ", "); // for debug
		//cout << endl << endl; // for debug
	}
	cerr << endl;
	//cout << "final:" << endl << std::flush;
	//print_vec(cout, theta, ", "); // for debug
	//cout << endl << endl << std::flush; // for debug
}

void readCounts_by_reads_posterior_probability_version_regular(Matrix_Double& q, double unit, Matrix_Double& readCounts)
{
	unsigned int ncol = q.get_column_num(); // T, number of known tissues
	unsigned int nrow = q.get_row_num(); // N, number of reads

	vector<int> uniq_marker_Ids;
	q.get_unique_row_labels(uniq_marker_Ids);
	int nMarker = uniq_marker_Ids.size(); // M, number of markers that all N reads cover. Some markers may not be covered by these N reads.

	double v;
	int curr_markerId=-1, prev_markerId=-1; // all marker index is 0-base
	vector<double> readCountsPerMarker(ncol, 0); // a vector of read counts for one marker (T elements) with default values 0
	for (int i=0; i<nrow; i++) {
		// for each row (or read)
		q.get_row_label(i, curr_markerId); // get the marker ID of read i
		if (curr_markerId != prev_markerId) {
			// this is the 1st read of a new marker
			// we need to summarize read counts of each tissue for the previous marker and deposit them to the matrix "readCounts"
			// then clear readCountsPerMarker for a new read counting of the next new marker
			if (prev_markerId!=-1) { // current marker is not the 1st marker
				multi_vec_by_number(readCountsPerMarker, unit); // normalize the raw read counts by unit.
				readCounts.append_row_vector(readCountsPerMarker, prev_markerId);
				assign_vector_zeros(readCountsPerMarker);
			}
			prev_markerId = curr_markerId;
		}
		for (int j=0; j<ncol; j++) {
			q.get_element(i,j,v);
			readCountsPerMarker[j] += v;
		}
	}
	multi_vec_by_number(readCountsPerMarker, unit); // normalize the raw read counts by unit.
	readCounts.append_row_vector(readCountsPerMarker, curr_markerId);
}

void readCounts_by_reads_posterior_probability_version_unknownclass(Matrix_Double& q, vector<double>& q_unknown, double unit, Matrix_Double& readCounts)
{
	unsigned int ncol = q.get_column_num(); // T, number of known tissues
	unsigned int nrow = q.get_row_num(); // N, number of reads
	if (nrow!=(unsigned int)q_unknown.size()) {
		cerr << "Error (readCounts_by_reads_posterior_probability_version_unknownclass): row number of q_unknown does not match with row number of q!" << endl;
		exit(EXIT_FAILURE);
	}

	vector<int> uniq_marker_Ids;
	q.get_unique_row_labels(uniq_marker_Ids);
	int nMarker = uniq_marker_Ids.size(); // M, number of markers that all N reads cover. Some markers may not be covered by these N reads.

	double v;
	int curr_markerId=-1, prev_markerId=-1; // all marker index is 0-base
	vector<double> readCountsPerMarker(ncol+1, 0); // a vector of read counts for one marker (T+1 elements) with default values 0. The last element is for the unknown class.
	for (int i=0; i<nrow; i++) {
		// for each row (or read)
		q.get_row_label(i, curr_markerId); // get the marker ID of read i
		if (curr_markerId != prev_markerId) {
			// this is the 1st read of a new marker
			// we need to summarize read counts of each tissue for the previous marker and deposit them to the matrix "readCounts"
			// then clear readCountsPerMarker for a new read counting of the next new marker
			if (prev_markerId!=-1) { // current marker is not the 1st marker
				multi_vec_by_number(readCountsPerMarker, unit); // normalize the raw read counts by unit.
				readCounts.append_row_vector(readCountsPerMarker, prev_markerId);
				assign_vector_zeros(readCountsPerMarker);
			}
			prev_markerId = curr_markerId;
		}
		for (int j=0; j<ncol; j++) {
			q.get_element(i,j,v);
			readCountsPerMarker[j] += v;
		}
		readCountsPerMarker[ncol] += q_unknown[i];
	}
	multi_vec_by_number(readCountsPerMarker, unit); // normalize the raw read counts by unit.
	readCounts.append_row_vector(readCountsPerMarker, curr_markerId);
}

