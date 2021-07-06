#ifndef MATRIX_H
#define MATRIX_H

#include <sstream>
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include "data_types.h"

using namespace std;

class Matrix_Double {
private:
	vector<vector<double> > mat; // a double matrix of size nrow X ncol
	vector<int> row_labels; // a int vector of size nrow X 1. The first column (integer) is assumed to be the label of a row. 
	unsigned int nrow, ncol;
	bool empty;
	void process_one_line_of_mat_file(string & line, int num_header_column);
	// putting constructors as private can make a class non-copyable
public:
	Matrix_Double() {nrow=0; ncol=0; empty=true;};
	Matrix_Double(const unsigned int ncol_) {nrow=0; ncol=ncol_; empty=true;};
	Matrix_Double(const unsigned int nrow_, const unsigned int ncol_, const double init_value);
	Matrix_Double(const string & mat_file, bool header_line, int num_header_column);
	void multipy_single_value(const double v) {
		if (!empty) {
			for (int i=0; i<nrow; i++) {
				for (int j=0; j<ncol; j++) {
					mat[i][j] *= v;
				}
			}
		}
	};
	long append_row_vector(vector<double>& v_, int row_label_);
	long append_row_vector_with_filter(vector<double>& v_, int row_label_, double min_threshold_maxv_minv=1.5);
	vector<double>* get_row_vector(int row_index) { // row_index is 0-based
		if (row_index>=nrow) return(NULL);
		else return(&mat[row_index]); };
	unsigned int get_row_num() {return nrow;};
	unsigned int get_column_num() {return ncol;};
	bool get_column_sum(const int j, double & v);
	bool get_row_sum(const int i, double & v);
	bool get_element(const int i, const int j, double & v);
	bool set_element(const int i, const int j, double v);
	bool isempty() {return empty;};
	bool exist_row_label(const int row_label) {
		if (find(row_labels.begin(), row_labels.end(), row_label) != row_labels.end()) {
			return(true);
		} else {
			return(false);
		}
	};
	int get_row_index(const int row_label) {
		if (empty) {
			return(-1);
		}
		int row_index = -1;
		vector<int>::iterator it = find(row_labels.begin(), row_labels.end(), row_label);
		if (it != row_labels.end()) {
			row_index = it - row_labels.begin();
		}
		return(row_index);
	};
	bool get_row_label(const int row_index, int & row_label) { // row_index is 0-base
		if (!empty) { row_label = row_labels[row_index]; return true; }
		else return false;};
	bool set_row_labels(vector<int>& row_labels_) {
		if (!empty) {
			if (nrow==(unsigned int)row_labels_.size()) {
				row_labels.clear();
				for (int i=0; i<nrow; i++) row_labels.push_back(row_labels_[i]);
				return true;
			} else return false;
		} else return false;
	};
	vector<int> & get_row_labels() {return row_labels;};
	void get_unique_row_labels(vector<int> & uniq_labels);
	void print_with_additional_column_of_Bins2Value(ostream & os, Bins2Value& additional_column_data);
};

ostream & operator<<(ostream & os, Matrix_Double & mat);
void em_supervise(Matrix_Double & p, int max_iter, vector<double> & theta, Matrix_Double& q);
void em_semisupervise(Matrix_Double & p, vector<int> & Rm, vector<int> & Rl,
	int max_iter, vector<double> & theta, Matrix_Double& q, vector<double>& q_unknown,
	vector<double> & m);
void readCounts_by_reads_posterior_probability_version_regular(Matrix_Double & q, double unit, Matrix_Double& readCounts);
void readCounts_by_reads_posterior_probability_version_unknownclass(Matrix_Double& q, vector<double>& q_unknown, double unit, Matrix_Double& readCounts);
#endif
