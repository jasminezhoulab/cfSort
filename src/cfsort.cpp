// Academic Software License: © 2021 UCLA (“Institution”). Academic or nonprofit researchers are
// permitted to use this Software (as defined below) subject to Paragraphs 1-3:
//
// Institution hereby grants to you free of charge, so long as you are an academic or nonprofit
// researcher, a nonexclusive license under Institution’s copyright ownership interest in this
// software and any derivative works made by you thereof (collectively, the “Software”) to use,
// copy, and make derivative works of the Software solely for educational or academic research
// purposes, and to distribute such Software free of charge to other academic or nonprofit
// researchers for their educational or academic research purposes, in all cases subject to the
// terms of this Academic Software License. Except as granted herein, all rights are reserved by
// Institution, including the right to pursue patent protection of the Software.
//
// Any distribution of copies of this Software -- including any derivative works made by you
// thereof -- must include a copy (including the copyright notice above), and be made subject to
// the terms, of this Academic Software License; failure by you to adhere to the requirements in
// Paragraphs 1 and 2 will result in immediate termination of the license granted to you pursuant
// to this Academic Software License effective as of the date you first used the Software.
//
// IN NO EVENT WILL INSTITUTION BE LIABLE TO ANY ENTITY OR PERSON FOR DIRECT, INDIRECT, SPECIAL,
// INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF INSTITUTION HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. INSTITUTION
// SPECIFICALLY DISCLAIMS ANY AND ALL WARRANTIES, EXPRESS AND IMPLIED, INCLUDING, BUT NOT LIMITED
// TO, ANY IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE
// IS PROVIDED “AS IS.” INSTITUTION HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
// ENHANCEMENTS, OR MODIFICATIONS OF THIS SOFTWARE.
//
// Commercial entities: please contact Prof. Xiangong Jasmine Zhou (xjzhou@mednet.ucla.edu) for
// licensing opportunities.
//

#include <cstdlib>
#include <vector>
#include <map>
#include <string>
#include <climits>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <boost/algorithm/string.hpp>
#include <boost/assign.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include "my_getopt-1.5/getopt.h"
#include "utils.h"
#include "data_types.h"
#include "matrix.h"

using namespace std;
using namespace boost;

#include "template_utils.cpp"

//#define DEBUG

// Global variables
string tissue_markers_file="";
int num_tissues=0;
string reads_binning_file="stdin"; // default is "stdin"
string output_file="stdout"; // default is "stdout"
string output_type="tissueFraction"; // default is "tissueFraction"
double min_likelihood_ratio_cutoff=-1.0; // default is -1.0, which means "No reads filtering with likelihood ratio cutoff"
string em_algorithm_type="em.global.unknown"; // default is "the EM algorithm with unknown class", whose short name: "em.global.unknown". Other option is "em_known"
int em_max_iterations=100; // default is 100 iterations

/*-----------------------------------
command line parser's data structure
-----------------------------------*/
struct option longopts[] =
{
	/* name,                      has_arg,          flag, val */ /* longind */
	{ "help",                           no_argument,       0,   'h' }, /*       0 */
	{ "reads_binning_file",             required_argument, 0,   'r' }, /*       2 */
	{ "tissue_markers_file",            required_argument, 0,   't' }, /*       1 */
	{ "num_tissues",                    required_argument, 0,   'T' }, /*       1 */
	{ "min_likelihood_ratio_cutoff",    required_argument, 0,   'p' }, /*       1 */
	{ "em_algorithm",                   required_argument, 0,   'a' }, /*       3 */
	{ "em_max_iterations",              required_argument, 0,   'x' }, /*       3 */
	{ "output_file",                    required_argument, 0,   'o' }, /*       3 */
	{ "output_type",                    required_argument, 0,   'O' }, /*       3 */
	/* end-of-list marker */
	{ 0, 0, 0, 0 }
};

char *shortopts = "hr:t:T:a:p:x:o:O:"; /* short options string */
int longind = 0; /* long option list index */

void print_usage() {
	cerr << "Calculate tissue-specific likelihood for each sequencing read and perform EM algorithm on all reads to obtain model parameters, then compute and output the tissue-specific read counts for each marker" << endl << endl;
	cerr << "USAGE: cfsort [options]" << endl << endl;
	cerr << "   Options:" << endl;
	cerr << "     -a [string]:  read-based tissue deconvolution EM algorithm type: em.global.unknown (default), em.global.known, em.local.unknown, em.local.known." << endl;
	cerr << "     -r [FILE]:  reads binning file (default: stdin)" << endl;
	cerr << "     -o [FILE]:  output file (default: stdout)" << endl;
	cerr << "     -O [FILE]:  output type tissueFraction (default), tissueFraction+readCountRaw, tissueFraction+readCountPerMillion, tissueFraction+readCountPerBillion" << endl;
    cerr << "     -T [positive integer]: number of tissue types." << endl;
    cerr << "     -p [positive real value]: all reads with likelihood ratio < cutoff (default: -1.0, meanin no reads filtering is used) will not be used. Likelihood ratio is the max(all tissues' likelihoods)/min(all tissues' likelihoods)" << endl;
	cerr << "     -t [FILE]:  tissue markers file: each marker for each tissue type has either a median beta value or a paired shape parameters of beta-distribution from the tissue population of tissue types." << endl;
	cerr << "     -x [positive integer]:  EM algorithm maximum iteration number: 100 (default)." << endl;
	cerr << endl;
    cerr << "version 1.0, June 15, 2021" << endl;
    cerr << "By Wenyuan Li" << endl;
	cerr << endl;
}

void parse_command_line(int argc, char * argv[])
{
	int opt; /* during argument parsing, opt contains the return value from getopt() */
	int longind = 0; /* long option list index */
	/* parse all options from the command line */
	while ((opt=getopt_long(argc, argv, shortopts, longopts, &longind)) != -1) {
		switch (opt) {
			case 'h': /* --help */
				print_usage(); /* 'parms_setting' is a global variable */
				exit(EXIT_FAILURE);
				break;
			case 'r': /* --reads_binning_file */
				reads_binning_file = optarg;
				break;
			case 'T': /* --num_tissues */
				num_tissues = (int)atoi(optarg);
				break;
			case 'p': /* --min_likelihood_ratio_cutoff */
				min_likelihood_ratio_cutoff = atof(optarg);
				break;
			case 't': /* --input_file_of_tissue_markers */
				tissue_markers_file = optarg;
				break;
			case 'x': /* --em_max_iterations */
				em_max_iterations = (int)atoi(optarg);
				break;
			case 'a': /* --em_algorithm */
				em_algorithm_type = optarg;
				break;
			case 'o': /* --output_file */
				output_file = optarg;
				break;
			case 'O': /* --output_type */
				output_type = optarg;
				break;
			default: /* something unexpected has happened */
				cerr << endl << "Error: wrong parameters!" << endl << endl;
				print_usage();
				exit(EXIT_FAILURE);
		}
	}
	if (output_type!="tissueFraction" && output_type!="tissueFraction+readCountRaw" && output_type!="tissueFraction+readCountPerMillion" && output_type!="tissueFraction+readCountPerBillion") {
		cerr << endl << "Error: output types should take one of four values: tissueFraction, tissueFraction+readCountRaw, tissueFraction+readCountPerMillion and tissueFraction+readCountPerBillion!" << endl << endl;
		print_usage();
		exit(EXIT_FAILURE);
	}
	if (em_algorithm_type.compare("em.global.known")!=0 && em_algorithm_type.compare("em.global.unknown")!=0 && em_algorithm_type.compare("em.local.unknown")!=0 && em_algorithm_type.compare("em.local.unknown")!=0) {
		cerr << endl << "Error: EM algorithm types should take one of three values: em.global.known, em.global.unknown, em.local.known, em.local.unknown!" << endl << endl;
		print_usage();
		exit(EXIT_FAILURE);
	}
	if (reads_binning_file.empty()) {
		cerr << endl << "Error: input reads binning file is required!" << endl << endl;
		print_usage();
		exit(EXIT_FAILURE);
	}
	if (num_tissues<=0) {
		cerr << endl << "Error: number of tissue types must be >1!" << endl << endl;
		print_usage();
	}
	if (tissue_markers_file.empty()) {
		cerr << endl << "Error: input file of tissue markers is required!" << endl << endl;
		print_usage();
		exit(EXIT_FAILURE);
	}
}

int main(int argc, char* argv[]) {
	parse_command_line(argc, argv);

	cerr << "reads methylation states file: " << reads_binning_file << endl;
	cerr << "tissue markers file: " << tissue_markers_file << endl;
	cerr << "number of tissue types: " << num_tissues << endl;
	if (min_likelihood_ratio_cutoff==-1.0)
		cerr << "No likilihood-ratio-based reads filtering" << endl;
	else
		cerr << "min likelihood ratio cutoff (for each read): " << min_likelihood_ratio_cutoff << endl;
	cerr << "tissue deconvolution algorithm: " << em_algorithm_type << endl;
	cerr << "EM max iterations: " << em_max_iterations << endl;
	cerr << "output type: " << output_type << endl;
	cerr << "output file: " << output_file << endl;
	cerr << endl;

	cerr << "reading " << tissue_markers_file << " ..." << endl;
	string fileext_gzip = ".gz";
	Bins2Values marker2beta;
	vector<string> tissue_names;
	if (str_ends_with(tissue_markers_file,fileext_gzip)) {
		read_tissue_markers_gz_file(tissue_markers_file, 2, num_tissues, marker2beta, tissue_names);
	} else {
		read_tissue_markers_txt_file(tissue_markers_file, 2, num_tissues, marker2beta, tissue_names);
	}
	//string debug_file="./debug.single_value.tissue_markers.txt";
	//write_Bins2Values(marker2beta, tissue_names, debug_file, false);
	//exit(EXIT_SUCCESS);

	cerr << "calculating reads likelihoods ..." << endl;
	Matrix_Double reads_likelihoods((unsigned int)num_tissues);
	Bins2UnsignedIntegers marker2rowindexes;
	Bins2Value marker2ambiguousreadcounts;
	vector<int> Rm, Rl;
	unsigned long num_total_reads;
	num_total_reads = calc_read_probability_by_marker2beta(reads_binning_file, marker2beta, reads_likelihoods, marker2rowindexes, marker2ambiguousreadcounts, Rm, Rl, min_likelihood_ratio_cutoff);
	int num_reads_non_ambiguous = (int)reads_likelihoods.get_row_num();
	int num_reads_ambiguous = 0;
	Bins2Value::iterator m2amb_iter;
	for (m2amb_iter=marker2ambiguousreadcounts.begin(); m2amb_iter!=marker2ambiguousreadcounts.end(); m2amb_iter++) {
		num_reads_ambiguous += (int)(m2amb_iter->second);
	}
	cerr << "#reads_total: " << num_total_reads << endl;
	cerr << "#reads_covering_markers (non-ambiguous): " << num_reads_non_ambiguous << endl;
	cerr << "#reads_covering_markers (ambiguous): " << num_reads_ambiguous << endl;

	//////////////////////////////
	////// begin for debug
#ifdef DEBUG
	cout << "reads likelihoods: " << endl;
	cout << reads_likelihoods << endl;
	cout << "marker2ambiguousreadcounts: " << endl;
	cout << marker2ambiguousreadcounts << endl;
	//for (int z=0; z<Rm.size(); z++) {
		//cerr << Rm[z] << "\t" << Rl[z] << endl;
	//}
	//cout << endl;
	cout << "marker2rowindexes: " << endl;
	print_Bins2UnsignedIntegers(marker2rowindexes);
	cout << "c++_done" << endl;
#endif
	//// end for debug
	//////////////////////////////

	cerr << "perform " << em_algorithm_type << " ..." << endl;
	vector<double> theta;
	//cout << "debug before em_type.compare" << endl << flush;
	if (em_algorithm_type.compare("em.global.known")==0) {
		// create and initialize q with the same size of p and with all elements initialized as 0
		Matrix_Double q(num_reads_non_ambiguous, num_tissues, 0);
		// estimate tissue fractions (theta) and tissue-specific posterior probability matrix (q)
		em_supervise(reads_likelihoods, em_max_iterations, theta, q);

		// output
		if (output_type=="tissueFraction") {
			ofstream fout(output_file.c_str());
			if (fout.is_open()) {
				print_vec(fout, tissue_names, "\t", "#");
				fout << endl;
				fout.precision(6);
				print_vec(fout, theta, "\t", "#");
				fout << endl << "#reads_total: " << num_total_reads << "\t#reads_covering_markers(non_ambiguous): " << num_reads_non_ambiguous << endl << "\t#reads_covering_markers(ambiguous): " << num_reads_ambiguous << endl;
				fout.close();
			}
		} else if (output_type.find("readCount")!=std::string::npos) {
			cerr << "Count reads of each marker by reads posterior probabilities: " << output_type << " ..." << endl;
			double unit = 1e6 / num_total_reads; // default unit is readCountPerMillion
			if (output_type.find("readCountPerMillion")!=std::string::npos) {
				unit = 1e6 / num_total_reads;
				cerr << "   unit: 1e6 / #total_reads = " << unit << endl;
			} else if (output_type.find("readCountPerBillion")!=std::string::npos) {
				unit = 1e9 / num_total_reads;
				cerr << "   unit: 1e9 / #total_reads = " << unit << endl;
			} else if (output_type.find("readCountRaw")!=std::string::npos) {
				unit = 1.0;
				cerr << "   unit: " << unit << endl;
			}
			q.set_row_labels(reads_likelihoods.get_row_labels());

			Matrix_Double readCounts(num_tissues);
			readCounts_by_reads_posterior_probability_version_regular(q, unit, readCounts);
			ofstream fout(output_file.c_str());
			if (fout.is_open()) {
				print_vec(fout, tissue_names, "\t", "#marker_index\t");
				fout << endl;
				fout.precision(6);
				print_vec(fout, theta, "\t", "#tissue_fractions\t");
				fout << endl << "#reads_total: " << num_total_reads << "\t#reads_covering_markers(non_ambiguous): " << num_reads_non_ambiguous << "\t#reads_covering_markers(ambiguous): " << num_reads_ambiguous << "\tunit: " << unit << endl;
				fout << readCounts;
				fout.close();
			}
		}
		//////////////////////////////
		////// begin for debug
		cerr << "tissue fractions: ";
		print_vec(cerr, theta, ", ");
		cerr << endl;
		////// end for debug
		//////////////////////////////
	}
	if (em_algorithm_type.compare("em.global.unknown")==0) {
		// create and initialize q with the same size of p and with all elements initialized as 0
		Matrix_Double q(num_reads_non_ambiguous, num_tissues, 0);
		// estimate tissue fractions (theta) and tissue-specific posterior probability matrix (q)
		em_supervise(reads_likelihoods, em_max_iterations, theta, q);
		// adjust 'theta' and 'q' by "num_reads_non_ambiguous/(num_reads_ambiguous + num_reads_non_ambiguous)", because we consider 'theta_unknown'
		double theta_unknown = (double)num_reads_ambiguous/(num_reads_ambiguous + num_reads_non_ambiguous);
		double tmp = (double)num_reads_non_ambiguous/(num_reads_ambiguous + num_reads_non_ambiguous);
		for (int i=0; i<theta.size(); i++) {
			theta[i] *= tmp;
		}
		q.multipy_single_value(tmp);
		// output
		if (output_type=="tissueFraction") {
			ofstream fout(output_file.c_str());
			if (fout.is_open()) {
				print_vec(fout, tissue_names, "\t", "#");
				fout << "\tunknown" << endl;
				fout.precision(6);
				print_vec(fout, theta, "\t", "#");
				fout << endl << "#total_reads: " << num_total_reads << "\t#reads_covering_markers(non_ambiguous): " << num_reads_non_ambiguous << endl << "\t#reads_covering_markers(ambiguous): " << num_reads_ambiguous << endl;
				fout.close();
			}
		} else if (output_type.find("readCount")!=std::string::npos) {
			cerr << "Count reads of each marker by reads posterior probabilities: " << output_type << " ..." << endl;
			double unit = 1e6 / num_total_reads; // default unit is readCountPerMillion
			if (output_type.find("readCountPerMillion")!=std::string::npos) {
				unit = 1e6 / num_total_reads;
				cerr << "   unit: 1e6 / #total_reads = " << unit << endl;
			} else if (output_type.find("readCountPerBillion")!=std::string::npos) {
				unit = 1e9 / num_total_reads;
				//unit = 1.0;  // for debug only
				cerr << "   unit: 1e9 / #total_reads = " << unit << endl;
			} else if (output_type.find("readCountRaw")!=std::string::npos) {
				unit = 1.0;
				cerr << "   unit: " << unit << endl;
			}
			q.set_row_labels(reads_likelihoods.get_row_labels());

			Matrix_Double readCounts(num_tissues); // the last column is for the unknown class.
			readCounts_by_reads_posterior_probability_version_regular(q, unit, readCounts); // collapse from reads to markers
			//tmp = (double)num_reads_ambiguous/(num_reads_ambiguous + num_reads_non_ambiguous);
			for (m2amb_iter=marker2ambiguousreadcounts.begin(); m2amb_iter!=marker2ambiguousreadcounts.end(); m2amb_iter++) {
				m2amb_iter->second *= theta_unknown * unit;
			}
			ofstream fout(output_file.c_str());
			if (fout.is_open()) {
				print_vec(fout, tissue_names, "\t", "#marker_index\t");
				fout << "\tunknown" << endl;
				fout.precision(6);
				print_vec(fout, theta, "\t", "#tissue_fractions\t");
				fout << "\t" << theta_unknown << endl;
				fout << "#reads_total: " << num_total_reads << "\t#reads_covering_markers(non_ambiguous): " << num_reads_non_ambiguous << "\t#reads_covering_markers(ambiguous): " << num_reads_ambiguous << "\tunit: " << unit << endl;
				readCounts.print_with_additional_column_of_Bins2Value(fout, marker2ambiguousreadcounts);
				//fout << readCounts;
				fout.close();
			}
		}
		//////////////////////////////
		////// begin for debug
		cerr << "tissue fractions: ";
		print_vec(cerr, theta, ", ");
		cerr << ", " << theta_unknown << endl;
		////// end for debug
		//////////////////////////////
	}
}

