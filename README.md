# cfSort

## INTRODUCTION

cfSort is a versatile toolset of deconvoluting the cell-free DNA fragments (or sequencing reads) into multiple tissue types or tumor types. cfSort can be freely used for educational and research purposes by non-profit institutions and U.S. government agencies only under the UCLA Academic Software License. For information on the use for a commercial purpose or by a commercial or for-profit entity, please contact Prof. Xiangong Jasmine Zhou (https://zhoulab.dgsom.ucla.edu/).

Given the cell-free DNA methylation sequencing data of a patient’s cfDNA sample, for each cancer methylation marker (or tissue methylation marker), we deconvolve the tumor-derived (or tissue-specific) sequencing reads from all reads within the marker's genomic region. This individual-sequencing-read-based deconvolution framework exploits the pervasiveness of DNA methylation for signal enhancement. Specifically, we improved upon our previous probabilistic read deconvolution algorithm, i.e., CancerDetector [1], by (1) adding an “unknown” class to represent reads that cannot be classified to any known class (tumor type or tissue type), and (2) expanding the 2-class likelihood model to a k-class (k>=2)  model, for classifying reads into k different classes, and (3) For a given set of markers, we construct a profile vector where the length of the vector is the number of markers and the value in each entry is the normalized counts of tumor-derived (tissue-derived) reads. 

Note that this code can deconvolute the different sources of cfDNA reads in two contexts: (1) separating reads into tumor-derived reads and background reads; and (2) separating the reads from different tissues. 

References:
[1] Li W, Li Q, Kang S, Same M, Zhou Y, Sun C, Liu CC, Matsuoka L, Sher L, Wong WH, Alber F, Zhou XJ. CancerDetector: ultrasensitive and non-invasive cancer detection at the resolution of individual reads using cell-free DNA methylation sequencing data. Nucleic Acids Res. 2018 Sep 6;46(15):e89. doi: 10.1093/nar/gky423. PMID: 29897492; PMCID: PMC6125664.

## PREREQUISITE PACKAGES

cfSort was developed on UCLA Hoffman2 cluster. All the analyses in this study were done there. The system information of UCLA Hoffman2 cluster is:
Linux n6015 2.6.32-754.14.2.el6.x86_64 #1 SMP Tue May 14 19:35:42 UTC 2019 x86_64 x86_64 x86_64 GNU/Linux
There is no non-standard hardware required.

1) Boost C++ Libraries, version 1.55.0
2) C++ Library zlib, version 1.2.8

## INSTALLATION

The code is programmed by C++. For running the code, please go to the directory "src/" and do the following steps:

1. Edit Makefile to provide the directories of include and library of three libraries: Boost and zlib.

2. Type the following command to build executable file "cfsort":

cd src
make -f ./Makefile

3. We also provided the executable file "cfsort", so that you can directly run this code in Linux system, and do not need to build the executable file by yourself.

## USAGE

Inputs:
1) input markers file, in plain text
2) input cfDNA methylation sequencing reads file, either in plain text or compressed form.
3) number of tissue types
4) likelihood ratio threshold: a positive float number. We suggest this number is 2.

Outputs:
5) output file format: "tissueFraction+readCountPerBillion"
6) output filename

If the input file of cfDNA methylation sequencing reads is gzip file (a compressed file), please use the following usage:

gzip -dc <input_file_of_cfDNA_reads> | cfsort -a em.global.unknown -p <likelihood_ratio_threshold> -O tissueFraction+readCountPerBillion -r stdin -t <input_file_of_markers> -o <output_file_of_read_counts_profile> -T <number_of_tissues>

If the input file of cfDNA methylation sequencing reads is the plain text file, please use the following usage:

cfsort -a em.global.unknown -p <likelihood_ratio_threshold> -O tissueFraction+readCountPerBillion -r <input_file_of_cfDNA_reads> -t <input_file_of_markers> -o <output_file_of_read_counts_profile> -T <number_of_tissues>

### INPUT FILE FORMAT

+++ input cfDNA methylation sequencing reads file
Each line is a read. All columns are delimited by TAB. There is one header line.
We used three column to represent each read's information as below:

Column 1: marker_index, an integer number that is a marker index.
Column 2: number of methylated CpG of the read that is mapped to this marker's genomic region.
Column 2: number of unmethylated CpG of the read that is mapped to this marker's genomic region.

For example:
marker_index	meth_count	unmeth_count
1	5	2
1	3	4
2	8	1
2	3	5
...

An example file see "demo/input/example.reads.txt.gz"

+++ markers file format
Each line is a marker including its genomic region and methylation pattern (two shape parameters of the beta distribution, which are delimited by :). All columns are delimited by TAB. There is one header line.
Column 1: marker_index
Column 2+: paired values (shape parameters of the beta distribution) for this marker (each column is a class). Each pair contains two values, delimited by ":". If you cannot estimate shape parameters of the beta distribution, you may also use a single value, which is the average methylation level of this marker across all samples in the same tissue type.

For example:
marker.index	class1	class2	class3	class4	class5	class6	class7
1	1.68:1.82	2:18	1.02:0.459	2.08:18	0.061	0.176	6.65:2.19
2	0.605	0.47	0.667	29.1:12.1	0.0528	0.542	7.31:3.47
3	12.8:2.68	9.35:1.75	19.4:4.06	6.28:7.06	8.44:20.4	21.1:3.48	6.69:4.01
...

An example file see "demo/input/example.markers.txt"


### OUTPUT FILE FORMAT

The output file is the normalized tissue-specific read counts in each marker. The first three lines are header lines, indicating the column names and the global fraction of tissue types and other information. Starting from the 4-th line, each line corresponds to a marker, consisting of the tissue-specific read counts of this marker.
Column 1: marker index
Column 2+: read count of each tissue type in the marker, and the read count of the unknown class is listed in the last column of each line

For example:
#marker_index	class1	class2	class3	class4	class5	class6	class7	unknown
#tissue_fractions	1.58246e-13	1.35298e-19	0.0442296	1.45249e-21	0.0157226	4.00821e-18	0.930494	0.00955414
#reads_total: 942	#reads_covering_markers(non_ambiguous): 933	#reads_covering_markers(ambiguous): 9	unit: 1.06157e+06
1	3.05614e-05	4.69734e-12	1.60738e+07	5.32786e-14	240091	3.10029e-10	4.33698e+08	0
2	3.39305e-06	2.33859e-12	1.11051e+06	4.09455e-14	3602.4	7.64251e-11	2.41202e+07	0
...

An example file see "demo/output/example.class.specific.read.counts.profile_reference"


## DEMO: EXAMPLE DATA

A demo is provided along with the scripts. Example data and required reference files are saved under demo folder.

This demo can be run with simple command lines. This demo was expected to run about 30 seconds. 
Please run the following commands:

### run cfsort
cd demo
./run_demo.sh

If you see the following message, then it runs well in your system:

cfsort Test run passed

