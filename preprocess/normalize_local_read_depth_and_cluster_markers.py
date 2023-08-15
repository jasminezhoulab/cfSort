
import sys
import gzip
import numpy as np
from collections import defaultdict

marker_file = sys.argv[1]
readcount_file = sys.argv[2]
depth_file = sys.argv[3]
output_file = sys.argv[4]

marker = defaultdict(list)
with gzip.open(marker_file, 'rt') as fmarker:
	line = fmarker.readline()
	sp = line.strip().split('\t')
	marker_cluster_index = sp.index("marker_cluster_index")
	row_id = 0
	for line in fmarker:
		sp = line.strip().split('\t')
		marker_cluster = sp[marker_cluster_index]
		marker[marker_cluster].append(row_id)
		row_id += 1

marker_key = marker.keys()


n_sample = 0
with gzip.open(output_file, 'wt') as foutput:
	with gzip.open(depth_file, 'rt') as fdepth:
		with gzip.open(readcount_file, 'rt') as freadcount:
			for readcount_line in freadcount:
				depth_line = fdepth.readline()
				readcount_sp = readcount_line.strip().split('\t')
				depth_sp = depth_line.strip().split('\t')
				assert depth_sp[0] == readcount_sp[0], "unmatched sample id"
				sample_id = depth_sp[0]
				readcount = np.array([float(readcount_sp[i]) if readcount_sp[i] != "NA" and readcount_sp[i] != '' else 0.0 for i in range(1, len(readcount_sp))])
				depth = np.array([float(depth_sp[i]) if depth_sp[i] != "NA" and depth_sp[i] != '' else 0.0 for i in range(1, len(depth_sp))])
				output = sample_id
				for key in marker_key:
					merge_depth = sum(depth[marker[key]])
					merge_readcount = sum(readcount[marker[key]])
					assert merge_readcount <= merge_depth, "wrong read count, more than depth" 
					if merge_depth == 0.0:
						frac = 0.0
					else:
						frac = merge_readcount/merge_depth
					output = output + '\t' + str(frac) 
				foutput.write(output + '\n')
				print(n_sample)
				n_sample += 1







