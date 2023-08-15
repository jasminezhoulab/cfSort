import re
import sys
import random
import gzip
import itertools
from itertools import islice
from itertools import izip_longest
from itertools import izip

fname=sys.argv[1]
outname=sys.argv[2]

with gzip.open(outname, 'w') as final:
	final.write("marker_index\tcpg_locs\tmeth_string\tmeth_count\tunmeth_count\tstrand\n")
	with open(fname, 'r') as f:
		for line in f:
			li=line.strip()
			#print(li)

			sp=re.split(r'\t', li)
			marker_index=int(sp[0])
			cpg_locs=str(sp[3])
			strand=str(sp[2])

			meth_string=str(sp[1])
			count_meth=meth_string.count("1")
			count_unmeth=meth_string.count("0")

			final.write(str(marker_index) + "\t" + str(cpg_locs) + "\t" + str(meth_string) + "\t"+ str(count_meth) + "\t" + str(count_unmeth) + "\t" + str(strand) +"\n")
