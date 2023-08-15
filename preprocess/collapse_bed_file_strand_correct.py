#SCRIPT TO COLLAPSE BED FILE (BEDTOOLS BAMTOBED OUTPUT) TO FRAGMENT LEVEL FOR PE READS
#input: bed file of paired end sequencing data produced by bamtools bamtobed -i ${bam_file} (unzipped or gzipped) (sort it first to avoid memory issues, but not necessary...sort -k1,1 -k2,2n $int_bed > $bed)
#output: bed file where each line is from a collapsed fragment: R1 and R2, strand of fragment is the strand R1 mapped to (unzipped)

#usage python collapse_bed_file_strand_correct.py $input_bed $output_bed
import re
import sys
import gzip
fname=sys.argv[1] #$input_bed, can be gzipped or not
outname=sys.argv[2] #$output_bed

read_id={}

if fname.endswith("gz"):
	f=gzip.open(fname)
else:
	f=open(fname)

with open(outname,'w') as final:
	for line in f:
		dup_line=re.split(r'\t',line)
		chrom=dup_line[0]
		start_frag=int(dup_line[1])
		end_frag=int(dup_line[2])
		input_strand=str(dup_line[5]).strip()
		read_name=re.split(r'/',dup_line[3])
		if len(read_name)>=2:
			read_pair=int(read_name[(len(read_name)-1)])
			if (read_pair==1):
				strand=input_strand.strip()
			else:
				if input_strand=="+":
					strand="-"
				else:
					strand="+"
			uniq_read_name=read_name[0]
			key=uniq_read_name
			if key in read_id: #already seen this read's mate, so collapse and write to file
				cur_start=read_id[key][1]
				cur_end=read_id[key][2]
				if start_frag<cur_start:
					read_id[key][1]=start_frag
				if end_frag>cur_end:
					read_id[key][2]=end_frag

				name=str(key)
				val=read_id[key]
				chrom=str(val[0])
				start=str(val[1])
				end=str(val[2])
				strand=str(val[3].strip())
				fragment_length=str(int(end)-int(start))
				final.write(chrom+"\t"+start+"\t"+end+"\t"+name+","+chrom+","+start+","+end+"\t"+fragment_length+"\t"+strand+"\n")
				read_id.pop(key, None) #pop key to avoid memory issues (if file is sorted R1 and R2 should be ~close together and won't have to store too many keys at once)
			else: #first time seeing this read, store info in dictionary
				read_id[key]=[chrom,start_frag,end_frag, strand]
		elif len(read_name)==1:
			next
		
print len(read_id)


