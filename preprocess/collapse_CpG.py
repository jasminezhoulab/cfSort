import operator
import re
import gzip
import sys

#usage: python collapse_CpG.py CpG_OT_file CpG_OB_file output_file_name

cpg_ot=sys.argv[1] #CpG_OT file (bismark output, still gzipped)
cpg_ob=sys.argv[2]  #CpG_OB file (bismark output, still gzipped)
outname=sys.argv[3] #output file

read_id={}

with open(outname,'a') as final:
	#Bismark's CpG_OT output file--all the fragments in this file mapped to the (+) strand, there is not R1/R2 in this file (if overall fragment mapped to (+) strand, both mates are present in this file, without any labeling. overlap (if any) is only counted once.
	with gzip.open(cpg_ot) as f_ot:
		next(f_ot) #skip header
		prev_read_name="first_read"
		for line in f_ot: #walk through CpG_OT file line by line, collapsing all lines corresponding to the same read onto one line in output (they are sorted by read name)
			dup_line=re.split(r'\t',line)
			read_name=dup_line[0]

			if read_name==prev_read_name: #same read as the one before it, continue collecting methylation info
				cpg_site=int(dup_line[3])
				meth_raw=str(dup_line[4]).strip()
				if meth_raw=="Z":
					meth_string.append(1)
				else:
					meth_string.append(0)
				cpg_string.append(cpg_site)
				num_cpgs=num_cpgs+1
			else: #new read, write the previous read to output file
				if (prev_read_name!="first_read"): #edge case for first line in file
					together=zip(cpg_string, meth_string)
					sorted_together=sorted(together) #the CpG_OT file is sorted WRT read name; so need to sort the methylation calls to be in numerical order (sort their position and call together so you don't mess up order_
					cpg_string_sorted=[x[0] for x in sorted_together]
					meth_string_sorted=[x[1] for x in sorted_together]
					final_cpg_string=",".join(map(str,cpg_string_sorted))
					final_meth_string="".join(map(str,meth_string_sorted))
					start=str(cpg_string_sorted[0])	#make it look like a bed file...but the coordinates are just the first and last CpG--NOT the coordinates of the actual read
					end=str(cpg_string_sorted[-1])
					#write to file
					final.write(chrom + "\t" + start + "\t" + end +"\t+\t" +str(num_cpgs)+ "\t" + str(final_cpg_string) + "\t" + str(final_meth_string) + "\t" + prev_read_name + "\n")
				
				#start reading in next read (or read in first read)
				prev_read_name=read_name[:]
				chrom=str(dup_line[2])
				cpg_site=int(dup_line[3])
				meth_raw=str(dup_line[4]).strip()
				if meth_raw=="Z":
					meth_string=[1]
				else:
					meth_string=[0]
				cpg_string=[int(cpg_site)]
				num_cpgs=1

		#edge case for last line in file
		together=zip(cpg_string, meth_string)
		sorted_together=sorted(together)
		cpg_string_sorted=[x[0] for x in sorted_together]
		meth_string_sorted=[x[1] for x in sorted_together]
		final_cpg_string=",".join(map(str,cpg_string_sorted))
		final_meth_string="".join(map(str,meth_string_sorted))
		start=str(cpg_string_sorted[0])
		end=str(cpg_string_sorted[-1])
		final.write(chrom + "\t" + start + "\t" + end +"\t+\t" +str(num_cpgs)+ "\t" + str(final_cpg_string) + "\t" + str(final_meth_string) + "\t" + prev_read_name + "\n")
		f_ot.close()

	#Bismark's CpG_OB output file--all the fragments in this file mapped to the (-) strand; same code as above but changed "+" to "-" in output file lines
	with gzip.open(cpg_ob) as f_ob:
		next(f_ob)
		prev_read_name="first_read"
		for line in f_ob:
			dup_line=re.split(r'\t',line)
			read_name=dup_line[0]

			if read_name==prev_read_name:
				cpg_site=int(dup_line[3])
				meth_raw=str(dup_line[4]).strip()
				if meth_raw=="Z":
					meth_string.append(1)
				else:
					meth_string.append(0)
				cpg_string.append(cpg_site)
				num_cpgs=num_cpgs+1
			else:
				if (prev_read_name!="first_read"):
					together=zip(cpg_string, meth_string)
					sorted_together=sorted(together)
					cpg_string_sorted=[x[0] for x in sorted_together]
					meth_string_sorted=[x[1] for x in sorted_together]
					final_cpg_string=",".join(map(str,cpg_string_sorted))
					final_meth_string="".join(map(str,meth_string_sorted))
					start=str(cpg_string_sorted[0])
					end=str(cpg_string_sorted[-1])
					final.write(chrom + "\t" + start + "\t" + end +"\t-\t" +str(num_cpgs)+ "\t" + str(final_cpg_string) + "\t" + str(final_meth_string) + "\t" + prev_read_name + "\n")

				prev_read_name=read_name[:]
				chrom=str(dup_line[2])
				cpg_site=int(dup_line[3])
				meth_raw=str(dup_line[4]).strip()
				if meth_raw=="Z":
					meth_string=[1]
				else:
					meth_string=[0]
				cpg_string=[int(cpg_site)]
				num_cpgs=1
		together=zip(cpg_string, meth_string)
		sorted_together=sorted(together)
		cpg_string_sorted=[x[0] for x in sorted_together]
		meth_string_sorted=[x[1] for x in sorted_together]
		final_cpg_string=",".join(map(str,cpg_string_sorted))
		final_meth_string="".join(map(str,meth_string_sorted))
		start=str(cpg_string_sorted[0])
		end=str(cpg_string_sorted[-1])
		final.write(chrom + "\t" + start + "\t" + end +"\t-\t" +str(num_cpgs)+ "\t" + str(final_cpg_string) + "\t" + str(final_meth_string) + "\t" + prev_read_name + "\n")
		f_ob.close()




