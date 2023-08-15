
###### REQUIRED INPUTS
unsorted_bam=WBC-006-G-N.bam   ## BISMARK-ALIGNED BAM FILE
genome_folder=hg19   ## PATH TO BOWTIE2-INDEXED HG19 GENOME
marker_bed=theoretical_RRBS_regions_cpg_nums_numbered_merged_index.bed.gz   ## INDEX OF THEORETICALLY SEQUENCED REGIONS IN RRBS ASSAY
marker_file=cfSort_markers.txt.gz   ## TISSUE MARKERS
tmp_dir=tmp     ## DIRECTORY TO SAVE TEMPORARY FILES (IT WILL BE REMOVED AT THE END)
output_dir_pipeline=.    ## DIRECTORY TO SAVE FINAL INPUT VECTOR
sample_id=WBC-006-G-N    ## SAMPLE ALIAS 
code_path=.    ## DIRECTORY THAT CONTAINS ALL THE PREPROCESSING SCRIPTS


###### REQUIRED TOOLS
bedtools=/PATH/TO/BEDTOOLS
samtools=/PATH/TO/SAMTOOLS
bismark_methylation_extractor=/PATH/TO/BISMARK_METHYLATION_EXTRACTOR
python2=/PATH/TO/PYTHON2
python3=/PATH/TO/PYTHON3
Rscript=/PATH/TO/Rscript






cd ${output_dir_pipeline}


### extract methylation

out_bed_sorted=${tmp_dir}/${sample_id}.sorted.bed.gz
out_bam_sortn=${tmp_dir}/${sample_id}.sortn.bam

${bedtools} bamtobed -i ${unsorted_bam} |sort -k1,1 -k2,2n | gzip > ${out_bed_sorted}
${samtools} sort -n ${unsorted_bam} -o ${out_bam_sortn}
${bismark_methylation_extractor} -p --bedGraph --gzip --counts --buffer_size 11G --genome_folder ${genome_folder} --output ${tmp_dir} ${out_bam_sortn}   
rm ${out_bam_sortn}


### get annotated bed

temp_dir=${tmp_dir}
output_dir_methyl=${tmp_dir}
output_dir_bed=${tmp_dir}
sample_name=${sample_id}
int_frag=${temp_dir}/${sample_name}.bed.fragment_level_TEMP
CpG_OT=$(\ls ${output_dir_methyl}/CpG_OT*${sample_name}.*gz)
CpG_OB=$(\ls ${output_dir_methyl}/CpG_OB*${sample_name}.*gz)
int_meth=${temp_dir}/${sample_name}.all_md.sorted.bed
refo_meth=${temp_dir}/${sample_name}_meth_refo.bed
refo_frag=${temp_dir}/${sample_name}_frag_refo.bed
out_bed_sorted=${output_dir_bed}/${sample_name}.sorted.bed.gz
temp_bed=${temp_dir}/${sample_name}.sorted.bed.fragment_level.meth.tmp
out_bed=${output_dir_bed}/${sample_name}.sorted.bed.fragment_level.meth.gz

${python2} ${code_path}/collapse_bed_file_strand_correct.py ${out_bed_sorted} ${int_frag} 
awk 'BEGIN {FS=OFS="\t"} {print $1, $2, $3, $5, $6, $4}' ${int_frag} | awk 'BEGIN {FS=","} {print $1}'  | sort -k6 > $refo_frag 
rm -f $int_frag
${python2} ${code_path}/collapse_CpG.py ${CpG_OT} ${CpG_OB} $int_meth
sort -k8 $int_meth > $refo_meth 
rm -f $int_meth
join -t $'\t' -1 8 -2 6 -a 2 -e'NA' -o '2.1 2.2 2.3 2.6 2.4 2.5 1.5 1.6 1.7' $refo_meth $refo_frag > $temp_bed 
sort -k1,1 -k2,2n $temp_bed |gzip > $out_bed
rm -f $temp_bed



### convert to binary meth

bed_dir=${tmp_dir}
temp_dir=${tmp_dir}
out_dir=${tmp_dir}
sample=${sample_id}
sample_bed=${bed_dir}/${sample}.sorted.bed.fragment_level.meth.gz
bed_one_hit_py_input=${temp_dir}/${sample}.RRBS_by_strand_for_py
py_out=${out_dir}/${sample}.binary_meth_values.txt.gz

${bedtools} intersect -sorted -a $sample_bed -b $marker_bed -f .5 -s -wo | awk 'BEGIN {FS=OFS="\t"} {print $20, $9, $6, $8}'  > $bed_one_hit_py_input
${python2} ${code_path}/theor_RRBS_binary_meth_values.py $bed_one_hit_py_input $py_out


### extract features

input_file=${py_out}
sample_alias=${sample_id}

cd ${output_dir_pipeline}

${python3} ${code_path}/convert_binary_meth_values_file_to_alpha_value_distribution_file.py ${input_file} ${sample_alias}.alpha_value_distr.txt.gz
echo ${sample_alias}.alpha_value_distr.txt.gz > TMP_SAMPLE
${python3} ${code_path}/generate_data_matrix_by_alpha_value_files_and_alpha_value_markers.py hypo_cfsort ${marker_file} TMP_SAMPLE ${sample_alias}.mat.txt.gz
rm TMP_SAMPLE
${Rscript} ${code_path}/collect_read_depth_at_markers.R ${sample_alias}.alpha_value_distr.txt.gz ${sample_alias}.depth.txt.gz ${marker_file}
${python3} ${code_path}/normalize_local_read_depth_and_cluster_markers.py ${marker_file} ${sample_alias}.mat.txt.gz ${sample_alias}.depth.txt.gz ${sample_alias}.features.txt.gz

rm ${sample_alias}.mat.txt.gz ${sample_alias}.depth.txt.gz ${sample_alias}.alpha_value_distr.txt.gz 
rm -rf ${tmp_dir}



