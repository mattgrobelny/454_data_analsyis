#!/usr/bin/bash -x
# requires seq_crumbs 
# requires SPades 

files="../mid_MID1.sff
../mid_MID2.sff
../mid_MID3.sff
../mid_MID4.sff
../mid_MID5.sff
../mid_MID6.sff"

time=$(date)
echo $time > mid1_6_cleaning_data.txt
echo "Cleaning of:" >> mid1_6_cleaning_data.txt
echo $files >> mid1_6_cleaning_data.txt

for file in $files
do
	echo ""
	echo "working on $file"

	# extract file name string 
	fastq=".fastq"
	
	filename=$(echo $file | grep -oE "mid_MID[0-9]")


	# extract sff to fastq
	echo "extracting sff"
	sff_extract $file > $filename$fastq 


	# trim adapters 
	echo "Trimming adapters..." 
	adapter_trim="_AT"
	filename_out1=$filename$vector_trim$adapter_trim
	adapter_array=($(tagcleaner -predict -fastq mid_MID1.fastq | sed -E "1d; s/tag.\t([ATGCN]+)\t.+\t.+/\1\n/"))
	tagcleaner -tag3 ${adapter_array[0]} -tag5 ${adapter_array[1]} -out_format 3 -out $filename_out1 -fastq $filename$fastq

	num_nuc_end2=$(sed -n 2~4p $filename$vector_trim$adapter_trim$fastq | tr -d '\n' | wc -m)
	num_nuc_start2=$(sed -n 2~4p $filename$fastq | tr -d '\n' | wc -m)
	delta2=$(bc <<< "scale = 2; (1-($num_nuc_end2/$num_nuc_start2))*100")
	echo "Tags used for $file:" >> mid1_6_cleaning_data.txt
	echo "-tag3 ${adapter_array[0]}" >> mid1_6_cleaning_data.txt
	echo "-tag5 ${adapter_array[1]}" >> mid1_6_cleaning_data.txt
	echo "Percent of nucleotides trimmed $delta2 %" >> mid1_6_cleaning_data.txt

	# trim edges 
	left_clip=10
	right_clip=10
	trim_edge="_TE"
	echo "Clipping edges..." 
	echo "Left clip: $left_clip"
	echo "Right clip: $right_clip"
	filename_out2=$filename$vector_trim$adapter_trim$trim_edge$fastq
	trim_edges -l $left_clip -r $right_clip -o $filename_out2 $filename_out1$fastq

	num_nuc_start3=$(sed -n 2~4p $filename_out1$fastq | tr -d '\n' | wc -m)
	num_nuc_end3=$(sed -n 2~4p $filename_out2 | tr -d '\n' | wc -m)
	delta3=$(bc <<< "scale = 2; (1-($num_nuc_end3/$num_nuc_start3))*100")
	echo "Percent of nucleotides trimmed using edge clipping: ~ $delta3 %" >> mid1_6_cleaning_data.txt
	
	# filter by blast 
	echo "filtering by blast "
	vector_trim="_VF"
	filtered_out_file="_filtered_out_seqs.fastq"
	filter_id=95
	filter_by_blast -b vectors_454.fasta -e $filename$filtered_out_file -s $filter_id -o $filename$vector_trim$fastq $filename_out2
	
	# Write how much was filtered
	echo "Filtering by blast stats for: $filename.sff" >> mid1_6_cleaning_data.txt
	start_num=$(grep -c "^@I" $filename$fastq)
	echo "Starting num of reads: $start_num" >> mid1_6_cleaning_data.txt
	after_filter=$(grep -c "^@I" $filename$vector_trim$fastq)
	echo "After filtering: $after_filter" >> mid1_6_cleaning_data.txt
	delta=$(bc <<< "scale = 2; $start_num-$after_filter")
	echo "Num of reads filtered out: $delta"

	
	# trim by quality
	echo "Trimming by quality..."
	quality_thresh=20
	quality_clip="_QC"
	echo "Quality threshold: $quality_thresh" >> mid1_6_cleaning_data.txt
	filename_out3=$filename$vector_trim$adapter_trim$trim_edge$quality_clip$fastq
	trim_quality -q $quality_thresh -o $filename_out3 $filename$vector_trim$fastq

	num_nuc_start4=$(sed -n 2~4p $filename$vector_trim$fastq | tr -d '\n' | wc -m)
	num_nuc_end4=$(sed -n 2~4p $filename_out3 | tr -d '\n' | wc -m)
	delta4=$(bc <<< "scale = 2; (1-($num_nuc_end4/$num_nuc_start4))*100")
	echo "Percent of nucleotides trimmed using quality clipping: $delta4 %" >> mid1_6_cleaning_data.txt

	# assembly with Spades
	echo "assembling with Spades"
	out_file="~/Data/mids_all/cleaned_mids/Spades_output_trimmed_careful_"
	if (("$filename" == "mid_MID6"));
	then	
	spades.py --trusted-contigs ~/Desktop/Spades_run/trust.fasta -m 10 -t 16 --careful --s1 $filename_out3 -o $out_file$filename

	else
	spades.py -m 10 -t 16 --careful --s1 $filename_out3 -o $out_file$filename
	fi 
	
	# Pull out assemby info 
	#cd $out_file$filename 
	contigs="/contigs.fasta"
	contigs_file=$out_file$filename$contigs
	contig_stats=($(cat $contigs_file |grep ">" | sed -E 's/>NODE_([0-9]+)_length_([0-9]+)_cov_([0-9]+)/\1\t\2\t\3/'))
	numb_contigs=($(cat $contigs_file |grep ">" | sed -E 's/>NODE_([0-9]+)_length_([0-9]+)_cov_([0-9]+)/\1\t\2\t\3/'| cut -f 1 | tail -n 1))
	largest_contig=($( cat $contigs_file |grep ">" | sed -E 's/>NODE_([0-9]+)_length_([0-9]+)_cov_([0-9]+)/\1\t\2\t\3/'| cut -f 2 | head -n 1))
		
	#write contig stats to file 
	echo "Contig stats for $filename_out3 assembly" >> mid1_6_cleaning_data.txt
	echo "Total contigs: $numb_contigs" >> mid1_6_cleaning_data.txt
	echo "Largest contig: $largest_contig" >> mid1_6_cleaning_data.txt
	echo
	echo "#-------------------------------------------------------------------------------#"

	assembly_data="_Asm_Data.tsv"
	echo "Contig_num	Contig_len	Contig_cov" > $filename$vector_trim$adapter_trim$trim_edge$quality_clip$assembly_data
	echo "$contig_stats >> $filename$vector_trim$adapter_trim$trim_edge$quality_clip$assembly_data"
	echo "done with $file.sff"
done




