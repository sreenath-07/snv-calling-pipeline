#!/bin/bash


default_path () {
	#Incase your software files are not on your path, use sudo cp bcftools /usr/local/bin or add to path using export PATH=$PATH:<file_path>
	gatk_file=GenomeAnalysisTK.jar
	bcf_file=bcftools
	sam_file=samtools
	bwa_file=bwa
	java_file=java
}

default_path

help_argument () {
	echo "
	PRE-REQUISITES
	To run this script please ensure that the following software/files are installed/downloaded and added to your PATH:

	SOFTWARE:
	1. bwa
	2. samtools
	3. gatk V3.7.0 (broad)
	4. bcftools
	5. java v1.8 (for compatibility with gatk)

	FILES:
	1. 1st Reads File
	2. 2nd Reads File 
	3. Reference genome
	4. Mills vcf file (For human genome - https://storage.cloud.google.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz?_ga=2.132323273.-613126285.1633368214)


	INPUT ARGUMENTS
	-a	Input reads file – pair 1
	-b	Input reads file – pair 2
	-r	Reference genome file
	-e	Perform read re-alignment
	-o	Output VCF file name
	-f	Mills file location
	-z	Output VCF file should be gunzipped (*.vcf.gz)
	-v	Verbose mode; print each instruction/command to tell the user what your script is doing right now
	-i	Index your output BAM file (using samtools index)
	-h	Print usage information (how to run your script and the arguments it takes in) and exit
	"
}

#Default values of arguments
realign=0
gunzip=0
v=0
index=0
# output="na"
# millsFile="na"

#Parsing input arguments
while getopts "a:b:r:eo:f:zvih" option
do
	case $option in
		a)reads1=$OPTARG;; #1st Read
		b)reads2=$OPTARG;; #2nd Read
		r)ref=$OPTARG;;     #Ref Genome Location
		e)realign=1;; # Perform read re-alignment
		o)output=$OPTARG;; #Output vcf file name
		f)millsFile=$OPTARG;; #Location of Mills Indel File
		z)gunzip=1;;
		v)set -x ;;#v=1;; #verbose=T
		i)index=1;; #Index your output BAM file (using samtools index)
		h)help_argument;;
		*)echo 'Error in command line parsing'
	esac
done

file_checks () {

	## Checking for Input Files
	if test ! -f "$reads1"; then
		echo "1st reads file is missing"
		exit
	fi

	if test ! -f "$reads2"; then
		echo "2nd reads file is missing"
	fi

	## Checking for reference genome
	if test ! -f "$ref"; then
		echo "Reference genome file is missing"
		exit
	fi

	## Checking if the output vcf exists already
	if (test -f "$output.vcf") || (test -f "$output.vcf.gz") ; then
		echo "Output VCF file already exists. Enter y to exit the program. To overwrite the existing file please enter anything except y." 
		read answer
		if [[ $answer == "y" ]]; then
			echo "Exiting"
			exit
		fi
		echo "Continuing"
	fi

	# ## Checking if output file name is specified, by comparing with default flag
	# if [[ $output == "na" ]]; then
	# 	echo "Output file argument is missing"
	# 	exit
	# fi

	# ## Checking if the millsfile is specified
	# if [[ $millsFile == "na" ]]; then
	# 	echo "Mills file argument is missing"
	# 	exit
	# fi

}

mapping () {
	#MAPPING
	#pre-requisites - bwa installed
	bwa index $ref
	bwa mem -R '@RG\tID:foo\tSM:bar\tLB:library1' $ref $reads1 $reads2 > lane.sam 
	samtools fixmate -O bam lane.sam lane_fixmate.bam
	samtools sort -O bam -o lane_sorted.bam -T /tmp/lane_temp lane_fixmate.bam
	samtools index lane_sorted.bam
}

generate_index_dict_files () {
	samtools faidx $ref #Creating an index file
	samtools dict $ref -o "${ref%.*}.dict" #Creating a dict file. # %.* removes suffix starting with .
}

realignment_improvement () {
	#REALIGNMENT
	#pre-requisites - java 8, gatk, mills file downloaded/installed & indexing done

	$java_file -Xmx2g -jar $gatk_file -T RealignerTargetCreator -R $ref -I lane_sorted.bam -o lane.intervals --known $millsFile 2> ssrikrishnan6.log #--log_to_file gatk.log
	$java_file -Xmx4g -jar $gatk_file -T IndelRealigner -R $ref -I lane_sorted.bam -targetIntervals lane.intervals -known $millsFile -o lane_realigned.bam 2>> ssrikrishnan6.log
}

index_improvement () {
	#Index bam using samtools
	samtools index lane_realigned.bam 	
}

variant_calling () {
	#VARIANT CALLING
	if [[ realign==1 ]]; then
		bcftools mpileup -Ou -f $ref lane_realigned.bam | bcftools call -vmO z -o $output.vcf.gz
	else
		bcftools mpileup -Ou -f $ref lane_sorted.bam | bcftools call -vmO z -o $output.vcf.gz #Using the lane_sorted file which is created before re-alignment
	fi 
}

vcf_to_bed (){
	#Converting the output vcf file into bed format and generating snps and indels 
	#Pre-requisite - output vcf file is unzipped to .vcf format
	sed '/^#/d' $output.vcf| awk '{print substr($1,4) "\t" $2 "\t" $2+length($5)-length($4) "\t" length($5)-length($4)}' > bedfile.txt
	awk '{if ($4=="0") print $0}' bedfile.txt > ssrikrishnan6_snps.txt 
	awk '{if ($4!="0") print $0}' bedfile.txt > ssrikrishnan6_indels.txt 
}

###### Pipeline Begins ########
file_checks #Performing file check

mapping #Executing the mapping function

generate_index_dict_files #Executing the index_dict function

if [[ realign==1 ]]; then
	realignment_improvement #Executing the realignment function
fi

if [[ index==1 ]]; then
	index_improvement #Executing the indexing function
fi

variant_calling #Executing the variant calling function

if [[ gunzip==0 ]]; then
	gzip -dk $output.vcf.gz #if gunzip flag isn't mention, the file is unzipped
fi

#vcf_to_bed #Executing the bed conversion function

###### Pipeline Ends ########

