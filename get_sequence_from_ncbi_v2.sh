#!/usr/bin/bash

Help()
{
   # Display Help
   echo "This script reads NCBI acc numbers from a file and retrieves the sequences from the specified NCBI db."
   echo
   echo "Syntax: ${0} [-i|-d|-o|-h]"
   echo "options:"
   echo "-i	input file with accession numbers."
   echo "-d	ncbi database."
   echo "-o	output file."
   echo "-h	print this help."
   echo
}

while getopts "hi:d:o:" option; do
   case $option in
	h) # display Help
        	Help
        	exit;;
	i) INPUT=${OPTARG};;
	d) DB=${OPTARG};;
	o) OUTPUT=${OPTARG};;
	\?) # incorrect option
		echo
        	echo "Error, Invalid option"
		echo
	 	Help
         	exit;;
   esac
done

# Job Name
#$ -N getseq

# Execute the script from the Current Working Directory
#$ -cwd

# Merge the output of the script, and any error messages generated to one file
#$ -j y

# Send the output of the script to a directory called 'UGE-output' uder current working directory (cwd)
if [ ! -d "getseq_output" ]; then #Create output directory in case it does NOT exist
    mkdir getseq_output
fi
#$ -o getseq_output/

# Tell the job your cpu and memory requirements
#$ -pe threaded 8
# -l mem_free=20G,h_vmem=24G

# Send mail when the job is submitted, and when the job completes
#$ -m be

#  Specify an email address to use
#$ -M lorenziha@nih.gov


pools=200
counter=0
ACC=''
while IFS= read -r acn 
do
	counter=$counter+1
	if (( counter ==  pools )); then
		echo ${ACC}
		wget -O getseq_tmp.fasta  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=${DB}&id=${ACC}&rettype=fasta"
		cat getseq_tmp.fasta >> ./getseq_output/${OUTPUT}
		counter=0
		ACC=''
		echo
	fi
	
	if [ "$ACC" != '' ]; then
		ACC=${ACC},${acn}
	else 
		ACC=${acn}
	fi
done < "${INPUT}"
rm getseq_tmp.fasta

# Run leftover ACCs
wget -O getseq_tmp.fasta  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=${DB}&id=${ACC}&rettype=fasta"
cat getseq_tmp.fasta >> ./getseq_output/${OUTPUT}
rm getseq_tmp.fasta



