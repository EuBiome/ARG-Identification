#!/bin/bash

#Processing pipeline:

echo Insert the paths of the read files

echo Forward:

read read_1
#read_1 = "$read_1"

echo Reverse:

read read_2
#read_2 = "$read_2"


basefolder=${read_1%/*}

echo Insert the path of the kraken2 database that you want to use
read db_path

echo Insert the path of the kraken2 plasmid database that you want to use
read plasmid_db_path

echo Do you want to precess the reads with Trimmomatic? type yes or no
read answer

if [ "$answer" = "yes" ]; then
	echo Insert the type of operation you want to make with Trimmomatic es:SLIDINGWINDOW:20:26
	read operation

	echo Trimmomatic just began

	trimmomatic PE -threads 4 "$read_1" "$read_2" "$basefolder/1P.fastq" "$basefolder/1U.fastq" "$basefolder/2P.fastq" "$basefolder/2U.fastq" "$operation"
else
	echo We proceed directly with the assembly
fi

echo We are now assemblying the reads with MEGAHIT

if [ "$answer" = "yes" ]; then
	megahit --k-list 29,39,59,79 --min-contig-len 500 -o "$basefolder/Megahit" -1 "$basefolder/1P.fastq" -2 "$basefolder/2P.fastq"
else
	megahit --k-list 29,39,59,79 --min-contig-len 500 -o "$basefolder/Megahit" -1 "$read_1" -2 "$read_2"
fi

echo MEGAHIT has terminated

echo Now we are going to run rgi over the contigs 

source ~/anaconda3/etc/profile.d/conda.sh
conda activate rgi

rgi main -i "$basefolder/Megahit/final.contigs.fa" -o "$basefolder/rgi" --clean --low_quality

conda deactivate

echo rgi has termined 
echo The classification is starting with kraken2


kraken2 --memory-mapping --db "$db_path" --use-names --report "$basefolder/k2_report" --output "$basefolder/k2_results" "$basefolder/Megahit/final.contigs.fa"
kraken2 --db "$plasmid_db_path" --use-names --report "$basefolder/k2_plasmid_report" --output "$basefolder/k2_plasmid_results" "$basefolder/Megahit/final.contigs.fa"

echo The pipeline is termined
