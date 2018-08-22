#!/bin/bash
#<add job submission details>

## navigate to directory
cd <path to genome assemblies>
mkdir B_cereus_de_novo_assembly/raw_data
cd B_cereus_de_novo_assembly/raw_data

## Download raw read files from NCBI-SRA
module load sratoolkit/2.8.1
fastq-dump <list of files to retrieve> --split-files # see Supp. Table 1 for SRR numbers

## Pre-process reads with Trim Galore: remove adapters, trim at low quality bases
module load TrimGalore/0.4.4
module load python

# run loop to process all R1/R2 files
files_trim=`ls | grep "_1.fastq"`
for read1 in $files_trim; do read2=`echo $read1 | cut -d '_' -f1`_2.fastq; trim_galore --paired -q 30 $read1 $read2; done
rm *txt
mkdir trimmed
mv *val*.fq trimmed

# zip SRA files
gzip *.fastq


### RUN SPADES ###

modules purge
mkdir spades

# run loop to process all trimmed R1/R2 files
files_spades=`ls trimmed | grep "val_1.fq"`
for read1 in $files_spades
do
  	read2=`echo $read1 | cut -d '_' -f1`_2_val_2.fq
                out=`echo $read1 | cut -d '_' -f1`_spades_out
                spades.py -1 trimmed/$read1 -2 trimmed/$read2 -o spades/$out
done

# rename contig files
files_rename=`ls spades/ | grep "_spades_out"`
for contig in $files_rename; do newname=`echo $contig | cut -d s -f1`spades_contigs.fa; mv spades/$contig/contigs.fasta spades/$contig/$newname; done
for scaffold in $files_rename; do newname=`echo $scaffold | cut -d s -f1`spades_scaffolds.fa; mv spades/$scaffold/scaffolds.fasta spades/$scaffold/$newname; done

# copy assemblies to B_cereus_spades directory
cd <path to genome assemblies>/B_cereus_de_novo_assembly
cp raw_data/spades/*spades_out/SRR* .

