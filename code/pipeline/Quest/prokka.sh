#!/bin/bash
#<add job submission details>

## navigate to directory and unzip assembly files
cd <path to genome assemblies>
mkdir checkm
gunzip */*.gz

## load Prokka (Seeman 2014)
module purge
module load perl/5.26 # program set to my path in Quest server; requires perl


### RUN PROKKA ###


## B. cereus GenBank assemblies ##

cd <path to genome assemblies>/B_cereus

# create file names for annotating files
files_annotate=`ls | grep "fna"`

# use loop to run prokka on all files
for fna in $files_annotate; do annotate=`echo $fna | cut -d "." -f1`_prok; prokka --kingdom Bacteria --outdir prokka_out/$annotate --genus Bacillus --locustag $annotate --prefix $annotate $fna; done

# zip assembly and output files
gzip *.fna
gzip prokka_out/*/*

# copy the output annotated files (.gff) to roary input directory
mkdir prokka_out/roary_inp
cd prokka_out/roary_inp
cp <path to genome assemblies>/B_cereus/prokka_out/*/*.gff.gz .


## S. aureus GenBank assemblies ##

cd <path to genome assemblies>/S_aureus

# create file names for annotating files
files_annotate=`ls | grep "fna"`

# use loop to run prokka on all files
for fna in $files_annotate; do annotate=`echo $fna | cut -d "." -f1`_prok; prokka --kingdom Bacteria --outdir prokka_out/$annotate --genus Staphylococcus --locustag $annotate --prefix $annotate $fna; done

# zip assembly and output files
gzip *.fna
gzip prokka_out/*/*

# copy the output annotated files (.gff) to roary input directory
mkdir prokka_out/roary_inp
cd prokka_out/roary_inp
cp <path to genome assemblies>/S_aureus/prokka_out/*/*.gff.gz .


## B. cereus genomes assembled de novo with Spades ##

cd <path to genome assemblies>/B_cereus_de_novo_assembly
gunzip *.gz

# create file names for annotating files
files_annotate=`ls | grep "fa"`

# use loop to run prokka on all files
for fa in $files_annotate; do annotate=`echo $fa | cut -d "." -f1`_prok; prokka --kingdom Bacteria --outdir prokka_out/$annotate --genus Bacillus --locustag $annotate --prefix $annotate $fa; done

# zip assembly and output files
gzip *.fa
gzip prokka_out/*/*

# copy the output annotated files (.gff) to roary input directory
mkdir prokka_out/roary_inp
cd prokka_out/roary_inp
cp <path to genome assemblies>/B_cereus_de_novo_assembly/prokka_out/*/*.gff.gz .

