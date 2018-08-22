#!/bin/bash
#<add job submission details>

## navigate to directory
cd <path to genome assemblies>

## load Roary (Page et al. 2015)
module purge
module load python/anaconda3 # program set to Quest server conda files

### RUN ROARY ###


## Annotated B. cereus GenBank assemblies ##

cd B_cereus/prokka_out/roary_inp

# retain only high-quality genomes for the pan-genome analysis
# i.e., those with >97% completeness and <3% contamination, according to CheckM (see Supp. Table 1)
# manually remove genomes that did not meet these criteria from the roary input folder before proceeding

gunzip *gff.gz
mkdir roary_out

roary -e -p 8 *.gff -i 90 -cd 95 -f roary_out


## Annotated S. aureus GenBank assemblies ##

cd S_aureus/prokka_out/roary_inp

# retain only high-quality genomes for the pan-genome analysis
# i.e., those with >97% completeness and <3% contamination, according to CheckM (see Supp. Table 1)
# manually remove genomes that did not meet these criteria from the roary input folder before proceeding

gunzip *gff.gz
mkdir roary_out

roary -e -p 8 *.gff -i 90 -cd 95 -f roary_out


## Annotated B. cereus genomes assembled de novo with Spades ##

cd B_cereus_de_novo_assembly/prokka_out/roary_inp

# retain only high-quality genomes for the pan-genome analysis
# i.e., those with >97% completeness and <3% contamination, according to CheckM (see Supp. Table 1)
# manually remove genomes that did not meet these criteria from the roary input folder before proceeding

gunzip *gff.gz
mkdir roary_out

roary -e -p 8 *.gff -i 90 -cd 95 -f roary_out

