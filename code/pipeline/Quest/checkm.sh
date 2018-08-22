#!/bin/bash
#<add job submission details>

## navigate to directory and unzip assembly files
cd <path to genome assemblies>
gunzip */*.gz

## load CheckM (Parks et al. 2015)
module load checkm
mkdir checkm_files

## B. cereus GenBank assemblies
mkdir checkm_files/out_Bac
checkm taxonomy_wf genus Bacillus B_cereus/ checkm_files/out_Bac -t 18 -x fna -f checkm_files/checkm_B_cereus.txt

## S. aureus GenBank assemblies
mkdir checkm_files/out_Staph
checkm taxonomy_wf genus Staphylococcus S_aureus/ checkm_files/out_Staph -t 18 -x fna -f checkm_files/checkm_S_aureus.txt

## B. cereus genomes assembled de novo with Spades
mkdir checkm_files/out_B_cereus_de_novo_assembly
checkm taxonomy_wf genus Bacillus B_cereus_de_novo_assembly/ checkm_files/out_B_cereus_de_novo_assembly -t 18 -x fa -f checkm_files/checkm_B_cereus_de_novo_assembly.txt

## zip assembly files
gzip */*.fna

