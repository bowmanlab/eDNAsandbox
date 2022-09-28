#!/bin/bash

#python3 parse_midori.py

#for GENE in CO1 CytB lrRNA srRNA;do
#	clustalo --auto --force -i MIDORI2_UNIQ_NUC_GB251_${GENE}_RAW.select.fasta -o MIDORI2_UNIQ_NUC_GB251_${GENE}_RAW.select.align.fasta
#done

#python3 concatenate_alignment.py

#modeltest-ng -i MIDORI2_UNIQ_NUC_GB251_CONCAT.select.align.fasta -d nt -p 8

## winning model is GTR+I+G4

#raxml-ng --check --msa MIDORI2_UNIQ_NUC_GB251_CONCAT.select.align.fasta --model GTR+I+G4

raxml-ng --all --msa MIDORI2_UNIQ_NUC_GB251_CONCAT.select.align.fasta --model GTR+I+G4