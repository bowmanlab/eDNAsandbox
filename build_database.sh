#!/bin/bash

#python3 parse_midori.py

for GENE in CO1 CytB lrRNA srRNA;do
	clustalo --auto --force -i MIDORI2_UNIQ_NUC_GB251_CO1_RAW.select.fasta -o MIDORI2_UNIQ_NUC_GB251_CO1_RAW.select.align.fasta
done
