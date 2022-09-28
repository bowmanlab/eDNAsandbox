# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 12:40:51 2022

@author: jeff
"""


## Data downloaded manually from http://www.reference-midori.info/download.php#
## Used RAW -> uniq data. Probably actually want RAW_sp -> uniq

## Scripted download possible from direct links like:
## wget http://www.reference-midori.info/download/Databases/GenBank250/SPINGO_sp/uniq/MIDORI2_UNIQ_SP_NUC_GB250_ND4L_SPINGO.fasta.gz

from Bio import SeqIO
import gzip
import pandas as pd

genes = ['A6', 'CO1', 'CO2', 'CO3', 'Cytb', 'lrRNA', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6', 'srRNA']

## Gather species for the purpose of initiating a dataframe quickly

taxa = []

for gene in genes:
    with gzip.open('MIDORI2_UNIQ_NUC_GB251_' + gene + '_RAW.fasta.gz', 'rt') as gzfile:
        for record in SeqIO.parse(gzfile, 'fasta'):
            record_genome = str(record.id).split('.')[0]
            taxa.append(record_genome)
            
taxa = set(taxa)

## Construct a dataframe of genes by genome
            
genes_present = pd.DataFrame(columns = genes, index = taxa)

for gene in genes:
    with gzip.open('MIDORI2_UNIQ_NUC_GB251_' + gene + '_RAW.fasta.gz', 'rt') as gzfile:
        for record in SeqIO.parse(gzfile, 'fasta'):
            record_genome = str(record.id).split('.')[0]
            genes_present.loc[record_genome, gene] = 1
            print(gene, record_genome)
            
## Which genomes have CO1, CytB, lrRNA, srRNA
            
genes_present_select = genes_present[~genes_present['CO1'].isnull()]
genes_present_select = genes_present_select[~genes_present_select['Cytb'].isnull()]
genes_present_select = genes_present_select[~genes_present_select['lrRNA'].isnull()]
genes_present_select = genes_present_select[~genes_present_select['srRNA'].isnull()]

## If you want a random subset for testing, do that here.

genes_present_select = genes_present_select.sample(1000)

## Iterate across input files again and identify those with all target genes.
            
target_genomes = set(genes_present_select.index)

for gene in ['CO1', 'CytB', 'lrRNA', 'srRNA']:
    found_records = []
    with gzip.open('MIDORI2_UNIQ_NUC_GB251_' + gene + '_RAW.fasta.gz', 'rt') as gzfile:
        for record in SeqIO.parse(gzfile, 'fasta'):
            record_genome = str(record.id).split('.')[0]
            if record_genome in target_genomes:
                found_records.append(record)
                print(gene, record_genome)
    with open('MIDORI2_UNIQ_NUC_GB251_' + gene + '_RAW.select.fasta', 'w') as fasta_out:
        SeqIO.write(found_records, fasta_out, 'fasta')
                