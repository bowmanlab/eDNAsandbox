# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 12:40:51 2022

@author: jeff
"""


## Data downloaded manually from http://www.reference-midori.info/download.php#
## Used RAW -> uniq data. Probably actually want RAW_sp -> uniq

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
            
genes_present.sum(axis = 1)