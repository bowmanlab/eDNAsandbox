# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 07:22:44 2022

@author: jeff
"""

## This script does not currently capture partition information, but this could
## be done in the future.

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

genes = ['CO1', 'CytB', 'lrRNA', 'srRNA']

## Get taxa from one of the sequence files, also map taxa to full taxonomy

genome_taxa = pd.DataFrame()

for record in SeqIO.parse('MIDORI2_UNIQ_NUC_GB251_CO1_RAW.select.align.fasta', 'fasta'):
    record_genome = str(record.id).split('.')[0]
    record_taxonomy = str(record.description).split(';')
    
    for rank in record_taxonomy:
        level = rank.split('_')[0]
        taxonomy = '_'.join(rank.split('_')[1:])
        genome_taxa.loc[record_genome, level] = taxonomy
        
## Get sequences from all files

all_seqs = pd.DataFrame(columns = genes, index = genome_taxa.index)

for gene in genes:
    for record in SeqIO.parse('MIDORI2_UNIQ_NUC_GB251_' + gene + '_RAW.select.align.fasta', 'fasta'):
        record_genome = str(record.id).split('.')[0]
        all_seqs.loc[record_genome, gene] = str(record.seq)
        
## Write out concatenated alignment
        
with open('MIDORI2_UNIQ_NUC_GB251_CONCAT.select.align.fasta', 'w') as fasta_out:
    for index, row in all_seqs.iterrows():
        seq = row.str.cat()
        record = SeqRecord(Seq(seq), id = index, description = '')
        SeqIO.write(record, fasta_out, 'fasta')

