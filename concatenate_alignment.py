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

#genes = ['CO1', 'CytB', 'lrRNA', 'srRNA']
genes = ['CO1', 'CytB']

## Read in taxonomy

genome_taxa = pd.read_csv('MIDORI2_UNIQ_NUC_GB251_RAW_genome_taxa_select.csv', index_col=0)
        
## Get sequences from all alignmen files

aligned_seqs = pd.DataFrame(columns = genes, index = genome_taxa.index)

for gene in genes:
    for record in SeqIO.parse('MIDORI2_UNIQ_NUC_GB251_' + gene + '_RAW.select.align.fasta', 'fasta'):
        record_genome = str(record.id).split('_')[-1] #clean this up later
        aligned_seqs.loc[record_genome, gene] = str(record.seq)
        
## Write phylum level trees and select phylum reps. Right now arthopoda and
## chordota have a large number of seqs and probably need to be subdivided.
## Taking care of this would also take care of entries with no phylum, as we
## could create a new column in the taxonomy dataframe of what subtree each
## member should belong to.
        
phylum_reps = []

for phylum in genome_taxa.phylum.unique():
    with open('MIDORI2_UNIQ_NUC_GB251_CONCAT.' + phylum + '.select.align.fasta', 'w') as fasta_out:
        try:
            temp = list(genome_taxa[genome_taxa.phylum == phylum].sample(20).index)
            phylum_reps = phylum_reps + temp
        except ValueError:
            temp = list(genome_taxa[genome_taxa.phylum == phylum].index)
            phylum_reps = phylum_reps + temp
        for index, row in aligned_seqs.reindex(genome_taxa.loc[genome_taxa['phylum'] == phylum].index).iterrows():
            seq = row.str.cat()
            try:
                record = SeqRecord(Seq(seq), id = index, description = genome_taxa.loc[index, 'class'])
            except TypeError:
                record = SeqRecord(Seq(seq), id = index, description = 'no_class')
                    
            SeqIO.write(record, fasta_out, 'fasta') 

## Write out phylum reps for guide tree
        
with open('MIDORI2_UNIQ_NUC_GB251_CONCAT.guide.select.align.fasta', 'w') as fasta_out:
    for index, row in aligned_seqs.loc[phylum_reps].iterrows():
        seq = row.str.cat()
        record = SeqRecord(Seq(seq), id = index, description = genome_taxa.loc[index, 'phylum'])
        SeqIO.write(record, fasta_out, 'fasta')        

