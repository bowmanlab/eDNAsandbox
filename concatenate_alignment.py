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
import os

genes = ['CO1', 'CytB', 'lrRNA', 'srRNA']

## Read in taxonomy

genome_taxa = pd.read_csv('MIDORI2_UNIQ_NUC_GB251_RAW_genome_taxa_select.csv', index_col=0)
        
## Get sequences from all alignment files

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
        
## Right now only 100 of each phyla being printed to make a faster test version
## of the database.
        
## RAxML needs >= 4 taxa to build a tree.  If you have few than that don't
## bother, and rely on placement in the guide tree.
        
phylum_reps = []

for phylum in genome_taxa.phylum.unique():
    if phylum != 'no_phylum':
        os.mkdir(phylum)
        with open(phylum + '/MIDORI2_UNIQ_NUC_GB251_CONCAT.' + phylum + '.select.align.fasta', 'w') as fasta_out:
            
            ## Select for guide tree.  If few than 20 reps in that phylum, put all
            ## on guide tree.
            
            try:
                temp = list(genome_taxa[genome_taxa.phylum == phylum].sample(20).index)
                phylum_reps = phylum_reps + temp
            except ValueError:
                temp = list(genome_taxa[genome_taxa.phylum == phylum].index)
                phylum_reps = phylum_reps + temp
                
            ## Outer try clause is just for producing small select alignments for testing.
                
            try:
                for index, row in aligned_seqs.reindex(genome_taxa.loc[genome_taxa['phylum'] == phylum].index).sample(100).iterrows():
                    seq = row.str.cat()
                    try:
                        record = SeqRecord(Seq(seq), id = index, description = '')
                    except TypeError:
                        record = SeqRecord(Seq(seq), id = index, description = '')
                    SeqIO.write(record, fasta_out, 'fasta') 
            except ValueError:
                for index, row in aligned_seqs.reindex(genome_taxa.loc[genome_taxa['phylum'] == phylum].index).iterrows():
                    seq = row.str.cat()
                    try:
                        record = SeqRecord(Seq(seq), id = index, description = '')
                    except TypeError:
                        record = SeqRecord(Seq(seq), id = index, description = '')                               
                    SeqIO.write(record, fasta_out, 'fasta') 
                    
    ## Don't create a reference directory for no_refs, all of these will be used
    ## on the guide tree.
                    
    else:
        temp = list(genome_taxa[genome_taxa.phylum == phylum].index)
        phylum_reps = phylum_reps + temp

## Write out phylum reps for guide tree
                
os.mkdir('guide')
        
with open('guide/MIDORI2_UNIQ_NUC_GB251_CONCAT.guide.select.align.fasta', 'w') as fasta_out:
    for index, row in aligned_seqs.loc[phylum_reps].iterrows():
        seq = row.str.cat()
        record = SeqRecord(Seq(seq), id = index + '|' + genome_taxa.loc[index, 'phylum'], description = '')
        SeqIO.write(record, fasta_out, 'fasta')        

