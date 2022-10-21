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
import shutil
import numpy as np

genes = ['CO1', 'CytB', 'lrRNA', 'srRNA']

## Read in taxonomy

genome_taxa = pd.read_csv('MIDORI2_UNIQ_NUC_GB251_RAW_genome_taxa_select.csv', index_col=0)
        
## Get sequences from all alignment files

aligned_seqs = pd.DataFrame(columns = genes, index = genome_taxa.index)

for gene in genes:
    for record in SeqIO.parse('MIDORI2_UNIQ_NUC_GB251_' + gene + '_RAW.select.align.fasta', 'fasta'):
        record_genome = str(record.id).split('_')[-1] #clean this up later
        aligned_seqs.loc[record_genome, gene] = str(record.seq)
        
## Assign appropriate subtrees
        
genome_taxa['subtree'] = genome_taxa['phylum']
genome_taxa['subtree_rank'] = 'phylum'
        
## Iterate across rows that don't have a subtree.  These either had no
## assigned phylum or the phylum was too large.
        
for index, row in genome_taxa.loc[pd.isnull(genome_taxa.subtree)].iterrows():
    if pd.isnull(row['class']) == False:
        genome_taxa.loc[index, 'subtree'] = row['class']
        genome_taxa.loc[index, 'subtree_rank'] = 'class'
    elif pd.isnull(row['order']) == False:
        genome_taxa.loc[index, 'subtree'] = row['order']
        genome_taxa.loc[index, 'subtree_rank'] = 'order'
    elif pd.isnull(row['family']) == False:
        genome_taxa.loc[index, 'subtree'] = row['family']
        genome_taxa.loc[index, 'subtree_rank'] = 'family'
    elif pd.isnull(row['genus']) == False:
        genome_taxa.loc[index, 'subtree'] = row['genus']
        genome_taxa.loc[index, 'subtree_rank'] = 'genus'
        
## Now iteratively reduce large subtrees until largest is 1000.  The best
## way to do this is probably not the most efficient.  You want to divide each
## subtree only once per iteration, so as to not unnecessarily impose a lower
## rank than necessary.

while genome_taxa.subtree.value_counts().max() > 1000 or genome_taxa[pd.isnull(genome_taxa.subtree)].shape[0] > 0:
    for subtree in genome_taxa.subtree.unique():
        print(subtree, genome_taxa[genome_taxa.subtree == subtree].shape[0])
        
        if genome_taxa[genome_taxa.subtree == subtree].shape[0] > 1000:
            
            ## Get next rank to right
            
            current_rank = genome_taxa.loc[genome_taxa.subtree == subtree, 'subtree_rank'][0]
            new_rank = genome_taxa.columns[list(genome_taxa.columns).index(current_rank) + 1]
            genome_taxa.loc[genome_taxa.subtree == subtree, 'subtree_rank'] = new_rank
            genome_taxa.loc[genome_taxa.subtree == subtree, 'subtree'] = genome_taxa.loc[genome_taxa.subtree == subtree, new_rank]
            
        elif pd.isnull(subtree) == True:
            
            ## Get next rank to right
            
            current_rank = genome_taxa.loc[pd.isnull(genome_taxa.subtree), 'subtree_rank'][0]
            new_rank = genome_taxa.columns[list(genome_taxa.columns).index(current_rank) + 1]
            genome_taxa.loc[pd.isnull(genome_taxa.subtree), 'subtree_rank'] = new_rank
            genome_taxa.loc[pd.isnull(genome_taxa.subtree), 'subtree'] = genome_taxa.loc[pd.isnull(genome_taxa.subtree), new_rank]
            
            #!!! This works, but is still suboptimal as it creates a lot of
            ## lower taxonomic ranks with just 1-2 members.  Maybe not an issue?
                  
## Write phylum level trees and select phylum reps. Right now arthopoda and
## chordota have a large number of seqs and probably need to be subdivided.
## Taking care of this would also take care of entries with no phylum, as we
## could create a new column in the taxonomy dataframe of what subtree each
## member should belong to.
        
## Right now only 100 of each phyla being printed to make a faster test version
## of the database.
        
## RAxML needs >= 4 taxa to build a tree.  If you have few than that don't
## bother, and rely on placement in the guide tree.
        
subtree_reps = []

for subtree in genome_taxa.subtree.unique():
    print(subtree, genome_taxa[genome_taxa.subtree == subtree].shape[0])
    if genome_taxa[genome_taxa.subtree == subtree].shape[0] >= 4:
        shutil.rmtree(subtree, ignore_errors = True)
        os.mkdir(subtree)
        with open(subtree + '/MIDORI2_UNIQ_NUC_GB251_CONCAT.' + subtree + '.select.align.fasta', 'w') as fasta_out:
            
            ## Select for guide tree.  If few than 20 reps in that phylum, put all
            ## on guide tree.
            
            try:
                temp = list(genome_taxa[genome_taxa.subtree == subtree].sample(10).index)
                subtree_reps = subtree_reps + temp
            except ValueError:
                temp = list(genome_taxa[genome_taxa.subtree == subtree].index)
                subtree_reps = subtree_reps + temp
                
            ## Outer try clause is just for producing small select alignments for testing.  To
            ## turn this off use arbitrarily large number of taxa.
                
            try:
                for index, row in aligned_seqs.reindex(genome_taxa.loc[genome_taxa['subtree'] == subtree].index).sample(9999).iterrows():
                    seq = row.str.cat()
                    try:
                        record = SeqRecord(Seq(seq), id = index, description = '')
                    except TypeError:
                        record = SeqRecord(Seq(seq), id = index, description = '')
                    SeqIO.write(record, fasta_out, 'fasta') 
            except ValueError:
                for index, row in aligned_seqs.reindex(genome_taxa.loc[genome_taxa['subtree'] == subtree].index).iterrows():
                    seq = row.str.cat()
                    try:
                        record = SeqRecord(Seq(seq), id = index, description = '')
                    except TypeError:
                        record = SeqRecord(Seq(seq), id = index, description = '')                               
                    SeqIO.write(record, fasta_out, 'fasta') 
                    
    ## If fewer than 4 in the subtree, don't create a tree.
                    
    else:
        temp = list(genome_taxa[genome_taxa.subtree == subtree].index)
        subtree_reps = subtree_reps + temp

## Write out phylum reps for guide tree

shutil.rmtree('guide', ignore_errors = True)                
os.mkdir('guide')
        
with open('guide/MIDORI2_UNIQ_NUC_GB251_CONCAT.guide.select.align.fasta', 'w') as fasta_out:
    for index, row in aligned_seqs.loc[subtree_reps].iterrows():
        seq = row.str.cat()
        record = SeqRecord(Seq(seq), id = index + '|' + genome_taxa.loc[index, 'subtree'], description = '')
        SeqIO.write(record, fasta_out, 'fasta')   
        
genome_taxa.to_csv('MIDORI2_UNIQ_NUC_GB251_RAW_genome_taxa_select.csv')

