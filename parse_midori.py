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
import numpy as np

#genes = ['A6', 'CO1', 'CO2', 'CO3', 'Cytb', 'lrRNA', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6', 'srRNA']
genes = ['CO1', 'Cytb', 'lrRNA', 'srRNA']
ranks = ['norank','superkingdom','kingdom','phylum','subphylum','class','order','family','genus','species']

## Gather species for the purpose of initiating a dataframe quickly

taxa = []

for gene in genes:
    with gzip.open('MIDORI2_UNIQ_NUC_GB251_' + gene + '_RAW.fasta.gz', 'rt') as gzfile:
        for record in SeqIO.parse(gzfile, 'fasta'):
            record_genome = str(record.id).split('.')[0]
            taxa.append(record_genome)
            
taxa = set(taxa)

## Map taxonomy to genome

genome_taxa = pd.DataFrame(index = taxa, columns = ranks)

for gene in genes:
    with gzip.open('MIDORI2_UNIQ_NUC_GB251_' + gene + '_RAW.fasta.gz', 'rt') as gzfile:
        for record in SeqIO.parse(gzfile, 'fasta'):
            record_genome = str(record.id).split('.')[0]
            record_taxonomy = str(record.description).split(';')
            
            try:
                np.isnan(genome_taxa.loc[record_genome, 'norank'])
                for rank in record_taxonomy:
                    level = rank.split('_')[0]
                    if level in ranks:
                        taxonomy = '_'.join(rank.split('_')[1:])
                        genome_taxa.loc[record_genome, level] = taxonomy
                print(gene, record_genome)
            except TypeError:
                continue

## Construct a dataframe of genes by genome
            
genes_present = pd.DataFrame(columns = genes, index = taxa)

for gene in genes:
    with gzip.open('MIDORI2_UNIQ_NUC_GB251_' + gene + '_RAW.fasta.gz', 'rt') as gzfile:
        for record in SeqIO.parse(gzfile, 'fasta'):
            record_genome = str(record.id).split('.')[0]
            genes_present.loc[record_genome, gene] = str(record.seq)
            print(gene, record_genome)
            
## Which genomes have CO1, CytB, lrRNA, srRNA
            
genes_present_select = genes_present[~genes_present['CO1'].isnull()]
genes_present_select = genes_present_select[~genes_present_select['Cytb'].isnull()]
genes_present_select = genes_present_select[~genes_present_select['lrRNA'].isnull()]
genes_present_select = genes_present_select[~genes_present_select['srRNA'].isnull()]

genome_taxa_select = genome_taxa.reindex(genes_present_select.index)

## It doesn't make sense to extract phylum refs and divide into guide and
## phylum fastas until after alignment.  So here, make a master fasta file for
## alignment.

## Some entries don't have phylum.  Going to cause problems later, here we'll
## just treat them with a separate analysis.

genome_taxa.loc[genome_taxa['phylum'].isnull(), 'phylum'] = 'no_phylum'
            
for gene in genes_present_select.columns:
    with open('MIDORI2_UNIQ_NUC_GB251_' + gene + '_RAW.select.fasta', 'w') as fasta_out:
        for index, value in genes_present_select.iterrows():
            phylum = genome_taxa.loc[index, 'phylum']
            print('>' + index + '\n' + value[gene], file = fasta_out)
            
genome_taxa_select.to_csv('MIDORI2_UNIQ_NUC_GB251_RAW_genome_taxa_select.csv')
genes_present_select.to_csv('MIDORI2_UNIQ_NUC_GB251_RAW_genes_present_select.csv')


## YOU ARE HERE, CHECK CSVS                