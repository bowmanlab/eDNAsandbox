# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 14:45:47 2022

@author: jeff
"""


## If you want a random subset for testing, do that here.

#genes_present_select = genes_present_select.sample(1000)

target_taxa = []

for phylum in genome_taxa_select.phylum.unique():
    try:
        temp = list(genome_taxa_select[genome_taxa_select.phylum == phylum].sample(20).index)
        target_taxa = target_taxa + temp
    except ValueError:
        temp = list(genome_taxa_select[genome_taxa_select.phylum == phylum].index)
        target_taxa = target_taxa + temp
        
temp = []
temp[1]


for gene in genes_present_select.columns:
    
    temp = genes_present_select.loc[target_taxa, gene]
    for index, value in temp.iterrow():
        phylum = genome_taxa_select.loc[index, 'Phylum']
        print(phylum + '_' + index + '\n' + value, fasta_out)