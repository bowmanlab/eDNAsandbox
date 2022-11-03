# -*- coding: utf-8 -*-
"""
Created on Fri Oct  7 12:57:02 2022

@author: jeff
"""

import pandas as pd
import subprocess

genome_taxa = pd.read_csv('MIDORI2_UNIQ_NUC_GB251_RAW_genome_taxa_select.csv', index_col=0)

#!!! it would be advantageous to run alignment check, then use the reduced.phy
#file for building the tree and phylogenetic placement.

#!!! actually this wouldn't work, because this would require aligning the query
## against each reduced alignment

def run_raxml(prefix, msa, subtree):
    
    ## Create trees with parsimony and random starts.
    
    subprocess.run('cd ' + subtree + ';raxml-ng \
                --redo \
                --model GTR+I+G4 \
                --search \
                --msa ' + msa + ' \
                --prefix ' + prefix + ' \
                --tree pars{9},rand{9}', shell = True, executable = '/bin/bash')
                
def run_raxml_light(prefix, msa, subtree):
    
    ## Create trees with parsimony and random starts.
    
    subprocess.run('cd ' + subtree + ';raxml-ng \
                --redo \
                --model GTR+I+G4 \
                --search \
                --msa ' + msa + ' \
                --prefix ' + prefix + ' \
                --tree pars{2},rand{2}', shell = True, executable = '/bin/bash')
                
for subtree in genome_taxa.subtree.unique():
    if subtree != 'guide':
        msa = 'MIDORI2_UNIQ_NUC_GB251_CONCAT.' + subtree + '.select.align.fasta'
        
        if genome_taxa.loc[genome_taxa.subtree == subtree].shape[0] < 400:
            run_raxml('MIDORI2_UNIQ_NUC_GB251_CONCAT.' + subtree + '.select', msa, subtree)
        else:
            run_raxml_light('MIDORI2_UNIQ_NUC_GB251_CONCAT.' + subtree + '.select', msa, subtree)
    
## Guide tree
    
msa = 'MIDORI2_UNIQ_NUC_GB251_CONCAT.guide.select.align.fasta'
run_raxml_light('MIDORI2_UNIQ_NUC_GB251_CONCAT.guide.select', msa, 'guide')   

## Write out genome_taxa in the taxonomy format that gappa wants.  This should
## probably be in concatenate_alignment.py

with open('MIDORI2_UNIQ_NUC_GB251_RAW_genome_taxa_select.taxonomy', 'w') as tax_out:
    for index, row in genome_taxa.iterrows():
        row[pd.isnull(row)] = 'NaN'
        temp = ';'.join(row[:-1])
        print(index + '\t' + temp, file = tax_out)
