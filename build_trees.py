# -*- coding: utf-8 -*-
"""
Created on Fri Oct  7 12:57:02 2022

@author: jeff
"""

import pandas as pd
import subprocess
import numpy as np

genome_taxa = pd.read_csv('MIDORI2_UNIQ_NUC_GB251_RAW_genome_taxa_select.csv', index_col=0)

def run_raxml(prefix, msa):
    
    ## Create trees with parsimony and random starts.
    
    subprocess.run('raxml-ng \
                --model GTR+I+G4 \
                --redo \
                --search \
                --msa ' + msa + ' \
                --prefix ' + prefix + ' \
                --workers 18 \
                --tree pars{9},rand{9} \
                --threads 36', shell = True, executable = '/bin/bash')
                
for phylum in genome_taxa.phylum.unique():
    msa = 'MIDORI2_UNIQ_NUC_GB251_CONCAT.' + phylum + '.select.align.fasta'
    run_raxml('MIDORI2_UNIQ_NUC_GB251_CONCAT.' + phylum + '.select', msa)
    
## Guide tree
    
msa = 'MIDORI2_UNIQ_NUC_GB251_CONCAT.guide.select.align.fasta'
run_raxml('MIDORI2_UNIQ_NUC_GB251_CONCAT.guide.select', msa)   

## Write out genome_taxa in the taxonomy format that gappa wants. This needs
## to be changed if you get rid of the "|.." trailing values, which you
## probably should.

with open('MIDORI2_UNIQ_NUC_GB251_RAW_genome_taxa_select.taxonomy', 'w') as tax_out:
    for index, row in genome_taxa.iterrows():
        row[pd.isnull(row)] = 'NaN'
        temp = ';'.join(row)
        print(index + '|' + genome_taxa.loc[index, 'phylum'] + '\t' + temp, file = tax_out)
