# -*- coding: utf-8 -*-
"""
Created on Sat Oct  8 07:19:49 2022

@author: jeff
"""

import subprocess
import pandas as pd
from Bio import SeqIO
import os

genome_taxa = pd.read_csv('MIDORI2_UNIQ_NUC_GB251_RAW_genome_taxa_select.csv', index_col=0)
taxonomy = 'MIDORI2_UNIQ_NUC_GB251_RAW_genome_taxa_select.taxonomy'
cwd = os.getcwd() + '/' # The current working directory.

#%% Define functions

def make_unique(query):
    
    ## Check for duplicate names.
    
    dup_names = False
    old_names = []
    
    for record in SeqIO.parse(query + '.fasta', 'fasta'):
        old_names.append(str(record.id))
    
    if len(old_names) != len(set(old_names)):
        dup_names = True

    seq_count = {}
    seq_names = {}
    name_seq = {}
    i = 0
        
    for record in SeqIO.parse(query + '.fasta', 'fasta'):
        i += 1
        name = str(record.id)
        if dup_names == True:
            name = str(i) + '|' + name
        seq = str(record.seq).upper()
        
        if seq in list(seq_count.keys()):
            seq_count[seq] = seq_count[seq] + 1
            temp_name_list = seq_names[seq]
            temp_name_list.append(name)
            seq_names[seq] = temp_name_list
        else:
            seq_count[seq] = 1
            seq_names[seq] = [name]
            
    with open(query + '.unique.fasta', 'w') as fasta_out, open(query + '.unique.count', 'w') as count_out:
        print('rep_name' + ',' + 'abundance', file = count_out)
        for seq in list(seq_names.keys()):
            print('>' + seq_names[seq][0], file = fasta_out)
            print(seq, file = fasta_out)
            
            print(seq_names[seq][0] + ',' + str(seq_count[seq]), file = count_out)
            name_seq[seq_names[seq][0]] = seq
            
    return seq_count, seq_names, name_seq

## Split alignment. For reasons that aren't clear, clustalo changes the alignment
## length by 1. So the new reference alignment is created, but is discarded after
## use by epa-ng.   
                
def split_alignment(combined_in, query_out, ref_out):
    with open(query_out, 'w') as fasta_out, open(ref_out, 'w') as ref_fasta_out:                       
        for record in SeqIO.parse(combined_in, 'fasta'):
            record_name = str(record.id).split('|')[0]
            if record_name not in genome_taxa.index:
                SeqIO.write(record, fasta_out, 'fasta')
            else:
                SeqIO.write(record, ref_fasta_out, 'fasta')
    
#%% Run program

## Align query to guide alignment. This results in a combined ref + query alignment.
## Clustalo is very slow, so it will be important to create a unique sequence
## file first.

ref = 'MIDORI2_UNIQ_NUC_GB251_CONCAT.guide.select'
query = '2111SR_20211107_40_16_40m_900mL_S8_L001.12SV5.exp'
phylum = 'guide'

## Make unique.  If program hangs here you probably didn't denoise your reads
## adequately!

seq_count, seq_names, name_seq = make_unique(cwd + query)

subprocess.run('clustalo --force --auto \
                -i ' + query + '.unique.fasta \
                --profile1 ' + ref + '.align.fasta \
                -o ' + query + '.' + phylum + '.align.fasta', shell = True)
                       
## Split alignment. For reasons that aren't clear, clustalo changes the alignment
## length by 1. So the new reference alignment is created, but is discarded after
## use by epa-ng.

split_alignment(query + '.' + phylum + '.align.fasta',
                query + '.align.fasta',
                'temp_ref.align.fasta')
                       
## Place to guide tree and create taxonomy.  This can be parsed for
## phylum tree.  The downside to this approach over how paprica does it, is
## that it might be less flexible if the downstream tree taxonomic level isn't
## fixed.
                       
subprocess.run('epa-ng --redo \
               -t ' + ref + '.raxml.bestTree \
               -q ' + query + '.align.fasta \
               -s temp_ref.align.fasta \
               -m ' + ref + '.raxml.bestModel', shell = True)

subprocess.run('mv epa_result.jplace ' + query + '.' + phylum + '.jplace', shell = True) 
subprocess.run('mv epa_info.log ' + query + '.' + phylum + '.epa_info.log', shell = True)
os.remove('temp_ref.align.fasta')


subprocess.run('gappa examine assign \
               --allow-file-overwriting \
               --taxon-file ' + taxonomy + ' \
               --file-prefix ' + query + '.' + phylum + '. \
               --jplace-path ' + query + '.' + phylum + '.jplace \
               --per-query-results', shell = True)

## Parse taxonomy.  The current solution is kludgy as it iterates across
## unneeded lines.
                           
qphyla = pd.Series()
                           
with open(query + '.' + phylum + '.per_query.tsv', 'r') as tax_in:
    for line in tax_in.readlines():
        line = line.strip()
        line = line.split('\t')
        qread = line[0]
        qtax = line[5].split(';')
        try:
            qphylum = qtax[3]
            print(qphylum)
            qphyla[qread] = qphylum
        except IndexError:
            continue
        
## Split query alignment according to guide tree phyla.
            
for phylum in qphyla.unique():
    with open(query + '.' + phylum + '.align.fasta', 'w') as fasta_out:                       
        for record in SeqIO.parse(query + '.align.fasta', 'fasta'):
            if qphyla[record.id] == phylum:
                SeqIO.write(record, fasta_out, 'fasta')   
                
## Now do the placement on each phylum tree.  It appears to be necessary
## to realign the reads to each reference alignment.

for phylum in qphyla.unique():                
    ref = 'MIDORI2_UNIQ_NUC_GB251_CONCAT.' + phylum + '.select'
                    
    subprocess.run('clustalo --force --auto \
                    --profile1 ' + query + '.' + phylum + '.align.fasta \
                    --profile2 ' + ref + '.align.fasta \
                    -o ' + query + '.' + phylum + '.align.fasta.tmp', shell = True)
                    
    split_alignment(query + '.' + phylum + '.align.fasta.tmp',
                    query + '.align.fasta.tmp',
                    ref + '.align.fasta.tmp')
    
    subprocess.run('epa-ng --redo \
               -t ' + ref + '.raxml.bestTree \
               -q ' + query + '.align.fasta.tmp \
               -s ' + ref + '.align.fasta.tmp \
               -m ' + ref + '.raxml.bestModel', shell = True)
               
    
    subprocess.run('mv epa_result.jplace ' + query + '.' + phylum + '.jplace', shell = True) 
    subprocess.run('mv epa_info.log ' + query + '.' + phylum + '.epa_info.log', shell = True)          
    subprocess.run('rm *tmp', shell = True)
                    
#stop_here = []
#stop_here[1]

## Not done yet: Need to classify ASVs, aggregate the classifications, and
## combine these data with the count data.
 



