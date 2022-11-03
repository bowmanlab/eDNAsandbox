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
            
    with open(query + '/' + query + '.unique.fasta', 'w') as fasta_out, open(query + '/' + query + '.unique.count', 'w') as count_out:
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
                
## Creat a file mapping every rank in taxonomy to subtree rank.  This 
## will allow you to look up the appropriate subtree regradless of level of
## taxonomic ranking.
                
#!!! Not currently used

taxonomy_subtree = {}

for name, row in genome_taxa.iterrows():
    subtreerank = row['subtree_rank']
    subtree = row[subtreerank]
    row = list(row)
    subtreei = row.index(subtree)
    for taxa in row[subtreei:-1]:
        taxonomy_subtree[taxa] = subtree

## Align query to guide alignment. This results in a combined ref + query alignment.
## Clustalo is very slow, so it will be important to create a unique sequence
## file first.
                
#!!! Why aren't you aligning once, to the full alignment?

subtree = 'guide'
ref = subtree + '/MIDORI2_UNIQ_NUC_GB251_CONCAT.' + subtree + '.select'
query = '2111SR_20211107_40_16_40m_900mL_S8_L001.12SV5.exp'

## Make a directory for results.

try:
    os.mkdir(query)
except FileExistsError:
    pass

## Make unique.  If program hangs here you probably didn't denoise your reads
## adequately!

seq_count, seq_names, name_seq = make_unique(query)

## From this point all work takes place in the query subdirectory.

query_path = query + '/' + query

subprocess.run('clustalo --force --auto \
                -i ' + query_path + '.unique.fasta \
                --profile1 ' + ref + '.align.fasta \
                -o ' + query_path + '.' + subtree + '.align.fasta', shell = True)
                       
## Split alignment. For reasons that aren't clear, clustalo changes the alignment
## length by 1. So the new reference alignment is created, but is discarded after
## use by epa-ng.  It may be possible to use epa-ng --split for this, but this
## works fine for now.

split_alignment(query_path + '.' + subtree + '.align.fasta',
                query_path + '.align.fasta',
                'temp_ref.align.fasta')
                       
## Place to guide tree and create taxonomy.  This can be parsed for
## phylum tree.  The downside to this approach over how paprica does it, is
## that it might be less flexible if the downstream tree taxonomic level isn't
## fixed.
                       
subprocess.run('epa-ng --redo \
               -t ' + ref + '.raxml.bestTree \
               -q ' + query_path + '.align.fasta \
               -s temp_ref.align.fasta \
               -m ' + ref + '.raxml.bestModel', shell = True)

subprocess.run('mv epa_result.jplace ' + query_path + '.' + subtree + '.jplace', shell = True) 
subprocess.run('mv epa_info.log ' + query_path + '.' + subtree + '.epa_info.log', shell = True)
os.remove('temp_ref.align.fasta')


subprocess.run('gappa examine assign \
               --allow-file-overwriting \
               --out-dir ' + query + ' \
               --taxon-file ' + taxonomy + ' \
               --file-prefix ' + query + '.' + subtree + '. \
               --jplace-path ' + query_path + '.' + subtree + '.jplace \
               --per-query-results \
               --best-hit', shell = True)

## Parse taxonomy.
               
#!!! need to rewrite parse-taxonomy
        
def parse_taxonomy(per_query_tsv, subtree_taxonomy):
    
    tax_in = pd.read_csv(per_query_tsv, index_col = 0, sep = '\t')
    for qread, row in tax_in.iterrows():
        qsubtree = row['taxopath'].split(';')[-1]
        if qsubtree == "NaN":
            qsubtree = 'guide'
        print(qread, qsubtree)
        subtree_taxonomy.loc[qread, 'taxonomy'] = row['taxopath']
        subtree_taxonomy.loc[qread, 'subtree'] = qsubtree
            
    return(subtree_taxonomy)

subtree_taxonomy = pd.DataFrame(columns = ['taxonomy', 'subtree'])
subtree_taxonomy = parse_taxonomy(query_path + '.guide.per_query.tsv', subtree_taxonomy)
        
## Split query alignment according to guide tree phyla.
            
for subtree in subtree_taxonomy.subtree.unique():
    if subtree != 'guide':
        with open(query_path + '.' + subtree + '.align.fasta', 'w') as fasta_out:                       
            for record in SeqIO.parse(query_path + '.align.fasta', 'fasta'):
                if subtree_taxonomy.loc[record.id, 'subtree'] == subtree:
                    SeqIO.write(record, fasta_out, 'fasta')   
                    
        ## Now do the placement on each phylum tree.  It appears to be necessary
        ## to realign the reads to each reference alignment.       
                    
        ref = subtree + '/MIDORI2_UNIQ_NUC_GB251_CONCAT.' + subtree + '.select' 
                        
        subprocess.run('clustalo --force --auto \
                        --profile1 ' + query_path + '.' + subtree + '.align.fasta \
                        --profile2 ' + ref + '.align.fasta \
                        -o ' + query_path + '.' + subtree + '.align.fasta.tmp', shell = True)
                        
        split_alignment(query_path + '.' + subtree + '.align.fasta.tmp',
                        query_path + '.align.fasta.tmp',
                        ref + '.align.fasta.tmp')
        
        subprocess.run('epa-ng --redo \
                   -t ' + ref + '.raxml.bestTree \
                   -q ' + query_path + '.align.fasta.tmp \
                   -s ' + ref + '.align.fasta.tmp \
                   -m ' + ref + '.raxml.bestModel', shell = True)
                       
        subprocess.run('mv epa_result.jplace ' + query_path + '.' + subtree + '.jplace', shell = True) 
        subprocess.run('mv epa_info.log ' + query_path + '.' + subtree + '.epa_info.log', shell = True)          
        subprocess.run('rm *tmp', shell = True)
        
        subprocess.run('gappa examine assign \
               --allow-file-overwriting \
               --out-dir ' + query + ' \
               --taxon-file ' + taxonomy + ' \
               --file-prefix ' + query + '.' + subtree + '. \
               --jplace-path ' + query_path + '.' + subtree + '.jplace \
               --per-query-results \
               --best-hit', shell = True)
               
        subprocess.run('gappa examine edpl \
               --allow-file-overwriting \
               --out-dir ' + query + ' \
               --file-prefix ' + query + '.' + subtree + '. \
               --jplace-path ' + query_path + '.' + subtree + '.jplace', shell = True)
               
        subtree_taxonomy = parse_taxonomy(query_path + '.' + subtree + '.per_query.tsv', subtree_taxonomy)

for key in name_seq.keys():               
    subtree_taxonomy.loc[key, 'asv'] = name_seq[key] 

subtree_taxonomy = subtree_taxonomy.set_index('asv')      
               
for asv in subtree_taxonomy.index():
    subtree_taxonomy.loc[asv, 'abundance'] = seq_count[asv]
    
subtree_taxonomy.to_csv(query_path + '.asv_taxcount.csv')
                        
#stop_here = []
#stop_here[1]

## Not done yet: Need to classify ASVs, aggregate the classifications, and
## combine these data with the count data.
 



