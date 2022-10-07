#!/bin/bash

REFERENCE=MIDORI2_UNIQ_NUC_GB251_CONCAT.select
QUERY=2111SR_20211107_40_16_40m_900mL_S8_L001.12SV5.uni

clustalo --auto -i $QUERY.fasta --profile1 $REFERENCE.align.fasta -o $QUERY.align.fasta

epa-ng --redo -t $REFERENCE.align.fasta.raxml.bestTree -q $QUERY.align.fasta -s $REFERENCE.align.fasta -m $REFERENCE.align.fasta.raxml.bestModel

mv epa_result.jplace $QUERY.jplace 

gappa examine edpl --allow-file-overwriting --file-prefix $QUERY --jplace-path $QUERY.jplace

gappa examine heat-tree --allow-file-overwriting --write-phyloxml-tree --file-prefix $QUERY --jplace-path $QUERY.jplace