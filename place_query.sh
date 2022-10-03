#!/bin/bash

REFERENCE=MIDORI2_UNIQ_NUC_GB251_CONCAT.select
QUERY=

clustalo -i $query.fasta --profile1 $REFERENCE.align.fasta -o $query.align.fasta

epa-ng --redo -t $REFERENCE.align.fasta.raxml.bestTree -s $REFERENCE.align.fasta -q $query.align.fasta -m $REFERENCE.select.align.fasta.raxml.bestModel