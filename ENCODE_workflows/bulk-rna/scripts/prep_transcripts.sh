#!/bin/bash -e

set -x

ref_fasta_gz=$1
annotation_gz=$2

gunzip -c ${ref_fasta_gz} > ref.fasta
gunzip -c ${annotation_gz} > anno.gtf
gffread -g ref.fasta -W -w transcripts.fasta anno.gtf
