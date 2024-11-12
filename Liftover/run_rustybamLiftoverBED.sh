#!/bin/bash
#SBATCH --mem=400G
#SBATCH -n 1
#SBATCH -N 1

BEDFILE=$1
PAFFILE=$2
OUTPUT=$3

echo "Bed file to lift over: "$BEDFILE
echo "PAF file: "$PAFFILE
echo "Output: "$OUTPUT

/ref/hllab/software/rustybam/rustybam liftover --bed $BEDFILE $PAFFILE >$OUTPUT

echo "Complete!"