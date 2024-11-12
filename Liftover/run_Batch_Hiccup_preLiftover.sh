#!/bin/bash
#SBATCH --mem=8G
#SBATCH -n 1
#SBATCH -N 1

#### the sbatch array parameters have to be passed in when invoking this script

# [Alignment parameters lookup file]

read BEDFILE PAFFILE OUTPUT < <( sed -n ${SLURM_ARRAY_TASK_ID}p $1 )


eval $( spack load --sh k8@0.2.4 )


echo "lift bed file: "$BEDFILE
echo "with paf file: "$PAFFILE
echo "store output in: "$OUTPUT

echo "Lines in input: "
wc -l $BEDFILE

k8 /ref/hllab/software/paftools/paftools.js liftover -q 5 -l 50000 -d 1 $PAFFILE $BEDFILE > $OUTPUT

echo "Complete!"

echo "Lines in lift over output: "
wc -l $OUTPUT

