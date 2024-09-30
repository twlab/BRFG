#!/bin/bash
#SBATCH --mem=16G
#SBATCH -n 1
#SBATCH -N 1

read INPUT OUTPUT < <( sed -n ${SLURM_ARRAY_TASK_ID}p $1 )

echo "Loading software..."
eval $( spack load --sh samtools@1.13 )

# print parameters to log
echo "INPUT: $INPUT"
echo "OUTPUT: $OUTPUT"

echo "Samtools stats..."
samtools stats $INPUT >$OUTPUT
echo "Done!"