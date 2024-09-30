#!/bin/bash
#SBATCH --mem=16G
#SBATCH -n 1
#SBATCH -N 1

read INPUT OUTPUT < <( sed -n ${SLURM_ARRAY_TASK_ID}p $1 )

echo "Loading software..."
eval $( spack load --sh python@3.8.12 )
eval $( spack load --sh py-pysam@0.21.0 )

# print parameters to log
echo "INPUT: $INPUT"
echo "OUTPUT: $OUTPUT"

echo "Running Calculate Cigar Composition..."
python3 /scratch/hllab/Juan/General_Code/Read_Handling/Calculate_Cigar_Composition_Bed.py $INPUT $OUTPUT
echo "Done!"