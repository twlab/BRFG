#!/bin/bash
#SBATCH --mem=8G
#SBATCH -n 1
#SBATCH -N 1

read BAM BED < <( sed -n ${SLURM_ARRAY_TASK_ID}p $1 )

eval $( spack load --sh bedtools2@2.30.0 )

echo "BAM: $BAM"
echo "BED: $BED"

bedtools bamtobed -i $BAM > $BED

echo "Done!"