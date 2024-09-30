#!/bin/bash
#SBATCH --mem=16G
#SBATCH --ntasks=1                   # Number of tasks (1 task)
#SBATCH --cpus-per-task=8            # Number of CPU cores per task (8 threads)

read INPUT OUTPUT < <( sed -n ${SLURM_ARRAY_TASK_ID}p $1 )

echo "Loading software..."
eval $( spack load --sh samtools@1.13 )

# print parameters to log
echo "INPUT: $INPUT"
echo "OUTPUT: $OUTPUT"

echo "Samtools sort..."
samtools sort --threads 8 $INPUT -o $OUTPUT -O BAM
echo "Samtools index..."
samtools index $OUTPUT
echo "Done!"