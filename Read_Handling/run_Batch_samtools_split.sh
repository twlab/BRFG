#!/bin/bash
#SBATCH --mem=16G
#SBATCH --ntasks=1                   # Number of tasks (1 task)
#SBATCH --cpus-per-task=4            # Number of CPU cores per task (8 threads)

# this script will take in a list of input files
# each is a bam file to be split by the chromosome field
# the output will be a set of bam files, one for each chromosome in the input file

read INPUT < <( sed -n ${SLURM_ARRAY_TASK_ID}p $1 )

echo "Loading software..."
eval $( spack load --sh samtools@1.13 )

# print parameters to log
echo "INPUT: $INPUT"

echo "Samtools split..."
base_name=$(basename "$INPUT" .bam)
for chr in $(samtools idxstats "$INPUT" | cut -f1 | grep -v '*'); do
  samtools view --threads 4 -b "$INPUT" "$chr" -o "${base_name}__${chr}.bam"
done
echo "Done!"
