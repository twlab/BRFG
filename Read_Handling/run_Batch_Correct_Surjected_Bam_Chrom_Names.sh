#!/bin/bash
#SBATCH --mem=16G
#SBATCH -n 1
#SBATCH -N 1

# Read the BAM file from the input list
read BAM_FILE < <(sed -n ${SLURM_ARRAY_TASK_ID}p $1)

# py-pysam@0.21.0
eval $( spack load --sh /usbkij7 )

BAM_ROOT=$(basename $BAM_FILE .bam)
OUTPUT_FIXED="${BAM_ROOT}_ChromNamefixed.bam"

echo "Processing BAM file: $BAM_FILE"
echo "Output will be saved to: $OUTPUT_FIXED"

# Run the Python script to fix chromosome names, specifying the prefix to remove
python3 /scratch/hllab/Juan/General_Code/Read_Handling/Correct_Surjected_Bam_Chrom_Names.py $BAM_FILE $OUTPUT_FIXED

echo "Complete!"
