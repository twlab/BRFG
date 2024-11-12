#!/bin/bash
#SBATCH --mem=32G
#SBATCH -n 1
#SBATCH -N 1


read LIFTOVERBED ORIGNALFILE OUTPUTFILE < <( sed -n ${SLURM_ARRAY_TASK_ID}p $1 )

echo "Loading software..."

eval $( spack load --sh python@3 )
eval $( spack load --sh /i7mcgz4 )

echo "Start conversion..."

python3 -u /scratch/hllab/Juan/General_Code/IncorporateLiftedOverPositions_BedGraph_PAFtools.py $LIFTOVERBED $ORIGNALFILE $OUTPUTFILE

echo "Complete!"



