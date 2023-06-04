#!/bin/bash

#SBATCH -n 48
#SBATCH --time=2:00:00
#SBATCH --job-name=nbody_run
#SBATCH --output=slurm/slurm_output.txt
#SBATCH --error=slurm/slurm_error.txt
#SBATCH --open-mode=truncate
#SBATCH --constraint=EPYC_7742

module load python
module load ffmpeg

lscpu

export OMP_NUM_THREADS=48
cd build/
make

echo ""; echo "starting simulation at "; date; echo ""

./nbody

echo ""; echo "ending simulation at "; date; echo ""

cd ../scripts/
python3 animation.py
python3 3Danimation.py
python3 energyplot.py