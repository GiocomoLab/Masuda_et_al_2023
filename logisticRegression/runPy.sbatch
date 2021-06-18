#!/bin/bash
#
#set the job name (output file and error file)
#SBATCH --job-name=positionDecoder_Ket
#SBATCH --output=run.out
#SBATCH --error=run.error
#################
#SBATCH --mail-type=ALL
#################
#set time needed
#SBATCH --time=10:00:00
#################
#SBATCH -p giocomo,owners,normal
#################
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3
#################
#request 3GB RAM on CPU
#SBATCH --mem=12G


module load python/3.6.1
module load py-pandas/0.23.0_py36
module load py-numpy/1.14.3_py36
module load py-scipy/1.1.0_py36
module load py-scikit-image/0.15.0_py36
module load py-scikit-learn/0.19.1_py36
module load py-matplotlib/2.2.2_py36

python3 positionDecoderRMSEplot.py
