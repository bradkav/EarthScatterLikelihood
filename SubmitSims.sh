#!/bin/bash

#SBATCH -N 1 --ntasks-per-node=16  
#SBATCH -t 12:00:00
#SBATCH -p normal

# #SBATCH -t 00:04:30
# #SBATCH -p short

#SBATCH -o slurm_output/slurm-%j.out # STDOUT                                                                                                          #SBATCH -e slurm_output/slurm-%j.err # STDERR


cd $HOME/EarthScatterLikelihood

#module load openmpi/gnu
#module load python/2.7.9

module load pre2019

module unload GCCcore
module load Python/2.7.12-intel-2016b
module load slurm-tools

export SLURM_CPU_BIND=none

#ID="Sapphire_Test"
#ID="Final_10eV_res3eV_Sapphire_clip2_refine_long"
ID="Final_60eV_res18eV_EDE_clip2_refine_long"
MASS="0.400"

EXE="calcContour_EDE_long"

#time mpirun -np 16 python2.7 RunMPI_Contours.py -m_x $MASS -data 1 -out_dir ${ID}1 -hemisphere N -exe $EXE
#time mpirun -np 16 python2.7 RunMPI_Contours.py -m_x $MASS -data 1 -out_dir ${ID}1 -hemisphere S -exe $EXE
#time mpirun -np 16 python2.7 RunMPI_Contours.py -m_x $MASS -data 3 -out_dir ${ID}3 -hemisphere N -exe $EXE
time mpirun -np 16 python2.7 RunMPI_Contours.py -m_x $MASS -data 3 -out_dir ${ID}3 -hemisphere S -exe $EXE
