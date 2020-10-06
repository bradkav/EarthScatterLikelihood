#!/bin/bash

#SBATCH -N 1 --ntasks-per-node=16  
#SBATCH -t 12:00:00                                                                                                                                       
#SBATCH -p normal                                                                                                                                          
#SBATCH -o slurm_output/slurm-%j.out # STDOUT                                                                                                          #SBATCH -e slurm_output/slurm-%j.err # STDERR


cd $HOME/EarthScatterLikelihood

#module load openmpi/gnu
#module load python/2.7.9

module load pre2019

module unload GCCcore
module load Python/2.7.12-intel-2016b
module load slurm-tools

export SLURM_CPU_BIND=none

ID="Final_75eV_res25eV_EdeBG_linB"

#time mpirun -np 16 python2.7 RunMPI_Contours.py -m_x 0.100 -data 1 -out_dir ${ID}1 -hemisphere N
#time mpirun -np 16 python2.7 RunMPI_Contours.py -m_x 0.100 -data 1 -out_dir ${ID}1 -hemisphere S
#time mpirun -np 16 python2.7 RunMPI_Contours.py -m_x 0.100 -data 3 -out_dir ${ID}3 -hemisphere N
time mpirun -np 16 python2.7 RunMPI_Contours.py -m_x 0.100 -data 3 -out_dir ${ID}3 -hemisphere S
