#!/bin/bash
 
#SBATCH -N 1 --ntasks-per-node=16  
#SBATCH -t 18:00:00                                                                                                                                        
#SBATCH -p normal                                                                                                                                          
#SBATCH -o slurm-%j.out # STDOUT                                                                                                                           
#SBATCH -e slurm-%j.err # STDERR



cd $HOME/EarthScatterLikelihood

#module load openmpi/gnu
#module load python/2.7.9

module load pre2019

module unload GCCcore
module load Python/2.7.12-intel-2016b
module load slurm-tools

#time mpirun -np 16 python2.7 RunMPI_Contours.py -m_x 0.200 -data 1 -out_dir Final1
time mpirun -np 16 python2.7 RunMPI_Contours.py -m_x 0.200 -data 3 -out_dir Final3