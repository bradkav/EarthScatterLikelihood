#!/usr/bin/env python
from mpi4py import MPI
from subprocess import call
import numpy as np

import sys

import argparse

comm = MPI.COMM_WORLD
#Get total number of MPI processes
nprocs = comm.Get_size()
#Get rank of current process
rank = comm.Get_rank()

#R_list = np.logspace(np.log10(args.R_ini), np.log10(args.R_fin), nprocs)

#Parse the arguments!
parser = argparse.ArgumentParser(description='...')
parser.add_argument('-m_x','--m_x', help='DM mass in GeV', type=float,default = 0.1)
#parser.add_argument('-sigma_p','--sigma_p', help='DM-nucleon cross section, sigma_p in cm^2', type=float, required=True)
#parser.add_argument('-sigma_SI','--sigma_SI', help='Cross section in cm^2', type = float, default=1e-34)
parser.add_argument('-data', '--data', help='Type of data (1, 2, 3)', type=int, default = 1)
#parser.add_argument('-lat', '--lat', help='Detector latitude', type=float, default=45.454)
parser.add_argument('-out_dir', '-out_dir', help='Output directory', type=str, required=True)
parser.add_argument('-hemisphere', '-hemisphere', help="N or S", type=str, required=True)
parser.add_argument('-exe', '-exe', type=str, default="calcContour")

args = parser.parse_args()


#MC_script.py -R 1 -N 20000

#Directory where the calc files are located
myDir = "/home/kavanagh/EarthScatterLikelihood/"
cmd = "cd "+myDir+" ; "
#cmd += "-R " + str(R_list[rank])
#cmd += " -N " + str(args.N_AMC)

sig_list = np.logspace(-37, -30, 15)

if (rank < 15):
    if (args.hemisphere == "N"):
        cmd += " ./" + args.exe + " "  + str(args.m_x) + " " + str(sig_list[rank]) + " " + str(args.data) + " 0 45.454 " + args.out_dir + "_N ; "
        cmd += " ./" + args.exe + " "  + str(args.m_x) + " " + str(sig_list[rank]) + " " + str(args.data) + " 1 45.454 " + args.out_dir + "_N ; "
    elif (args.hemisphere == "S"):
        cmd += " ./" + args.exe + " "  + str(args.m_x) + " " + str(sig_list[rank]) + " " + str(args.data) + " 0 -37.07 " + args.out_dir + "_S; " 
        cmd += " ./" + args.exe + " "  + str(args.m_x) + " " + str(sig_list[rank]) + " " + str(args.data) + " 1 -37.07 " + args.out_dir + "_S" 

#cmd += " >> script_output"
    
sts = call(cmd,shell=True)
comm.Barrier()
