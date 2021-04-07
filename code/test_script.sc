#!/bin/bash -l
#SBATCH -J 870_6_6
#SBATCH -o output_%j.txt
#SBATCH -e errors_%j.txt
#SBATCH -t 40:00:00
#SBATCH -n 6
#SBATCH -c 2
#SBATCH -p batch
#
hostname
slurmnodes
OMP_NUM_THREADS=2 mpiexec -hostfile hostfile -np 6 python3 run_test.py "/home/sc20-790674/pdb/pdb_870_samples" 6 6
