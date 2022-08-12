#!/bin/bash
#PBS -q single
#PBS -l nodes=1:ppn=1
#PBS -l walltime=10:00:00
#PBS -o output
#PBS -N setup_phylocap
#PBS -A hpc_epic_frog

cd /home/chutter/

bash bash_manual_install.sh
