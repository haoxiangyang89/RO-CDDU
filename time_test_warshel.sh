#!/bin/sh
#SBATCH -J time_test
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH -p warshel-cpu
##SBATCH --nodelist=warshel-cpu02
##SBATCH --exclusive

module load gurobi/912 metis/4.0 openblas/0.3.15 julia/1.6.2

julia $1 $2 $3 $4 $5 $6
