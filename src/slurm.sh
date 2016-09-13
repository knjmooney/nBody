#!/bin/bash

#SBATCH -n 32         # 8 cores = 1 node on lonsdale
#SBATCH -p compute    # dubug application compute?
#SBATCH -t 2-00:00:00   # Time limit
#SBATCH -U mschpc     # Group
#SBATCH -J n-body     # Job name

# source the module commands
source /etc/profile.d/modules.sh

# load the modules used to build the xhpl binary
module load apps libs cports openmpi/1.8.6-gnu gcc/4.9.3-gnu

# run it
mpirun -n 1 ./n-body  -n 1536 -s 100000 -w 100 -h 100 -t 0.00000001
mpirun -n 2 ./n-body  -n 1536 -s 100000 -w 100 -h 100 -t 0.00000001
mpirun -n 4 ./n-body  -n 1536 -s 100000 -w 100 -h 100 -t 0.00000001
mpirun -n 8 ./n-body  -n 1536 -s 100000 -w 100 -h 100 -t 0.00000001
mpirun -n 16 ./n-body -n 1536 -s 100000 -w 100 -h 100 -t 0.00000001
mpirun -n 24 ./n-body -n 1536 -s 100000 -w 100 -h 100 -t 0.00000001
mpirun -n 32 ./n-body -n 1536 -s 100000 -w 100 -h 100 -t 0.00000001

mpirun -n 1 ./n-body  -n 2304 -s 100000 -w 100 -h 100 -t 0.00000001
mpirun -n 2 ./n-body  -n 2304 -s 100000 -w 100 -h 100 -t 0.00000001
mpirun -n 4 ./n-body  -n 2304 -s 100000 -w 100 -h 100 -t 0.00000001
mpirun -n 8 ./n-body  -n 2304 -s 100000 -w 100 -h 100 -t 0.00000001
mpirun -n 16 ./n-body -n 2304 -s 100000 -w 100 -h 100 -t 0.00000001
mpirun -n 24 ./n-body -n 2304 -s 100000 -w 100 -h 100 -t 0.00000001
mpirun -n 32 ./n-body -n 2304 -s 100000 -w 100 -h 100 -t 0.00000001

