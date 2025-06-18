#!/bin/bash -l
#SBATCH --mem=100G
#SBATCH --cpus-per-task 75
#SBATCH --time=1:00:00
#SBATCH --partition=gpu-h100-80g
module load mamba
module load gcc
module load rust
mamba activate sage
export LD_LIBRARY_PATH=./cpp-bindings/hexl/build/hexl/lib:$(pwd)

srun  ./target/release/main

