#!/bin/bash -l
#SBATCH --mem=380G
#SBATCH --cpus-per-task 75
#SBATCH --time=00:10:00
#SBATCH --partition=gpu-h100-80g
module load mamba
module load gcc
module load rust
export LD_LIBRARY_PATH=./cpp-bindings/hexl/build/hexl/lib:$(pwd)
mamba activate sage

srun  ./target/release/main
