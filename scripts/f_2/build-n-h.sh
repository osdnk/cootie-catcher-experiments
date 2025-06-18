#!/bin/bash

export LD_LIBRARY_PATH=./cpp-bindings/hexl/build/hexl/lib:$(pwd)
module load gcc mamba cmake

mamba deactivate
cargo build --release --bin main --features F_NEW_2,use-hardware
mv target/release/main target/release/F_NEW_2_H
mamba activate sage

