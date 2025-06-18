#!/bin/bash

export LD_LIBRARY_PATH=./cpp-bindings/hexl/build/hexl/lib:$(pwd)
module load gcc mamba cmake

mamba deactivate
cargo build --release --bin main --features D_NEW,use-hardware
mv target/release/main target/release/D_NEW_H
mamba activate sage

