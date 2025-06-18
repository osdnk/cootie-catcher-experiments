#!/bin/bash

export LD_LIBRARY_PATH=./cpp-bindings/hexl/build/hexl/lib:$(pwd)
module load gcc mamba cmake

mamba deactivate
cargo build --release --bin main --features B_OLD,use-hardware
mv target/release/main target/release/B_OLD_H
mamba activate sage

