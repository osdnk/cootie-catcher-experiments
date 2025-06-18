#!/bin/bash

export LD_LIBRARY_PATH=./cpp-bindings/hexl/build/hexl/lib:$(pwd)
module load gcc mamba cmake

mamba deactivate
cargo build --release --bin main --features F_OLD_3,use-hardware
mv target/release/main target/release/F_OLD_3_H
mamba activate sage

