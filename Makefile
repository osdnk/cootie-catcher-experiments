CC = g++
INCLUDE_DIRS = ./cpp-bindings/ ./cpp-bindings/hexl/hexl/include/
FLAGS = -Wall -Wno-unused-function -Wno-unused-result -funroll-all-loops -march=native -lm -Wno-sign-compare -Wno-write-strings
LD_LIBS = hexl
LIB_DIRS = ./cpp-bindings/hexl/build/hexl/lib ./cpp-bindings/hexl/build/hexl/lib64
LIBS += $(addprefix -L, $(LIB_DIRS)) $(addprefix -l, $(LD_LIBS))
INCLUDE_FLAGS = $(addprefix -I, $(INCLUDE_DIRS))
FLAGS += $(INCLUDE_FLAGS)
OPT_FLAGS = -O3 -flto  $(FLAGS)
LIB_FLAGS = -O3 $(FLAGS)

ifeq ($(ENABLE_VAES), true)
	FLAGS += -DVAES_OPT
endif

ifeq ($(ENABLE_AVX512), true)
	FLAGS += -DAVX512_OPT
endif

SRC=polynomial.cpp misc.cpp

ALL_SRC = $(addprefix ./, $(SRC))

all: lib lib/libarith

wrapper:
	$(CC) -std=c++17 -fPIC -shared -o libhexl_wrapper.so cpp-bindings/hexl_wrapper.cpp $(OPT_FLAGS) $(LIBS)

wrapper_dummy:
	$(CC) -std=c++11 -fPIC -shared -o libhexl_wrapper.so cpp-bindings/empty.cpp $(OPT_FLAGS)

hexl: hexl/build
	cmake --build ./cpp-bindings/hexl/build -DCMAKE_C_COMPILER=gcc

hexl-triton: hexl/build-triton
	cmake --build ./cpp-bindings/hexl/build

hexl/build:
	cmake -S ./cpp-bindings/hexl/ -B ./cpp-bindings/hexl/build

hexl/build-triton:
	cmake -S ./cpp-bindings/hexl/ -B ./cpp-bindings/hexl/build -DCMAKE_C_COMPILER=/appl/scibuilder-spack/aalto-rhel9-prod/2024-01-compilers/software/linux-rhel9-haswell/gcc-11.4.1/gcc-12.3.0-xh5vv5d/bin/gcc

rust:
	export LD_LIBRARY_PATH=./cpp-bindings/hexl/build/hexl/lib:$(pwd)
	cargo bench
