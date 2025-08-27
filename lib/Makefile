CC = gcc
CXX = g++
INCLUDE_DIRS = ./include/ ./src/third-party/hexl/hexl/include/
FLAGS = -Wall -Wno-unused-function -Wno-unused-result -funroll-all-loops -march=native -lm -Wno-sign-compare -Wno-write-strings 
LD_LIBS = hexl
LIB_DIRS = ./src/third-party/hexl/build/hexl/lib 
LIBS += $(addprefix -L, $(LIB_DIRS)) $(addprefix -l, $(LD_LIBS)) 
INCLUDE_FLAGS = $(addprefix -I, $(INCLUDE_DIRS))
FLAGS += $(INCLUDE_FLAGS) -DAVX512_OPT
OPT_FLAGS = -O3 -fwhole-program -flto  $(FLAGS)
LIB_FLAGS = -O3 $(FLAGS)

ifeq ($(ENABLE_VAES), true)
	FLAGS += -DVAES_OPT
endif

ifeq ($(ENABLE_AVX512), true)
	FLAGS += -DAVX512_OPT
endif

SRC= third-party/misc_tp.cpp sumcheck.cpp
BLAKE_SRC=blake3.c blake3_dispatch.c blake3_portable.c blake3_sse2_x86-64_unix.S blake3_sse41_x86-64_unix.S blake3_avx2_x86-64_unix.S blake3_avx512_x86-64_unix.S
RNS_SRC=polynomial.cpp rlwe.cpp misc.cpp vector.cpp 
MP_SRC=mp_polynomial.cpp

ALL_SRC = $(addprefix ./src/, $(SRC))  $(addprefix ./src/third-party/blake3/, $(BLAKE_SRC)) $(addprefix ./src/rns-rlwe/, $(RNS_SRC))

all: lib lib/librings

main: $(ALL_SRC) main_benchmark.cpp
	$(CXX) -g -o main $^ $(OPT_FLAGS) $(LIBS)

lib/librings: $(ALL_SRC) 
	$(CXX) -g -fPIC -shared -o lib/librings.so $^ $(LIB_FLAGS) $(LIBS)

hexl: hexl/build
	cmake --build ./src/third-party/hexl/build

hexl/build:
	cmake -D CMAKE_C_COMPILER=$(CC) CMAKE_CXX_COMPILER=$(CXX) -S ./src/third-party/hexl/ -B ./src/third-party/hexl/build

lib:
	mkdir -p lib

clean:
	touch lib; touch src/third-party/hexl/build/; rm -rf lib/ src/third-party/hexl/build/

