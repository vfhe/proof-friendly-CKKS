# Proof-Friendly CKKS



### Paper

[eprint 2025/286](https://ia.cr/2025/286)

```
@InProceedings{verifiable_ckks,
  author="Cascudo, Ignacio and Costache, Anamaria and Cozzo, Daniele and Fiore, Dario and Guimar{\~a}es, Antonio and Soria-Vazquez, Eduardo",
  editor="Tauman Kalai, Yael and Kamara, Seny F.",
  title="Verifiable Computation for Approximate Homomorphic Encryption Schemes",
  booktitle="Advances in Cryptology -- CRYPTO 2025",
  year="2025",
  publisher="Springer Nature Switzerland",
  address="Cham",
  pages="643--677",
  doi="10.1007/978-3-032-01907-3_21",
  isbn="978-3-032-01907-3"
}
```

> [!NOTE]
> This commit stores the reference code for the results presented in the [paper](https://ia.cr/2025/286).\
> Please note this is just a proof-of-concept implementation. \
> [Click here to see the latest version of our proof-friendly CKKS library](https://github.com/vfhe/proof-friendly-CKKS).

## Requirements

- Basic development packages (On Debian-based systems: `sudo apt install make gcc cmake`)
- Python 3 (Tested on Python 3.10.12)
- SageMath (tested with version 9.5)
- GCC/G++ (Tested on GCC 11.4.0)
- Strongly recommended: A processor with [AVX-512 IFMA](https://en.wikipedia.org/wiki/Advanced_Vector_Extensions#AVX-IFMA) (Tested on Intel i7-1165G7 and Intel Xeon 8488C).

## Compiling and running

### Clone this repository

```
git clone https://github.com/vfhe/proof-friendly-CKKS.git
```

### Running the Python code


The Python code will compile the C library during the first execution. You can recompile the C code using the following flags **at the end** of any python command:
```
  --full-recompile
    Deletes all building files and recompiles the entire code, including HEXL.
```

```
  --recompile
    Recompiles our C library, but not HEXL. Use it if you make any modifications in our C code.
```

```
  --cc=X
  --cxx=X
    Defines the compiler for C and CPP code, respectively. For example `--cxx=g++-15`.
```

> [!WARNING]
> The optimized version of this implementation requires [AVX-512 IFMA](https://en.wikipedia.org/wiki/Advanced_Vector_Extensions#AVX-IFMA). It should be significantly slower and possibly unstable without it (but it should still work).

#### Proof-friendly CKKS multiplication benchmark:

```
python3 -m pyckks.ckks
```

- Notice that the `-m pyckks.ckks` is only required because you are running a module inside our library.
- You can change parameters in lines 184 to 187 of [pyckks/ckks.py](./pyckks/ckks.py). 
- Defaults are `(N,d,L) = (2^14, 4, 6)`.

#### SumCheck benchmark:
```
python3 -m pyckks.bench_sumcheck 10
```
- Replace 10 with the number of variables.
- You can change parameters in line 46 of [pyckks/bench_sumcheck.py](./pyckks/bench_sumcheck.py). 
- Please note that this benchmark uses 48 threads. You can change this in line 47 of [lib/src/sumcheck.cpp](./lib/src/sumcheck.cpp). Remember to recompile the code with:
`
python3 pyckks.bench_sumcheck 10 --recompile
`

### Running CPP code
We also provide a benchmark for piecewise Reed-Solomon encoding.
```
cd lib
make main
./main
```
- You can change parameters in lines 26 to 30 of [lib/main_benchmark.cpp](./lib/main_benchmark.cpp). For default parameters, `MATRIX_y` is $5\sqrt{n}$, for a commitment with $\log_2(n)$ variables. `RS_cw_sz` must be adjusted accordingly, such that `MATRIX_x/RS_cw_sz` is the RS code rate.
- You can change the number of threads in line 24 of [lib/main_benchmark.cpp](./lib/main_benchmark.cpp). Notice that this will run the same code (but with different inputs of the same size) multiple times in parallel (since we always perform many commitments at once, there is no reason for internal parallelism). 

### Instructions for MAC

The Python library should also work on MAC, but performance will be **significantly** worse. We suggest following these steps:
1. Install a GNU C++ compiler, if you don't have one.
  ```
brew install gcc
```
2. Compile our code with the following command, replacing `X` with the version of g++ installed by brew. 
  ```
python3 -m pyckks.ckks --cxx=g++-X
```
3. If, for some reason, Python does not have access to the SageMath libraries, you can also run:
  ```
sage --python -m pyckks.ckks --cxx=g++-X
```
Note that:
- [Sage can also be installed on Mac OS with brew](https://formulae.brew.sh/cask/sage).
- The option `--cxx=g++-X` is only necessary when compiling for the first time, or when recompiling the C library. It is not needed afterwards. 

# License

[Apache License Version 2.0](./LICENSE)

Additionally, this repository contains copies or snippets of code from:

- [MOSFHET](https://github.com/antoniocgj/MOSFHET): [Apache License Version 2.0](https://github.com/antoniocgj/MOSFHET/blob/main/LICENSE) - Copyright Antonio Guimarães et al. - See their [detailed copyright information](https://github.com/antoniocgj/MOSFHET/tree/main?tab=readme-ov-file#license).
- [HELIOPOLIS](https://github.com/antoniocgj/HELIOPOLIS): [Apache License Version 2.0](https://github.com/antoniocgj/HELIOPOLIS/blob/main/LICENSE)
- [Intel HEXL](https://github.com/intel/hexl): [Apache License 2.0](https://github.com/intel/hexl/blob/development/LICENSE) - Copyright 2020 Intel Corporation
- [BLAKE3](https://github.com/BLAKE3-team/BLAKE3): [Apache License 2.0](https://github.com/BLAKE3-team/BLAKE3/blob/master/LICENSE_A2) - Copyright 2019 Jack O'Connor and Samuel Neves
