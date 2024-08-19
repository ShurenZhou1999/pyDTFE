
# Delaunay Tessellation Field Interpolation (DTFE) with Python-wrapper

## Introduction

In this project, we provide easy Python operation of **Delaunay Tessellation Field Interpolation (DTFE)** method. The source code is rewritten from the public C++ [DTFE](https://github.com/MariusCautun/DTFE/) code, while the functions are wrapped by `Pybind11`. One can easily access to the DTFE method in friendly Python while keep the high efficiency of C++.


## Compilation
To use this tool, the C++ codes and its dependencies have to be previously installed and compiled. 

The wrapper `Pybind11` is header-only library, but the rest of libraries should be pre-compiled. Modify the bash file [compile.sh](./src/compile.sh) and make sure that you have loaded the necessary environments, then run it. 

* [CMake](https://cmake.org/)
* [Boost library](https://www.boost.org/)
* [CGAL](https://www.cgal.org/)
* [Pybind11](https://github.com/pybind/pybind11)
