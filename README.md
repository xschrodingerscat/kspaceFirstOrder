## Overview

k-Wave website (http://www.k-wave.org) for further details.


## Repository structure

    .
    +--Containers    - Matrix and output stream containers
    +--Data          - Small test data
    +--GetoptWin64   - Windows version of the getopt routine
    +--Hdf5          - HDF5 classes (file access)
    +--KSpaceSolver  - Solver classes with all the kernels
    +--Logger        - Logger class for reporting progress and errors
    +--MatrixClasses - Matrix classes holding simulation data
    +--OutputStreams - Output streams for sampling data
    +--Parameters    - Parameters of the simulation
    +--Utils         - Utility routines
    +--Docs					- Fftw and hdf5 source code, user's manuel for k-wave and papers 
    Changelog.md     - Change log
    License.md       - License file
    Makefile         - GNU Makefile
    Readme.md        - Read me
    Doxyfile         - Doxygen documentation file
    header_bg.png    - Doxygen logo
    main.cpp         - Main file of the project


## Compilation

The source codes of `kspaceFirstOrder-OMP`

1. Cmake version update

   ```bash
   wget https://cmake.org/files/v3.13/cmake-3.13.0-Linux-x86_64.tar.gz
   sudo apt autoremove cmake
   tar -zxvf cmake-3.13.0-Linux-x86_64.tar.gz
   mv cmake-3.13.0-Linux-x86_64 /opt/cmake-3.13.0
   ln -sf /opt/cmake-3.13.0/bin/* /usr/bin
   cmake --version
   ```

2. Ffwt lib installation

   ```bash
   tar -zxvf fftw-3.3.9.tar.gz
   cd fftw-3.3.9
   ./configure --enable-single --enable-avx --enable-openmp --enable-share
   make -j && make install
   ```

3. Hdf5 lib installation

   ```bash
   tar -zxvf CMake-hdf5-1.10.7.tar.gz
   cd CMake-hdf5-1.10.7
   ./build-unix.sh
   ./HDF5-1.10.7-Linux.sh
   cp -rf HDF_Group/HDF5/1.10.7/* /usr/
   ldconfig
   ```
   
4. kspaceFirstOrder app compilation

   ```bash
   make
   ```

   

## Usage

The C++ codes offers a lot of parameters and output flags to be used. For more
information, please type:

```bash
./kspaceFirstOrder-OMP --help
```

