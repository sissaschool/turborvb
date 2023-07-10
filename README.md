<img src="logo/turborvb_logo.png" width="70%">

# SISSA Quantum Monte Carlo package

[TurboRVB](https://turborvb.sissa.it) is an open-source computational package for **ab initio Quantum Monte Carlo (QMC) simulations** of both molecular and bulk electronic systems. The code was initially launched by **Sandro Sorella** and **Michele Casula** in SISSA and has been continuously developed by many contributors for over 20 years. The code implements two types of well established QMC algorithms: Variational Monte Carlo (VMC), and Diffusion Monte Carlo in its robust and efficient lattice regularized variant (LRDMC).

TurboRVB is distinguishable from other QMC codes in the following features:

- The code employs a resonating valence bond (RVB)-type wave function, such as the Jastrow Geminal/Jastrow Pfaffian. This wave function includes the correlation effect beyond the Jastrow-Slater wave function, which is commonly used in other QMC codes.

- Implemented state-of-art optimization algorithms, such as the stochastic reconfiguration and the linear method, realize a stable optimization of the amplitude and nodal surface of many-body wave functions at the variational quantum Monte Carlo level.

- The code implements the so-called lattice regularized diffusion Monte Carlo method, which provides a numerically stable diffusion quantum Monte Carlo calculation.

- The implementation of an adjoint algorithmic differentiation allows us to differentiate many-body wave functions efficiently and to perform structural optimizations and calculate molecular dynamics, thanks to an efficient and accurate evaluation of ionic forces.

- The possibility of performing different flavors of molecular dynamics, from classical to path integral, seamlessly integrated in the package, with original integration algorithms.

# Developers

- Kosuke Nakano (NIMS/Japan), kousuke_1123_at_icloud.com
- Otto Kohulák (CNRS/France), pravod_at_gmail.com
- Michele Casula (CNRS/France), michele.casula_at_gmail.com

# Installation guide

## CMake build

TurboRVB can be compiled with `CMake` build system providing convenient way to compile TurboRVB package assuming minimum version of `CMake` is 3.20.0 and compilers and math libraries (BLAS and LAPACK) are loaded on the system. One can run the cmake build by executing:

`cmake -S . -B build`

Option `-S` specifies path to the source directory (root directory of the repository) and option `-B` specifies build directory. This will execute the "configure" state of the compilation. It is recomended to run out-of-the-source compilation all the time. `CMake` will find the best suitable compilers and libraries on the system. At the end of the run a summary with compilation details will be shown such as:

```
 ######################### SUMMARY #########################

   C compiler:         NVHPC (nvc ver. 22.5.0)
   Fortran compiler:   NVHPC (nvfortran ver. 22.5.0)

   BLAS used:
   /opt/addman/spack/var/spack/environments/nvidia-openmpi/.spack-env/view/lib/libblas.so

   Compiling QMC code:          ON
   Compiling DFT code:          ON
   Compiling tools:             ON

   Compiling serial version:    ON
   Compiling parallel version:  ON
   Compiling GPU:               OFF

   Optimizations:               ON
   LTO optimizations:           OFF

   Base preprocessor flags are set to:
     ::  _SIMD
     ::  __FFTW
     ::  __PORT
     ::  __USE_INTERNAL_FFTW

   Parallel specific preprocessor flags are set to:
     ::  PARALLEL
     ::  __SCALAPACK

   Base Fortran flags:

   Agressive flags:         -O3  (55 files)
   Non-agressive flags:     -O1/-O0  (184 files)

 ######################### xxxxxxx #########################
```

One can, however, change compilers, adjust parameters or add flags by `CMake` options. Informations about these option can be found in the `CMake` [documentation](https://cmake.org/documentation/). The most important are:

`-DCMAKE_C_COMPILER=` (STRING) T C compiler to be used (e.g. `gcc`, `/opt/intel/bin/icc`, ...)

`-DCMAKE_Fortran_COMPILER=` (STRING) This specifies Fortran compiler to be used

`-DCMAKE_INSTALL_PREFIX=` (STRING) Install prefix for `make install`

Besides `CMake` native options there are other options specific to TurboRVB build, all starting with `EXT_` prefix:

`-DEXT_SERIAL=` (BOOL, default = ON) Build serial code

`-DEXT_PARALLEL=` (BOOL, default = ON) Build parallel code

`-DEXT_QMC=` (BOOL, default = ON) Build QMC code

`-DEXT_DFT=` (BOOL, default = ON) Build DFT code

`-DEXT_TOOLS=` (BOOL, default = ON) Build TurboRVB tools

`-DEXT_GPU=` (BOOL, default = ON) Build GPU accelerated code

`-DEXT_OPT=` (BOOL, default = ON) Turn on optimizations

`-DEXT_LTO=` (BOOL, default = OFF) Turn on link time optimizations. Beware, this can make compilation very long (even hours), also sometimes creates instable code. But can provide some extra optimizations.

`-DEXT_SPEEDTEST=` (BOOL, default = OFF) Enable speed test targets

`-DEXT_DEBUG=` (BOOL, default = OFF) Turn on debug build

`-DEXT_TIME=` (BOOL, default = OFF) Turn on internal timer

`-DEXT_BLAS_LIB=` (STRING, default = "") Specify BLAS library, path to BLAS library that should be used

`-DEXT_LAPACK_LIB=` (STRING, default = "") Specify LAPACK library, path to LAPACK library that should be used

`-DEXT_OTHER_LIB=` (LIST of STRINGs, default = "") Other libraries that should be linked.

`-DEXT_FLAGS=` (LIST of STRINGs, default = "") Add preprocessor flags

`-DEXT_BLKL_FLAGS=` (LIST of STRINGs, default = "") Black list preprocessor flags (these will not be used)

`-DEXT_AGRS=` (LIST of STRINGs, default = "") Add compiler option for agressive optimization

`-DEXT_BLKL_AGRS=` (LIST of STRINGs, default = "") Black list compiler flags for aggresive optimization (these will not be used)

`-DEXT_GPUTYPE=` (STRING, default = "") Specify gpu type if necessary (compute capability, only for Nvidia GPUs, e.g. cc75, cc80 ...)

Here we provide an example `CMake` command. Let's say we would like to compile only parallel version of TurboRVB without build in DFT code. We would like to use our own OpenBLAS library, turning off internal `_SIMD` but adding flags `_NVTX`, `_DUNREL` and internal timers. There are multiple compilers installed on our system we would like to use nvidia one. The `CMake` command then should look like:

```
cmake -S. -Bbuild                                   \
          -DCMAKE_C_COMPILER=nvc                    \
          -DCMAKE_Fortran_COMPILER=nvfortran        \
          -DEXT_TIME=ON                             \
          -DEXT_DFT=OFF                             \
          -DEXT_SERIAL=OFF                          \
          -DEXT_BLAS_LIB=/opt/lib/libopenblas.a     \
          -DEXT_LAPACK_LIB=/opt/lib/libopenblas.a   \
          -DEXT_BLKL_FLAGS=_SIMD                    \
          -DEXT_FLAGS="_NVTX;_DUNREL"
```

These setting can be changed even after the cmake command was executed. The easiest way to do so is via TUI interface called `ccmake`. One only has to specify build directory:

`ccmake build`

Next we need to compile the code. We recomend to use some implementation of `Make` as a generator. We can report `Ninja` does not work in this case. One can invoke compilation by:

`cmake --build build`

or manually by executing make:

`cd build && make VERBOSE=1`

It is important to note, the parallel make does not work. Successfull compilation can checked by executing `ctest` inside the build directory. If all of the tests pass, one can safely use TurboRVB package.

By addition one can run `make speed_tests` evaluating speed of TurboRVB package on her/his machine.

## Notes & known compilation issues

1) If present, intel `ifort` compiler with `MKL` should be preferred to gnu `gfortran` compiler with `BLAS` and `LAPACK`, because it has been tested more and yields consistently more reliable and faster binaries.

2) There is an issue with the BLAS `gfortran` with `Mac`, in particular in zaxpy for the complex case. (see also Quantum Espresso documentation)

3) There is an issue with some version of BLAS provided by Accelerate Framework library. If test zdotc will not pas try to add -DEXT_FLAGS="_FIXBUG_ZDOTC", it should fix the problem.

# Running a docker container

Before you start, make sure Docker is installed and set up on your machine.

You can run TurboRVB from a Docker container by pulling the following Docker image:

`docker pull addman151/turborvb:latest`

The Docker container has all the required executables in its PATH. You can run them directly, just ensure the working directory is properly mounted and environment variables are set as needed.

In the commands below:

The `-i` flag starts the container in interactive mode.
The `-e` flag sets an environment variable inside the container. Here `OMP_NUM_THREADS=4` specifies the number of threads that OpenMP should use.
The `-v` flag mounts the current directory (as returned by $(pwd)) to `/app` inside the container. This allows the container to read and write files from your current directory.
The `-w` flag sets the working directory inside the Docker container. Here, it's set to `/app`.

Note: Ensure that the datasvmc.input and prep.input files are present in your current directory before running the commands.

Here are a few example commands to run TurboRVB:

```
docker run -i -e OMP_NUM_THREADS=4 -v "$(pwd):/app" -w /app addman151/turborvb:latest turborvb-serial.x < datasvmc.input
docker run -i -e OMP_NUM_THREADS=4 -v "$(pwd):/app" -w /app addman151/turborvb:latest mpirun -np 2 --oversubscribe turborvb-mpi.x < datasvmc.input
docker run -i -e OMP_NUM_THREADS=4 -v "$(pwd):/app" -w /app addman151/turborvb:latest prep-serial.x < prep.input
docker run -i -e OMP_NUM_THREADS=4 -v "$(pwd):/app" -w /app addman151/turborvb:latest mpirun -np 2 --oversubscribe prep-mpi.x < prep.input
```

# Environmental variable

In order to run correctly all the tools, you should put in your path the `turborvb/bin` directory of `TurboRVB`;  e.g. by using bash shell you should edit the file `.bashrc` in your home directory and you should add the line:

# Tests of compilation

We recommend you should run `jobcheck_serial.sh` in `test` directory every time you compile `TurboRVB` on your machine. If you compile `TurboRVB` with the modern `CMake`, you can do a more comprehensive test with the `ctest` command (See the CMake section).

# Manuals for users

There is a Read the Docs in the `doc` directory. Documented by `Sphinx`. The manual is deployed from [GitHub Pages](https://sissaschool.github.io/turborvb/).

# Tutorials for users

Tutorials are provided in another project. Visit the [turbotutorials repository](https://github.com/kousuke-nakano/turbotutorials) on the GitHub.

# Reference(s)

When you publish a paper using `TurboRVB`, please cite the following paper:

K. Nakano et al., "TurboRVB: A many-body toolkit for ab initio electronic simulations by quantum Monte Carlo" [**J. Chem. Phys.** 152, 204121 (2020)](https://aip.scitation.org/doi/10.1063/5.0005037).

# Manuals for developers

`TurboRVB` developer manual is automatically generated using [ford](https://forddocs.readthedocs.io/en/latest/index.html). To generate the developer manual, you should install `ford` according to the instruction written in the [website](https://forddocs.readthedocs.io/en/latest/index.html), go to `devel_tools/` directory, and make the document by the command `ford ford_config.md`. All developers should follow the rules written in `devel_tools/rules.md`

# License

All the material included in this distribution is free software;
you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation;
either version 3 of the License, or (at your option) any later version.

These programs are distributed in the hope that they will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
675 Mass Ave, Cambridge, MA 02139, USA.
