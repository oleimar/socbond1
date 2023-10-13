
# socbond1: C++ code for evolutionary simulation of helping with social bonds


## Overview

This repository contains C++ code and example data.
The executable program `Evo`, built from this code, will run evolutionary simulations of a population of several groups of individuals that build up social bonds through helping interactions.
The program was used to produce results for the paper "Social bond dynamics and the evolution of helping" by Olof Leimar and Redouan Bshary.


## System requirements

The program has been compiled and run on a Linux server with Ubuntu 22.04 LTS.
The C++ compiler was g++ version 11.4.0, provided by Ubuntu, with compiler flags for c++17, and `cmake` (<https://cmake.org/>) was used to build the program.
It can be run multithreaded using OpenMP, which speeds up execution times.
Most likely the instructions below will work for many Linux distributions.
The program has also been compiled and run on macOS, using the Apple supplied Clang version of g++.

The program reads input parameters from TOML files (<https://github.com/toml-lang/toml>), using the open source `cpptoml.h` header file (<https://github.com/skystrife/cpptoml>), which is included in this repository.

The program stores evolving populations in HDF5 files (<https://www.hdfgroup.org/>), which is an open source binary file format.
The program uses the open source HighFive library (<https://github.com/BlueBrain/HighFive>) to read and write to such files.
These pieces of software need to be installed in order for `cmake` to successfully build the program.


## Installation guide

Install the repository from Github to a local computer.
There is a single directory `socbond1` for source code and executable, a subdirectory `Data` where input data and data files containing simulated populations are kept, and a subdirectory `build` used by `cmake` for files generated during building, including the executable `Evo`.


## Building the program

The CMake build system is used.
If it does not exist, create a build subdirectory in the project folder (`mkdir build`) and make it the current directory (`cd build`).
If desired, for a build from scratch, delete any previous content (`rm -rf *`).
Run CMake from the build directory. For a release build:
```
cmake -D CMAKE_BUILD_TYPE=Release ../
```
and for a debug build replace Release with Debug.
If this succeeds, i.e. if the `CMakeLists.txt` file in the project folder is processed without problems, build the program:
```
cmake --build . --config Release
```
This should produce an executable in the `build` directory.


## Running

Make the Data directory current.
Assuming that the executable is called `Evo` and with an input file called `Run16.toml`, corresponding to case 1 in Table S1, run the program as
```
../build/Evo Run16.toml
```
Alternatively, using an R script file `Run16_run.R`, run the script as
```
Rscript Run16_run.R
```
where `Rscript` is the app for running R scripts.
You need to have `R` installed for this to work.


## Description of the evolutionary simulations

There is an input file, for instance `Run16.toml`, for a case, which typically simulates 10,000 periods, each consisting of `T = 20` time steps, inputting the population from, e.g., the HDF5 file `Run16.h5` and outputting to the same file.
Without an existing `Run16.h5` data file, the program can start by constructing individuals with genotypes from the allelic values given by `all0`in the input file.
To make this happen, use `read_from_file = false` in the input file.
Once you have a HDF5 file with a simulated population, there is an R script, e.g. `Run16_run.R`, which repeats runs a number of times, for instance 50, and for each run computes statistics on the evolving traits and adds a row to a TSV data file, e.g. `Run16_data.tsv`.

Each case is first run for many generations. When a seeming evolutionary equilibrium has been reached, 100 runs are kept in the summary data file, for instance `Run16_data.tsv`. These then represent evolution over `100*10,000*20 = 20,000,000` time steps. The final population file, e.g. `Run16.h5`, is then used for a simulation of helping using another input file, e.g. `Run16d.toml`, which collects and writes more data to file (such simulations are used for figures, e.g. for Fig. 2, 3 , 4, in the paper).


## License

The `Evo` program runs evolutionary simulations of helping in social groups.

Copyright (C) 2023  Olof Leimar

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
