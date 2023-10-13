cmake_minimum_required(VERSION 3.10)

# The Evo program runs evolutionary simulations
# Copyright (C) 2023  Olof Leimar
# See Readme.md for copyright notice

# set the project name
project(Evo VERSION 1.0 LANGUAGES CXX C)

# set c++ standard
set(CMAKE_CXX_STANDARD 17)

# we need HighFive to handle hdf5 files
find_package(HighFive REQUIRED)

find_package(OpenMP)

# add the executable
add_executable(Evo Evo.cpp EvoReci.cpp hdf5code.cpp Utils.cpp)

# link target
target_link_libraries(Evo PUBLIC HighFive OpenMP::OpenMP_CXX)
