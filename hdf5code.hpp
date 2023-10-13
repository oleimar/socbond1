#ifndef HDF5CODE_HPP
#define HDF5CODE_HPP

#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <string>

// The Evo program runs evolutionary simulations
// Copyright (C) 2023  Olof Leimar
// See Readme.md for copyright notice

//************************* Class hdf5R **********************************

// This class is used to read a population of individuals from a HDF5 file. The
// file is assumed to have number of datasets in its root group, with names and
// content (data spaces) corresponding to the data fields in the individuals.
// The first dimension of each dataset is assumed to be equal to the number of
// individuals in the population.

class h5R {
public:
    using flt = double;
    using f_type = std::vector<flt>;
    using vf_type = std::vector<f_type>;
    using vvf_type = std::vector<vf_type>;
    using ui_type = std::vector<unsigned>;
    using i_type = std::vector<int>;
    using vi_type = std::vector<i_type>;
    // constructor opens file for reading
    h5R(std::string in_name) : file(in_name, HighFive::File::ReadOnly) {}
    void read_flt(std::string ds_name, f_type& dat);
    void read_flt_arr(std::string ds_name, vf_type& dat);
    void read_flt_mat(std::string ds_name, vvf_type& dat);
    void read_uint(std::string ds_name, ui_type& dat);
    void read_int(std::string ds_name, i_type& dat);
    void read_int_arr(std::string ds_name, vi_type& dat);
private:
    // H5::H5File file;  // file (closed when object goes out of scope)
    HighFive::File file;
};

//*************************** Class h5W **********************************

// This class is used to write a population of individuals to a HDF5 file,
// overwriting any content if the file already exists. The data is written as a
// number of datasets in the files root group, with names and content
// (data spaces) corresponding to the data fields in the individuals. The first
// dimension of each dataset is equal to the number of individuals in the
// population.

class h5W {
public:
    using flt = double;
    using f_type = std::vector<flt>;
    using vf_type = std::vector<f_type>;
    using vvf_type = std::vector<vf_type>;
    using ui_type = std::vector<unsigned>;
    using i_type = std::vector<int>;
    using vi_type = std::vector<i_type>;
    // constructor opens file for writing, truncating any previous file/content
    h5W(std::string out_name) :
        file(out_name,
             HighFive::File::ReadWrite |
             HighFive::File::Create |
             HighFive::File::Truncate) {}
    void write_flt(std::string ds_name, const f_type& dat);
    void write_flt_arr(std::string ds_name, const vf_type& dat);
    void write_flt_mat(std::string ds_name, const vvf_type& dat);
    void write_uint(std::string ds_name, const ui_type& dat);
    void write_int(std::string ds_name, const i_type& dat);
    void write_int_arr(std::string ds_name, const vi_type& dat);
private:
    HighFive::File file;
};

#endif // HDF5CODE_HPP
