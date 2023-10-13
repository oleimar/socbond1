#include "hdf5code.hpp"
#include <iostream>
#include <string>

// The Evo program runs evolutionary simulations
// Copyright (C) 2023  Olof Leimar
// See Readme.md for copyright notice

//*************************** Class h5R **********************************

void h5R::read_flt(std::string ds_name, f_type& dat)
{
    HighFive::DataSet ds = file.getDataSet(ds_name);
    ds.read(dat);
}

void h5R::read_flt_arr(std::string ds_name, vf_type& dat)
{
    HighFive::DataSet ds = file.getDataSet(ds_name);
    ds.read(dat);
}

void h5R::read_flt_mat(std::string ds_name, vvf_type& dat)
{
    HighFive::DataSet ds = file.getDataSet(ds_name);
    ds.read(dat);
}

void h5R::read_uint(std::string ds_name, ui_type& dat)
{
    HighFive::DataSet ds = file.getDataSet(ds_name);
    ds.read(dat);
}

void h5R::read_int(std::string ds_name, i_type& dat)
{
    HighFive::DataSet ds = file.getDataSet(ds_name);
    ds.read(dat);
}

void h5R::read_int_arr(std::string ds_name, vi_type& dat)
{
    HighFive::DataSet ds = file.getDataSet(ds_name);
    ds.read(dat);
}


//*************************** Class h5W **********************************

void h5W::write_flt(std::string ds_name, const f_type& dat)
{
    HighFive::DataSet ds =
        file.createDataSet<flt>(ds_name, HighFive::DataSpace::From(dat));
    ds.write(dat);
}

void h5W::write_flt_arr(std::string ds_name, const vf_type& dat)
{
    HighFive::DataSet ds =
        file.createDataSet<flt>(ds_name, HighFive::DataSpace::From(dat));
    ds.write(dat);
}

void h5W::write_flt_mat(std::string ds_name, const vvf_type& dat)
{
    HighFive::DataSet ds =
        file.createDataSet<flt>(ds_name, HighFive::DataSpace::From(dat));
    ds.write(dat);
}

void h5W::write_uint(std::string ds_name, const ui_type& dat)
{
    HighFive::DataSet ds =
        file.createDataSet<unsigned>(ds_name, HighFive::DataSpace::From(dat));
    ds.write(dat);
}

void h5W::write_int(std::string ds_name, const i_type& dat)
{
    HighFive::DataSet ds =
        file.createDataSet<int>(ds_name, HighFive::DataSpace::From(dat));
    ds.write(dat);
}

void h5W::write_int_arr(std::string ds_name, const vi_type& dat)
{
    HighFive::DataSet ds =
        file.createDataSet<int>(ds_name, HighFive::DataSpace::From(dat));
    ds.write(dat);
}
