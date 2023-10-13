#ifndef EVOCODE_HPP
#define EVOCODE_HPP

// NOTE: for easier debugging one can comment out the following
#ifdef _OPENMP
#define PARA_RUN
#endif

#include "Genotype.hpp"
#include "Phenotype.hpp"
#include "Individual.hpp"
#include "Group.hpp"
#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <random>
#include <algorithm>

// The Evo program runs evolutionary simulations
// Copyright (C) 2023  Olof Leimar
// See Readme.md for copyright notice

// An individual has num_loci genetically determined traits (see
// Phenotype.hpp), and there is one locus for each trait

//************************* Class EvoInpData ***************************

// This class is used to 'package' input data in a single place; the
// constructor extracts data from an input file

class EvoInpData {
public:
    using flt = double;         // the TOML class requires doubles
    using ftype = std::vector<flt>;
    std::size_t max_num_thrds;  // Max number of threads to use
    int num_loci;               // Number of loci in individual's genotype
    int ng;                     // Number of groups
    int N;                      // Group size
    int K;                      // Number of places
    int T;                      // Number of steps in a reproductive interval
    int nrepi;                  // Number of reproductive intervals
    flt rho0;                   // initial value of association estimate
    flt xs;                     // initial value of bond strength
    flt xa;                     // asymptotic (maximum) value of bond strength
    flt shp;                    // shape parameter for beta distribution of q
    flt pd;                     // probability detect success in foraging 
    flt ps0;                    // probability of foraging success for q = 0
    flt g0;                     // parameter for helping asymptote
    flt alphp;                  // learning rate for place and association
    flt betd;                   // parameter for place preference function
    flt epsp;                   // parameter for random place choice
    flt ntld0;                  // place carrying capacity
    flt y0;                     // midpoint for bond strength in help function 
    flt mu0;                    // background mortality
    flt a;                      // parameter for survival function
    flt b;                      // parameter for survival function
    flt c;                      // parameter for survival function
    flt d;                      // parameter for helping function
    flt u0;                     // minimum amount of help
    int mxtries;                // max number of tries to find donor of help
    ftype mut_rate;             // Probability of mutation at each locus
    ftype SD;                   // SD of mutational increments at each locus
    ftype max_val;              // Maximal allelic value at each locus
    ftype min_val;              // Minimal allelic value at each locus
    ftype rho;                  // Recombination rates
    ftype all0;                 // Starting allelic values (if not from file)
    bool ind_reco;              // Whether to have individual recognition
    bool dhst;                 // Whether to save first group detailed history
    bool hist;                  // Whether to compute and save history
    bool read_from_file;        // Whether to read population from file
    std::string h5InName;       // File name for input of population
    std::string h5OutName;      // File name for output of population
    std::string h5dHistName;    // File name for output of detailed history
    std::string h5HistName;     // File name for output of history
    std::string h5IndsName;     // File name for output of individuals

    std::string InpName;  // Name of input data file
    bool OK;              // Whether input data has been successfully read

    EvoInpData(const char* filename);
};


//***************************** Class Evo ******************************

class Evo {
public:
    // types needed to define individual
    using mut_rec_type = MutRec<MutIncrBiExp<>>;
    // using mut_rec_type = MutRec<MutIncrNorm<>>;
    using gam_type = Gamete<mut_rec_type>;
    using gen_type = Haplotype<gam_type>;
    using dip_type = Diplotype<gam_type>;
    using phen_type = Phenotype<gen_type>;
    using ind_type = Individual<gen_type, phen_type>;
    using dhstat_type = DetHistStat<phen_type>;
    using hstat_type = HistStat<phen_type>;
    using istat_type = IndStat<phen_type>;
    using ss_type = typename phen_type::ss_type;
    using fVec = typename ss_type::fVec;
    using iVec = typename ss_type::iVec;
    // use std::vector containers for (sub)populations
    using vind_type = std::vector<ind_type>;
    using vph_type = std::vector<phen_type>;
    using grp_type = Group<phen_type>;
    using flt = double;
    using ftype = std::vector<flt>;
    using vftype = std::vector<ftype>;
    using vvftype = std::vector<vftype>;
    using uitype = std::vector<unsigned>;
    using itype = std::vector<int>;
    using vitype = std::vector<itype>;
    using vdhs_type = std::vector<dhstat_type>;
    using vhs_type = std::vector<hstat_type>;
    using vis_type = std::vector<istat_type>;
    using rand_eng = std::mt19937;
    using vre_type = std::vector<rand_eng>;
    using rand_int = std::uniform_int_distribution<int>;
    using rand_uni = std::uniform_real_distribution<flt>;
    using rand_norm = std::normal_distribution<flt>;
    using rand_discr = std::discrete_distribution<int>;
    using rand_gam = std::gamma_distribution<flt>;
    Evo(const EvoInpData& eid);
    void Run();
    void h5_read_pop(const std::string& infilename);
    void h5_write_pop(const std::string& outfilename) const;
    void h5_write_dhst(const std::string& dhstfilename) const;
    void h5_write_hist(const std::string& histfilename) const;
    void h5_write_inds(const std::string& indsfilename) const;
private:
    vind_type SelectReproduce(const vind_type& ppop, mut_rec_type& mr);

    EvoInpData id;
    int num_loci;
    int ng;
    int N;
    int Ntot;
    int K;
    int T;
    int nrepi;
    flt rho0;
    flt xs;
    flt xa;
    flt shp;
    flt pd;
    flt ps0;
    flt g0;
    flt alphp;
    flt betd;
    flt epsp;
    flt ntld0;
    flt y0;
    flt mu0;
    flt a;
    flt b;
    flt c;
    flt d;
    flt u0;
    int mxtries;
    bool ind_reco;
    bool dhst;
    bool hist;
    std::size_t num_thrds;
    uitype sds;
    vre_type vre;
    vind_type pop;
    vdhs_type dhstat;
    vhs_type hstat;
    vis_type istat;
};

#endif // EVOCODE_HPP
