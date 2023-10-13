#include "cpptoml.h"   // to read input parameters from TOML file
#include "EvoReci.hpp"
#include "hdf5code.hpp"
#include "Utils.hpp"
#include <algorithm>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>

#ifdef PARA_RUN
#include <omp.h>
#endif

// The Evo program runs evolutionary simulations
// Copyright (C) 2023  Olof Leimar
// See Readme.md for copyright notice

//************************** Read and ReadArr ****************************

// convenience functions to read from TOML input file

// this template function can be used for any type of single value
template<typename T>
void Get(std::shared_ptr<cpptoml::table> infile,
         T& value, const std::string& name)
{
    auto val = infile->get_as<T>(name);
    if (val) {
        value = *val;
    } else {
        std::cerr << "Read failed for identifier " << name << "\n";
    }
}

// this template function can be used for a vector or array (but there is no
// checking how many elements are read)
template<typename It>
void GetArr(std::shared_ptr<cpptoml::table> infile,
            It beg, const std::string& name)
{
    using valtype = typename std::iterator_traits<It>::value_type;
    auto vp = infile->get_array_of<valtype>(name);
    if (vp) {
        std::copy(vp->begin(), vp->end(), beg);
    } else {
        std::cerr << "Read failed for identifier " << name << "\n";
    }
}


//************************** class EvoInpData ****************************

EvoInpData::EvoInpData(const char* filename) :
      OK(false)
{
    auto idat = cpptoml::parse_file(filename);
    Get(idat, max_num_thrds, "max_num_thrds");
    Get(idat, num_loci, "num_loci");
    Get(idat, ng, "ng");
    Get(idat, N, "N");
    Get(idat, K, "K");
    Get(idat, T, "T");
    Get(idat, nrepi, "nrepi");
    Get(idat, rho0, "rho0");
    Get(idat, xs, "xs");
    Get(idat, xa, "xa");
    Get(idat, shp, "shp");
    Get(idat, pd, "pd");
    Get(idat, ps0, "ps0");
    Get(idat, g0, "g0");
    Get(idat, alphp, "alphp");
    Get(idat, betd, "betd");
    Get(idat, epsp, "epsp");
    Get(idat, ntld0, "ntld0");
    Get(idat, y0, "y0");
    Get(idat, mu0, "mu0");
    Get(idat, a, "a");
    Get(idat, b, "b");
    Get(idat, c, "c");
    Get(idat, d, "d");
    Get(idat, u0, "u0");
    Get(idat, mxtries, "mxtries");
    mut_rate.resize(num_loci);
    GetArr(idat, mut_rate.begin(), "mut_rate");
    SD.resize(num_loci);
    GetArr(idat, SD.begin(), "SD");
    max_val.resize(num_loci);
    GetArr(idat, max_val.begin(), "max_val");
    min_val.resize(num_loci);
    GetArr(idat, min_val.begin(), "min_val");
    rho.resize(num_loci);
    GetArr(idat, rho.begin(), "rho");
    Get(idat, ind_reco, "ind_reco");
    Get(idat, dhst, "dhst");
    Get(idat, hist, "hist");
    Get(idat, read_from_file, "read_from_file");
    if (read_from_file) {
        Get(idat, h5InName, "h5InName");
    } else {
        all0.resize(num_loci);
        GetArr(idat, all0.begin(), "all0");
    }
    Get(idat, h5OutName, "h5OutName");
    if (dhst) {
        Get(idat, h5dHistName, "h5dHistName");
    }
    if (hist) {
        Get(idat, h5HistName, "h5HistName");
        Get(idat, h5IndsName, "h5IndsName");
    }
    InpName = std::string(filename);
    OK = true;
}


//****************************** Class Evo *****************************

Evo::Evo(const EvoInpData& eid) :
    id{eid},
    num_loci{id.num_loci},
    ng{id.ng},
    N{id.N},
    Ntot{ng*N},
    K{id.K},
    T{id.T},
    nrepi{id.nrepi},
    rho0{static_cast<flt>(id.rho0)},
    xs{static_cast<flt>(id.xs)},
    xa{static_cast<flt>(id.xa)},
    shp{static_cast<flt>(id.shp)},
    pd{static_cast<flt>(id.pd)},
    ps0{static_cast<flt>(id.ps0)},
    g0{static_cast<flt>(id.g0)},
    alphp{static_cast<flt>(id.alphp)},
    betd{static_cast<flt>(id.betd)},
    epsp{static_cast<flt>(id.epsp)},
    ntld0{static_cast<flt>(id.ntld0)},
    y0{static_cast<flt>(id.y0)},
    mu0{static_cast<flt>(id.mu0)},
    a{static_cast<flt>(id.a)},
    b{static_cast<flt>(id.b)},
    c{static_cast<flt>(id.c)},
    d{static_cast<flt>(id.d)},
    u0{static_cast<flt>(id.u0)},
    mxtries{id.mxtries},
    ind_reco{id.ind_reco},
    dhst{id.dhst},
    hist{id.hist},
    num_thrds{1}
{
    // decide on number of threads for parallel processing
#ifdef PARA_RUN
    num_thrds = omp_get_max_threads();
    if (num_thrds > id.max_num_thrds) num_thrds = id.max_num_thrds;
    std::cout << "Number of threads: " << num_thrds << '\n';
#endif
    // generate one seed for each thread
    sds.resize(num_thrds);
    std::random_device rd;
    for (unsigned i = 0; i < num_thrds; ++i) {
        sds[i] = rd();
        // sds[i] = 1234;
        // sds[i] = 5678;
        vre.push_back(rand_eng(sds[i]));
    }

    // Note concerning thread safety: in order to avoid possible problems with
    // multiple threads, the std::vector container pop is allocated once and
    // for all here, and thread-local data are then copied into position in
    // pop (thus avoiding potentially unsafe push_back and insert).

    // create Ntot "placeholder individuals" in population
    gam_type gam(num_loci);
    ind_type indi(N, K, rho0, xs, ntld0, gam);
    pop.resize(Ntot, indi);
    // learning history stats
    if (hist) {
        if (ind_reco) {
            hstat.reserve(ng*T);
        }
        istat.reserve(ng*N);
    }

    // check if population data should be read from file
    if (id.read_from_file) {
        h5_read_pop(id.h5InName);
    } else {
        // random number engine
        rand_eng& eng = vre[0];
        // rand_uni uni(0, 1);
        rand_int rint(1, 10*Ntot); // for 'unique' individual IDs
        rand_gam rgam(shp, 1);
        // construct all individuals as essentially the same
        gam_type gam(num_loci); // starting gamete
        for (int l = 0; l < num_loci; ++l) {
            gam.gamdat[l] = static_cast<flt>(id.all0[l]);
        }
        int j = 0;
        for (int gn = 0; gn < ng; ++gn) { // groups
            for (int i = 0; i < N; ++i) { // inds in group
                ind_type ind(N, K, rho0, xs, ntld0, gam);
                phen_type& ph = ind.phenotype;
                // initialize some variables
                ph.gnum = gn;   // set group number
                ph.Set_inum(i); // set individual number
                ph.ID = rint(eng);
                // set random 'individual quality', i.e. q
                flt gm1 = rgam(eng);
                flt gm2 = rgam(eng);
                // q assumed to have a Beta(shp, shp) distribution
                ph.q = gm1/(gm1 + gm2);
                ph.ps = ps0 + (1 - ps0)*ph.q;
                ph.hx = ph.ha*(g0 + (1 - g0)*ph.q);
                ss_type& ss = ph.ss;
                if (ind_reco) {
                    ss.rho[i] = 0;
                    ss.x[i] = 0;
                }
                pop[j++] = ind;
            }
        }
    }
}

void Evo::Run()
{
    Timer timer(std::cout);
    timer.Start();
    ProgressBar PrBar(std::cout, nrepi);
    // set up  "global" mutation record, with engine and parameters controlling
    // mutation, segregation and recombination
    mut_rec_type mr0(vre[0], num_loci);
    for (int l = 0; l < num_loci; ++l) {
        mr0.mut_rate[l] = static_cast<flt>(id.mut_rate[l]);
        mr0.SD[l] = static_cast<flt>(id.SD[l]);
        mr0.max_val[l] = static_cast<flt>(id.max_val[l]);
        mr0.min_val[l] = static_cast<flt>(id.min_val[l]);
        mr0.rho[l] = static_cast<flt>(id.rho[l]);
    }
    // run through reproductive intervals
    for (int repi = 0; repi < nrepi; ++repi) {
        
        // start a reproductive interval by copying living individuals from the
        // previous pop (including the pop read from file or newly created), and
        // then replace dead individuals with new ones
        vind_type nextpop = SelectReproduce(pop, mr0);
        // copy individuals in nextpop to pop
        for (int gn = 0; gn < ng; ++gn) { // groups
            for (int i = 0; i < N; ++i) { // inds in group
                int j = gn*N + i;
                pop[j] = nextpop[j];
            }
        }
        
        // use parallel for processing over the interactions in groups
#pragma omp parallel for num_threads(num_thrds)
        for (int gn = 0; gn < ng; ++gn) {
#ifdef PARA_RUN
            int threadn = omp_get_thread_num();
#else
            int threadn = 0;
#endif
            // thread-local random number engine
            rand_eng& eng = vre[threadn];
            // rand_uni uni(0, 1);

            // thread-local container for group phenotypes
            vph_type gph;
            gph.reserve(N);
            for (int i = 0; i < N; ++i) {
                // copy phenotypes from population to group
                gph.push_back(pop[gn*N + i].phenotype);
            }

            // get history only for single threaded
            bool hst = false;
            if (num_thrds == 1 && hist) {
                hst = true;
            }

            // set up group interaction
            grp_type grp(N, K, T, repi, rho0, xs, xa, 
                pd, alphp, betd, epsp, ntld0, y0,
                mu0, a, b, c, d, u0, mxtries,
                gph, eng, dhst, hst);
            if (ind_reco) {
                // group interaction with individual recognition
                grp.Interact();
            } else {
                // group interaction without individual recognition
                grp.InteractNIR();
            }
            // get resulting phenotype of group members
            const vph_type& memb = grp.Get_memb();
            // copy individuals from group back into (global) container
            for (int i = 0; i < N; ++i) {
                pop[gn*N + i].phenotype = memb[i];
            }
            if (dhst) {
                // append detailed history from the group
                const vdhs_type& st = grp.Get_dhstat();
                dhstat.insert(dhstat.end(), st.begin(), st.end());
            }
            if (hst) {
                if (ind_reco) {
                    // append history from the group
                    const vhs_type& st = grp.Get_hstat();
                    hstat.insert(hstat.end(), st.begin(), st.end());
                }
                // append individuals from the group that died
                const vis_type& ist = grp.Get_istat();
                istat.insert(istat.end(), ist.begin(), ist.end());
            }
        }  // end of parallel for (over groups)

        // all set to start next reproductive interval (or to save the final
        // reproductive interval population)
        ++PrBar;
    }
    PrBar.Final();
    timer.Stop();
    timer.Display();
    h5_write_pop(id.h5OutName);
    if (dhst) {
        // for detailed history, typically only the first group is used to
        // collect data
        h5_write_dhst(id.h5dHistName);
    }
    if (hist) {
        h5_write_hist(id.h5HistName);
        h5_write_inds(id.h5IndsName);
    }
}

// return vector of Ntot individuals from the parent population in ppop, with
// living individuals being equally likely to deliver gametes to new individuals
// that replace dead ones, using mutation and recombination parameters from mr
Evo::vind_type Evo::SelectReproduce(const vind_type& ppop, mut_rec_type& mr)
{
    // create Ntot "placeholder individuals" in nextpop
    gam_type gam(num_loci);
    ind_type ind(N, K, rho0, xs, ntld0, gam);
    vind_type nextpop(Ntot, ind);
    int np = ppop.size(); // we ought to have np == Ntot
    if (np > 0) {
        // get discrete distribution with parental survival as weights
        ftype wei(np);
        for (int i = 0; i < np; ++i) {
            const phen_type& ph = ppop[i].phenotype;
            wei[i] = static_cast<flt>(ph.alive);
            if (wei[i] < 0.0) {
                wei[i] = 0.0;
            }
        }
        rand_uni uni(0, 1);
        rand_discr dscr(wei.begin(), wei.end());
        rand_int rint(1, 10*Ntot); // for 'unique' individual IDs
        rand_gam rgam(shp, 1);
        // get nextpop; copy living individuals in ppop to nextpop, and replace
        // dead ones through reproduction
        for (int gn = 0; gn < ng; ++gn) { // groups
            for (int i = 0; i < N; ++i) { // inds in group
                int ii = gn*N + i;
                const ind_type& indii = ppop[ii];
                const phen_type& phi = indii.phenotype;
                if (phi.alive) {
                    // copy this individual to its correct position in nextpop
                    // (i.e. same as previous position)
                    nextpop[ii] = indii;
                } else {
                    // replace dead individual with a new one
                    // find mother for individual to be constructed
                    int imat = dscr(mr.eng);
                    const ind_type& matind = ppop[imat];
                    // copy of maternal gamete
                    gam_type matgam = matind.genotype.gam;
                    // find father for individual to be constructed
                    int ipat = dscr(mr.eng);
                    const ind_type& patind = ppop[ipat];
                    // copy of paternal gamete
                    gam_type patgam = patind.genotype.gam;
                    // construct intermediate diploid
                    dip_type dip(matgam, patgam);
                    // construct new individual with gamete (with mutation and
                    // recombination) from dip
                    gam_type gam = dip.GetGamete(mr);
                    // construct the new individual
                    ind_type indi(N, K, rho0, xs, ntld0, gam);
                    phen_type& ph = indi.phenotype;
                    // initialize some variables
                    ph.gnum = gn;   // set group number
                    ph.Set_inum(i); // set individual number
                    ph.ID = rint(mr.eng);
                    // set random individual quality, i.e. q
                    flt gm1 = rgam(mr.eng);
                    flt gm2 = rgam(mr.eng);
                    // q assumed to have a Beta(shp, shp) distribution
                    ph.q = gm1/(gm1 + gm2);
                    // note that mean q is 0.5
                    ph.ps = ps0 + (1 - ps0)*ph.q;
                    ph.hx = ph.ha*(g0 + (1 - g0)*ph.q);
                    ss_type& ss = ph.ss;
                    if (ind_reco) {
                        ss.rho[i] = 0;
                        ss.x[i] = 0;
                    }
                    // place new individual in position k of nextpop
                    nextpop[ii] = indi;
                }
            }
        }
    }
    return nextpop;
}

void Evo::h5_read_pop(const std::string& infilename)
{
    // read data and put in pop
    h5R h5(infilename);
    std::vector<ftype> gams(Ntot, ftype(num_loci));
    // read gametes
    h5.read_flt_arr("Gam", gams);
    for (int i = 0; i < Ntot; ++i) {
        gam_type& gam = pop[i].genotype.gam;
        for (int l = 0; l < num_loci; ++l) {
            gam[l] = gams[i][l];
        }
    }
    ftype fval(Ntot);
    // bet
    h5.read_flt("bet", fval);
    for (int i = 0; i < Ntot; ++i) {
        pop[i].phenotype.bet = fval[i];
    }
    // pn0
    h5.read_flt("pn0", fval);
    for (int i = 0; i < Ntot; ++i) {
        pop[i].phenotype.pn0 = fval[i];
    }
    // alph
    h5.read_flt("alph", fval);
    for (int i = 0; i < Ntot; ++i) {
        pop[i].phenotype.alph = fval[i];
    }
    // ha
    h5.read_flt("ha", fval);
    for (int i = 0; i < Ntot; ++i) {
        pop[i].phenotype.ha = fval[i];
    }
    // hs
    h5.read_flt("hs", fval);
    for (int i = 0; i < Ntot; ++i) {
        pop[i].phenotype.hs = fval[i];
    }
    // std::vector to hold fVec member
    vftype fvec(Ntot, ftype(N));
    // rho
    h5.read_flt_arr("rho", fvec);
    for (int i = 0; i < Ntot; ++i) {
        fVec& rho = pop[i].phenotype.ss.rho;
        for (int j = 0; j < N; ++j) {
            rho[j] = fvec[i][j];
        }
    }
    // x
    h5.read_flt_arr("x", fvec);
    for (int i = 0; i < Ntot; ++i) {
        fVec& x = pop[i].phenotype.ss.x;
        for (int j = 0; j < N; ++j) {
            x[j] = fvec[i][j];
        }
    }
    // don
    h5.read_flt_arr("don", fvec);
    for (int i = 0; i < Ntot; ++i) {
        fVec& don = pop[i].phenotype.ss.don;
        for (int j = 0; j < N; ++j) {
            don[j] = fvec[i][j];
        }
    }
    // rec
    h5.read_flt_arr("rec", fvec);
    for (int i = 0; i < Ntot; ++i) {
        fVec& rec = pop[i].phenotype.ss.rec;
        for (int j = 0; j < N; ++j) {
            rec[j] = fvec[i][j];
        }
    }
    // std::vector to hold iVec member
    vitype ivec(Ntot, itype(N));
    // nd
    h5.read_int_arr("nd", ivec);
    for (int i = 0; i < Ntot; ++i) {
        iVec& nd = pop[i].phenotype.ss.nd;
        for (int j = 0; j < N; ++j) {
            nd[j] = ivec[i][j];
        }
    }
    // nr
    h5.read_int_arr("nr", ivec);
    for (int i = 0; i < Ntot; ++i) {
        iVec& nr = pop[i].phenotype.ss.nr;
        for (int j = 0; j < N; ++j) {
            nr[j] = ivec[i][j];
        }
    }
    // t0
    h5.read_int_arr("t0", ivec);
    for (int i = 0; i < Ntot; ++i) {
        iVec& t0 = pop[i].phenotype.ss.t0;
        for (int j = 0; j < N; ++j) {
            t0[j] = ivec[i][j];
        }
    }
    // std::vector to hold ftype member
    vftype pvec(Ntot, ftype(K));
    // nhat
    h5.read_flt_arr("nhat", pvec);
    for (int i = 0; i < Ntot; ++i) {
        ftype& nhat = pop[i].phenotype.nhat;
        for (int k = 0; k < K; ++k) {
            nhat[k] = pvec[i][k];
        }
    }
    // yhat
    h5.read_flt_arr("yhat", pvec);
    for (int i = 0; i < Ntot; ++i) {
        ftype& yhat = pop[i].phenotype.yhat;
        for (int k = 0; k < K; ++k) {
            yhat[k] = pvec[i][k];
        }
    }
    // q
    h5.read_flt("q", fval);
    for (int i = 0; i < Ntot; ++i) {
        pop[i].phenotype.q = fval[i];
    }
    // ps
    h5.read_flt("ps", fval);
    for (int i = 0; i < Ntot; ++i) {
        pop[i].phenotype.ps = fval[i];
    }
    // ha
    h5.read_flt("ha", fval);
    for (int i = 0; i < Ntot; ++i) {
        pop[i].phenotype.hx = fval[i];
    }
    // u
    h5.read_flt("u", fval);
    for (int i = 0; i < Ntot; ++i) {
        pop[i].phenotype.u = fval[i];
    }
    // std::vector to hold int (or bool) member
    itype ival(Ntot);
    // z
    h5.read_int("z", ival);
    for (int i = 0; i < Ntot; ++i) {
        pop[i].phenotype.z = ival[i];
    }
    // zet
    h5.read_flt("zet", fval);
    for (int i = 0; i < Ntot; ++i) {
        pop[i].phenotype.zet = fval[i];
    }
    // dont
    h5.read_flt("dont", fval);
    for (int i = 0; i < Ntot; ++i) {
        pop[i].phenotype.dont = fval[i];
    }
    // rect
    h5.read_flt("rect", fval);
    for (int i = 0; i < Ntot; ++i) {
        pop[i].phenotype.rect = fval[i];
    }
    // kp
    h5.read_int("kp", ival);
    for (int i = 0; i < Ntot; ++i) {
        pop[i].phenotype.kp = ival[i];
    }
    // age
    h5.read_int("age", ival);
    for (int i = 0; i < Ntot; ++i) {
        pop[i].phenotype.age = ival[i];
    }
    // ndh
    h5.read_int("ndh", ival);
    for (int i = 0; i < Ntot; ++i) {
        pop[i].phenotype.ndh = ival[i];
    }
    // nz0
    h5.read_int("nz0", ival);
    for (int i = 0; i < Ntot; ++i) {
        pop[i].phenotype.nz0 = ival[i];
    }
    // nrh
    h5.read_int("nrh", ival);
    for (int i = 0; i < Ntot; ++i) {
        pop[i].phenotype.nrh = ival[i];
    }
    // nc0
    h5.read_int("nc0", ival);
    for (int i = 0; i < Ntot; ++i) {
        pop[i].phenotype.nc0 = ival[i];
    }
    // nc1
    h5.read_int("nc1", ival);
    for (int i = 0; i < Ntot; ++i) {
        pop[i].phenotype.nc1 = ival[i];
    }
    // ID
    h5.read_int("ID", ival);
    for (int i = 0; i < Ntot; ++i) {
        pop[i].phenotype.ID = ival[i];
    }
    // inum
    h5.read_int("inum", ival);
    for (int i = 0; i < Ntot; ++i) {
        pop[i].phenotype.inum = ival[i];
    }
    // gnum
    h5.read_int("gnum", ival);
    for (int i = 0; i < Ntot; ++i) {
        pop[i].phenotype.gnum = ival[i];
    }
    // female
    h5.read_int("female", ival);
    for (int i = 0; i < Ntot; ++i) {
        pop[i].phenotype.female = ival[i];
    }
    // alive
    h5.read_int("alive", ival);
    for (int i = 0; i < Ntot; ++i) {
        pop[i].phenotype.alive = ival[i];
    }

}

void Evo::h5_write_pop(const std::string& outfilename) const
{
    h5W h5(outfilename);
    std::vector<ftype> gams(Ntot, ftype(num_loci));
    // write gametes
    for (int i = 0; i < Ntot; ++i) {
        const gam_type& gam = pop[i].genotype.gam;
        for (int l = 0; l < num_loci; ++l) {
            gams[i][l] = gam[l];
        }
    }
    h5.write_flt_arr("Gam", gams);
    // write members of phenotypes
    // std::vector to hold flt member
    ftype fval(Ntot);
    // bet
    std::transform(pop.begin(), pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.bet; });
    h5.write_flt("bet", fval);
    // pn0
    std::transform(pop.begin(), pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.pn0; });
    h5.write_flt("pn0", fval);
    // alph
    std::transform(pop.begin(), pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.alph; });
    h5.write_flt("alph", fval);
    // ha
    std::transform(pop.begin(), pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.ha; });
    h5.write_flt("ha", fval);
    // hs
    std::transform(pop.begin(), pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.hs; });
    h5.write_flt("hs", fval);
    // std::vector to hold fVec member
    vftype fvec(Ntot, ftype(N));
    // rho
    for (unsigned i = 0; i < Ntot; ++i) {
        const fVec& rho = pop[i].phenotype.ss.rho;
        for (int j = 0; j < N; ++j) {
            fvec[i][j] = rho[j];
        }
    }
    h5.write_flt_arr("rho", fvec);
    // x
    for (unsigned i = 0; i < Ntot; ++i) {
        const fVec& x = pop[i].phenotype.ss.x;
        for (int j = 0; j < N; ++j) {
            fvec[i][j] = x[j];
        }
    }
    h5.write_flt_arr("x", fvec);
    // don
    for (int i = 0; i < Ntot; ++i) {
        const fVec& don = pop[i].phenotype.ss.don;
        for (int j = 0; j < N; ++j) {
            fvec[i][j] = don[j];
        }
    }
    h5.write_flt_arr("don", fvec);
    // rec
    for (int i = 0; i < Ntot; ++i) {
        const fVec& rec = pop[i].phenotype.ss.rec;
        for (int j = 0; j < N; ++j) {
            fvec[i][j] = rec[j];
        }
    }
    h5.write_flt_arr("rec", fvec);
    // std::vector to hold iVec member
    vitype ivec(Ntot, itype(N));
    // nd
    for (int i = 0; i < Ntot; ++i) {
        const iVec& nd = pop[i].phenotype.ss.nd;
        for (int j = 0; j < N; ++j) {
            ivec[i][j] = nd[j];
        }
    }
    h5.write_int_arr("nd", ivec);
    // nr
    for (int i = 0; i < Ntot; ++i) {
        const iVec& nr = pop[i].phenotype.ss.nr;
        for (int j = 0; j < N; ++j) {
            ivec[i][j] = nr[j];
        }
    }
    h5.write_int_arr("nr", ivec);
    // t0
    for (int i = 0; i < Ntot; ++i) {
        const iVec& t0 = pop[i].phenotype.ss.t0;
        for (int j = 0; j < N; ++j) {
            ivec[i][j] = t0[j];
        }
    }
    h5.write_int_arr("t0", ivec);
    // std::vector to hold ftype member
    vftype pvec(Ntot, ftype(K));
    // nhat
    for (unsigned i = 0; i < Ntot; ++i) {
        const ftype& nhat = pop[i].phenotype.nhat;
        for (int k = 0; k < K; ++k) {
            pvec[i][k] = nhat[k];
        }
    }
    h5.write_flt_arr("nhat", pvec);
    // yhat
    for (unsigned i = 0; i < Ntot; ++i) {
        const ftype& yhat = pop[i].phenotype.yhat;
        for (int k = 0; k < K; ++k) {
            pvec[i][k] = yhat[k];
        }
    }
    h5.write_flt_arr("yhat", pvec);
    // q
    std::transform(pop.begin(), pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.q; });
    h5.write_flt("q", fval);
    // ps
    std::transform(pop.begin(), pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.ps; });
    h5.write_flt("ps", fval);
    // ha
    std::transform(pop.begin(), pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.hx; });
    h5.write_flt("ha", fval);
    // u
    std::transform(pop.begin(), pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.u; });
    h5.write_flt("u", fval);
    // std::vector to hold int member
    itype ival(Ntot);
    // z
    std::transform(pop.begin(), pop.end(), ival.begin(),
                   [](const ind_type& i) -> int
                   { return i.phenotype.z; });
    h5.write_int("z", ival);
    // zet
    std::transform(pop.begin(), pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.zet; });
    h5.write_flt("zet", fval);
    // dont
    std::transform(pop.begin(), pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.dont; });
    h5.write_flt("dont", fval);
    // rect
    std::transform(pop.begin(), pop.end(), fval.begin(),
                   [](const ind_type& i) -> flt
                   { return i.phenotype.rect; });
    h5.write_flt("rect", fval);
    // kp
    std::transform(pop.begin(), pop.end(), ival.begin(),
                   [](const ind_type& i) -> int
                   { return i.phenotype.kp; });
    h5.write_int("kp", ival);
    // age
    std::transform(pop.begin(), pop.end(), ival.begin(),
                   [](const ind_type& i) -> int
                   { return i.phenotype.age; });
    h5.write_int("age", ival);
    // ndh
    std::transform(pop.begin(), pop.end(), ival.begin(),
                   [](const ind_type& i) -> int
                   { return i.phenotype.ndh; });
    h5.write_int("ndh", ival);
    // nz0
    std::transform(pop.begin(), pop.end(), ival.begin(),
                   [](const ind_type& i) -> int
                   { return i.phenotype.nz0; });
    h5.write_int("nz0", ival);
    // nrh
    std::transform(pop.begin(), pop.end(), ival.begin(),
                   [](const ind_type& i) -> int
                   { return i.phenotype.nrh; });
    h5.write_int("nrh", ival);
    // nc0
    std::transform(pop.begin(), pop.end(), ival.begin(),
                   [](const ind_type& i) -> int
                   { return i.phenotype.nc0; });
    h5.write_int("nc0", ival);
    // nc1
    std::transform(pop.begin(), pop.end(), ival.begin(),
                   [](const ind_type& i) -> int
                   { return i.phenotype.nc1; });
    h5.write_int("nc1", ival);
    // ID
    std::transform(pop.begin(), pop.end(), ival.begin(),
                   [](const ind_type& i) -> int
                   { return i.phenotype.ID; });
    h5.write_int("ID", ival);
    // inum
    std::transform(pop.begin(), pop.end(), ival.begin(),
                   [](const ind_type& i) -> int
                   { return i.phenotype.inum; });
    h5.write_int("inum", ival);
    // gnum
    std::transform(pop.begin(), pop.end(), ival.begin(),
                   [](const ind_type& i) -> int
                   { return i.phenotype.gnum; });
    h5.write_int("gnum", ival);
    // female
    std::transform(pop.begin(), pop.end(), ival.begin(),
                   [](const ind_type& i) -> int
                   { return i.phenotype.female; });
    h5.write_int("female", ival);
    // alive
    std::transform(pop.begin(), pop.end(), ival.begin(),
                   [](const ind_type& i) -> int
                   { return i.phenotype.alive; });
    h5.write_int("alive", ival);
}

void Evo::h5_write_dhst(const std::string& dhstfilename) const
{
    // h5 file to save history in
    h5W lh5(dhstfilename);
    int hlen = dhstat.size();
    // std::vector to hold int int member
    itype hival(hlen);
    // gnum
    std::transform(dhstat.begin(), dhstat.end(), hival.begin(),
                   [](const dhstat_type& st) -> int
                   { return st.gnum; });
    lh5.write_int("gnum", hival);
    // kp
    std::transform(dhstat.begin(), dhstat.end(), hival.begin(),
                   [](const dhstat_type& st) -> int
                   { return st.kp; });
    lh5.write_int("kp", hival);
    // IDi
    std::transform(dhstat.begin(), dhstat.end(), hival.begin(),
                   [](const dhstat_type& st) -> int
                   { return st.IDi; });
    lh5.write_int("IDi", hival);
    // i
    std::transform(dhstat.begin(), dhstat.end(), hival.begin(),
                   [](const dhstat_type& st) -> int
                   { return st.i; });
    lh5.write_int("i", hival);
    // IDj
    std::transform(dhstat.begin(), dhstat.end(), hival.begin(),
                   [](const dhstat_type& st) -> int
                   { return st.IDj; });
    lh5.write_int("IDj", hival);
    // j
    std::transform(dhstat.begin(), dhstat.end(), hival.begin(),
                   [](const dhstat_type& st) -> int
                   { return st.j; });
    lh5.write_int("j", hival);
    // agei
    std::transform(dhstat.begin(), dhstat.end(), hival.begin(),
                   [](const dhstat_type& st) -> int
                   { return st.agei; });
    lh5.write_int("agei", hival);
    // agej
    std::transform(dhstat.begin(), dhstat.end(), hival.begin(),
                   [](const dhstat_type& st) -> int
                   { return st.agej; });
    lh5.write_int("agej", hival);
    // zi
    std::transform(dhstat.begin(), dhstat.end(), hival.begin(),
                   [](const dhstat_type& st) -> int
                   { return st.zi; });
    lh5.write_int("zi", hival);
    // zj
    std::transform(dhstat.begin(), dhstat.end(), hival.begin(),
                   [](const dhstat_type& st) -> int
                   { return st.zj; });
    lh5.write_int("zj", hival);
    // tstep
    std::transform(dhstat.begin(), dhstat.end(), hival.begin(),
                   [](const dhstat_type& st) -> int
                   { return st.tstep; });
    lh5.write_int("tstep", hival);
    // repi
    std::transform(dhstat.begin(), dhstat.end(), hival.begin(),
                   [](const dhstat_type& st) -> int
                   { return st.repi; });
    lh5.write_int("repi", hival);
    // std::vector to hold flt member
    ftype hfval(hlen);
    // qi
    std::transform(dhstat.begin(), dhstat.end(), hfval.begin(),
                   [](const dhstat_type& st) -> flt
                   { return st.qi; });
    lh5.write_flt("qi", hfval);
    // qj
    std::transform(dhstat.begin(), dhstat.end(), hfval.begin(),
                   [](const dhstat_type& st) -> flt
                   { return st.qj; });
    lh5.write_flt("qj", hfval);
    // uji
    std::transform(dhstat.begin(), dhstat.end(), hfval.begin(),
                   [](const dhstat_type& st) -> flt
                   { return st.uji; });
    lh5.write_flt("uji", hfval);
    // zeti
    std::transform(dhstat.begin(), dhstat.end(), hfval.begin(),
                   [](const dhstat_type& st) -> flt
                   { return st.zeti; });
    lh5.write_flt("zeti", hfval);
    // zetj
    std::transform(dhstat.begin(), dhstat.end(), hfval.begin(),
                   [](const dhstat_type& st) -> flt
                   { return st.zetj; });
    lh5.write_flt("zetj", hfval);
    // rhoij
    std::transform(dhstat.begin(), dhstat.end(), hfval.begin(),
                   [](const dhstat_type& st) -> flt
                   { return st.rhoij; });
    lh5.write_flt("rhoij", hfval);
    // xij
    std::transform(dhstat.begin(), dhstat.end(), hfval.begin(),
                   [](const dhstat_type& st) -> flt
                   { return st.xij; });
    lh5.write_flt("xij", hfval);
    // xji
    std::transform(dhstat.begin(), dhstat.end(), hfval.begin(),
                   [](const dhstat_type& st) -> flt
                   { return st.xji; });
    lh5.write_flt("xji", hfval);
    // recij
    std::transform(dhstat.begin(), dhstat.end(), hfval.begin(),
                   [](const dhstat_type& st) -> flt
                   { return st.recij; });
    lh5.write_flt("recij", hfval);
    // donij
    std::transform(dhstat.begin(), dhstat.end(), hfval.begin(),
                   [](const dhstat_type& st) -> flt
                   { return st.donij; });
    lh5.write_flt("donij", hfval);
    // nrij
    std::transform(dhstat.begin(), dhstat.end(), hival.begin(),
                   [](const dhstat_type& st) -> int
                   { return st.nrij; });
    lh5.write_int("nrij", hival);
    // ndij
    std::transform(dhstat.begin(), dhstat.end(), hival.begin(),
                   [](const dhstat_type& st) -> int
                   { return st.ndij; });
    lh5.write_int("ndij", hival);
}

void Evo::h5_write_hist(const std::string& histfilename) const
{
    // h5 file to save history in
    h5W lh5(histfilename);
    int hlen = hstat.size();
    // std::vector to hold int int member
    itype hival(hlen);
    // gnum
    std::transform(hstat.begin(), hstat.end(), hival.begin(),
                   [](const hstat_type& st) -> int
                   { return st.gnum; });
    lh5.write_int("gnum", hival);
    // IDi
    std::transform(hstat.begin(), hstat.end(), hival.begin(),
                   [](const hstat_type& st) -> int
                   { return st.IDi; });
    lh5.write_int("IDi", hival);
    // i
    std::transform(hstat.begin(), hstat.end(), hival.begin(),
                   [](const hstat_type& st) -> int
                   { return st.i; });
    lh5.write_int("i", hival);
    // IDj
    std::transform(hstat.begin(), hstat.end(), hival.begin(),
                   [](const hstat_type& st) -> int
                   { return st.IDj; });
    lh5.write_int("IDj", hival);
    // j
    std::transform(hstat.begin(), hstat.end(), hival.begin(),
                   [](const hstat_type& st) -> int
                   { return st.j; });
    lh5.write_int("j", hival);
    // t0
    std::transform(hstat.begin(), hstat.end(), hival.begin(),
                   [](const hstat_type& st) -> int
                   { return st.t0; });
    lh5.write_int("t0", hival);
    // t1
    std::transform(hstat.begin(), hstat.end(), hival.begin(),
                   [](const hstat_type& st) -> int
                   { return st.t1; });
    lh5.write_int("t1", hival);
    // End
    std::transform(hstat.begin(), hstat.end(), hival.begin(),
                   [](const hstat_type& st) -> int
                   { return st.End; });
    lh5.write_int("End", hival);
    // std::vector to hold flt member
    ftype hfval(hlen);
    // qi
    std::transform(hstat.begin(), hstat.end(), hfval.begin(),
                   [](const hstat_type& st) -> flt
                   { return st.qi; });
    lh5.write_flt("qi", hfval);
    // qj
    std::transform(hstat.begin(), hstat.end(), hfval.begin(),
                   [](const hstat_type& st) -> flt
                   { return st.qj; });
    lh5.write_flt("qj", hfval);
    // rhoij
    std::transform(hstat.begin(), hstat.end(), hfval.begin(),
                   [](const hstat_type& st) -> flt
                   { return st.rhoij; });
    lh5.write_flt("rhoij", hfval);
    // xij
    std::transform(hstat.begin(), hstat.end(), hfval.begin(),
                   [](const hstat_type& st) -> flt
                   { return st.xij; });
    lh5.write_flt("xij", hfval);
    // donij
    std::transform(hstat.begin(), hstat.end(), hfval.begin(),
                   [](const hstat_type& st) -> flt
                   { return st.donij; });
    lh5.write_flt("donij", hfval);
    // recij
    std::transform(hstat.begin(), hstat.end(), hfval.begin(),
                   [](const hstat_type& st) -> flt
                   { return st.recij; });
    lh5.write_flt("recij", hfval);
    // ndij
    std::transform(hstat.begin(), hstat.end(), hival.begin(),
                   [](const hstat_type& st) -> int
                   { return st.ndij; });
    lh5.write_int("ndij", hival);
    // nrij
    std::transform(hstat.begin(), hstat.end(), hival.begin(),
                   [](const hstat_type& st) -> int
                   { return st.nrij; });
    lh5.write_int("nrij", hival);
    // kpi
    std::transform(hstat.begin(), hstat.end(), hival.begin(),
                   [](const hstat_type& st) -> int
                   { return st.kpi; });
    lh5.write_int("kpi", hival);
    // kpj
    std::transform(hstat.begin(), hstat.end(), hival.begin(),
                   [](const hstat_type& st) -> int
                   { return st.kpj; });
    lh5.write_int("kpj", hival);
}

void Evo::h5_write_inds(const std::string& indsfilename) const
{
    // h5 file to save individuals in
    h5W h5(indsfilename);
    int ilen = istat.size();
    // std::vector to hold flt member
    ftype fval(ilen);
    // bet
    std::transform(istat.begin(), istat.end(), fval.begin(),
                   [](const istat_type& i) -> flt
                   { return i.bet; });
    h5.write_flt("bet", fval);
    // pn0
    std::transform(istat.begin(), istat.end(), fval.begin(),
                   [](const istat_type& i) -> flt
                   { return i.pn0; });
    h5.write_flt("pn0", fval);
    // alph
    std::transform(istat.begin(), istat.end(), fval.begin(),
                   [](const istat_type& i) -> flt
                   { return i.alph; });
    h5.write_flt("alph", fval);
    // ha
    std::transform(istat.begin(), istat.end(), fval.begin(),
                   [](const istat_type& i) -> flt
                   { return i.ha; });
    h5.write_flt("ha", fval);
    // hs
    std::transform(istat.begin(), istat.end(), fval.begin(),
                   [](const istat_type& i) -> flt
                   { return i.hs; });
    h5.write_flt("hs", fval);
    // q
    std::transform(istat.begin(), istat.end(), fval.begin(),
                   [](const istat_type& i) -> flt
                   { return i.q; });
    h5.write_flt("q", fval);
    // ps
    std::transform(istat.begin(), istat.end(), fval.begin(),
                   [](const istat_type& i) -> flt
                   { return i.ps; });
    h5.write_flt("ps", fval);
    // ha
    std::transform(istat.begin(), istat.end(), fval.begin(),
                   [](const istat_type& i) -> flt
                   { return i.hx; });
    h5.write_flt("ha", fval);
    // std::vector to hold int int member
    itype ival(ilen);
    // z
    std::transform(istat.begin(), istat.end(), ival.begin(),
                   [](const istat_type& i) -> int
                   { return i.z; });
    h5.write_int("z", ival);
    // zet
    std::transform(istat.begin(), istat.end(), fval.begin(),
                   [](const istat_type& i) -> flt
                   { return i.zet; });
    h5.write_flt("zet", fval);
    // dont
    std::transform(istat.begin(), istat.end(), fval.begin(),
                   [](const istat_type& i) -> flt
                   { return i.dont; });
    h5.write_flt("dont", fval);
    // rect
    std::transform(istat.begin(), istat.end(), fval.begin(),
                   [](const istat_type& i) -> flt
                   { return i.rect; });
    h5.write_flt("rect", fval);
    // kp
    std::transform(istat.begin(), istat.end(), ival.begin(),
                   [](const istat_type& i) -> int
                   { return i.kp; });
    h5.write_int("kp", ival);
    // age
    std::transform(istat.begin(), istat.end(), ival.begin(),
                   [](const istat_type& i) -> int
                   { return i.age; });
    h5.write_int("age", ival);
    // ndh
    std::transform(istat.begin(), istat.end(), ival.begin(),
                   [](const istat_type& i) -> int
                   { return i.ndh; });
    h5.write_int("ndh", ival);
    // nz0
    std::transform(istat.begin(), istat.end(), ival.begin(),
                   [](const istat_type& i) -> int
                   { return i.nz0; });
    h5.write_int("nz0", ival);
    // nrh
    std::transform(istat.begin(), istat.end(), ival.begin(),
                   [](const istat_type& i) -> int
                   { return i.nrh; });
    h5.write_int("nrh", ival);
    // nc0
    std::transform(istat.begin(), istat.end(), ival.begin(),
                   [](const istat_type& i) -> int
                   { return i.nc0; });
    h5.write_int("nc0", ival);
    // nc1
    std::transform(istat.begin(), istat.end(), ival.begin(),
                   [](const istat_type& i) -> int
                   { return i.nc1; });
    h5.write_int("nc1", ival);
    // ID
    std::transform(istat.begin(), istat.end(), ival.begin(),
                   [](const istat_type& i) -> int
                   { return i.ID; });
    h5.write_int("ID", ival);
    // inum
    std::transform(istat.begin(), istat.end(), ival.begin(),
                   [](const istat_type& i) -> int
                   { return i.inum; });
    h5.write_int("inum", ival);
    // gnum
    std::transform(istat.begin(), istat.end(), ival.begin(),
                   [](const istat_type& i) -> int
                   { return i.gnum; });
    h5.write_int("gnum", ival);
    // tstep
    std::transform(istat.begin(), istat.end(), ival.begin(),
                   [](const istat_type& i) -> int
                   { return i.tstep; });
    h5.write_int("tstep", ival);
    // repi
    std::transform(istat.begin(), istat.end(), ival.begin(),
                   [](const istat_type& i) -> int
                   { return i.repi; });
    h5.write_int("repi", ival);
}
