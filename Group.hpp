#ifndef GROUP_HPP
#define GROUP_HPP

#include <vector>
#include <random>
#include <cmath>
#include <algorithm>

// The Evo program runs evolutionary simulations
// Copyright (C) 2023  Olof Leimar
// See Readme.md for copyright notice

//************************* struct DetHistStat ***************************

// This struct stores detailed data on a single helping interaction between
// two group members (i and j); a sequence of structs can be saved as a record
// of the interaction history in the group (typically the first group)

template<typename PhenType>
struct DetHistStat {
    using phen_type = PhenType;
    using flt = typename phen_type::flt;
    int gnum;      // group number
    int kp;        // place for i and j
    int IDi;       // ID for recipient individual
    int i;         // inum for recipient individual
    int IDj;       // ID for donor individual 
    int j;         // inum for donor individual 
    int agei;      // age for i
    int agej;      // age for j
    int zi;        // resource state for i
    int zj;        // resource state for j
    int tstep;     // time in current reproductive interval
    int repi;      // current reproductive interval
    flt qi;        // quality for i
    flt qj;        // quality for j
    flt uji;       // amount transferred
    flt zeti;      // current resource balance for i
    flt zetj;      // current resource balance for i
    flt rhoij;     // association between i and j
    flt xij;       // current social state for i
    flt xji;       // current social state for i
    flt recij;     // previous total received by i from j
    flt donij;     // previous total donated by i to j
    int nrij;      // previous times received by i from j
    int ndij;      // previous times donated by i to j
};

//************************* struct HistStat *******************************

// This struct stores data on the relationship history between two group
// members (i and j); a sequence of structs can be saved as a record of the
// relationship history in the group

// Ending types: 0 is i dies, 1 is j dies

template<typename PhenType>
struct HistStat {
    using phen_type = PhenType;
    using flt = typename phen_type::flt;
    int gnum;      // group number
    int IDi;       // ID for focal individual
    int i;         // inum for focal individual
    int IDj;       // ID for other individual 
    int j;         // inum for other individual 
    int t0;        // age for i when partnership starts
    int t1;        // age for i when partnership ends
    int End;       // ending type
    flt qi;        // quality for i
    flt qj;        // quality for j
    flt rhoij;     // social state for i at end
    flt xij;       // social state for i at end
    flt donij;     // don for i to j
    flt recij;     // rec for i from j
    int ndij;      // number of times donated for i to j
    int nrij;      // number of times received for i from j
    int kpi;       // place for i at end
    int kpj;       // place for j at end
};

//************************** struct IndStat *******************************

// This struct stores data on the life of an individual, up to death; it is
// mostly the same as the data for the class Phenotype, except that the social
// state data (ss) are not included.

template<typename PhenType>
struct IndStat {
    using phen_type = PhenType;
    using flt = typename phen_type::flt;
    flt bet;     // place choice parameter (genetically determined)
    flt pn0;     // bond plus for new partner (genetically determined)
    flt alph;    // learning rate (genetically determined)
    flt ha;      // helping asymptote
    flt hs;      // helping sensitivity to helping balance
    flt q;       // individual quality, 0 < q < 1
    flt ps;      // probability succeed in foraging
    flt hx;      // max helping amount (when z = 1)
    int z;       // resource state (0/1)
    flt zet;     // net help received final period
    flt dont;    // total help donated
    flt rect;    // total help received
    int kp;      // current place
    int age;     // age (number of periods)
    int ndh;     // number of donations of help
    int nz0;     // number of interactions as recipient
    int nrh;     // number of times receiving help
    int nc0;     // number of partnerships started
    int nc1;     // number of partnerships ended by death
    int ID;      // 'unique' identifying number
    int inum;    // individual number
    int gnum;    // group number
    int tstep;   // time step of interaction in group
    int repi;    // reproductive interval number
};

//**************************** class Group **********************************

// This class deals with the interactions in one group

template<typename PhenType>
class Group {
public:
    using phen_type = PhenType;
    using ss_type = typename phen_type::ss_type;
    using fVec = typename ss_type::fVec;
    using iVec = typename ss_type::iVec;
    using vph_type = std::vector<phen_type>;
    using flt = typename phen_type::flt;
    using ftype = std::vector<flt>;
    using itype = std::vector<int>;
    using dhstat_type = DetHistStat<phen_type>;
    using vdhs_type = std::vector<dhstat_type>;
    using hstat_type = HistStat<phen_type>;
    using vhs_type = std::vector<hstat_type>;
    using istat_type = IndStat<phen_type>;
    using vis_type = std::vector<istat_type>;
    using rand_eng = std::mt19937;
    using rand_uni = std::uniform_real_distribution<flt>;
    using rand_int = std::uniform_int_distribution<int>;
    using rand_discr = std::discrete_distribution<int>;
    Group(int a_N,
        int a_K,
        int a_T,
        int a_repi,
        flt a_rho0,
        flt a_xs,
        flt a_xa,
        flt a_pd,
        flt a_alphp,
        flt a_betd,
        flt a_epsp,
        flt ntld0,
        flt a_y0,
        flt a_mu0,
        flt a_a,
        flt a_b,
        flt a_c,
        flt a_d,
        flt a_u0,
        int a_mxtries,
        const vph_type& a_memb,
        rand_eng& a_eng,
        bool a_dhst = false,
        bool a_hist = false);
    const vph_type& Get_memb() const { return memb; }
    const vdhs_type& Get_dhstat() const { return dhstat; }
    const vhs_type& Get_hstat() const { return hstat; }
    const vis_type& Get_istat() const { return istat; }
    void Interact();
    void InteractNIR(); // without individual recognition

private:
    // amount of help as function of 'bond strength' 
    flt hlp(flt rhox, int z, flt zet, flt hx, flt hsw) {
        flt hz = hx/(1 + std::exp(-b*z - c*zet - hsw)); 
        return hz/(1 + std::exp(- d*(rhox - y0))); 
    }
    // mortality as a function of states z and zet
    flt mu(int z, flt zet) { 
        return (1 - mu0)/(1 +  std::exp(a + b*z + c*zet)) + mu0;
    }
    void Add_dhstat(int tstep, int i, int j);
    void Add_hstat(int i, int j, int End);
    void Add_istat(int tstep, int i);
    int N;           // group size
    int K;           // number of places
    int T;           // number of periods in group
    int repi;        // reproductive interval number
    flt rho0;        // initial value of association estimate
    flt xs;          // initial value of bond strength
    flt xa;          // asymptotic (maximum) value of bond strength
    flt pd;          // probability detection of success in foraging
    flt alphp;       // learning rate for place and association
    flt betd;        // parameter for place preference function
    flt epsp;        // parameter for random place choice
    ftype ntld;      // carrying capacities of places
    flt y0;          // midpoint for bond strength in helping function
    flt mu0;         // background mortality
    flt a;           // parameter for survival function
    flt b;           // parameter for survival function
    flt c;           // parameter for survival function
    flt d;           // parameter for helping function
    flt u0;          // cut-off for no giving
    int mxtries;     // max number of tries to find donor
    vph_type memb;   // phenotypes of members of the group
    rand_eng& eng;   // reference to random number engine
    bool dhst;       // whether to collect detailed history stats
    bool hist;       // whether to collect history stats
    vdhs_type dhstat;  // detailed history statistics
    vhs_type hstat;  // history statistics
    vis_type istat;  // individual statistics
};

template<typename PhenType>
Group<PhenType>::Group(int a_N,
    int a_K,
    int a_T,
    int a_repi,
    flt a_rho0,
    flt a_xs,
    flt a_xa,
    flt a_pd,
    flt a_alphp,
    flt a_betd,
    flt a_epsp,
    flt ntld0,
    flt a_y0,
    flt a_mu0,
    flt a_a,
    flt a_b,
    flt a_c,
    flt a_d,
    flt a_u0,
    int a_mxtries,
    const vph_type& a_memb,
    rand_eng& a_eng,
    bool a_dhst,
    bool a_hist) :
    N{a_N},
    K{a_K},
    T{a_T},
    repi{a_repi},
    rho0{a_rho0},
    xs{a_xs},
    xa{a_xa},
    pd{a_pd},
    alphp{a_alphp},
    betd{a_betd},
    epsp{a_epsp},
    ntld(K, ntld0),
    y0{a_y0},
    mu0{a_mu0},
    a{a_a},
    b{a_b},
    c{a_c},
    d{a_d},
    u0{a_u0},
    mxtries{a_mxtries},
    memb{a_memb},
    eng{a_eng},
    dhst{a_dhst},
    hist{a_hist}
{
    if (dhst) {
        dhstat.reserve(N*T);
    }
    if (hist) {
        hstat.reserve(T);
        istat.reserve(N);
    }
}

template<typename PhenType>
void Group<PhenType>::Interact()
{
    rand_uni uni(0, 1);
    rand_int rint(0, K - 1);
    // run through periods (days) in this interval of periods
    for (int tstep = 0; tstep < T; ++tstep) {
        // initialize, set resource states for current period, and store
        // indices of living group members in inds
        itype inds;        
        for (int i = 0; i < N; ++i) {
            phen_type& phi = memb[i];
            if (phi.alive) {
                inds.push_back(i);
                phi.u = 0;
                // resource state
                phi.z = (uni(eng) < phi.ps) ? 1 : 0;
                // starting resource balance
                phi.zet = 0.5*phi.q; // assume starting zet is 50% of q
            }
        }
        int Ncurr = inds.size();     
        // distribute individuals over places
        itype np(K, 0); // count per place
        for (int ii = 0; ii < Ncurr; ++ii) {
            int i = inds[ii];
            phen_type& phi = memb[i];
            int kp = 0;
            if (uni(eng) < epsp) {
                // select random place
                kp = rint(eng);
            } else {
                // select according to preference
                ftype wpl(K);
                for (int k = 0; k < K; ++k) {
                    wpl[k] = std::exp(phi.bet*phi.yhat[k])/
                        (1 + std::exp(betd*(phi.nhat[k] - ntld[k])));
                }
                rand_discr dpl(wpl.begin(), wpl.end());
                kp = dpl(eng);
            }
            ++np[kp];
            phi.kp = kp;
        }
        // run through all pairs of alive individuals and update associations
        for (int ii = 0; ii < Ncurr; ++ii) {
            int i = inds[ii];
            phen_type& phi = memb[i];
            // update count for individual's chosen place
            int kp = phi.kp;
            phi.nhat[kp] += alphp*(np[kp] - phi.nhat[kp]);
            ss_type& ssi = phi.ss;
            for (int jj = ii + 1; jj < Ncurr; ++jj) {
                int j = inds[jj];
                phen_type& phj = memb[j];
                flt chij = (phj.kp == kp) ? 1 : 0;
                ssi.rho[j] += alphp*(chij - ssi.rho[j]);
                // put small association values equal to zero
                if (ssi.rho[j] < 0.01) { 
                    ssi.rho[j] = 0;
                }
                ss_type& ssj = phj.ss;
                ssj.rho[i] += alphp*(chij - ssj.rho[i]);
                if (ssj.rho[i] < 0.01) {
                    ssj.rho[i] = 0;
                }
            }
            // update i's bond strength yhat for place k
            for (int jj = 0; jj < Ncurr; ++jj) {
                int j = inds[jj];
                phen_type& phj = memb[j];
                if (j != i && phj.kp == kp) {
                    // j is another individual in the same place as i
                    phi.yhat[kp] += alphp*
                        (ssi.rho[j]*ssi.x[j] - phi.yhat[kp]);
                }
            }
        }
        // randomly shuffle indices (to avoid order effects)
        std::shuffle(inds.begin(), inds.end(), eng);
        // run through individuals and check requests for help
        if (Ncurr > 1) {
            for (int ii = 0; ii < Ncurr; ++ii) {
                int i = inds[ii];
                phen_type& phi = memb[i];
                int zi = phi.z;
                // individual i only requests help when zi = 0 
                if (zi == 0) {
                    ++phi.nz0; // count times when z = 0 for i
                    ss_type& ssi = phi.ss;
                    // select action j (request help from j)
                    int kp = phi.kp;
                    itype indsci; // available individuals in place kp
                    ftype xvec;
                    for (int jj = 0; jj < Ncurr; ++jj) {
                        int j = inds[jj];
                        phen_type& phj = memb[j];
                        // check if j is available at kp
                        if (j != i && phj.kp == kp) {
                            // i's estimate of state z of j
                            int zij = (uni(eng) < pd) ? phj.z : 1 - phj.z;
                            // only consider j with zij = 1 as donors
                            if (zij == 1) {
                                // add individual to indsci
                                indsci.push_back(j);
                                // if this is first time asking help from j,
                                // add the value newpi to the rho*x bond
                                // strength
                                flt newpi = (ssi.nr[j] == 0) ? phi.pn0 : 0;
                                xvec.push_back(std::exp(ssi.rho[j]*ssi.x[j] + 
                                    newpi));
                            }
                        }
                    }
                    int nchoices = indsci.size();
                    if (nchoices > 0) {
                        int ntries = std::min(mxtries, nchoices);
                        for (int tr = 0; tr < ntries; ++ tr) {
                            rand_discr dhlp(xvec.begin(), xvec.end());
                            int jj = dhlp(eng);
                            int j = indsci[jj]; // ask help from j
                            // individual i requests help from j and
                            // individual j provides help to i
                            phen_type& phj = memb[j];
                            ss_type& ssj = phj.ss;
                            // values for social states at beginning of
                            // interaction
                            flt rhoji = ssj.rho[i];
                            flt xji = ssj.x[i];
                            flt hxj = phj.hx;
                            int zj = phj.z;
                            flt zetj = phj.zet;
                            flt hswji = phj.hs*(ssj.rec[i] - ssj.don[i]);
                            flt rhoxji = rhoji*xji;
                            // if this is first time donating help to i, add
                            // the value newpj to the rho*x bond strength
                            flt newpj = (ssj.nd[i] == 0) ? (phj.pn0) : 0;
                            flt uji = hlp(rhoxji + newpj, zj, zetj, hxj, hswji);
                            if (uji < u0) {
                                // let small amounts correspond to not giving
                                // any help
                                uji = 0;
                            }
                            phj.u = uji;       // help provided by j to i
                            if (dhst && phi.gnum == 0) {
                                // only collect detailed data for one group
                                // (assume, the first group)
                                Add_dhstat(tstep, i, j);
                            }
                            phi.zet += uji;    // adjust helping balance
                            phi.rect += uji;   // increase total received
                            ssi.rec[j] += uji; // increase amount received
                            // check if this is a new partnership
                            if (ssi.nr[j] == 0 && ssi.nd[j] == 0) {
                                ++phi.nc0;
                                ssi.t0[j] = phi.age;
                            }
                            ++ssi.nr[j];  // increase count of help from j
                            ++phi.nrh;
                            // adjust helping balance
                            phj.zet -= uji;
                            phj.dont += uji;    // increase total donated
                            ssj.don[i] += uji;  // increase amount donated
                            // check if this is a new partnership
                            if (ssj.nr[i] == 0 && ssj.nd[i] == 0) {
                                ++phj.nc0;
                                ssj.t0[i] = phj.age;
                            }
                            ++ssj.nd[i];  // increase count of help to i
                            ++phj.ndh;
                            // perform updates of x (inspired by
                            // Rescorla-Wagner)
                            flt dij = (xa - ssi.x[j])/(xa - xs);
                            ssi.x[j] += phi.alph*dij*uji;
                            flt dji = (xa - ssj.x[i])/(xa - xs);
                            ssj.x[i] += phj.alph*dji*uji;
                            // prepare to find next and different j
                            xvec[jj] = 0;
                        }
                    }
                }
            }
        }
        // take into account survival to next period
        for (int i = 0; i < N; ++i) {
            phen_type& phi = memb[i];
            if (phi.alive) {
                ++phi.age; // increase age
                if (uni(eng) < mu(phi.z, phi.zet)) {
                    // the individual dies
                    phi.alive = false;
                    // update other group members
                    for (int j = 0; j < N; ++j) {
                        phen_type& phj = memb[j];
                        ss_type& ssj = phj.ss;
                        if (phj.alive) {
                            // update current partners
                            if (ssj.nr[i] + ssj.nd[i] > 0) {
                                if (hist) {
                                    Add_hstat(i, j, 0);
                                    Add_hstat(j, i, 1);
                                }
                                ssj.x[i] = xs;
                                ssj.don[i] = 0;
                                ssj.rec[i] = 0;
                                ssj.nd[i] = 0;
                                ssj.nr[i] = 0;
                                ssj.t0[i] = 0;
                                ++phj.nc1;
                                ++phi.nc1;
                            }
                        }
                        ssj.rho[i] = rho0;
                    }
                    if (hist) {
                        Add_istat(tstep, i);
                    }
                }
            }
        }
        // all set for next period
    }
}

// Interaction in group without individual recognition: this is implemented by
// using elements like, for individual i, phi.ss.x[i] to represent all other
// group members j (these elements are zero in the version with individual
// recognition), and assuming that all rho are equal to 1
template<typename PhenType>
void Group<PhenType>::InteractNIR()
{
    rand_uni uni(0, 1);
    rand_int rint(0, K - 1);
    // run through periods (days) in this interval of periods
    for (int tstep = 0; tstep < T; ++tstep) {
        // initialize, set resource states for current period, and store
        // indices of living group members in inds
        itype inds;        
        for (int i = 0; i < N; ++i) {
            phen_type& phi = memb[i];
            if (phi.alive) {
                inds.push_back(i);
                phi.u = 0;
                // resource state
                phi.z = (uni(eng) < phi.ps) ? 1 : 0;
                // starting resource balance
                phi.zet = 0.5*phi.q; // assume starting zet is 50% of q
            }
        }
        int Ncurr = inds.size();     
        // distribute individuals over places
        itype np(K, 0); // count per place
        for (int ii = 0; ii < Ncurr; ++ii) {
            int i = inds[ii];
            phen_type& phi = memb[i];
            int kp = 0;
            if (uni(eng) < epsp) {
                // select random place
                kp = rint(eng);
            } else {
                // select according to preference
                ftype wpl(K);
                for (int k = 0; k < K; ++k) {
                    wpl[k] = std::exp(phi.bet*phi.yhat[k])/
                        (1 + std::exp(betd*(phi.nhat[k] - ntld[k])));
                }
                rand_discr dpl(wpl.begin(), wpl.end());
                kp = dpl(eng);
            }
            ++np[kp];
            phi.kp = kp;
        }
        // Note: without individual recognition, there is no update of
        // associations (we can assume all rho equal to one)

        // randomly shuffle indices (to avoid order effects)
        std::shuffle(inds.begin(), inds.end(), eng);
        // run through individuals and check requests for help
        if (Ncurr > 1) {
            for (int ii = 0; ii < Ncurr; ++ii) {
                int i = inds[ii];
                phen_type& phi = memb[i];
                int zi = phi.z;
                // individual i only requests help when zi = 0 
                if (zi == 0) {
                    ++phi.nz0; // count times when z = 0 for i
                    ss_type& ssi = phi.ss;
                    // select action j (request help from j)
                    int kp = phi.kp;
                    itype indsci; // available individuals in place kp
                    ftype xvec;
                    for (int jj = 0; jj < Ncurr; ++jj) {
                        int j = inds[jj];
                        phen_type& phj = memb[j];
                        // check if j is available at kp
                        if (j != i && phj.kp == kp) {
                            // i's estimate of state z of j
                            int zij = (uni(eng) < pd) ? phj.z : 1 - phj.z;
                            // only consider j with zij = 1 as donors
                            if (zij == 1) {
                                // add individual to indsci
                                indsci.push_back(j);
                                // Note: without individual recognition, there
                                // is no taking into account if i has previously
                                // asked help from j
                            }
                        }
                    }
                    int nchoices = indsci.size();
                    if (nchoices > 0) {
                        rand_int dhlp(0, nchoices - 1);
                        int ntries = std::min(mxtries, nchoices);
                        for (int tr = 0; tr < ntries; ++ tr) {
                            int jj = dhlp(eng);
                            int j = indsci[jj]; // ask help from j
                            // individual i requests help from j and
                            // individual j provides help to i
                            phen_type& phj = memb[j];
                            ss_type& ssj = phj.ss;
                            // values for social states at beginning of
                            // interaction
                            flt xji = ssj.x[j];     // no ind rec
                            flt hxj = phj.hx;
                            int zj = phj.z;
                            flt zetj = phj.zet;
                            flt hswji = phj.hs*(ssj.rec[j] - ssj.don[j]);
                            // Note: without individual recognition, there is no
                            // taking into account if j has previously donated
                            // to j, and rho assumed to be 1
                            flt uji = hlp(xji, zj, zetj, hxj, hswji);
                            if (uji < u0) {
                                // let small amounts correspond to not giving
                                // any help
                                uji = 0;
                            }
                            phj.u = uji;       // help provided by j to i
                            if (dhst) {
                                Add_dhstat(tstep, i, j);
                            }
                            phi.zet += uji;    // adjust helping balance
                            phi.rect += uji;   // increase total received
                            ssi.rec[i] += uji; // increase amount received
                            ++ssi.nr[i];  // increase count of help 
                            ++phi.nrh;
                            // adjust helping balance
                            phj.zet -= uji;
                            phj.dont += uji;    // increase total donated
                            ssj.don[j] += uji;  // increase amount donated
                            ++ssj.nd[j];  // increase count of help
                            ++phj.ndh;
                            // perform updates of x (inspired by
                            // Rescorla-Wagner)
                            flt dij = (xa - ssi.x[i])/(xa - xs);
                            ssi.x[i] += phi.alph*dij*uji;
                            flt dji = (xa - ssj.x[j])/(xa - xs);
                            ssj.x[j] += phj.alph*dji*uji;
                        }
                    }
                }
            }
        }
        // take into account survival to next period
        for (int i = 0; i < N; ++i) {
            phen_type& phi = memb[i];
            if (phi.alive) {
                ++phi.age; // increase age
                if (uni(eng) < mu(phi.z, phi.zet)) {
                    // the individual dies
                    phi.alive = false;
                    if (hist) {
                        Add_istat(tstep, i);
                    }
                }
            }
        }
        // all set for next period
    }
}

template<typename PhenType>
void Group<PhenType>::Add_dhstat(int tstep, int i, int j)
{
    phen_type& phi = memb[i];
    phen_type& phj = memb[j];
    ss_type& ssi = phi.ss;
    ss_type& ssj = phj.ss;
    dhstat_type st;
    st.gnum = phi.gnum;
    st.kp = phi.kp;
    st.IDi = phi.ID;
    st.i = phi.inum;
    st.IDj = phj.ID;
    st.j = phj.inum;
    st.agei = phi.age;
    st.agej = phj.age;
    st.zi = phi.z;
    st.zj = phj.z;
    st.tstep = tstep;
    st.repi = repi;
    st.qi = phi.q;
    st.qj = phj.q;
    st.uji = phj.u;
    st.zeti = phi.zet;
    st.zetj = phj.zet;
    st.rhoij = ssi.rho[j];
    st.xij = ssi.x[j];
    st.xji = ssj.x[i];
    st.donij = ssi.don[j];
    st.recij = ssi.rec[j];
    st.ndij = ssi.nd[j];
    st.nrij = ssi.nr[j];
    dhstat.push_back(st);
}

template<typename PhenType>
void Group<PhenType>::Add_hstat(int i, int j, int End)
{
    phen_type& phi = memb[i];
    phen_type& phj = memb[j];
    ss_type& ssi = phi.ss;
    ss_type& ssj = phj.ss;
    hstat_type st;
    st.gnum = phi.gnum;
    st.IDi = phi.ID;
    st.i = phi.inum;
    st.IDj = phj.ID;
    st.j = phj.inum;
    st.t0 = ssi.t0[j];
    st.t1 = phi.age;
    st.End = End;
    st.qi = phi.q;
    st.qj = phj.q;
    st.rhoij = ssi.rho[j];
    st.xij = ssi.x[j];
    st.donij = ssi.don[j];
    st.recij = ssi.rec[j];
    st.ndij = ssi.nd[j];
    st.nrij = ssi.nr[j];
    st.kpi = phi.kp;
    st.kpj = phj.kp;
    hstat.push_back(st);
}

template<typename PhenType>
void Group<PhenType>::Add_istat(int tstep, int i)
{
    phen_type& phi = memb[i];
    istat_type st;
    st.bet = phi.bet;
    st.pn0 = phi.pn0;
    st.alph = phi.alph;
    st.ha = phi.ha;
    st.hs = phi.hs;
    st.q = phi.q;
    st.ps = phi.ps;
    st.hx = phi.hx;
    st.z = phi.z;
    st.zet = phi.zet;
    st.dont = phi.dont;
    st.rect = phi.rect;
    st.kp = phi.kp;
    st.age = phi.age;
    st.ndh = phi.ndh;
    st.nz0 = phi.nz0;
    st.nrh = phi.nrh;
    st.nc0 = phi.nc0;
    st.nc1 = phi.nc1;
    st.ID = phi.ID;
    st.inum = phi.inum;
    st.gnum = phi.gnum;
    st.tstep = tstep;
    st.repi = repi;
    istat.push_back(st);
}

#endif // GROUP_HPP
