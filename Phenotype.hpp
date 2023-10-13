#ifndef PHENOTYPE_HPP
#define PHENOTYPE_HPP

#include <string>
#include <vector>
#include <ostream>
#include <istream>
#include <cmath>

// The Evo program runs evolutionary simulations
// Copyright (C) 2023  Olof Leimar
// See Readme.md for copyright notice

//************************** struct SocState ******************************

// This struct stores an individual's current associations and bond strengths
// to other group members, together with data on the amounts and number of
// times of donating and receiving help

template<typename flt>
struct SocState {
    using fVec = std::vector<flt>;
    using iVec = std::vector<int>;
    SocState(int N, flt rho0 = 0, flt xs = 0) :
        rho(N, rho0),
        x(N, xs),
        don(N, 0),
        rec(N, 0),
        nd(N, 0),
        nr(N, 0),
        t0(N, 0) {}
    fVec rho;    // association estimates
    fVec x;      // bond strength values
    fVec don;    // amount donated
    fVec rec;    // amount received
    iVec nd;     // number of times donated
    iVec nr;     // number of times received
    iVec t0;     // age for i when helping starts
};


//************************* struct Phenotype ******************************

template<typename GenType>
struct Phenotype {
// public:
    using flt = double;
    using ss_type = SocState<flt>;
    using vp_type = std::vector<flt>;
    using gen_type = GenType;
    using val_type = typename gen_type::val_type;
    Phenotype(int a_N = 0, int a_K = 0, 
        flt rho0 = 0, flt xs = 0, flt nhat0 = 0, 
        const gen_type& gt = gen_type()) :
        N{a_N},
        K{a_K},
        ss(N, rho0, xs),
        nhat(K, nhat0),
        yhat(K, rho0*xs) { Assign(gt); }
    void Assign(const gen_type& gt);
    void Set_inum(int a_inum);
    bool Female() const { return female; }
    // public data members
    int N;        // group size
    int K;        // number of places
    flt bet;      // place choice parameter (genetically determined)
    flt pn0;      // bond plus for new partner (genetically determined)
    flt alph;     // learning rate (genetically determined)
    flt ha;       // helping asymptote (genetically determined)
    flt hs;       // helping sensitivity (genetically determined)
    ss_type ss;   // social state
    vp_type nhat; // estimated density of places
    vp_type yhat; // estimated bond strengths of places
    flt q;        // individual quality, 0 < q < 1
    flt ps;       // probability succeed in foraging
    flt hx;       // max helping amount
    flt u;        // action used (amount of help)
    int z;        // resource state (0/1)
    flt zet;      // net help received current period
    flt dont;     // total help donated
    flt rect;     // total help received
    int kp;       // current place
    int age;      // age (number of periods)
    int ndh;      // number of donations of help
    int nz0;      // number of interactions as recipient (i.e., z == 0)
    int nrh;      // number of times receiving help
    int nc0;      // number of partnerships started
    int nc1;      // number of partnerships ended by death
    int ID;       // 'unique' identifying number
    int inum;     // individual number
    int gnum;     // group number
    bool female;
    bool alive;
};

template<typename GenType>
void Phenotype<GenType>::Assign(const gen_type& gt)
{
    val_type val = gt.Value();
    // assume val is a vector with components corresponding to the genetically
    // determined traits
    bet = val[0];
    pn0 = val[1];
    alph = val[2];
    ha = val[3];
    hs = val[4];
    q = 0.5;
    ps = 1;
    hx = 1;
    u = 0;
    z = 1;
    zet = 0;
    dont = 0;
    rect = 0;
    kp = 0;
    age = 0;
    ndh = 0;
    nz0 = 0;
    nrh = 0;
    nc0 = 0;
    nc1 = 0;
    ID = 0;
    inum = 0;
    gnum = 0;
    female = true;
    alive = true;
}

template<typename GenType>
void Phenotype<GenType>::Set_inum(int a_inum)
{
    inum = a_inum;
}

#endif // PHENOTYPE_HPP
