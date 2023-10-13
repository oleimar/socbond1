# functions to read data from hdf5 file created by simulation program

# return data table for phenotype data
h5_dt <- function(hf_name) {
    require(hdf5r)
    require(data.table)
    f.h5 <- H5File$new(hf_name, mode = "r")
    bet <- f.h5[["bet"]][]
    pn0 <- f.h5[["pn0"]][]
    alph <- f.h5[["alph"]][]
    ha <- f.h5[["ha"]][]
    hs <- f.h5[["hs"]][]
    q <- f.h5[["q"]][]
    ps <- f.h5[["ps"]][]
    hx <- f.h5[["hx"]][]
    u <- f.h5[["u"]][]
    z <- f.h5[["z"]][]
    zet <- f.h5[["zet"]][]
    dont <- f.h5[["dont"]][]
    rect <- f.h5[["rect"]][]
    kp <- f.h5[["kp"]][] + 1
    age <- f.h5[["age"]][]
    ndh <- f.h5[["ndh"]][]
    nz0 <- f.h5[["nz0"]][]
    nrh <- f.h5[["nrh"]][]
    nc0 <- f.h5[["nc0"]][]
    nc1 <- f.h5[["nc1"]][]
    ID <- f.h5[["ID"]][]
    inum <- f.h5[["inum"]][] + 1
    gnum <- f.h5[["gnum"]][] + 1
    female <- f.h5[["female"]][]
    alive <- f.h5[["alive"]][]
    f.h5$close_all()
    data.table(bet = bet, pn0 = pn0, alph = alph, ha = ha, hs = hs, 
               q = q, ps = ps, hx = hx, u = u, z = z, zet = zet, 
               dont = dont, rect = rect, age = age, kp = kp, 
               ndh = ndh, nz0 = nz0, nrh = nrh, nc0 = nc0, nc1 = nc1,
               inum = inum, ID = ID, gnum = gnum, 
               female = female, alive = alive)
}

# return matrix where each row is a gamete value
h5_gam <- function(hf_name) {
    require(hdf5r)
    f.h5 <- H5File$new(hf_name, mode = "r")
    gam <- t(f.h5[["Gam"]][,])
    f.h5$close_all()
    gam
}

# return two-dimensional array (first index is individual) 
h5_rho <- function(hf_name) {
    require(hdf5r)
    f.h5 <- H5File$new(hf_name, mode = "r")
    rho <- t(f.h5[["rho"]][,])
    f.h5$close_all()
    rho
}

# return two-dimensional array (first index is individual) 
h5_x <- function(hf_name) {
    require(hdf5r)
    f.h5 <- H5File$new(hf_name, mode = "r")
    x <- t(f.h5[["x"]][,])
    f.h5$close_all()
    x
}

# return two-dimensional array (first index is individual) 
h5_don <- function(hf_name) {
    require(hdf5r)
    f.h5 <- H5File$new(hf_name, mode = "r")
    don <- t(f.h5[["don"]][,])
    f.h5$close_all()
    don
}

# return two-dimensional array (first index is individual) 
h5_rec <- function(hf_name) {
    require(hdf5r)
    f.h5 <- H5File$new(hf_name, mode = "r")
    rec <- t(f.h5[["rec"]][,])
    f.h5$close_all()
    rec
}

# return two-dimensional array (first index is individual) 
h5_nd <- function(hf_name) {
    require(hdf5r)
    f.h5 <- H5File$new(hf_name, mode = "r")
    nd <- t(f.h5[["nd"]][,])
    f.h5$close_all()
    nd
}

# return two-dimensional array (first index is individual) 
h5_nr <- function(hf_name) {
    require(hdf5r)
    f.h5 <- H5File$new(hf_name, mode = "r")
    nr <- t(f.h5[["nr"]][,])
    f.h5$close_all()
    nr
}

# return two-dimensional array (first index is individual) 
h5_t0 <- function(hf_name) {
    require(hdf5r)
    f.h5 <- H5File$new(hf_name, mode = "r")
    t0 <- t(f.h5[["t0"]][,])
    f.h5$close_all()
    t0
}

# return two-dimensional array (first index is individual) 
h5_nhat <- function(hf_name) {
    require(hdf5r)
    f.h5 <- H5File$new(hf_name, mode = "r")
    nhat <- t(f.h5[["nhat"]][,])
    f.h5$close_all()
    nhat
}

# return two-dimensional array (first index is individual) 
h5_yhat <- function(hf_name) {
    require(hdf5r)
    f.h5 <- H5File$new(hf_name, mode = "r")
    yhat <- t(f.h5[["yhat"]][,])
    f.h5$close_all()
    yhat
}

# return data table for detailed history data
h5_ddt <- function(hf_name) {
    require(hdf5r)
    require(data.table)
    f.h5 <- H5File$new(hf_name, mode = "r")
    gnum <- f.h5[["gnum"]][] + 1
    kp <- f.h5[["kp"]][] + 1
    IDi <- f.h5[["IDi"]][]
    i <- f.h5[["i"]][] + 1
    IDj <- f.h5[["IDj"]][]
    j <- f.h5[["j"]][] + 1
    agei <- f.h5[["agei"]][]
    agej <- f.h5[["agej"]][]
    zi <- f.h5[["zi"]][]
    zj <- f.h5[["zj"]][]
    tstep <- f.h5[["tstep"]][] + 1
    repi <- f.h5[["repi"]][] + 1
    qi <- f.h5[["qi"]][]
    qj <- f.h5[["qj"]][]
    uji <- f.h5[["uji"]][]
    zeti <- f.h5[["zeti"]][]
    zetj <- f.h5[["zetj"]][]
    rhoij <- f.h5[["rhoij"]][]
    xij <- f.h5[["xij"]][]
    xji <- f.h5[["xji"]][]
    recij <- f.h5[["recij"]][]
    donij <- f.h5[["donij"]][]
    nrij <- f.h5[["nrij"]][]
    ndij <- f.h5[["ndij"]][]

    f.h5$close_all()
    data.table(gnum = gnum, kp = kp, IDi = IDi, i = i, IDj = IDj, j = j, 
               agei = agei, agej = agej, zi = zi, zj = zj, 
               tstep = tstep, repi = repi, qi = qi, qj = qj, uji = uji,
               zeti = zeti, zetj = zetj, rhoij = rhoij, xij = xij, xji = xji,
               recij = recij, donij = donij, nrij = nrij, ndij = ndij)
}

# return data table for partnership history data
h5_hdt <- function(hf_name) {
    require(hdf5r)
    require(data.table)
    f.h5 <- H5File$new(hf_name, mode = "r")
    gnum <- f.h5[["gnum"]][] + 1
    IDi <- f.h5[["IDi"]][]
    i <- f.h5[["i"]][] + 1
    IDj <- f.h5[["IDj"]][]
    j <- f.h5[["j"]][] + 1
    t0 <- f.h5[["t0"]][]
    t1 <- f.h5[["t1"]][]
    End <- f.h5[["End"]][]
    qi <- f.h5[["qi"]][]
    qj <- f.h5[["qj"]][]
    rhoij <- f.h5[["rhoij"]][]
    xij <- f.h5[["xij"]][]
    donij <- f.h5[["donij"]][]
    recij <- f.h5[["recij"]][]
    ndij <- f.h5[["ndij"]][]
    nrij <- f.h5[["nrij"]][]
    kpi <- f.h5[["kpi"]][] + 1
    kpj <- f.h5[["kpj"]][] + 1

    f.h5$close_all()
    data.table(gnum = gnum, IDi = IDi, i = i, IDj = IDj, j = j, 
               t0 = t0, t1 = t1, dur = t1 - t0, 
               End = End, qi = qi, qj = qj, rhoij = rhoij, xij = xij, 
               donij = donij, recij = recij, ndij = ndij, nrij = nrij, 
               kpi = kpi, kpj = kpj)
}

# return data table for individual data
h5_idt <- function(hf_name) {
    require(hdf5r)
    require(data.table)
    f.h5 <- H5File$new(hf_name, mode = "r")
    bet <- f.h5[["bet"]][]
    pn0 <- f.h5[["pn0"]][]
    alph <- f.h5[["alph"]][]
    ha <- f.h5[["ha"]][]
    hs <- f.h5[["hs"]][]
    q <- f.h5[["q"]][]
    ps <- f.h5[["ps"]][]
    hx <- f.h5[["hx"]][]
    z <- f.h5[["z"]][]
    zet <- f.h5[["zet"]][]
    dont <- f.h5[["dont"]][]
    rect <- f.h5[["rect"]][]
    kp <- f.h5[["kp"]][] + 1
    age <- f.h5[["age"]][]
    ndh <- f.h5[["ndh"]][]
    nz0 <- f.h5[["nz0"]][]
    nrh <- f.h5[["nrh"]][]
    nc0 <- f.h5[["nc0"]][]
    nc1 <- f.h5[["nc1"]][]
    ID <- f.h5[["ID"]][]
    inum <- f.h5[["inum"]][] + 1
    gnum <- f.h5[["gnum"]][] + 1
    tstep <- f.h5[["tstep"]][] + 1
    repi <- f.h5[["repi"]][] + 1
    f.h5$close_all()
    data.table(bet = bet, pn0 = pn0, alph = alph, ha = ha, hs = hs, q = q, 
               ps = ps, hx = hx, z = z, zet = zet, dont = dont, rect = rect, 
               kp = kp, age = age, ndh = ndh, nz0 = nz0, nrh = nrh, 
               nc0 = nc0, nc1 = nc1,
               inum = inum, ID = ID, gnum = gnum, tstep = tstep, repi = repi)
}
