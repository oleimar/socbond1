# Run Evo program and add some results to file,

source("h5_funcs.R")

reps <- 50

for (rep in 1:reps) {

system("../build/Evo Run16.toml")

d1 <- h5_dt("Run16.h5")
av_alive <- with(d1, mean(alive))
av_bet0 <- with(d1, mean(bet0))
sd_bet0 <- with(d1, sd(bet0))
av_pn0 <- with(d1, mean(pn0))
sd_pn0 <- with(d1, sd(pn0))
av_alph <- with(d1, mean(alph))
sd_alph <- with(d1, sd(alph))
av_h0 <- with(d1, mean(h0))
sd_h0 <- with(d1, sd(h0))
av_h1 <- with(d1, mean(h1))
sd_h1 <- with(d1, sd(h1))
av_q <- with(d1, mean(q))
sd_q <- with(d1, sd(q))
av_z <- with(d1, mean(z))
sd_z <- with(d1, sd(z))
av_age <- with(d1, mean(age))
sd_age <- with(d1, sd(age))
av_ndh <- with(d1, mean(ndh))
sd_ndh <- with(d1, sd(ndh))
av_nz0 <- with(d1, mean(nz0))
sd_nz0 <- with(d1, sd(nz0))
av_nrh <- with(d1, mean(nrh))
sd_nrh <- with(d1, sd(nrh))
av_nc0 <- with(d1, mean(nc0))
sd_nc0 <- with(d1, sd(nc0))
av_nc1 <- with(d1, mean(nc1))
sd_nc1 <- with(d1, sd(nc1))
dr <- data.frame(av_alive = av_alive,
                 av_bet0 = av_bet0,
                 av_pn0 = av_pn0,
                 av_alph = av_alph,
                 av_h0 = av_h0,
                 av_h1 = av_h1,
                 av_q = av_q,
                 av_z = av_z,
                 av_age = av_age,
                 av_ndh = av_ndh,
                 av_nz0 = av_nz0,
                 av_nrh = av_nrh,
                 av_nc0 = av_nc0,
                 av_nc1 = av_nc1,
                 sd_bet0 = sd_bet0,
                 sd_pn0 = sd_pn0,
                 sd_alph = sd_alph,
                 sd_h0 = sd_h0,
                 sd_h1 = sd_h1,
                 sd_q = sd_q,
                 sd_z = sd_z,
                 sd_age = sd_age,
                 sd_ndh = sd_ndh,
                 sd_nz0 = sd_nz0,
                 sd_nrh = sd_nrh,
                 sd_nc0 = sd_nc0,
                 sd_nc1 = sd_nc1)

dr1 <- read.delim("Run16_data.tsv")
dr1 <- rbind(dr1, dr)
write.table(dr1, "Run16_data.tsv", quote = FALSE,
            sep = "\t", row.names = FALSE)
}
