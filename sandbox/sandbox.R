# ###############################################
# Some test with log-rnorm dist
# ###############################################
rm(list=ls())

set.seed(124)
r0 <- 1; s2 <- 1
r <- rlnorm(1e5, meanlog=log(r0) - s2/2, sdlog=sqrt(s2))
mean(r) # [1] 1.000487 Nice! Very close to one!

# ###############################################
# Test tree simulation with strict clock
# ###############################################
rm(list=ls())

require(ape)
data(pri10s)

set.seed(12377)
relaxed.tree(pri10s, r=1, model="hi") # failed as it should! :-)
tt <- relaxed.tree(pri10s, model="clk", r=0.1)
plot(tt)
plot(tt$edge.length, pri10s$edge.length); abline(0, 1)

# ###############################################
# Test tree simulation with iln clock
# ###############################################
rm(list=ls())

require(ape)
data(pri10s)

set.seed(23451)
tt <- relaxed.tree(pri10s, model="iln", r=.1, s2=1)
plot(tt)
plot(tt$edge.length, pri10s$edge.length); abline(0, 1)

# ###############################################
# Test tree simulation with gbm_RY07 clock
# ###############################################
rm(list=ls())

require(ape)
data(pri10s)

set.seed(777)
tt <- relaxed.tree(pri10s, model="gbm_RY07", r=.1, s2=1)
plot(tt)
plot(tt$edge.length, pri10s$edge.length); abline(0, 1)

