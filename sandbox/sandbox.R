# ###############################################
# Tests with multi-furcating trees and gbm
# ###############################################
rm(list=ls())

set.seed(157)
tt <- ape::read.tree(text = "((a:0.5,b:0.5):0.5,c:0.8);")
rtt <- simclock::relaxed.tree(tt, model="gbm", r=1, s2=1.2)
sim.blens <- c(0.3243708, 0.1583958, 0.7181730, 0.6346086)
all.equal(rtt$edge.length, sim.blens) # should be < 1e-7

mtt <- ape::read.tree(text = "((a:0.5,b:0.5,d:0.3):0.5,c:0.8);")
rmtt <- simclock::relaxed.tree(mtt, model="gbm", r=1, s2=1.2)
sim.blens2 <- c(0.4975896, 0.2720016, 0.6193712, 0.4945066, 2.0349189)
all.equal(rmtt$edge.length, sim.blens2) # should be < 1e-7
# above should fail with old code

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

# ###############################################
# Spit tree simulation with clk, iln, and gbm
# ###############################################
rm(list=ls())

require(ape)
data(pri10s)

r <- .1; s2 <- .01
set.seed(777)
tc <- relaxed.tree(pri10s, model="clk", r=r)
ti <- relaxed.tree(pri10s, model="iln", r=r, s2=s2)
tg <- relaxed.tree(pri10s, model="gbm_RY07", r=r, s2=s2)

par(mfrow=c(1,3))
plot(pri10s$edge.length, tc$edge.length); abline(0, r)
plot(pri10s$edge.length, ti$edge.length); abline(0, r)
plot(pri10s$edge.length, tg$edge.length); abline(0, r)

# ################################################
# Test rate path with 26s ladder tree
# ################################################
rm(list=ls())

require(ape)
require(simclock)

tt <- "(((((((((((((((((((((((((a,b),c),d),e),f),g),h),i),j),k),l),m),n),o),p),q),r),s),t),u),v),w),x),y),z);"
t26s = read.tree(text=tt)
t26s$edge.length <- rep(1, 50)
#plot(t16s)

#set.seed(3121)
plot(1:24, ty='n', ylim=c(-5, 5))

ra <- 3; s2 <- .1

for (i in 1:200) {
  # simulate tree:
  tr <- relaxed.tree(t26s, "gbm_RY07", ra, s2)
  # extract simulated rates:
  rs <- c(ra, tr$edge.length / t26s$edge.length)

  # plot internal rate path from root to tip number one:
  lines(log(rs[1:25]), col="gray")
}

abline(h=log(ra), lty=2, col="red")
lines(gbm_RY07q(.025, ra, s2, 0:24, log=TRUE), lty=2)
lines(gbm_RY07q(.5, ra, s2, 0:24, log=TRUE), lty=2)
lines(gbm_RY07q(.975, ra, s2, 0:24, log=TRUE), lty=2)
