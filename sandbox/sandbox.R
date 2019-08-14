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
