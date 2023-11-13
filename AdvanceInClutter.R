rm(list = ls())

library(spatstat)
library(segmented)
source("AdvanceInClutterFunctions.R")

## Figure 1

#clutter
set.seed(22)
pp <- rpoispp(2, win=owin(c(0,10),c(0,10)))
plot(pp, main = "")
pp
hist(nndist(pp, k = 10), 
     breaks = seq(min(nndist(pp, k = 10)), max(nndist(pp, k = 10)), l = 30),
     main = "", xlab = "Distances")

#feature
set.seed(2)
pp2 <- rpoispp(100, win=owin(c(3,4),c(3,4)))
plot(pp2, main = "")

# superimposition
pp3 <- superimpose(clutter= pp, feature = pp2)
plot(pp3, main = "", pch = c(1, 20))
hist(nndist(pp3, k = 10), 
     breaks = seq(min(nndist(pp3, k = 10)), max(nndist(pp3, k = 10)), l = 30),
     main = "", xlab = "Distances")



## Figure 2


set.seed(22)
pp <- rpoispp(200, win = owin(c(0, 1), c(0, 1)))
set.seed(2)
pp2 <- rpoispp(1200, win=owin(c(0,0.5),c(0,0.5)))
pp3 <- superimpose(clutter = pp, feature = pp2)
plot(pp3, main = "", pch = c(1, 20))
deltas <- vector()
x0 <- 5:35
for(i in x0){
  print(i)
  deltas[i] <- class_net_Engine(X = pp3, K = i)$entropy
}
deltas <- deltas[x0]
z0 <- - x0
out.lm <- lm(deltas ~ 1)
o <- segmented(out.lm, ~ z0)
K <- - round(o$psi[2])
K    
plot(x0, deltas, type = "l", xlab = "K", ylab = "Entropy")
lines(x0, round(o$fitted.values), lty = 2)
abline(v = K, lty = 2, lwd = 2)


## Figures 3 and 4

cc <- class_net(pp3, n_feat = 4, verbose = F)
cc$K_track

par(mfrow = c(2, 2))
plot(cc$class_track[[1]], which.marks = 1, cols = 5:6, main = "Iter 1")
p0 <- ppp(cc$class_track[[2]]$x, cc$class_track[[2]]$y, window = cc$X$window)
p0$marks <- cc$class_track[[2]]$marks
plot(p0, which.marks = 1, cols = 5:6, main = "Iter 2")
p1 <- ppp(cc$class_track[[3]]$x, cc$class_track[[3]]$y, window = cc$X$window)
p1$marks <- cc$class_track[[3]]$marks
plot(p1, which.marks = 1, cols = 5:6, main = "Iter 3")
p2 <- ppp(cc$class_track[[4]]$x, cc$class_track[[4]]$y, window = cc$X$window)
p2$marks <- cc$class_track[[4]]$marks
plot(p2, which.marks = 1, cols = 5:6, main = "Iter 4")


par(mfrow = c(1, 1))
plot(cc$n_entr, cc$deltas_track[[1]], type = "b", ylim = range(cc$deltas_track[[1]],
                                                               cc$deltas_track[[2]],
                                                               cc$deltas_track[[3]],
                                                               cc$deltas_track[[4]]),
     col = 1,
     xlab = "K", ylab = "Entropy")
lines(cc$n_entr, cc$deltas_track[[2]], type = "b", col = 2)
lines(cc$n_entr, cc$deltas_track[[3]], type = "b", col = 3)
lines(cc$n_entr, cc$deltas_track[[4]], type = "b", col = 4)
legend("topright", legend=c("Iter 1", "Iter 2", "Iter 3", "Iter 4"),
       col=1:4, lwd = 1)

c(sum(cc$deltas_track[[1]]), sum(cc$deltas_track[[2]]), sum(cc$deltas_track[[3]]), 
  sum(cc$deltas_track[[4]]))


cc <- class_net(pp3, n_feat = 1, verbose = F)
rates.localdetect(cc)
cc <- class_net(pp3, n_feat = 2, verbose = F)
rates.localdetect(cc)
cc <- class_net(pp3, n_feat = 3, verbose = F)
rates.localdetect(cc)
cc <- class_net(pp3, n_feat = 4, verbose = F)
rates.localdetect(cc)


## Figures 5 and 6 + Simulation scenarios

n_sim <- 100 

sim <- list(l = n_sim)
sim_c <- list(l = n_sim)
sim_f <- list(l = n_sim)
n <- vector(l = n_sim)

nclust <- function(x0, y0, radius, n) {
  return(runifdisc(n, radius, centre=c(x0, y0)))
}

set.seed(22)
for(i in 1:n_sim){
  progressreport(i, n_sim)
  sim_c[[i]] <- rpoispp(200, win=owin(c(0,1),c(0,1))) #sim
  sim_f[[i]] <- rpoispp(800, win=owin(c(0.25,0.5),c(0.25,0.5))) # sim
  # sim_f[[i]] <- rPoissonCluster(7.5, 0.2, nclust, radius=0.2, n=20) # scenario 1
  # sim_f[[i]] <- rPoissonCluster(15, 0.2, nclust, radius=0.2, n=10) # scenario 2
  # sim_f[[i]] <- rpoispp(600, win=owin(c(0,0.5),c(0,0.5))) # scenario 3
  # sim_f[[i]] <- rpoispp(300, win=owin(c(0.25,0.5),c(0.25,0.5))) # scenario 4
  n[i] <- npoints(sim_f[[i]])
  sim[[i]] <- superimpose(noise = sim_c[[i]], feature = sim_f[[i]])
}

mean(n)


cc <- class_net(sim[[1]], n_feat = 4, verbose = F)
cc$K_track

par(mfrow = c(2, 2))
plot(cc$class_track[[1]], which.marks = 1, cols = 5:6, main = "Iter 1")
p0 <- ppp(cc$class_track[[2]]$x, cc$class_track[[2]]$y, window = cc$X$window)
p0$marks <- cc$class_track[[2]]$marks
plot(p0, which.marks = 1, cols = 5:6, main = "Iter 2")
p1 <- ppp(cc$class_track[[3]]$x, cc$class_track[[3]]$y, window = cc$X$window)
p1$marks <- cc$class_track[[3]]$marks
plot(p1, which.marks = 1, cols = 5:6, main = "Iter 3")
p2 <- ppp(cc$class_track[[4]]$x, cc$class_track[[4]]$y, window = cc$X$window)
p2$marks <- cc$class_track[[4]]$marks
plot(p2, which.marks = 1, cols = 5:6, main = "Iter 4")

par(mfrow = c(1, 1))
plot(cc$n_entr, cc$deltas_track[[1]], type = "b", ylim = range(cc$deltas_track[[1]],
                                                               cc$deltas_track[[2]],
                                                               cc$deltas_track[[3]],
                                                               cc$deltas_track[[4]]),
     col = 1,
     xlab = "K", ylab = "Entropy")
lines(cc$n_entr, cc$deltas_track[[2]], type = "b", col = 2)
lines(cc$n_entr, cc$deltas_track[[3]], type = "b", col = 3)
lines(cc$n_entr, cc$deltas_track[[4]], type = "b", col = 4)
legend("topright", legend=c("Iter 1", "Iter 2", "Iter 3", "Iter 4"),
       col=1:4, lwd = 1)

c(sum(cc$deltas_track[[1]]), sum(cc$deltas_track[[2]]), sum(cc$deltas_track[[3]]), 
  sum(cc$deltas_track[[4]]))






plot(sim[[1]], cols = 5:6, main = " ")


res <- matrix(NA, ncol = 39, nrow = n_sim)

for(i in 1:n_sim){
  progressreport(i, n_sim)
  res[i, 1:3] <- as.numeric(rates.localdetect(class_net(sim[[i]], K = 10,
                                                        n_feat = 1, verbose = F)))
  res[i, 4:6] <- as.numeric(rates.localdetect(class_net(sim[[i]], K = 10,
                                                        n_feat = 2, verbose = F)))
  res[i, 7:9] <- as.numeric(rates.localdetect(class_net(sim[[i]], K = 10,
                                                        n_feat = 3, verbose = F)))
  res[i, 10:12] <- as.numeric(rates.localdetect(class_net(sim[[i]], K = 20,
                                                          n_feat = 1, verbose = F)))
  res[i, 13:15] <- as.numeric(rates.localdetect(class_net(sim[[i]], K = 20,
                                                          n_feat = 2, verbose = F)))
  res[i, 16:18] <- as.numeric(rates.localdetect(class_net(sim[[i]], K = 20,
                                                          n_feat = 3, verbose = F)))
  res[i, 19:21] <- as.numeric(rates.localdetect(class_net(sim[[i]], K = 30,
                                                          n_feat = 1, verbose = F)))
  res[i, 22:24] <- as.numeric(rates.localdetect(class_net(sim[[i]], K = 30,
                                                          n_feat = 2, verbose = F)))
  res[i, 25:27] <- as.numeric(rates.localdetect(class_net(sim[[i]], K = 30,
                                                          n_feat = 3, verbose = F)))
  cc <- class_net(sim[[i]], n_feat = 1, verbose = F)
  res[i, 28:31] <- c(as.numeric(rates.localdetect(cc)), cc$K_track)
  cc <- class_net(sim[[i]], n_feat = 2, verbose = F)
  res[i, 32:35] <- c(as.numeric(rates.localdetect(cc)), cc$K_track[2])
  cc <- class_net(sim[[i]], n_feat = 3, verbose = F)
  res[i, 36:39] <- c(as.numeric(rates.localdetect(cc)), cc$K_track[3])
}


round(apply(res[, 1:3], 2, function(x) mean(x, na.rm = T)), 2)
round(apply(res[, 4:6], 2, function(x) mean(x, na.rm = T)), 2)
round(apply(res[, 7:9], 2, function(x) mean(x, na.rm = T)), 2)

round(apply(res[, 10:12], 2, function(x) mean(x, na.rm = T)), 2)
round(apply(res[, 13:15], 2, function(x) mean(x, na.rm = T)), 2)
round(apply(res[, 16:18], 2, function(x) mean(x, na.rm = T)), 2)

round(apply(res[, 19:21], 2, function(x) mean(x, na.rm = T)), 2)
round(apply(res[, 22:24], 2, function(x) mean(x, na.rm = T)), 2)
round(apply(res[, 25:27], 2, function(x) mean(x, na.rm = T)), 2)

round(apply(res[, 28:30], 2, function(x) mean(x, na.rm = T)), 2)
round(apply(res[, 32:34], 2, function(x) mean(x, na.rm = T)), 2)
round(apply(res[, 36:38], 2, function(x) mean(x, na.rm = T)), 2)


## Case studies 


## Murchison gold data

cal <- read.table("california.txt", sep = ",", row.names = NULL, header = T)
str(cal)
range(cal$DateTime)

cal <- cal[cal$DateTime < "1971/10/02 19:38:54.96", ]
str(cal)

ppcal <- ppp(cal$Longitude, cal$Latitude,
             owin(range(cal$Longitude), range(cal$Latitude)))
par(mfrow = c(1, 1))
ppcal
plot(ppcal)
axis(1); axis(2)
w <- owin(poly = list(x = c(-123.5, - 121 , - 119.5, - 121  ),
                      y = c(39, 34.92233, 34.92233, 39)))
plot(w, add = TRUE)
ppcal2 <- ppcal[w]
ppcal2
range(cal$DateTime)

plot(ppcal2)
axis(1); axis(2)
ppcal2



cc <- class_net(ppcal2, n_feat = 4, verbose = F)
cc

cc$K_track

c(sum(cc$deltas_track[[1]]), sum(cc$deltas_track[[2]]), sum(cc$deltas_track[[3]]), 
  sum(cc$deltas_track[[4]]))

par(mfrow = c(1, 1))
plot(cc$class_track[[1]], which.marks = 1, cols = 5:6, main = "")
p0 <- ppp(cc$class_track[[2]]$x, cc$class_track[[2]]$y, window = cc$X$window)
p0$marks <- cc$class_track[[2]]$marks


## Detecting seismic faults

cc <- class_net(murchison[[1]], n_feat = 4, verbose = F)

plot(murchison[[2]], main = "", col = "grey")
plot(murchison[[1]], add = TRUE)

cc

cc$K_track

c(sum(cc$deltas_track[[1]]), sum(cc$deltas_track[[2]]), sum(cc$deltas_track[[3]]), 
  sum(cc$deltas_track[[4]]))


par(mfrow = c(1, 2))
plot(murchison[[2]], main = "Iter 1", col = "grey")
plot(cc$class_track[[1]], which.marks = 1, cols = 5:6, add = TRUE)
p0 <- ppp(cc$class_track[[2]]$x, cc$class_track[[2]]$y, window = cc$X$window)
p0$marks <- cc$class_track[[2]]$marks
plot(murchison[[2]], main = "Iter 2", col = "grey")
plot(p0, which.marks = 1, cols = 5:6, add = TRUE)

