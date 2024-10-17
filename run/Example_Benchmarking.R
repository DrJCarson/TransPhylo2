# Corresponds to "Benchmarking using multiple simulations" in the manuscript. Seeds used for simulation (line 31) were 1-50. In all cases inference was continued using the inferTTreecont function for 5 lots of 50 000 MCMC iterations (total chain length 250 000).

library(ape)
library(TransPhylo2)
library(lhs)

starttime <- Sys.time()

# Number of locations
ndemes <- 2

nSampled <- 250

mcmcIterations <- 50000
mcmcthinning <- 10

set.seed(50)

samples <- create_oalhs(n = 50, k = 8, T, F)

samples[, 1] <- 1 + 5 * samples[, 1]
samples[, 2] <- 1 + 5 * samples[, 2]
samples[, 3] <- 0.1 + 0.9 * samples[, 3]
samples[, 4] <- 0.1 + 0.9 * samples[, 4]
samples[, 5] <- 0.5 + 0.4 * samples[, 5]
samples[, 6] <- 0.5 + 0.4 * samples[, 6]
samples[, 7] <- samples[, 7]
samples[, 8] <- samples[, 8]

# Set seed and parameter values for simulation
seed <- 1
set.seed(seed)

off.r <- samples[seed, 1:2]
off.p <- c(0.5, 0.5)
pi <- samples[seed, 3:4]
rho <- samples[seed, 5:6]
kappa <- samples[seed, 7]
lambda <- samples[seed, 8]

true.parms <- list(off.r = off.r,
                   off.p = off.p,
                   pi = pi,
                   rho = rho,
                   kappa = kappa,
                   lambda = lambda)

sim_ctree = simulateOutbreakmulti(nSampled = nSampled, sec.n = 0, sec.t = 0.25, ndemes = ndemes, rho = rho, off.r = off.r, pi = pi, kappa = kappa, lambda = lambda)

ptree <- extractPTree(sim_ctree)

init.r <- rexp(2, rate = 0.5)
init.pi <- rbeta(2, shape1 = 0.5, shape2 = 0.5)
init.rho <- rbeta(2, shape1 = 0.5, shape2 = 0.5)
init.kappa <- rexp(1, rate = 0.5)
init.lambda <- rexp(1, rate = 0.5)

tp_res_1 <- inferTTreemulti(ptree = ptree, mcmcIterations = mcmcIterations, ndemes = ndemes, dateT = sim_ctree$dateT, update.ctree = T, init.r = init.r, init.p = off.p, init.pi = init.pi, init.rho = init.rho, init.kappa = init.kappa, init.lambda = init.lambda, update.r = T, update.p = F, update.pi = T, update.rho = T, update.kappa = T, update.lambda = T, thinning = mcmcthinning)

endtime <- Sys.time()
