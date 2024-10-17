# Corresponds to "Exemplary analysis of a simulated dataset where demes have different offspring distributions and sampling probabilities" in the manuscript. Seeds used for inference (line 31) were 1-4.

library(ape)
library(TransPhylo2)

starttime <- Sys.time()

# Number of locations
ndemes <- 2

nSampled <- 250

mcmcIterations <- 100000
mcmcthinning <- 10

set.seed(1)

off.r <- c(1.2, 2.2)
off.p <- c(0.5, 0.5)
pi <- c(0.4, 0.9)
rho <- c(0.9, 0.7)
kappa <- 0.1
lambda <- 0.2

sim_ctree = simulateOutbreakmulti(nSampled = nSampled, sec.n = 0, sec.t = 0.25, ndemes = ndemes, rho = rho, off.r = off.r, pi = pi, kappa = kappa, lambda = lambda)

ptree <- extractPTree(sim_ctree)
ptree$demes <- ptree$demes[1:nSampled]

# Seed for inference
seed <- 1

set.seed(seed)

init.r <- rexp(2, rate = 0.5)
init.pi <- rbeta(2, shape1 = 0.5, shape2 = 0.5)
init.rho <- rbeta(2, shape1 = 0.5, shape2 = 0.5)
init.kappa <- rexp(1, rate = 0.5)
init.lambda <- rexp(1, rate = 0.5)

tp_res_1 <- inferTTreemulti(ptree = ptree, mcmcIterations = mcmcIterations, ndemes = ndemes, dateT = sim_ctree$dateT, update.ctree = T, init.r = init.r, init.p = off.p, init.pi = init.pi, init.rho = init.rho, init.kappa = init.kappa, init.lambda = init.lambda, update.r = T, update.p = F, update.pi = T, update.rho = T, update.kappa = T, update.lambda = T, thinning = mcmcthinning)

endtime <- Sys.time()
