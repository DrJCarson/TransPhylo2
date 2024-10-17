# Corresponds to "Exemplary analysis of a simulated dataset where all demes have the same offspring distribution and sampling probabilities" in the manuscript. Seeds used for inference (line 31) were 1-4.

library(ape)
library(TransPhylo2)

starttime <- Sys.time()

# Number of locations
ndemes <- 5

nSampled <- 250

mcmcIterations <- 100000
mcmcthinning <- 10

set.seed(1)

off.r <- 2
off.p <- 0.5
pi <- 0.8
rho <- 0.8
kappa <- 0.1
lambda <- 0.2

sim_ctree = simulateOutbreak(nSampled = nSampled, sec.n = 0, sec.t = 0.25, ndemes = ndemes, rho = rho, off.r = off.r, pi = pi, kappa = kappa, lambda = lambda)

ptree <- extractPTree(sim_ctree)
ptree$demes <- ptree$demes[1:nSampled]

# Seed for inference
seed <- 1

set.seed(seed)

init.r <- rexp(1, rate = 0.5)
init.pi <- rbeta(1, shape1 = 0.5, shape2 = 0.5)
init.rho <- rbeta(1, shape1 = 0.5, shape2 = 0.5)
init.kappa <- rexp(1, rate = 0.5)
init.lambda <- rexp(1, rate = 0.5)

ptree$demes <- NULL

tp_res_2 <- inferTTree(ptree = ptree, mcmcIterations = mcmcIterations, ndemes = 1, dateT = sim_ctree$dateT, update.ctree = T, init.r = init.r, init.p = off.p, init.pi = init.pi, init.rho = 1, init.kappa = init.kappa, init.lambda = init.lambda, update.r = T, update.p = F, update.pi = T, update.rho = F, update.kappa = T, update.lambda = T, thinning = mcmcthinning)

endtime <- Sys.time()
