# Corresponds to "Application to an outbreak of tuberculosis in Argentina" in the manuscript. Seeds used for simulation (line 55) were 1-4. In all cases inference was continued using the inferTTreecont function for 2 lots of 60 000 MCMC iterations (total chain length 200 000).

library(ape)
library(TransPhylo2)

starttime <- Sys.time()

raw_tree <- read.nexus("eldholm2016.nex")

split_list <- (strsplit(raw_tree$tip.label, split = "[_]"))

demes <- numeric(length(raw_tree$tip.label))

dLS <- 0

for (i in 1:length(demes)) {

  if (split_list[[i]][4] == "neg") {

    demes[i] <- 1

  } else if (split_list[[i]][4] == "pos") {

    demes[i] <- 2

  }

  ds <- as.numeric(split_list[[i]][3])

  if (ds > dLS) {

    dLS <- ds

  }

}

ptree <- ptreeFromPhylo(raw_tree, dateLastSample = dLS)

ptree$demes <- demes


ndemes <- 2

mcmcIterations <- 80000
mcmcthinning <- 8

w.shape <- 1.3
w.scale <- 1 / 0.3

ws.shape <- 1.1
ws.scale <- 1 / 0.4

#Seed for inference
seed <- 1

set.seed(seed)

init.r <- rexp(ndemes, rate = 0.5)
init.pi <- rbeta(ndemes, shape1 = 0.5, shape2 = 0.5)
init.rho <- rbeta(ndemes, shape1 = 0.5, shape2 = 0.5)
init.kappa <- rexp(1, rate = 0.5)
init.lambda <- rexp(1, rate = 0.5)

tp_res <- inferTTreemulti(ptree = ptree, mcmcIterations = mcmcIterations, ndemes = ndemes, w.shape = w.shape, w.scale = w.scale, ws.shape = ws.shape, ws.scale = ws.scale, dateS = 1996.749, dateT = 2009.915, update.ctree = T, init.r = init.r, init.p = 0.5, init.pi = init.pi, init.rho = init.rho, init.kappa = init.kappa, init.lambda = init.lambda, update.r = T, update.p = F, update.pi = T, update.rho = T, update.kappa = T, update.lambda = T, thinning = mcmcthinning)

endtime <- Sys.time()

