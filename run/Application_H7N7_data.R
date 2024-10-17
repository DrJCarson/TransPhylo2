# Corresponds to "Application to an outbreak of avian flu H7N7 outbreak in the Netherlands" in the manuscript. Seeds used for simulation (line 38) were 1-4. In all cases inference was continued using the inferTTreecont function for 50 000 + 40 000 + 40 000 + 20 000 MCMC iterations (total chain length 200 000).

library(ape)
library(TransPhylo2)

starttime <- Sys.time()

raw_tree <- read.nexus("beast-H7N7.nex")
split_list <- (strsplit(raw_tree$tip.label, split = "[@]"))
demes <- numeric(length(raw_tree$tip.label))
dLS <- 0

for (i in 1:length(demes)) {
  if (split_list[[i]][3] == "G") demes[i] <- 1
  else if (split_list[[i]][3] == "L") demes[i] <- 2
  else if (split_list[[i]][3] == "C") demes[i] <- 3
  else if (split_list[[i]][3] == "S") demes[i] <- 4
  ds <- as.numeric(split_list[[i]][2])
  if (ds > dLS) dLS <- ds
}

ptree <- ptreeFromPhylo(raw_tree, dateLastSample = dLS)
ptree$demes <- demes
ndemes <- 4

mcmcIterations <- 50000
mcmcthinning <- 6

w.shape <- 3
w.scale <- 5/w.shape
ws.shape <- 10
ws.scale <- 7/ws.shape

dateS <- 50
dateT <- 125

# Seed for inference
seed <- 1

set.seed(seed)

init.r <- rexp(ndemes, rate = 0.5)
init.pi <- rbeta(ndemes, shape1 = 0.5, shape2 = 0.5)
init.rho <- rbeta(ndemes, shape1 = 0.5, shape2 = 0.5)
init.kappa <- rexp(1, rate = 0.5)
init.lambda <- rexp(1, rate = 0.5)

#Case 2
tp_res <- inferTTreemulti(ptree = ptree, mcmcIterations = mcmcIterations, ndemes = ndemes, w.shape = w.shape, w.scale = w.scale, ws.shape = ws.shape, ws.scale = ws.scale, dateS = dateS, dateT = dateT, update.ctree = T, init.r = init.r, init.p = 0.5, init.pi = init.pi, init.rho = init.rho, init.kappa = init.kappa, init.lambda = init.lambda, update.r = T, update.p = F, update.pi = T, update.rho = T, update.kappa = T, update.lambda = T, thinning = mcmcthinning)

endtime <- Sys.time()

