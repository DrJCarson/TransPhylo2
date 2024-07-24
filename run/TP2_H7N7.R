library(ape)
library(TransPhylo2)
rm(list=ls())
set.seed(0)

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
mcmcIterations <- 30000

w.shape <- 3
w.scale <- 5/w.shape
ws.shape <- 10
ws.scale <- 7/ws.shape
dateS <- 50
dateT <- 125

tp_time <- system.time({
#Case 1
#  tp_res <- inferTTree(ptree = ptree, mcmcIterations = mcmcIterations, ndemes = ndemes, w.shape = w.shape, w.scale = w.scale, ws.shape = ws.shape, ws.scale = ws.scale, dateS = dateS, dateT = dateT, update.ctree = T, update.r = T, update.p = F, update.pi = T, update.rho = T, update.kappa = T, update.lambda = T)
#Case 2
  tp_res <- inferTTreemulti(ptree = ptree, mcmcIterations = mcmcIterations, ndemes = ndemes, w.shape = w.shape, w.scale = w.scale, ws.shape = ws.shape, ws.scale = ws.scale, dateS = dateS, dateT = dateT, update.ctree = T, update.r = T, update.p = F, update.pi = T, update.rho = T, update.kappa = T, update.lambda = T)
})

save.image(file = "H7N7_Res_30k_1.rdata")

plot(tp_res)

par(mfrow=c(2,2))
res <- sapply(tp_res, function(x) x[["rho"]])
for (i in 1:4) plot(res[i,],type='l', xlab = "Iteration", ylab = sprintf("rho%d",i))

par(mfrow=c(2,2))
res <- sapply(tp_res, function(x) x[["off.r"]])
for (i in 1:4) plot(res[i,],type='l', xlab = "Iteration", ylab = sprintf("R%d",i))

par(mfrow=c(2,2))
res <- sapply(tp_res, function(x) x[["pi"]])
for (i in 1:4) plot(res[i,],type='l', xlab = "Iteration", ylab = sprintf("pi%d",i))
