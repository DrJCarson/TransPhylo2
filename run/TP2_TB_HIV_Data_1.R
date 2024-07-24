library(ape)
library(TransPhylo2)

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

mcmcIterations <- 30000

w.shape <- 1.3
w.scale <- 1 / 0.3

ws.shape <- 1.1
ws.scale <- 1 / 0.4

tp_time <- system.time({

  tp_res <- inferTTreemulti(ptree = ptree, mcmcIterations = mcmcIterations, ndemes = ndemes, w.shape = w.shape, w.scale = w.scale, ws.shape = ws.shape, ws.scale = ws.scale, dateS = lubridate::decimal_date(as.Date('1996-10-01')), dateT = lubridate::decimal_date(as.Date('2009-12-01')), update.ctree = T, update.r = T, update.p = F, update.pi = T, update.rho = T, update.kappa = T, update.lambda = T)

})

save.image(file = "TB_HIV_Res_30k_1.rdata")





if (0) {

  pi_res <- sapply(tp_res, function(x) x[["pi"]])

  pi1 <- pi_res[1, ]
  pi2 <- pi_res[2, ]

  par(mfrow = c(2, 1), bty = "n", mar = c(5, 5, 1, 1))

  plot(pi1, type = "l", ylim = c(0, 1), xlab = "Iteration", ylab = "pi negative")
  plot(pi2, type = "l", ylim = c(0, 1), xlab = "Iteration", ylab = "pi positive")

  dev.copy2pdf(file = "pi.pdf")




  r_res <- sapply(tp_res, function(x) x[["off.r"]])

  r1 <- r_res[1, ]
  r2 <- r_res[2, ]

  par(mfrow = c(2, 1), bty = "n", mar = c(5, 5, 1, 1))

  plot(r1, type = "l", ylim = c(0.5, 3.1), xlab = "Iteration", ylab = "R negative")
  plot(r2, type = "l", ylim = c(0.5, 3.1), xlab = "Iteration", ylab = "R positive")

  dev.copy2pdf(file = "r.pdf")




  rho_res <- sapply(tp_res, function(x) x[["rho"]])

  rho1 <- rho_res[1, ]
  rho2 <- rho_res[2, ]

  par(mfrow = c(2, 1), bty = "n", mar = c(5, 5, 1, 1))

  plot(rho1, type = "l", ylim = c(0.4, 1), xlab = "Iteration", ylab = "rho negative")
  plot(rho2, type = "l", ylim = c(0.4, 1), xlab = "Iteration", ylab = "rho positive")

  dev.copy2pdf(file = "rho.pdf")






  kappa <- sapply(tp_res, function(x) x[["kappa"]])

  lambda <- sapply(tp_res, function(x) x[["lambda"]])


  par(mfrow = c(2, 1), bty = "n", mar = c(5, 5, 1, 1))

  plot(kappa, type = "l", ylim = c(0, 10), xlab = "Iteration", ylab = "kappa")
  plot(lambda, type = "l", ylim = c(0, 2), xlab = "Iteration", ylab = "lambda")

  dev.copy2pdf(file = "coa.pdf")


}



