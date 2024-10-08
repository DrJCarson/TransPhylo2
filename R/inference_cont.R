#' Continues an MCMC chain from a previous set of results
#'
#' @param prevrun resTransPhylo object
#' @param mcmcIterations Number of new MCMC iterations to perform
#' @param thinning MCMC thinning interval between two sampled iterations
#' @param verbose Whether or not to use verbose mode (default is false)
#' @export
inferTTreecont <- function(prevrun,
                           mcmcIterations = 12000,
                           thinning = 1,
                           verbose = F) {

  fst <- prevrun[[1]]$finalstate

  ptree <- fst$ptree

  # Determine vector of primary observation times
  prim_obs_times <- calc_prim_obs(ptree)

  const.pi <- 3.14159265359

  dateS <- fst$dateS
  dateT <- fst$dateT

  r.shape <- fst$r.shape
  r.scale <- fst$r.scale
  p.shape1 <- fst$p.shape1
  p.shape2 <- fst$p.shape2
  pi.shape1 <- fst$pi.shape1
  pi.shape2 <- fst$pi.shape2
  kappa.shape <- fst$kappa.shape
  kappa.scale <- fst$kappa.scale
  lambda.shape <- fst$lambda.shape
  lambda.scale <- fst$lambda.scale
  rho.shape1 <- fst$rho.shape1
  rho.shape2 <- fst$rho.shape2

  w.shape <- fst$w.shape
  w.scale <- fst$w.scale
  ws.shape <- fst$ws.shape
  ws.scale <- fst$ws.scale

  mcmc.tree.updates <- fst$mcmc.tree.updates

  ctree <- fst$ctree

  grid.delta <- fst$grid.delta

  ndemes <- fst$ndemes
  demes.prior <- fst$demes.prior

  ss.a <- 0.234

  update.r <- fst$update.r
  update.p <- fst$update.p
  update.pi <- fst$update.pi
  update.rho <- fst$update.rho
  update.kappa <- fst$update.kappa
  update.lambda <- fst$update.lambda
  update.ctree <- fst$update.ctree

  mcmcn1 <- fst$mcmcIterations
  mcmcn2 <- mcmcIterations

  mcmcntot <- mcmcn1 + mcmcn2

  ss.f <- floor(0.5 * (1:mcmcntot))

  ss.min <- 0.1

  ttree <- extractTTree(ctree)

  thin0 <- floor(mcmcn1 / thinning)
  thin1 <- floor(mcmcntot / thinning)

  record <- vector('list', thin1 - thin0)

  grid <- fst$grid

  if (fst$case == 2) {

    tr.demes <- fst$tr.demes
    pi.demes <- fst$pi.demes
    rho.demes <- fst$rho.demes

    tr.dim <- length(unique(tr.demes))
    pi.dim <- length(unique(pi.demes))
    rho.dim <- length(unique(rho.demes))
    max.dim <- max(c(tr.dim, pi.dim, rho.dim))

    parms.init.tr <- fst$parms.init.tr
    parms.init.coa <- fst$parms.init.coa

    parms.curr.tr <- fst$parms.curr.tr
    parms.curr.coa <- fst$parms.curr.coa

    within_bounds <- numeric(4)

    mcmc.mu.tr <- fst$mcmc.mu.tr
    mcmc.mu.coa <- fst$mcmc.mu.coa
    mcmc.cov.tr <- fst$mcmc.cov.tr
    mcmc.cov.coa <- fst$mcmc.cov.coa

    update.tr <- c(update.r, update.p, update.pi, update.rho)
    update.coa <- c(update.kappa, update.lambda)

    ss.d.tr <- numeric(max.dim)

    for (i in 1:max.dim) {

      ss.d.tr[i] <- as.numeric(update.r) * as.numeric(!is.na(parms.curr.tr["r", i])) +
        as.numeric(update.p) * as.numeric(!is.na(parms.curr.tr["p", i])) +
        as.numeric(update.pi) * as.numeric(!is.na(parms.curr.tr["pi", i])) +
        as.numeric(update.rho) * as.numeric(!is.na(parms.curr.tr["rho", i]))

    }

    ss.d.coa <- sum(as.numeric(update.coa))

    ss.v0.tr <- ss.d.tr
    ss.v0.coa <- ss.d.coa

    ss.c.tr <- numeric(max.dim)

    for (i in 1:max.dim) {

      ss.c.tr[i] <- 2.38 ^ 2 / ss.d.tr[i]

    }

    ss.c.coa <- 2.38 ^ 2 / ss.d.coa

    ss.lamstart <- 1

    ss.lam.tr <- fst$ss.lam.tr
    ss.lam.coa <- fst$ss.lam.coa

    ss.nstart <- 5 / (ss.a * (1 - ss.a))

    ss.A <- -qnorm(ss.a / 2)

    ss.del.tr <- numeric(max.dim)

    for (i in 1:max.dim) {

      ss.del.tr[i] <- (1 - (1 / ss.d.tr[i])) * ((sqrt(2 * const.pi) * exp(ss.A ^ 2 / 2)) / (2 * ss.A)) + (1 / (ss.d.tr[i] * ss.a * (1 - ss.a)))

    }

    ss.del.coa <- (1 - (1 / ss.d.coa)) * ((sqrt(2 * const.pi) * exp(ss.A ^ 2 / 2)) / (2 * ss.A)) + (1 / (ss.d.coa * ss.a * (1 - ss.a)))

    trace.r <- array(dim = c(mcmcntot, tr.dim))
    trace.p <- array(dim = c(mcmcntot, tr.dim))
    trace.pi <- array(dim = c(mcmcntot, pi.dim))
    trace.rho <- array(dim = c(mcmcntot, rho.dim))

    trace.kappa <- numeric(mcmcntot)
    trace.lambda <- numeric(mcmcntot)

    ftr <- prevrun[[1]]$fulltrace

    d <- 1

    if (update.r) {
      for (de in  1:tr.dim) {
        trace.r[1:mcmcn1, de] <- ftr[, d]
        d <- d + 1
      }
    }

    if (update.p) {
      for (de in  1:tr.dim) {
        trace.p[1:mcmcn1, de] <- ftr[, d]
        d <- d + 1
      }
    }

    if (update.pi) {
      for (de in  1:pi.dim) {
        trace.pi[1:mcmcn1, de] <- ftr[, d]
        d <- d + 1
      }
    }

    if (update.rho) {
      for (de in  1:rho.dim) {
        trace.rho[1:mcmcn1, de] <- ftr[, d]
        d <- d + 1
      }
    }

    if (update.kappa) {
      trace.kappa[1:mcmcn1] <- ftr[, d]
      d <- d + 1
    }

    if (update.lambda) {
      trace.lambda[1:mcmcn1] <- ftr[, d]
      d <- d + 1
    }

    int.r <- parms.curr.tr["r", ]
    int.p <- parms.curr.tr["p", ]
    int.pi <- parms.curr.tr["pi", ]
    int.rho <- parms.curr.tr["rho", ]

    ext.r <- int.r[tr.demes]
    ext.p <- int.p[tr.demes]
    ext.pi <- int.pi[pi.demes]

    ext.rho <- int.rho[rho.demes]

    pm <- matrix(0, nrow = ndemes, ncol = ndemes)

    for (i in 1:ndemes) {

      pm[i, ] <- (1 - ext.rho[i]) / (ndemes - 1)

      pm[i, i] <- ext.rho[i]

    }

    fn_list <- num_approx_disc_multi(grid, ext.r, ext.p, ext.pi, w.shape, w.scale,
                                     ws.shape, ws.scale, dateS, dateT, ndemes, pm)

    llt_out <- log_lik_ttree_multiparm(ttree, grid, fn_list, ext.r, ext.p, ext.pi, w.shape, w.scale, ws.shape,
                                       ws.scale, dateS, dateT, grid.delta, ndemes, pm, demes.prior)

    pTTree <- llt_out$loglik
    dyn_L <- llt_out$dyn_L

    pPTree <- log_lik_ptree_given_ctree(ctree, parms.curr.coa[["kappa"]], parms.curr.coa[["lambda"]], NA)


    if (verbose == F) {

      pb <- utils::txtProgressBar(min = 0, max = mcmcIterations, style = 3)

    }

    tree_prop_count <- c("add" = 0,
                         "remove" = 0,
                         "move" = 0)

    tree_acc_count <- c("add" = 0,
                        "remove" = 0,
                        "move" = 0)

    tree_acc_rates <- c("add" = 0,
                        "remove" = 0,
                        "move" = 0)

    .Random.seed <- fst$Random.seed

    for (i in 1:mcmcIterations) {

      i2 <- mcmcn1 + i

      if (update.ctree) {

        if (verbose) {

          message("Proposing ttree update")

        }

        tree_prop_count["add"] <- 0
        tree_prop_count["remove"] <- 0
        tree_prop_count["move"] <- 0

        tree_acc_count["add"] <- 0
        tree_acc_count["remove"] <- 0
        tree_acc_count["move"] <- 0

        for (j in 1:mcmc.tree.updates) {

          u <- runif(1)
          if (u < 1 / 3) {

            proptree <- add_transmission(ctree = ctree)
            prop_type <- 1
            tree_prop_count["add"] <- tree_prop_count["add"] + 1

            host_map <- proptree$host_map

            inv_host_map <- rep(NA, length(host_map) + 1)
            inv_host_map[(1:(length(host_map) + 1)) %in% host_map] <- order(host_map, na.last = NA)

          } else if (u < 2 / 3) {

            proptree <- remove_transmission(ctree = ctree)
            prop_type <- 2
            tree_prop_count["remove"] <- tree_prop_count["remove"] + 1

            host_map <- proptree$host_map

            inv_host_map <- order(host_map, na.last = NA)

          } else {

            proptree <- remove_add(ctree = ctree)
            prop_type <- 3
            tree_prop_count["move"] <- tree_prop_count["move"] + 1

            host_map <- proptree$host_map

            inv_host_map <- order(host_map, na.last = NA)

          }

          if (proptree$is_possible == 1) {

            ctree2 <- proptree$ctree

            ttree2 <- extractTTree(ctree2)


            dyn_L2 <- dyn_L[inv_host_map, , drop = F]

            llt_out2 <- log_lik_ttree_multiparm_part(ttree2, grid, fn_list, ext.r, ext.p, ext.pi, w.shape, w.scale, ws.shape,
                                                     ws.scale, dateS, dateT, grid.delta, ndemes, pm, demes.prior,
                                                     dyn_L2, proptree$prop_hosts)

            pTTree2 <- llt_out2$loglik

            dyn_L2 <- llt_out2$dyn_L

            pPTree_part <- log_lik_ptree_given_ctree(ctree, parms.curr.coa[["kappa"]], parms.curr.coa[["lambda"]], proptree$curr_hosts)

            pPTree_part2 <- log_lik_ptree_given_ctree(ctree2, parms.curr.coa[["kappa"]], parms.curr.coa[["lambda"]], proptree$prop_hosts)



            if (log(runif(1)) < (pTTree2 + pPTree_part2 + proptree$rev_density -
                                 pTTree - pPTree_part - proptree$prop_density)) {

              ctree <- ctree2
              ttree <- ttree2
              dyn_L <- dyn_L2

              pTTree <- pTTree2
              pPTree <- pPTree + pPTree_part2 - pPTree_part

              if (prop_type == 1) {

                tree_acc_count["add"] <- tree_acc_count["add"] + 1

              } else if (prop_type == 2) {

                tree_acc_count["remove"] <- tree_acc_count["remove"] + 1

              } else {

                tree_acc_count["move"] <- tree_acc_count["move"] + 1

              }

              if (grid[length(grid)] > (min(ttree$ttree[, 1]) - 0.5 * 1)) {

                grid <- seq(dateT, min(ttree$ttree[, 1]) - 0.5 * 1 - grid.delta, by = - grid.delta)

                fn_list <- num_approx_disc_multi(grid, ext.r, ext.p, ext.pi, w.shape, w.scale,
                                                 ws.shape, ws.scale, dateS, dateT, ndemes, pm)

              }

            }

          }

        }

        if (i2 %% thinning == 0) {

          tree_acc_rates["add"] <- tree_acc_count["add"] / tree_prop_count["add"]
          tree_acc_rates["remove"] <- tree_acc_count["remove"] / tree_prop_count["remove"]
          tree_acc_rates["move"] <- tree_acc_count["move"] / tree_prop_count["move"]

          record[[(i2 / thinning) - thin0]]$tree.prop.counts <- tree_prop_count
          record[[(i2 / thinning) - thin0]]$tree.acc.rates <- tree_acc_rates

        }

      }

      if (max(ss.d.tr) > 0) {

        for (k in 1:max.dim) {

          if (ss.d.tr[k] > 0) {

            updt_idx <- which(!is.na(parms.curr.tr[, k]) &
                                c(update.r, update.p, update.pi, update.rho))

            parms.prop.tr <- MASS::mvrnorm(1, mu = parms.curr.tr[, k], Sigma = ss.lam.tr[k] * ss.c.tr[k] * mcmc.cov.tr[[k]])

            ss.u.tr <- log(runif(1))

            within_bounds[1] <- parms.prop.tr["r"] > 0
            within_bounds[2] <- parms.prop.tr["p"] > 0 & parms.prop.tr["p"] <= 1
            within_bounds[3] <- parms.prop.tr["pi"] > 0 & parms.prop.tr["pi"] <= 1
            within_bounds[4] <- parms.prop.tr["rho"] > 0 & parms.prop.tr["rho"] <= 1

            if (prod(within_bounds[updt_idx])){

              int.r <- parms.curr.tr["r", ]
              int.r[k] <- parms.prop.tr["r"]

              int.p <- parms.curr.tr["p", ]
              int.p[k] <- parms.prop.tr["p"]

              int.pi <- parms.curr.tr["pi", ]
              int.pi[k] <- parms.prop.tr["pi"]

              int.rho <- parms.curr.tr["rho", ]
              int.rho[k] <- parms.prop.tr["rho"]

              ext.r2 <- int.r[tr.demes]
              ext.p2 <- int.p[tr.demes]
              ext.pi2 <- int.pi[pi.demes]
              ext.rho2 <- int.rho[rho.demes]

              pm2 <- matrix(0, nrow = ndemes, ncol = ndemes)

              for (l in 1:ndemes) {

                pm2[l, ] <- (1 - ext.rho2[l]) / (ndemes - 1)

                pm2[l, l] <- ext.rho2[l]

              }

              fn_list2 <- num_approx_disc_multi(grid, ext.r2, ext.p2, ext.pi2, w.shape, w.scale,
                                                ws.shape, ws.scale, dateS, dateT, ndemes, pm2)

              llt_out2 <- log_lik_ttree_multiparm(ttree, grid, fn_list2, ext.r2, ext.p2, ext.pi2, w.shape, w.scale, ws.shape,
                                                  ws.scale, dateS, dateT, grid.delta, ndemes, pm2, demes.prior)

              pTTree2 <- llt_out2$loglik
              dyn_L2 <- llt_out2$dyn_L

              ss.alpha.tr <- (pTTree2 - pTTree)

              if (k <= tr.dim) {

                ss.alpha.tr <- ss.alpha.tr +
                  (dgamma(parms.prop.tr[["r"]], shape = r.shape[k], scale = r.scale[k], log = T) -
                     dgamma(parms.curr.tr[["r", k]], shape = r.shape[k], scale = r.scale[k], log = T)) +
                  (dbeta(parms.prop.tr[["p"]], shape1 = p.shape1[k], shape2 = p.shape2[k], log = T) -
                     dbeta(parms.curr.tr[["p", k]], shape1 = p.shape1[k], shape2 = p.shape2[k], log = T))

              }

              if (k <= pi.dim) {

                ss.alpha.tr <- ss.alpha.tr +
                  (dbeta(parms.prop.tr[["pi"]], shape1 = pi.shape1[k], shape2 = pi.shape2[k], log = T) -
                     dbeta(parms.curr.tr[["pi", k]], shape1 = pi.shape1[k], shape2 = pi.shape2[k], log = T))

              }

              if (k <= rho.dim) {

                ss.alpha.tr <- ss.alpha.tr +
                  (dbeta(parms.prop.tr[["rho"]], shape1 = rho.shape1[k], shape2 = rho.shape2[k], log = T) -
                     dbeta(parms.curr.tr[["rho", k]], shape1 = rho.shape1[k], shape2 = rho.shape2[k], log = T))

              }

            } else {

              ss.alpha.tr <- -Inf

            }

            if (ss.u.tr < ss.alpha.tr) {

              parms.curr.tr[updt_idx, k] <- parms.prop.tr[updt_idx]

              ext.r <- ext.r2
              ext.p <- ext.p2
              ext.pi <- ext.pi2
              ext.rho <- ext.rho2

              pm <- pm2

              fn_list <- fn_list2
              pTTree <- pTTree2
              dyn_L <- dyn_L2

            }

            if (k <= tr.dim) {

              trace.r[i2, k] <- parms.curr.tr[["r", k]]
              trace.p[i2, k] <- parms.curr.tr[["p", k]]

            }

            if (k <= pi.dim) {

              trace.pi[i2, k] <- parms.curr.tr[["pi", k]]

            }

            if (k <= rho.dim) {

              trace.rho[i2, k] <- parms.curr.tr[["rho", k]]

            }

            if (i2 == 1) {

              mcmc.mu.tr[[k]] <- 0.5 * (parms.init.tr[, k] + parms.curr.tr[, k])

              mcmc.cov.tr[[k]][updt_idx, updt_idx] <- ((1 / (ss.v0.tr[k] + ss.d.tr[k] + 3)) * (parms.init.tr[, k] %*% t(parms.init.tr[, k]) + parms.curr.tr[, k] %*% t(parms.curr.tr[, k]) + (ss.v0.tr[k] + ss.d.tr[k] + 1) * mcmc.cov.tr[[k]] - 2 * mcmc.mu.tr[[k]] %*% t(mcmc.mu.tr[[k]])))[updt_idx, updt_idx]

            } else if (ss.f[i2] == ss.f[i2 - 1]) {

              mcmc.mu.new.tr <- ((i2 - ss.f[i2]) / (i2 - ss.f[i2] + 1)) * mcmc.mu.tr[[k]] +
                (1 / (i2 - ss.f[i2] + 1)) * parms.curr.tr[, k]

              mcmc.cov.tr[[k]][updt_idx, updt_idx] <- ((1 / (i2 - ss.f[i2] + ss.v0.tr[k] + ss.d.tr[k] + 2)) * ((i2 - ss.f[i2] + ss.v0.tr[k] + ss.d.tr[k] + 1) * mcmc.cov.tr[[k]] + parms.curr.tr[, k] %*% t(parms.curr.tr[, k]) + (i2 - ss.f[i2]) * mcmc.mu.tr[[k]] %*% t(mcmc.mu.tr[[k]]) - (i2 - ss.f[i2] + 1) * mcmc.mu.new.tr %*% t(mcmc.mu.new.tr)))[updt_idx, updt_idx]

              mcmc.mu.tr[[k]] <- mcmc.mu.new.tr

            } else {

              rem.el <- ss.f[i2] - 1

              if (rem.el == 0) {

                parms.rem.tr <- parms.init.tr[, k]

              } else {

                parms.rem.tr <- rep(NA, 4)

                if (k <= tr.dim) {

                  parms.rem.tr[1:2] <-  c(r = trace.r[rem.el, k],
                                          p = trace.p[rem.el, k])

                }

                if (k <= pi.dim) {

                  parms.rem.tr[3] <- c(pi = trace.pi[rem.el, k])

                }

                if (k <= rho.dim) {

                  parms.rem.tr[4] <- c(rho = trace.rho[rem.el, k])

                }

              }

              mcmc.mu.new.tr <- mcmc.mu.tr[[k]] + (1 / (i2 - ss.f[i2] + 1)) * (parms.curr.tr[, k] - parms.rem.tr)

              mcmc.cov.tr[[k]][updt_idx, updt_idx] <- (mcmc.cov.tr[[k]] + (1 / (i2 - ss.f[i2] + ss.v0.tr[k] + ss.d.tr[k] + 2)) * (parms.curr.tr[, k] %*% t(parms.curr.tr[, k]) - parms.rem.tr %*% t(parms.rem.tr) + (i2 - ss.f[i2] + 1) * (mcmc.mu.tr[[k]] %*% t(mcmc.mu.tr[[k]]) - mcmc.mu.new.tr %*% t(mcmc.mu.new.tr))))[updt_idx, updt_idx]

              mcmc.mu.tr[[k]] <- mcmc.mu.new.tr

            }

            ss.lam.tr[k] <- max(c(ss.min, ss.lam.tr[k] * exp((ss.del.tr[k] / (ss.nstart + i2)) * (min(c(1, exp(ss.alpha.tr))) - ss.a))))

            update.idx.tr <- updt_idx

            zero.rc <- which(update.idx.tr == 0)

            if (length(zero.rc) > 0) {

              mcmc.cov.tr[[k]][zero.rc, ] <- 0
              mcmc.cov.tr[[k]][, zero.rc] <- 0

            }

          }

        }

      }

      if (ss.d.coa > 0) {

        parms.prop.coa <- MASS::mvrnorm(1, mu = parms.curr.coa, Sigma = ss.lam.coa * ss.c.coa * mcmc.cov.coa)

        ss.u.coa <- log(runif(1))

        if (parms.prop.coa["kappa"] >= 0 & parms.prop.coa["lambda"] >= 0) {

          pPTree2 <- log_lik_ptree_given_ctree(ctree, parms.prop.coa[["kappa"]],
                                               parms.prop.coa[["lambda"]], NA)

          ss.alpha.coa <- (pPTree2 - pPTree) +
            (dgamma(parms.prop.coa[["kappa"]], shape = kappa.shape, scale = kappa.scale, log = T) -
               dgamma(parms.curr.coa[["kappa"]], shape = kappa.shape, scale = kappa.scale, log = T)) +
            (dgamma(parms.prop.coa[["lambda"]], shape = lambda.shape, scale = lambda.scale, log = T) -
               dgamma(parms.curr.coa[["lambda"]], shape = lambda.shape, scale = lambda.scale, log = T))

        } else {

          ss.alpha.coa <- -Inf

        }

        if (ss.u.coa < ss.alpha.coa) {

          parms.curr.coa[which(update.coa)] <- parms.prop.coa[which(update.coa)]

          pPTree <- pPTree2

        }

        trace.kappa[i2] <- parms.curr.coa[["kappa"]]
        trace.lambda[i2] <- parms.curr.coa[["lambda"]]

        if (i2 == 1) {

          mcmc.mu.coa <- 0.5 * (parms.init.coa + parms.curr.coa)

          mcmc.cov.coa <- (1 / (ss.v0.coa + ss.d.coa + 3)) * (parms.init.coa %*% t(parms.init.coa) +
                                                                parms.curr.coa %*% t(parms.curr.coa) +
                                                                (ss.v0.coa + ss.d.coa + 1) * mcmc.cov.coa -
                                                                2 * mcmc.mu.coa %*% t(mcmc.mu.coa))

        } else if (ss.f[i2] == ss.f[i2 - 1]) {

          mcmc.mu.new.coa <- ((i2 - ss.f[i2]) / (i2 - ss.f[i2] + 1)) * mcmc.mu.coa +
            (1 / (i2 - ss.f[i2] + 1)) * parms.curr.coa

          mcmc.cov.coa <- (1 / (i2 - ss.f[i2] + ss.v0.coa + ss.d.coa + 2)) *
            ((i2 - ss.f[i2] + ss.v0.coa + ss.d.coa + 1) * mcmc.cov.coa +
               parms.curr.coa %*% t(parms.curr.coa) +
               (i2 - ss.f[i2]) * mcmc.mu.coa %*% t(mcmc.mu.coa) -
               (i2 - ss.f[i2] + 1) * mcmc.mu.new.coa %*% t(mcmc.mu.new.coa))

          mcmc.mu.coa <- mcmc.mu.new.coa

        } else {

          rem.el <- ss.f[i2] - 1

          if (rem.el == 0) {

            parms.rem.coa <- parms.init.coa

          } else {

            parms.rem.coa <- c(kappa = trace.kappa[rem.el],
                               lambda = trace.lambda[rem.el])

          }

          mcmc.mu.new.coa <- mcmc.mu.coa + (1 / (i2 - ss.f[i2] + 1)) * (parms.curr.coa - parms.rem.coa)

          mcmc.cov.coa <- mcmc.cov.coa + (1 / (i2 - ss.f[i2] + ss.v0.coa + ss.d.coa + 2)) *
            (parms.curr.coa %*% t(parms.curr.coa) - parms.rem.coa %*% t(parms.rem.coa) +
               (i2 - ss.f[i2] + 1) * (mcmc.mu.coa %*% t(mcmc.mu.coa) - mcmc.mu.new.coa %*% t(mcmc.mu.new.coa)))

          mcmc.mu.coa <- mcmc.mu.new.coa

        }

        ss.lam.coa <- max(c(ss.min, ss.lam.coa * exp((ss.del.coa / (ss.nstart + i2)) * (min(c(1, exp(ss.alpha.coa))) - ss.a))))

        update.idx.coa <- c(update.kappa, update.lambda)

        zero.rc <- which(update.idx.coa == 0)

        if (length(zero.rc) > 0) {

          mcmc.cov.coa[zero.rc, ] <- 0
          mcmc.cov.coa[, zero.rc] <- 0

        }

      }


      if (i2 %% thinning == 0) {

        if (verbose == F) {

          utils::setTxtProgressBar(pb, i)

        }

        if (verbose==T) {

          message(sprintf('it = %d, r = %f, off.p = %f, pi = %f, kappa = %f, lambda = %f, prior = %e, likelihood = %e, nind = %d', i, parms.curr.tr["r"], parms.curr.tr["p"], parms.curr.tr["pi"], parms.curr.coa["kappa"], parms.curr.coa["lambda"], pTTree, pPTree, nrow(ttree$ttree)))

        }

        rec_ctree <- trim_root(ctree)

        record[[(i2 / thinning) - thin0]]$ctree <- rec_ctree
        record[[(i2 / thinning) - thin0]]$pTTree <- pTTree
        record[[(i2 / thinning) - thin0]]$pPTree <- pPTree
        record[[(i2 / thinning) - thin0]]$kappa <- parms.curr.coa[["kappa"]]
        record[[(i2 / thinning) - thin0]]$lambda <- parms.curr.coa[["lambda"]]

        record[[(i2 / thinning) - thin0]]$off.r <- parms.curr.tr["r", ]
        record[[(i2 / thinning) - thin0]]$off.p <- parms.curr.tr["p", ]
        record[[(i2 / thinning) - thin0]]$pi <- parms.curr.tr["pi", ]
        record[[(i2 / thinning) - thin0]]$rho <- parms.curr.tr["rho", ]

        record[[(i2 / thinning) - thin0]]$w.shape <- w.shape
        record[[(i2 / thinning) - thin0]]$w.scale <- w.scale
        record[[(i2 / thinning) - thin0]]$ws.shape <- ws.shape
        record[[(i2 / thinning) - thin0]]$ws.scale <- ws.scale
        record[[(i2 / thinning) - thin0]]$mcmc.cov.tr <- mcmc.cov.tr
        record[[(i2 / thinning) - thin0]]$ss.lam.tr <- ss.lam.tr
        record[[(i2 / thinning) - thin0]]$mcmc.cov.coa <- mcmc.cov.coa
        record[[(i2 / thinning) - thin0]]$ss.lam.coa <- ss.lam.coa

        record[[(i2 / thinning) - thin0]]$source <- rec_ctree$ctree[rec_ctree$ctree[which(rec_ctree$ctree[, 4] == 0), 2], 4]
        if (record[[(i2 / thinning) - thin0]]$source <= length(rec_ctree$nam)) {

          record[[(i2 / thinning) - thin0]]$source <- rec_ctree$nam[record[[(i2 / thinning) - thin0]]$source]

        } else {

          record[[(i2 / thinning) - thin0]]$source <- 'Unsampled'

        }

      }

    }#End of main MCMC loop

    mcmcdim <- tr.dim * as.numeric(update.r) + tr.dim * as.numeric(update.p) +
      pi.dim * as.numeric(update.pi) + rho.dim * as.numeric(update.rho) +
      as.numeric(update.kappa) + as.numeric(update.lambda)

    cnames <- character(mcmcdim)

    d <- 1

    if (update.r) {

      for (de in  1:tr.dim) {

        cnames[d] <- paste("r", de, sep = "")

        d <- d + 1

      }

    }

    if (update.p) {

      for (de in  1:tr.dim) {

        cnames[d] <- paste("p", de, sep = "")

        d <- d + 1

      }

    }

    if (update.pi) {

      for (de in  1:pi.dim) {

        cnames[d] <- paste("pi", de, sep = "")

        d <- d + 1

      }

    }

    if (update.rho) {

      for (de in  1:rho.dim) {

        cnames[d] <- paste("rho", de, sep = "")

        d <- d + 1

      }

    }

    if (update.kappa) {

      cnames[d] <- "kappa"

      d <- d + 1

    }

    if (update.lambda) {

      cnames[d] <- "lambda"

      d <- d + 1

    }

    fulltrace <- array(dim = c(mcmcntot, mcmcdim), dimnames = list(NULL, cnames))

    d <- 1

    if (update.r) {

      for (de in  1:tr.dim) {

        fulltrace[, d] <- trace.r[, de]

        d <- d + 1

      }

    }

    if (update.p) {

      for (de in  1:tr.dim) {

        fulltrace[, d] <- trace.p[, de]

        d <- d + 1

      }

    }

    if (update.pi) {

      for (de in  1:pi.dim) {

        fulltrace[, d] <- trace.pi[, de]

        d <- d + 1

      }

    }

    if (update.rho) {

      for (de in  1:rho.dim) {

        fulltrace[, d] <- trace.rho[, de]

        d <- d + 1

      }

    }

    if (update.kappa) {

      fulltrace[, d] <- trace.kappa

      d <- d + 1

    }

    if (update.lambda) {

      fulltrace[, d] <- trace.lambda

      d <- d + 1

    }

    fulltrace <- coda::as.mcmc(fulltrace)

    record[[1]]$fulltrace <- fulltrace

    finalstate <- list()
    finalstate$ptree <- ptree
    finalstate$dateS <- dateS
    finalstate$dateT <- dateT
    finalstate$r.shape <- r.shape
    finalstate$r.scale <- r.scale
    finalstate$p.shape1 <- p.shape1
    finalstate$p.shape2 <- p.shape2
    finalstate$pi.shape1 <- pi.shape1
    finalstate$pi.shape2 <- pi.shape2
    finalstate$kappa.shape <- kappa.shape
    finalstate$kappa.scale <- kappa.scale
    finalstate$lambda.shape <- lambda.shape
    finalstate$lambda.scale <- lambda.scale
    finalstate$rho.shape1 <- rho.shape1
    finalstate$rho.shape2 <- rho.shape2
    finalstate$w.shape <- w.shape
    finalstate$w.scale <- w.scale
    finalstate$ws.shape <- ws.shape
    finalstate$ws.scale <- ws.scale
    finalstate$mcmc.tree.updates <- mcmc.tree.updates
    finalstate$ctree <- ctree
    finalstate$grid.delta <- grid.delta
    finalstate$ndemes <- ndemes
    finalstate$demes.prior <- demes.prior
    finalstate$tr.demes <- tr.demes
    finalstate$pi.demes <- pi.demes
    finalstate$rho.demes <- rho.demes
    finalstate$parms.init.tr <- parms.init.tr
    finalstate$parms.init.coa <- parms.init.coa
    finalstate$parms.curr.tr <- parms.curr.tr
    finalstate$parms.curr.coa <- parms.curr.coa
    finalstate$mcmc.mu.tr <- mcmc.mu.tr
    finalstate$mcmc.mu.coa <- mcmc.mu.coa
    finalstate$mcmc.cov.tr <- mcmc.cov.tr
    finalstate$mcmc.cov.coa <- mcmc.cov.coa
    finalstate$update.r <- update.r
    finalstate$update.p <- update.p
    finalstate$update.pi <- update.pi
    finalstate$update.rho <- update.rho
    finalstate$update.kappa <- update.kappa
    finalstate$update.lambda <- update.lambda
    finalstate$update.ctree <- update.ctree
    finalstate$mcmcIterations <- mcmcntot
    finalstate$ss.lam.tr <- ss.lam.tr
    finalstate$ss.lam.coa <- ss.lam.coa
    finalstate$grid <- grid
    finalstate$Random.seed <- .Random.seed
    finalstate$case <- 2

    record[[1]]$finalstate <- finalstate

    class(record)<-'resTransPhylo'
    return(record)

  } else {

    parms.init.tr <- fst$parms.init.tr
    parms.init.coa <- fst$parms.init.coa
    parms.init.rho <- fst$parms.init.rho

    parms.curr.tr <- fst$parms.curr.tr
    parms.curr.coa <- fst$parms.curr.coa
    parms.curr.rho <- fst$parms.curr.rho

    mcmc.mu.tr <- fst$mcmc.mu.tr
    mcmc.mu.coa <- fst$mcmc.mu.coa
    mcmc.mu.rho <- fst$mcmc.mu.rho
    mcmc.cov.tr <- fst$mcmc.cov.tr
    mcmc.cov.coa <- fst$mcmc.cov.coa
    mcmc.cov.rho <- fst$mcmc.cov.rho

    update.tr <- c(update.r, update.p, update.pi)
    update.coa <- c(update.kappa, update.lambda)

    ss.d.tr <- sum(as.numeric(update.tr))
    ss.d.coa <- sum(as.numeric(update.coa))
    ss.d.rho <- sum(as.numeric(update.rho))

    ss.v0.tr <- ss.d.tr
    ss.v0.coa <- ss.d.coa
    ss.v0.rho <- ss.d.rho

    ss.c.tr <- 2.38 ^ 2 / ss.d.tr
    ss.c.coa <- 2.38 ^ 2 / ss.d.coa
    ss.c.rho <- 2.38 ^ 2 / ss.d.rho

    ss.lamstart <- 1

    ss.lam.tr <- fst$ss.lam.tr
    ss.lam.coa <- fst$ss.lam.coa
    ss.lam.rho <- fst$ss.lam.rho

    ss.nstart <- 5 / (ss.a * (1 - ss.a))

    ss.A <- -qnorm(ss.a / 2)

    ss.del.tr <- (1 - (1 / ss.d.tr)) * ((sqrt(2 * const.pi) * exp(ss.A ^ 2 / 2)) / (2 * ss.A)) + (1 / (ss.d.tr * ss.a * (1 - ss.a)))
    ss.del.coa <- (1 - (1 / ss.d.coa)) * ((sqrt(2 * const.pi) * exp(ss.A ^ 2 / 2)) / (2 * ss.A)) + (1 / (ss.d.coa * ss.a * (1 - ss.a)))
    ss.del.rho <- (1 - (1 / ss.d.rho)) * ((sqrt(2 * const.pi) * exp(ss.A ^ 2 / 2)) / (2 * ss.A)) + (1 / (ss.d.rho * ss.a * (1 - ss.a)))

    trace.r <- numeric(mcmcntot)
    trace.p <- numeric(mcmcntot)
    trace.pi <- numeric(mcmcntot)
    trace.rho <- numeric(mcmcntot)

    trace.kappa <- numeric(mcmcntot)
    trace.lambda <- numeric(mcmcntot)

    ftr <- prevrun[[1]]$fulltrace

    d <- 1

    if (update.r) {
      trace.r[1:mcmcn1] <- ftr[, d]
      d <- d + 1
    }

    if (update.p) {
      trace.p[1:mcmcn1] <- ftr[, d]
      d <- d + 1
    }

    if (update.pi) {
      trace.pi[1:mcmcn1] <- ftr[, d]
      d <- d + 1
    }

    if (update.rho) {
      trace.rho[1:mcmcn1] <- ftr[, d]
      d <- d + 1
    }

    if (update.kappa) {
      trace.kappa[1:mcmcn1] <- ftr[, d]
      d <- d + 1
    }

    if (update.lambda) {
      trace.lambda[1:mcmcn1] <- ftr[, d]
      d <- d + 1
    }

    fn_list <- num_approx_disc(grid, parms.curr.tr[["r"]], parms.curr.tr[["p"]], parms.curr.tr[["pi"]], w.shape, w.scale,
                               ws.shape, ws.scale, dateS, dateT)

    pTTree <- log_lik_ttree(ttree, grid, fn_list, parms.curr.tr[["r"]], parms.curr.tr[["p"]], parms.curr.tr[["pi"]], w.shape, w.scale, ws.shape, ws.scale, dateS, dateT, grid.delta, NA)

    pPTree <- log_lik_ptree_given_ctree(ctree, parms.curr.coa[["kappa"]], parms.curr.coa[["lambda"]], NA)

    if (update.rho) {

      pm <- matrix((1 - parms.curr.rho) / (ndemes - 1), nrow = ndemes, ncol = ndemes)
      diag(pm) <- parms.curr.rho

      ll_out <- log_lik_locs_felsenstein(ttree, pm, demes.prior)

      pLocs <- ll_out$loglik
      dyn_L <- ll_out$dyn_L

    }

    if (verbose == F) {

      pb <- utils::txtProgressBar(min = 0, max = mcmcIterations, style = 3)

    }

    tree_prop_count <- c("add" = 0,
                         "remove" = 0,
                         "move" = 0)

    tree_acc_count <- c("add" = 0,
                        "remove" = 0,
                        "move" = 0)

    tree_acc_rates <- c("add" = 0,
                        "remove" = 0,
                        "move" = 0)

    .Random.seed <- fst$Random.seed

    for (i in 1:mcmcIterations) {

      i2 <- mcmcn1 + i

      if (update.ctree) {

        if (verbose) {

          message("Proposing ttree update")

        }

        tree_prop_count["add"] <- 0
        tree_prop_count["remove"] <- 0
        tree_prop_count["move"] <- 0

        tree_acc_count["add"] <- 0
        tree_acc_count["remove"] <- 0
        tree_acc_count["move"] <- 0

        for (j in 1:mcmc.tree.updates) {

          u <- runif(1)
          if (u < 1 / 3) {

            proptree <- add_transmission(ctree = ctree)
            prop_type <- 1
            tree_prop_count["add"] <- tree_prop_count["add"] + 1

            host_map <- proptree$host_map

            inv_host_map <- rep(NA, length(host_map) + 1)
            inv_host_map[(1:(length(host_map) + 1)) %in% host_map] <- order(host_map, na.last = NA)

          } else if (u < 2 / 3) {

            proptree <- remove_transmission(ctree = ctree)
            prop_type <- 2
            tree_prop_count["remove"] <- tree_prop_count["remove"] + 1

            host_map <- proptree$host_map

            inv_host_map <- order(host_map, na.last = NA)

          } else {

            proptree <- remove_add(ctree = ctree)
            prop_type <- 3
            tree_prop_count["move"] <- tree_prop_count["move"] + 1

            host_map <- proptree$host_map

            inv_host_map <- order(host_map, na.last = NA)

          }

          if (proptree$is_possible == 1) {

            ctree2 <- proptree$ctree

            ttree2 <- extractTTree(ctree2)

            pTTree_part <- log_lik_ttree(ttree, grid, fn_list, parms.curr.tr[["r"]], parms.curr.tr[["p"]], parms.curr.tr[["pi"]], w.shape, w.scale, ws.shape,
                                         ws.scale, dateS, dateT, grid.delta, proptree$curr_hosts)

            pTTree_part2 <- log_lik_ttree(ttree2, grid, fn_list, parms.curr.tr[["r"]], parms.curr.tr[["p"]], parms.curr.tr[["pi"]], w.shape, w.scale, ws.shape,
                                          ws.scale, dateS, dateT, grid.delta, proptree$prop_hosts)

            pPTree_part <- log_lik_ptree_given_ctree(ctree, parms.curr.coa[["kappa"]], parms.curr.coa[["lambda"]], proptree$curr_hosts)

            pPTree_part2 <- log_lik_ptree_given_ctree(ctree2, parms.curr.coa[["kappa"]], parms.curr.coa[["lambda"]], proptree$prop_hosts)

            if (update.rho) {

              dyn_L2 <- dyn_L[inv_host_map, , drop = F]

              ll_out <- log_lik_locs_felsenstein_part(ttree, pm, demes.prior, dyn_L, proptree$curr_hosts)
              ll_out2 <- log_lik_locs_felsenstein_part(ttree2, pm, demes.prior, dyn_L2, proptree$prop_hosts)

              dyn_L2 <- ll_out2$dyn_L

              #            pLocs2 <- pLocs + ll_out2$loglik - ll_out$loglik

              pLocs_diff <- ll_out2$loglik - ll_out$loglik

            } else {

              pLocs_diff <- 0

            }


            if (log(runif(1)) < (pTTree_part2 + pPTree_part2 + proptree$rev_density -
                                 pTTree_part - pPTree_part - proptree$prop_density +
                                 pLocs_diff)) {

              ctree <- ctree2
              ttree <- ttree2

              if (update.rho) {

                pLocs <- pLocs + pLocs_diff

                dyn_L <- dyn_L2

              }

              pTTree <- pTTree + pTTree_part2 - pTTree_part
              pPTree <- pPTree + pPTree_part2 - pPTree_part

              if (prop_type == 1) {

                tree_acc_count["add"] <- tree_acc_count["add"] + 1

              } else if (prop_type == 2) {

                tree_acc_count["remove"] <- tree_acc_count["remove"] + 1

              } else {

                tree_acc_count["move"] <- tree_acc_count["move"] + 1

              }

              if (grid[length(grid)] > (min(ttree$ttree[, 1]) - 0.5 * 1)) {

                grid <- seq(dateT, min(ttree$ttree[, 1]) - 0.5 * 1 - grid.delta, by = - grid.delta)

                fn_list <- num_approx_disc(grid, parms.curr.tr[["r"]], parms.curr.tr[["p"]], parms.curr.tr[["pi"]], w.shape, w.scale,
                                           ws.shape, ws.scale, dateS, dateT)

              }

            }

          }

        }

        if (i2 %% thinning == 0) {

          tree_acc_rates["add"] <- tree_acc_count["add"] / tree_prop_count["add"]
          tree_acc_rates["remove"] <- tree_acc_count["remove"] / tree_prop_count["remove"]
          tree_acc_rates["move"] <- tree_acc_count["move"] / tree_prop_count["move"]

          record[[(i2 / thinning) - thin0]]$tree.prop.counts <- tree_prop_count
          record[[(i2 / thinning) - thin0]]$tree.acc.rates <- tree_acc_rates

        }

      }

      if (ss.d.tr > 0) {

        parms.prop.tr <- MASS::mvrnorm(1, mu = parms.curr.tr, Sigma = ss.lam.tr * ss.c.tr * mcmc.cov.tr)

        ss.u.tr <- log(runif(1))

        if (parms.prop.tr["r"] > 0 & parms.prop.tr["p"] > 0 & parms.prop.tr["p"] <= 1 &
            parms.prop.tr["pi"] > 0 & parms.prop.tr["pi"] <= 1) {

          fn_list2 <- num_approx_disc(grid, parms.prop.tr[["r"]], parms.prop.tr[["p"]], parms.prop.tr[["pi"]], w.shape, w.scale,
                                      ws.shape, ws.scale, dateS, dateT)

          pTTree2 <- log_lik_ttree(ttree, grid, fn_list2, parms.prop.tr[["r"]], parms.prop.tr[["p"]],
                                   parms.prop.tr[["pi"]], w.shape, w.scale, ws.shape,
                                   ws.scale, dateS, dateT, grid.delta, NA)

          ss.alpha.tr <- (pTTree2 - pTTree) +
            (dgamma(parms.prop.tr[["r"]], shape = r.shape, scale = r.scale, log = T) -
               dgamma(parms.curr.tr[["r"]], shape = r.shape, scale = r.scale, log = T)) +
            (dbeta(parms.prop.tr[["p"]], shape1 = p.shape1, shape2 = p.shape2, log = T) -
               dbeta(parms.curr.tr[["p"]], shape1 = p.shape1, shape2 = p.shape2, log = T)) +
            (dbeta(parms.prop.tr[["pi"]], shape1 = pi.shape1, shape2 = pi.shape2, log = T) -
               dbeta(parms.curr.tr[["pi"]], shape1 = pi.shape1, shape2 = pi.shape2, log = T))

        } else {

          ss.alpha.tr <- -Inf

        }

        if (ss.u.tr < ss.alpha.tr) {

          parms.curr.tr[which(update.tr)] <- parms.prop.tr[which(update.tr)]

          fn_list <- fn_list2
          pTTree <- pTTree2

        }

        trace.r[i2] <- parms.curr.tr[["r"]]
        trace.p[i2] <- parms.curr.tr[["p"]]
        trace.pi[i2] <- parms.curr.tr[["pi"]]

        if (i2 == 1) {

          mcmc.mu.tr <- 0.5 * (parms.init.tr + parms.curr.tr)

          mcmc.cov.tr <- (1 / (ss.v0.tr + ss.d.tr + 3)) * (parms.init.tr %*% t(parms.init.tr) +
                                                             parms.curr.tr %*% t(parms.curr.tr) +
                                                             (ss.v0.tr + ss.d.tr + 1) * mcmc.cov.tr -
                                                             2 * mcmc.mu.tr %*% t(mcmc.mu.tr))

        } else if (ss.f[i2] == ss.f[i2 - 1]) {

          mcmc.mu.new.tr <- ((i2 - ss.f[i2]) / (i2 - ss.f[i2] + 1)) * mcmc.mu.tr +
            (1 / (i2 - ss.f[i2] + 1)) * parms.curr.tr

          mcmc.cov.tr <- (1 / (i2 - ss.f[i2] + ss.v0.tr + ss.d.tr + 2)) *
            ((i2 - ss.f[i2] + ss.v0.tr + ss.d.tr + 1) * mcmc.cov.tr +
               parms.curr.tr %*% t(parms.curr.tr) +
               (i2 - ss.f[i2]) * mcmc.mu.tr %*% t(mcmc.mu.tr) -
               (i2 - ss.f[i2] + 1) * mcmc.mu.new.tr %*% t(mcmc.mu.new.tr))

          mcmc.mu.tr <- mcmc.mu.new.tr

        } else {

          rem.el <- ss.f[i2] - 1

          if (rem.el == 0) {

            parms.rem.tr <- parms.init.tr

          } else {

            parms.rem.tr <- c(r = trace.r[rem.el],
                              p = trace.p[rem.el],
                              pi = trace.pi[rem.el])

          }

          mcmc.mu.new.tr <- mcmc.mu.tr + (1 / (i2 - ss.f[i2] + 1)) * (parms.curr.tr - parms.rem.tr)

          mcmc.cov.tr <- mcmc.cov.tr + (1 / (i2 - ss.f[i2] + ss.v0.tr + ss.d.tr + 2)) *
            (parms.curr.tr %*% t(parms.curr.tr) - parms.rem.tr %*% t(parms.rem.tr) +
               (i2 - ss.f[i2] + 1) * (mcmc.mu.tr %*% t(mcmc.mu.tr) - mcmc.mu.new.tr %*% t(mcmc.mu.new.tr)))

          mcmc.mu.tr <- mcmc.mu.new.tr

        }

        ss.lam.tr <- max(c(ss.min, ss.lam.tr * exp((ss.del.tr / (ss.nstart + i2)) * (min(c(1, exp(ss.alpha.tr))) - ss.a))))

        update.idx <- c(update.r, update.p, update.pi)

        zero.rc <- which(update.idx == 0)

        if (length(zero.rc) > 0) {

          mcmc.cov.tr[zero.rc, ] <- 0
          mcmc.cov.tr[, zero.rc] <- 0

        }

      }

      if (ss.d.coa > 0) {

        parms.prop.coa <- MASS::mvrnorm(1, mu = parms.curr.coa, Sigma = ss.lam.coa * ss.c.coa * mcmc.cov.coa)

        ss.u.coa <- log(runif(1))

        if (parms.prop.coa["kappa"] >= 0 & parms.prop.coa["lambda"] >= 0) {

          pPTree2 <- log_lik_ptree_given_ctree(ctree, parms.prop.coa[["kappa"]],
                                               parms.prop.coa[["lambda"]], NA)

          ss.alpha.coa <- (pPTree2 - pPTree) +
            (dgamma(parms.prop.coa[["kappa"]], shape = kappa.shape, scale = kappa.scale, log = T) -
               dgamma(parms.curr.coa[["kappa"]], shape = kappa.shape, scale = kappa.scale, log = T)) +
            (dgamma(parms.prop.coa[["lambda"]], shape = lambda.shape, scale = lambda.scale, log = T) -
               dgamma(parms.curr.coa[["lambda"]], shape = lambda.shape, scale = lambda.scale, log = T))

        } else {

          ss.alpha.coa <- -Inf

        }

        if (ss.u.coa < ss.alpha.coa) {

          parms.curr.coa[which(update.coa)] <- parms.prop.coa[which(update.coa)]

          pPTree <- pPTree2

        }

        trace.kappa[i2] <- parms.curr.coa[["kappa"]]
        trace.lambda[i2] <- parms.curr.coa[["lambda"]]

        if (i2 == 1) {

          mcmc.mu.coa <- 0.5 * (parms.init.coa + parms.curr.coa)

          mcmc.cov.coa <- (1 / (ss.v0.coa + ss.d.coa + 3)) * (parms.init.coa %*% t(parms.init.coa) +
                                                                parms.curr.coa %*% t(parms.curr.coa) +
                                                                (ss.v0.coa + ss.d.coa + 1) * mcmc.cov.coa -
                                                                2 * mcmc.mu.coa %*% t(mcmc.mu.coa))

        } else if (ss.f[i2] == ss.f[i2 - 1]) {

          mcmc.mu.new.coa <- ((i2 - ss.f[i2]) / (i2 - ss.f[i2] + 1)) * mcmc.mu.coa +
            (1 / (i2 - ss.f[i2] + 1)) * parms.curr.coa

          mcmc.cov.coa <- (1 / (i2 - ss.f[i2] + ss.v0.coa + ss.d.coa + 2)) *
            ((i2 - ss.f[i2] + ss.v0.coa + ss.d.coa + 1) * mcmc.cov.coa +
               parms.curr.coa %*% t(parms.curr.coa) +
               (i2 - ss.f[i2]) * mcmc.mu.coa %*% t(mcmc.mu.coa) -
               (i2 - ss.f[i2] + 1) * mcmc.mu.new.coa %*% t(mcmc.mu.new.coa))

          mcmc.mu.coa <- mcmc.mu.new.coa

        } else {

          rem.el <- ss.f[i2] - 1

          if (rem.el == 0) {

            parms.rem.coa <- parms.init.coa

          } else {

            parms.rem.coa <- c(kappa = trace.kappa[rem.el],
                               lambda = trace.lambda[rem.el])

          }

          mcmc.mu.new.coa <- mcmc.mu.coa + (1 / (i2 - ss.f[i2] + 1)) * (parms.curr.coa - parms.rem.coa)

          mcmc.cov.coa <- mcmc.cov.coa + (1 / (i2 - ss.f[i2] + ss.v0.coa + ss.d.coa + 2)) *
            (parms.curr.coa %*% t(parms.curr.coa) - parms.rem.coa %*% t(parms.rem.coa) +
               (i2 - ss.f[i2] + 1) * (mcmc.mu.coa %*% t(mcmc.mu.coa) - mcmc.mu.new.coa %*% t(mcmc.mu.new.coa)))

          mcmc.mu.coa <- mcmc.mu.new.coa

        }

        ss.lam.coa <- max(c(ss.min, ss.lam.coa * exp((ss.del.coa / (ss.nstart + i2)) * (min(c(1, exp(ss.alpha.coa))) - ss.a))))

        update.idx <- c(update.kappa, update.lambda)

        zero.rc <- which(update.idx == 0)

        if (length(zero.rc) > 0) {

          mcmc.cov.coa[zero.rc, ] <- 0
          mcmc.cov.coa[, zero.rc] <- 0

        }

      }

      if (update.rho) {

        parms.prop.rho <- rnorm(1, mean = parms.curr.rho, sd = sqrt(ss.lam.rho * ss.c.rho * mcmc.cov.rho))

        ss.u.rho <- log(runif(1))

        if (parms.prop.rho >= 0 & parms.prop.rho <= 1) {

          pm2 <- matrix((1 - parms.prop.rho) / (ndemes - 1), nrow = ndemes, ncol = ndemes)
          diag(pm2) <- parms.prop.rho

          ll_out2 <- log_lik_locs_felsenstein(ttree, pm2, demes.prior)

          pLocs2 <- ll_out2$loglik
          dyn_L2 <- ll_out2$dyn_L

          ss.alpha.rho <- (pLocs2 - pLocs) +
            (dbeta(parms.prop.rho, shape1 = rho.shape1, shape2 = rho.shape2, log = T) -
               dbeta(parms.curr.rho, shape1 = rho.shape1, shape2 = rho.shape2, log = T))

        } else {

          ss.alpha.rho <- -Inf

        }

        if (ss.u.rho < ss.alpha.rho) {

          parms.curr.rho <- parms.prop.rho

          pm <- pm2
          pLocs <- pLocs2
          dyn_L <- dyn_L2

        }

        trace.rho[i2] <- parms.curr.rho

        if (i2 == 1) {

          mcmc.mu.rho <- 0.5 * (parms.init.rho + parms.curr.rho)

          mcmc.cov.rho <- (1 / (ss.v0.rho + ss.d.rho + 3)) * (parms.init.rho * parms.init.rho +
                                                                parms.curr.rho * parms.curr.rho +
                                                                (ss.v0.rho + ss.d.rho + 1) * mcmc.cov.rho -
                                                                2 * mcmc.mu.rho * mcmc.mu.rho)

        } else if (ss.f[i2] == ss.f[i2 - 1]) {

          mcmc.mu.new.rho <- ((i2 - ss.f[i2]) / (i2 - ss.f[i2] + 1)) * mcmc.mu.rho +
            (1 / (i2 - ss.f[i2] + 1)) * parms.curr.rho

          mcmc.cov.rho <- (1 / (i2 - ss.f[i2] + ss.v0.rho + ss.d.rho + 2)) *
            ((i2 - ss.f[i2] + ss.v0.rho + ss.d.rho + 1) * mcmc.cov.rho +
               parms.curr.rho * parms.curr.rho +
               (i2 - ss.f[i2]) * mcmc.mu.rho * mcmc.mu.rho -
               (i2 - ss.f[i2] + 1) * mcmc.mu.new.rho * mcmc.mu.new.rho)

          mcmc.mu.rho <- mcmc.mu.new.rho

        } else {

          rem.el <- ss.f[i2] - 1

          if (rem.el == 0) {

            parms.rem.rho <- parms.init.rho

          } else {

            parms.rem.rho <- trace.rho[rem.el]

          }

          mcmc.mu.new.rho <- mcmc.mu.rho + (1 / (i2 - ss.f[i2] + 1)) * (parms.curr.rho - parms.rem.rho)

          mcmc.cov.rho <- mcmc.cov.rho + (1 / (i2 - ss.f[i2] + ss.v0.rho + ss.d.rho + 2)) *
            (parms.curr.rho * parms.curr.rho - parms.rem.rho * parms.rem.rho +
               (i2 - ss.f[i2] + 1) * (mcmc.mu.rho * mcmc.mu.rho - mcmc.mu.new.rho * mcmc.mu.new.rho))

          mcmc.mu.rho <- mcmc.mu.new.rho

        }

        ss.lam.rho <- max(c(ss.min, ss.lam.rho * exp((ss.del.rho / (ss.nstart + i2)) * (min(c(1, exp(ss.alpha.rho))) - ss.a))))

      }

      if (i2 %% thinning == 0) {

        if (verbose == F) {

          utils::setTxtProgressBar(pb, i)

        }

        if (verbose==T) {

          message(sprintf('it = %d, r = %f, off.p = %f, pi = %f, kappa = %f, lambda = %f, prior = %e, likelihood = %e, nind = %d', i, parms.curr.tr["r"], parms.curr.tr["p"], parms.curr.tr["pi"], parms.curr.coa["kappa"], parms.curr.coa["lambda"], pTTree, pPTree, nrow(ttree$ttree)))

        }

        rec_ctree <- trim_root(ctree)

        record[[(i2 / thinning) - thin0]]$ctree <- rec_ctree
        record[[(i2 / thinning) - thin0]]$pTTree <- pTTree
        record[[(i2 / thinning) - thin0]]$pPTree <- pPTree
        record[[(i2 / thinning) - thin0]]$kappa <- parms.curr.coa[["kappa"]]
        record[[(i2 / thinning) - thin0]]$lambda <- parms.curr.coa[["lambda"]]
        record[[(i2 / thinning) - thin0]]$rho <- parms.curr.rho
        record[[(i2 / thinning) - thin0]]$off.r <- parms.curr.tr[["r"]]
        record[[(i2 / thinning) - thin0]]$off.p <- parms.curr.tr[["p"]]
        record[[(i2 / thinning) - thin0]]$pi <- parms.curr.tr[["pi"]]
        record[[(i2 / thinning) - thin0]]$w.shape <- w.shape
        record[[(i2 / thinning) - thin0]]$w.scale <- w.scale
        record[[(i2 / thinning) - thin0]]$ws.shape <- ws.shape
        record[[(i2 / thinning) - thin0]]$ws.scale <- ws.scale
        record[[(i2 / thinning) - thin0]]$mcmc.cov.tr <- mcmc.cov.tr
        record[[(i2 / thinning) - thin0]]$ss.lam.tr <- ss.lam.tr
        record[[(i2 / thinning) - thin0]]$mcmc.cov.coa <- mcmc.cov.coa
        record[[(i2 / thinning) - thin0]]$ss.lam.coa <- ss.lam.coa
        record[[(i2 / thinning) - thin0]]$mcmc.cov.rho <- mcmc.cov.rho
        record[[(i2 / thinning) - thin0]]$ss.lam.rho <- ss.lam.rho

        record[[(i2 / thinning) - thin0]]$source <- rec_ctree$ctree[rec_ctree$ctree[which(rec_ctree$ctree[, 4] == 0), 2], 4]
        if (record[[(i2 / thinning) - thin0]]$source <= length(rec_ctree$nam)) {

          record[[(i2 / thinning) - thin0]]$source <- rec_ctree$nam[record[[(i2 / thinning) - thin0]]$source]

        } else {

          record[[(i2 / thinning) - thin0]]$source <- 'Unsampled'

        }

      }

    }#End of main MCMC loop

    update.vec <- c(update.r, update.p, update.pi, update.rho, update.kappa, update.lambda)

    mcmcdim <- sum(as.numeric(update.vec))

    cnames <- c("r", "p", "pi", "rho", "kappa", "lambda")[which(update.vec)]

    fulltrace <- array(dim = c(mcmcntot, mcmcdim), dimnames = list(NULL, cnames))

    d <- 1

    if (update.r) {

      fulltrace[, d] <- trace.r
      d <- d + 1

    }

    if (update.p) {

      fulltrace[, d] <- trace.p
      d <- d + 1

    }

    if (update.pi) {

      fulltrace[, d] <- trace.pi
      d <- d + 1

    }

    if (update.rho) {

      fulltrace[, d] <- trace.rho
      d <- d + 1

    }

    if (update.kappa) {

      fulltrace[, d] <- trace.kappa
      d <- d + 1

    }

    if (update.lambda) {

      fulltrace[, d] <- trace.lambda
      d <- d + 1

    }

    fulltrace <- coda::as.mcmc(fulltrace)

    record[[1]]$fulltrace <- fulltrace

    finalstate <- list()
    finalstate$ptree <- ptree
    finalstate$dateS <- dateS
    finalstate$dateT <- dateT
    finalstate$r.shape <- r.shape
    finalstate$r.scale <- r.scale
    finalstate$p.shape1 <- p.shape1
    finalstate$p.shape2 <- p.shape2
    finalstate$pi.shape1 <- pi.shape1
    finalstate$pi.shape2 <- pi.shape2
    finalstate$kappa.shape <- kappa.shape
    finalstate$kappa.scale <- kappa.scale
    finalstate$lambda.shape <- lambda.shape
    finalstate$lambda.scale <- lambda.scale
    finalstate$rho.shape1 <- rho.shape1
    finalstate$rho.shape2 <- rho.shape2
    finalstate$w.shape <- w.shape
    finalstate$w.scale <- w.scale
    finalstate$ws.shape <- ws.shape
    finalstate$ws.scale <- ws.scale
    finalstate$mcmc.tree.updates <- mcmc.tree.updates
    finalstate$ctree <- ctree
    finalstate$grid.delta <- grid.delta
    finalstate$ndemes <- ndemes
    finalstate$demes.prior <- demes.prior
    finalstate$parms.init.tr <- parms.init.tr
    finalstate$parms.init.coa <- parms.init.coa
    finalstate$parms.init.rho <- parms.init.rho
    finalstate$parms.curr.tr <- parms.curr.tr
    finalstate$parms.curr.coa <- parms.curr.coa
    finalstate$parms.curr.rho <- parms.curr.rho
    finalstate$mcmc.mu.tr <- mcmc.mu.tr
    finalstate$mcmc.mu.coa <- mcmc.mu.coa
    finalstate$mcmc.mu.rho <- mcmc.mu.rho
    finalstate$mcmc.cov.tr <- mcmc.cov.tr
    finalstate$mcmc.cov.coa <- mcmc.cov.coa
    finalstate$mcmc.cov.rho <- mcmc.cov.rho
    finalstate$update.r <- update.r
    finalstate$update.p <- update.p
    finalstate$update.pi <- update.pi
    finalstate$update.rho <- update.rho
    finalstate$update.kappa <- update.kappa
    finalstate$update.lambda <- update.lambda
    finalstate$update.ctree <- update.ctree
    finalstate$mcmcIterations <- mcmcntot
    finalstate$ss.lam.tr <- ss.lam.tr
    finalstate$ss.lam.coa <- ss.lam.coa
    finalstate$ss.lam.rho <- ss.lam.rho
    finalstate$grid <- grid
    finalstate$Random.seed <- .Random.seed
    finalstate$case <- 1

    record[[1]]$finalstate <- finalstate

    class(record)<-'resTransPhylo'
    return(record)

  }

}
