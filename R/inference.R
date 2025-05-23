#' Infer transmission tree given a phylogenetic tree
#' @param ptree Phylogenetic tree
#' @param w.shape Shape parameter of the Gamma distribution representing the generation time
#' @param w.scale Scale parameter of the Gamma distribution representing the generation time
#' @param ws.shape Shape parameter of the Gamma distribution representing the sampling time
#' @param ws.scale Scale parameter of the Gamma distribution representing the sampling time
#' @param w.mean Mean of the Gamma distribution representing the generation time
#' @param w.std Std of the Gamma distribution representing the generation time
#' @param ws.mean Mean of the Gamma distribution representing the sampling time
#' @param ws.std Std of the Gamma distribution representing the sampling time
#' @param mcmcIterations Number of MCMC iterations to run the algorithm for
#' @param thinning MCMC thinning interval between two sampled iterations
#' @param mcmc.tree.updates Number of transmission tree updates per parameter update
#' @param init.r Starting value of offspring distribution parameter r
#' @param init.p Starting value of offspring distribution parameter p
#' @param init.pi Starting value of sampling proportion pi
#' @param init.kappa Starting value for the initial within-host population size
#' @param init.lambda Starting value for the within-host population growth rate
#' @param init.rho Starting value for the probability of transmission to the same location
#' @param init.ctree Starting value for the combined tree
#' @param update.r Whether to update the parameter r
#' @param update.p Whether to update the parameter p
#' @param update.pi Whether to update the parameter pi
#' @param update.kappa Whether to update the parameter kappa
#' @param update.lambda Whether to update the parameter lambda
#' @param update.rho Whether to update parameter rho
#' @param update.ctree Whether to update the transmission tree
#' @param ndemes Number of possible demes
#' @param r.shape Shape parameter for the Gamma prior of parameter r
#' @param r.scale Scale parameter for the Gamma prior of parameter r
#' @param p.shape1 Shape1 parameter for the Beta prior of parameter p
#' @param p.shape2 Shape2 parameter for the Beta prior of parameter p
#' @param pi.shape1 Shape1 parameter for the Beta prior of parameter pi
#' @param pi.shape2 Shape2 parameter for the Beta prior of parameter pi
#' @param kappa.shape Shape parameter for the Gamma prior of parameter kappa
#' @param kappa.scale Scale parameter for the Gamma prior of parameter kappa
#' @param lambda.shape Shape parameter for the Gamma prior of parameter lambda
#' @param lambda.scale Scale parameter for the Gamma prior of parameter lambda
#' @param rho.shape1 Shape1 parameter for the Beta prior of parameter rho
#' @param rho.shape2 Shape2 parameter for the Beta prior of parameter rho
#' @param demes.prior Prior probability for the deme of the root host
#' @param dateS Start date for observations
#' @param dateT End date for observations
#' @param grid.delta Grid resolution for approximating exclusion probabilities
#' @param verbose Whether or not to use verbose mode (default is false)
#' @export
inferTTree <- function(ptree,
                        w.shape = 2,
                        w.scale = 1,
                        ws.shape = NA,
                        ws.scale = NA,
                        w.mean = NA,
                        w.std = NA,
                        ws.mean = NA,
                        ws.std = NA,
                        mcmcIterations = 12000,
                        thinning = 1,
                        mcmc.tree.updates = NA,
                        init.r = 2,
                        init.p = 0.5,
                        init.pi = 0.5,
                        init.kappa = 0.5,
                        init.lambda = 0.5,
                        init.rho = 0.8,
                        init.ctree = NA,
                        update.r = T,
                        update.p = F,
                        update.pi = T,
                        update.kappa = T,
                        update.lambda = T,
                        update.rho = F,
                        update.ctree = T,
                        ndemes = NA,
                        r.shape = 1,
                        r.scale = 1,
                        p.shape1 = 1,
                        p.shape2 = 1,
                        pi.shape1 = 1,
                        pi.shape2 = 1,
                        kappa.shape = 1,
                        kappa.scale = 1,
                        lambda.shape = 1,
                        lambda.scale = 1,
                        rho.shape1 = 1,
                        rho.shape2 = 1,
                        demes.prior = NA,
                        dateS = -Inf,
                        dateT = NA,
                        grid.delta = NA,
                        verbose = F) {

  # Ensure that all leaves have unique times
  ptree$ptree[, 1] <- ptree$ptree[, 1] + runif(nrow(ptree$ptree)) * 1e-10

  # Remove excess deme information (i.e. from simulated data)
  ptree$demes <- ptree$demes[1:length(ptree$nam)]

  # Ensure branch lengths of ptree are positive
  for (i in (ceiling(nrow(ptree$ptree) / 2) + 1):nrow(ptree$ptree)) {

    for (j in 2:3) {

      if (ptree$ptree[ptree$ptree[i, j], 1] - ptree$ptree[i, 1] < 0) {

        stop("The phylogenetic tree contains negative branch lengths!")

      }

    }

  }

  # Determine vector of primary observation times
  prim_obs_times <- calc_prim_obs(ptree)

  # Check that demes are given if updating rho
  if (length(ptree$demes) == 0 & update.rho) {

    stop('Demes are needed in ptree to update parameter rho')

  }

  if (length(ptree$demes) == 0) {

    ptree$demes[1:length(ptree$nam)] <- 1

  }

  # Ensure observation start time is consistent with observation times
  if (dateS > min(ptree$ptree[which(ptree$ptree[, 2] == 0), 1])) {

    stop('The parameter dateS cannot be later than any observation dates')

  }

  # If no observation end date is provided, set to approximately 00:00 today
  if (is.na(dateT)) {

    dateT <- 1970 + as.numeric(Sys.Date()) / 365.25

  }

  # Ensure observation end time is consistent with primary observation times
  if (dateT < max(prim_obs_times)) {

    stop('The parameter dateT cannot be earlier than any primary observation dates')

  }

  if (!is.na(w.mean) && !is.na(w.std)) {

    w.shape <- w.mean ^ 2 /  w.std ^ 2
    w.scale <- w.std ^ 2 / w.mean

  }

  if (!is.na(ws.mean)&&!is.na(ws.std)) {

    ws.shape <- ws.mean ^ 2 / ws.std ^ 2
    ws.scale <- ws.std ^ 2 / ws.mean

  }

  if (is.na(ws.shape)) {

    ws.shape <- w.shape

  }

  if (is.na(ws.scale)) {

    ws.scale <- w.scale

  }

  if (is.na(mcmc.tree.updates)) {

    mcmc.tree.updates <- length(prim_obs_times)

  }

#  if (is.na(sum(mcmc.cov.tr))) {

    mcmc.cov.tr <- diag(c(0.5 ^ 2 * as.numeric(update.r),
                          0.25 ^ 2 * as.numeric(update.p),
                          0.25 ^ 2 * as.numeric(update.pi)))

#  }

#  if (is.na(sum(mcmc.cov.coa))) {

    mcmc.cov.coa <- diag(c(0.1 ^ 2 * as.numeric(update.kappa),
                           0.1 ^ 2 * as.numeric(update.lambda)))

#  }

  mcmc.cov.rho <- 0.1 ^ 2 * as.numeric(update.rho)

  if (sum(is.na(init.ctree))) {

    ctree <- init_ctree(ptree)

  } else {

    if (length(init.ctree$demes) == 0) {

      init.ctree$demes[1:length(ptree$nam)] <- 1

    }

    ctree <- init.ctree

  }

  if (is.na(grid.delta)) {

    grid.delta <- 0.001 * (max(ptree$ptree[, 1]) - min(ptree$ptree[, 1]))

  }

  const.pi <- 3.14159265359

  parms.init.tr <- c(r = init.r,
                     p = init.p,
                     pi = init.pi)
  parms.init.coa <- c(kappa = init.kappa,
                      lambda = init.lambda)
  parms.init.rho <- init.rho

  parms.curr.tr <- parms.init.tr
  parms.curr.coa <- parms.init.coa
  parms.curr.rho <- parms.init.rho

  ss.a <- 0.234

  update.tr <- c(update.r, update.p, update.pi)
  update.coa <- c(update.kappa, update.lambda)

  mcmc.mu.tr <- NA
  mcmc.mu.coa <- NA
  mcmc.mu.rho <- NA

  ss.d.tr <- sum(as.numeric(update.tr))
  ss.d.coa <- sum(as.numeric(update.coa))
  ss.d.rho <- sum(as.numeric(update.rho))

  ss.v0.tr <- ss.d.tr
  ss.v0.coa <- ss.d.coa
  ss.v0.rho <- ss.d.rho

  ss.f <- floor(0.5 * (1:mcmcIterations))

  ss.min <- 0.1

  ss.c.tr <- 2.38 ^ 2 / ss.d.tr
  ss.c.coa <- 2.38 ^ 2 / ss.d.coa
  ss.c.rho <- 2.38 ^ 2 / ss.d.rho

  ss.lamstart <- 1

  ss.lam.tr <- ss.lamstart
  ss.lam.coa <- ss.lamstart
  ss.lam.rho <- ss.lamstart

  ss.nstart <- 5 / (ss.a * (1 - ss.a))

  ss.A <- -qnorm(ss.a / 2)

  ss.del.tr <- (1 - (1 / ss.d.tr)) * ((sqrt(2 * const.pi) * exp(ss.A ^ 2 / 2)) / (2 * ss.A)) + (1 / (ss.d.tr * ss.a * (1 - ss.a)))
  ss.del.coa <- (1 - (1 / ss.d.coa)) * ((sqrt(2 * const.pi) * exp(ss.A ^ 2 / 2)) / (2 * ss.A)) + (1 / (ss.d.coa * ss.a * (1 - ss.a)))
  ss.del.rho <- (1 - (1 / ss.d.rho)) * ((sqrt(2 * const.pi) * exp(ss.A ^ 2 / 2)) / (2 * ss.A)) + (1 / (ss.d.rho * ss.a * (1 - ss.a)))

  ttree <- extractTTree(ctree)

  record <- vector('list', mcmcIterations / thinning)

  trace.r <- numeric(mcmcIterations)
  trace.p <- numeric(mcmcIterations)
  trace.pi <- numeric(mcmcIterations)
  trace.kappa <- numeric(mcmcIterations)
  trace.lambda <- numeric(mcmcIterations)
  trace.rho <- numeric(mcmcIterations)

  grid <- seq(dateT, min(ttree$ttree[, 1]) - 0.5 * 1 - grid.delta, by = - grid.delta)

  fn_list <- num_approx_disc(grid, init.r, init.p, init.pi, w.shape, w.scale,
                             ws.shape, ws.scale, dateS, dateT)

  pTTree <- log_lik_ttree(ttree, grid, fn_list, init.r, init.p, init.pi, w.shape, w.scale,
                          ws.shape, ws.scale, dateS, dateT, grid.delta, NA)

  pPTree <- log_lik_ptree_given_ctree(ctree, init.kappa, init.lambda, NA)

  if (update.rho) {

    if (is.na(ndemes)) {

      ndemes <- max(ptree$demes)

    }

    if (is.na(demes.prior)) {

      demes.prior <- rep(1 / ndemes, ndemes)

    }

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

  for (i in 1:mcmcIterations) {

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

      if (i %% thinning == 0) {

        tree_acc_rates["add"] <- tree_acc_count["add"] / tree_prop_count["add"]
        tree_acc_rates["remove"] <- tree_acc_count["remove"] / tree_prop_count["remove"]
        tree_acc_rates["move"] <- tree_acc_count["move"] / tree_prop_count["move"]

        record[[i / thinning]]$tree.prop.counts <- tree_prop_count
        record[[i / thinning]]$tree.acc.rates <- tree_acc_rates

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

      trace.r[i] <- parms.curr.tr[["r"]]
      trace.p[i] <- parms.curr.tr[["p"]]
      trace.pi[i] <- parms.curr.tr[["pi"]]

      if (i == 1) {

        mcmc.mu.tr <- 0.5 * (parms.init.tr + parms.curr.tr)

        mcmc.cov.tr <- (1 / (ss.v0.tr + ss.d.tr + 3)) * (parms.init.tr %*% t(parms.init.tr) +
                                                  parms.curr.tr %*% t(parms.curr.tr) +
                                                  (ss.v0.tr + ss.d.tr + 1) * mcmc.cov.tr -
                                                  2 * mcmc.mu.tr %*% t(mcmc.mu.tr))

      } else if (ss.f[i] == ss.f[i - 1]) {

        mcmc.mu.new.tr <- ((i - ss.f[i]) / (i - ss.f[i] + 1)) * mcmc.mu.tr +
          (1 / (i - ss.f[i] + 1)) * parms.curr.tr

        mcmc.cov.tr <- (1 / (i - ss.f[i] + ss.v0.tr + ss.d.tr + 2)) *
          ((i - ss.f[i] + ss.v0.tr + ss.d.tr + 1) * mcmc.cov.tr +
             parms.curr.tr %*% t(parms.curr.tr) +
             (i - ss.f[i]) * mcmc.mu.tr %*% t(mcmc.mu.tr) -
             (i - ss.f[i] + 1) * mcmc.mu.new.tr %*% t(mcmc.mu.new.tr))

        mcmc.mu.tr <- mcmc.mu.new.tr

      } else {

        rem.el <- ss.f[i] - 1

        if (rem.el == 0) {

          parms.rem.tr <- parms.init.tr

        } else {

          parms.rem.tr <- c(r = trace.r[rem.el],
                            p = trace.p[rem.el],
                            pi = trace.pi[rem.el])

        }

        mcmc.mu.new.tr <- mcmc.mu.tr + (1 / (i - ss.f[i] + 1)) * (parms.curr.tr - parms.rem.tr)

        mcmc.cov.tr <- mcmc.cov.tr + (1 / (i - ss.f[i] + ss.v0.tr + ss.d.tr + 2)) *
          (parms.curr.tr %*% t(parms.curr.tr) - parms.rem.tr %*% t(parms.rem.tr) +
             (i - ss.f[i] + 1) * (mcmc.mu.tr %*% t(mcmc.mu.tr) - mcmc.mu.new.tr %*% t(mcmc.mu.new.tr)))

        mcmc.mu.tr <- mcmc.mu.new.tr

      }

      ss.lam.tr <- max(c(ss.min, ss.lam.tr * exp((ss.del.tr / (ss.nstart + i)) * (min(c(1, exp(ss.alpha.tr))) - ss.a))))

      update.idx.coa <- c(update.r, update.p, update.pi)

      zero.rc <- which(update.idx.coa == 0)

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

      trace.kappa[i] <- parms.curr.coa[["kappa"]]
      trace.lambda[i] <- parms.curr.coa[["lambda"]]

      if (i == 1) {

        mcmc.mu.coa <- 0.5 * (parms.init.coa + parms.curr.coa)

        mcmc.cov.coa <- (1 / (ss.v0.coa + ss.d.coa + 3)) * (parms.init.coa %*% t(parms.init.coa) +
                                                           parms.curr.coa %*% t(parms.curr.coa) +
                                                           (ss.v0.coa + ss.d.coa + 1) * mcmc.cov.coa -
                                                           2 * mcmc.mu.coa %*% t(mcmc.mu.coa))

      } else if (ss.f[i] == ss.f[i - 1]) {

        mcmc.mu.new.coa <- ((i - ss.f[i]) / (i - ss.f[i] + 1)) * mcmc.mu.coa +
          (1 / (i - ss.f[i] + 1)) * parms.curr.coa

        mcmc.cov.coa <- (1 / (i - ss.f[i] + ss.v0.coa + ss.d.coa + 2)) *
          ((i - ss.f[i] + ss.v0.coa + ss.d.coa + 1) * mcmc.cov.coa +
             parms.curr.coa %*% t(parms.curr.coa) +
             (i - ss.f[i]) * mcmc.mu.coa %*% t(mcmc.mu.coa) -
             (i - ss.f[i] + 1) * mcmc.mu.new.coa %*% t(mcmc.mu.new.coa))

        mcmc.mu.coa <- mcmc.mu.new.coa

      } else {

        rem.el <- ss.f[i] - 1

        if (rem.el == 0) {

          parms.rem.coa <- parms.init.coa

        } else {

          parms.rem.coa <- c(kappa = trace.kappa[rem.el],
                             lambda = trace.lambda[rem.el])

        }

        mcmc.mu.new.coa <- mcmc.mu.coa + (1 / (i - ss.f[i] + 1)) * (parms.curr.coa - parms.rem.coa)

        mcmc.cov.coa <- mcmc.cov.coa + (1 / (i - ss.f[i] + ss.v0.coa + ss.d.coa + 2)) *
          (parms.curr.coa %*% t(parms.curr.coa) - parms.rem.coa %*% t(parms.rem.coa) +
             (i - ss.f[i] + 1) * (mcmc.mu.coa %*% t(mcmc.mu.coa) - mcmc.mu.new.coa %*% t(mcmc.mu.new.coa)))

        mcmc.mu.coa <- mcmc.mu.new.coa

      }

      ss.lam.coa <- max(c(ss.min, ss.lam.coa * exp((ss.del.coa / (ss.nstart + i)) * (min(c(1, exp(ss.alpha.coa))) - ss.a))))

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

      trace.rho[i] <- parms.curr.rho

      if (i == 1) {

        mcmc.mu.rho <- 0.5 * (parms.init.rho + parms.curr.rho)

        mcmc.cov.rho <- (1 / (ss.v0.rho + ss.d.rho + 3)) * (parms.init.rho * parms.init.rho +
                                                              parms.curr.rho * parms.curr.rho +
                                                              (ss.v0.rho + ss.d.rho + 1) * mcmc.cov.rho -
                                                              2 * mcmc.mu.rho * mcmc.mu.rho)

      } else if (ss.f[i] == ss.f[i - 1]) {

        mcmc.mu.new.rho <- ((i - ss.f[i]) / (i - ss.f[i] + 1)) * mcmc.mu.rho +
          (1 / (i - ss.f[i] + 1)) * parms.curr.rho

        mcmc.cov.rho <- (1 / (i - ss.f[i] + ss.v0.rho + ss.d.rho + 2)) *
          ((i - ss.f[i] + ss.v0.rho + ss.d.rho + 1) * mcmc.cov.rho +
             parms.curr.rho * parms.curr.rho +
             (i - ss.f[i]) * mcmc.mu.rho * mcmc.mu.rho -
             (i - ss.f[i] + 1) * mcmc.mu.new.rho * mcmc.mu.new.rho)

        mcmc.mu.rho <- mcmc.mu.new.rho

      } else {

        rem.el <- ss.f[i] - 1

        if (rem.el == 0) {

          parms.rem.rho <- parms.init.rho

        } else {

          parms.rem.rho <- trace.rho[rem.el]

        }

        mcmc.mu.new.rho <- mcmc.mu.rho + (1 / (i - ss.f[i] + 1)) * (parms.curr.rho - parms.rem.rho)

        mcmc.cov.rho <- mcmc.cov.rho + (1 / (i - ss.f[i] + ss.v0.rho + ss.d.rho + 2)) *
          (parms.curr.rho * parms.curr.rho - parms.rem.rho * parms.rem.rho +
             (i - ss.f[i] + 1) * (mcmc.mu.rho * mcmc.mu.rho - mcmc.mu.new.rho * mcmc.mu.new.rho))

        mcmc.mu.rho <- mcmc.mu.new.rho

      }

      ss.lam.rho <- max(c(ss.min, ss.lam.rho * exp((ss.del.rho / (ss.nstart + i)) * (min(c(1, exp(ss.alpha.rho))) - ss.a))))

    }

    if (i %% thinning == 0) {

      if (verbose == F) {

        utils::setTxtProgressBar(pb, i)

      }

      if (verbose==T) {

        message(sprintf('it = %d, r = %f, off.p = %f, pi = %f, kappa = %f, lambda = %f, prior = %e, likelihood = %e, nind = %d', i, parms.curr.tr["r"], parms.curr.tr["p"], parms.curr.tr["pi"], parms.curr.coa["kappa"], parms.curr.coa["lambda"], pTTree, pPTree, nrow(ttree$ttree)))

      }

      rec_ctree <- trim_root(ctree)

      record[[i / thinning]]$ctree <- rec_ctree
      record[[i / thinning]]$pTTree <- pTTree
      record[[i / thinning]]$pPTree <- pPTree
      record[[i / thinning]]$kappa <- parms.curr.coa[["kappa"]]
      record[[i / thinning]]$lambda <- parms.curr.coa[["lambda"]]
      record[[i / thinning]]$rho <- parms.curr.rho
      record[[i / thinning]]$off.r <- parms.curr.tr[["r"]]
      record[[i / thinning]]$off.p <- parms.curr.tr[["p"]]
      record[[i / thinning]]$pi <- parms.curr.tr[["pi"]]
      record[[i / thinning]]$w.shape <- w.shape
      record[[i / thinning]]$w.scale <- w.scale
      record[[i / thinning]]$ws.shape <- ws.shape
      record[[i / thinning]]$ws.scale <- ws.scale
      record[[i / thinning]]$mcmc.cov.tr <- mcmc.cov.tr
      record[[i / thinning]]$ss.lam.tr <- ss.lam.tr
      record[[i / thinning]]$mcmc.cov.coa <- mcmc.cov.coa
      record[[i / thinning]]$ss.lam.coa <- ss.lam.coa
      record[[i / thinning]]$mcmc.cov.rho <- mcmc.cov.rho
      record[[i / thinning]]$ss.lam.rho <- ss.lam.rho

      record[[i / thinning]]$source <- rec_ctree$ctree[rec_ctree$ctree[which(rec_ctree$ctree[, 4] == 0), 2], 4]
      if (record[[i / thinning]]$source <= length(rec_ctree$nam)) {

        record[[i / thinning]]$source <- rec_ctree$nam[record[[i / thinning]]$source]

      } else {

        record[[i / thinning]]$source <- 'Unsampled'

      }

    }

  }#End of main MCMC loop

  update.vec <- c(update.r, update.p, update.pi, update.rho, update.kappa, update.lambda)

  mcmcdim <- sum(as.numeric(update.vec))

  cnames <- c("r", "p", "pi", "rho", "kappa", "lambda")[which(update.vec)]

  fulltrace <- array(dim = c(mcmcIterations, mcmcdim), dimnames = list(NULL, cnames))

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
  finalstate$mcmcIterations <- mcmcIterations
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
