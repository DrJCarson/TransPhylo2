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
#' @param mcmc.cov.tr Initial proposal covariance for transmission parameters
#' @param mcmc.cov.coa Initial proposal covariance for coalescent parameters
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
#' @param ndemes Number of possible locations
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
#' @param demes.prior Prior probability for the location of the root host
#' @param dateS Start date for observations
#' @param dateT End date for observations
#' @param grid.delta Grid resolution for approximating exclusion probabilities
#' @param verbose Whether or not to use verbose mode (default is false)
#' @export
inferTTreemulti <- function(ptree,
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
                            mcmc.cov.tr = NA,
                            mcmc.cov.coa = NA,
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
                            update.rho = T,
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
                            tr.demes = NA,
                            pi.demes = NA,
                            rho.demes = NA,
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

  # Check that locations are given if updating rho
  if (length(ptree$demes) == 0 & update.rho) {

    stop('Locations are needed in ptree to update parameter rho')

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

  if (is.na(ndemes)) {

    ndemes <- max(ptree$demes)

  }

  if (is.na(demes.prior)) {

    demes.prior <- rep(1 / ndemes, ndemes)

  }

  if (is.na(tr.demes[1])) {

    tr.demes <- 1:ndemes

  }

  if (is.na(pi.demes[1])) {

    pi.demes <- 1:ndemes

  }

  if (is.na(rho.demes[1])) {

    rho.demes <- 1:ndemes

  }

  tr.dim <- length(unique(tr.demes))

  pi.dim <- length(unique(pi.demes))

  rho.dim <- length(unique(rho.demes))

  max.dim <- max(c(tr.dim, pi.dim, rho.dim))

  if (length(init.r) == 1) {

    init.r <- rep(init.r, tr.dim)

  }

  if (length(init.p) == 1) {

    init.p <- rep(init.p, tr.dim)

  }

  if (length(init.pi) == 1) {

    init.pi <- rep(init.pi, pi.dim)

  }

  if (length(init.rho) == 1) {

    init.rho <- rep(init.rho, rho.dim)

  }


#  if (is.na(sum(mcmc.cov.tr))) {

#    mcmc.cov.tr <- diag(c(0.5 ^ 2 * as.numeric(update.r),
#                          0.25 ^ 2 * as.numeric(update.p),
#                          0.25 ^ 2 * as.numeric(update.pi)),
#                          0.25 ^ 2 * as.numeric(update.rho))

#  }

#  if (is.na(sum(mcmc.cov.coa))) {

#    mcmc.cov.coa <- diag(c(0.1 ^ 2 * as.numeric(update.kappa),
#                           0.1 ^ 2 * as.numeric(update.lambda)))

#  }



  parms.init.tr <- array(NA, dim = c(4, max.dim),
                         dimnames = list(c("r", "p", "pi", "rho"), NULL))

  parms.init.tr["r", 1:tr.dim] <- init.r
  parms.init.tr["p", 1:tr.dim] <- init.p
  parms.init.tr["pi", 1:pi.dim] <- init.pi
  parms.init.tr["rho", 1:rho.dim] <- init.rho

  parms.init.coa <- c(kappa = init.kappa,
                      lambda = init.lambda)

  parms.curr.tr <- parms.init.tr
  parms.curr.coa <- parms.init.coa

  within_bounds <- numeric(4)

  mcmc.mu.tr <- list()

  if (is.na(sum(mcmc.cov.tr))) {

    mcmc.cov.tr <- list()

    for (i in 1:max.dim) {

      mcmc.cov.tr[[i]] <- diag(c(0.5 ^ 2 * as.numeric(update.r) * as.numeric(!is.na(parms.init.tr["r", i])),
                                 0.25 ^ 2 * as.numeric(update.p) * as.numeric(!is.na(parms.init.tr["p", i])),
                                 0.25 ^ 2 * as.numeric(update.pi) * as.numeric(!is.na(parms.init.tr["pi", i])),
                               0.25 ^ 2 * as.numeric(update.rho) * as.numeric(!is.na(parms.init.tr["rho", i]))))

    }

  }

  if (is.na(sum(mcmc.cov.coa))) {

    mcmc.cov.coa <- diag(c(0.1 ^ 2 * as.numeric(update.kappa),
                           0.1 ^ 2 * as.numeric(update.lambda)))

  }

  ss.a <- 0.234

  update.tr <- c(update.r, update.p, update.pi, update.rho)
  update.coa <- c(update.kappa, update.lambda)

  ss.d.tr <- numeric(max.dim)

  for (i in 1:max.dim) {

    ss.d.tr[i] <- as.numeric(update.r) * as.numeric(!is.na(parms.init.tr["r", i])) +
      as.numeric(update.p) * as.numeric(!is.na(parms.init.tr["p", i])) +
      as.numeric(update.pi) * as.numeric(!is.na(parms.init.tr["pi", i])) +
      as.numeric(update.rho) * as.numeric(!is.na(parms.init.tr["rho", i]))

  }

  ss.d.coa <- sum(as.numeric(update.coa))

  ss.v0.tr <- ss.d.tr
  ss.v0.coa <- ss.d.coa

  ss.f <- floor(0.5 * (1:mcmcIterations))

  ss.min <- 0.1

  ss.c.tr <- numeric(max.dim)

  for (i in 1:max.dim) {

    ss.c.tr[i] <- 2.38 ^ 2 / ss.d.tr[i]

  }

  ss.c.coa <- 2.38 ^ 2 / ss.d.coa

  ss.lamstart <- 1

  ss.lam.tr <- rep(ss.lamstart, max.dim)
  ss.lam.coa <- ss.lamstart

  ss.nstart <- 5 / (ss.a * (1 - ss.a))

  ss.A <- -qnorm(ss.a / 2)

  ss.del.tr <- numeric(max.dim)

  for (i in 1:max.dim) {

    ss.del.tr[i] <- (1 - (1 / ss.d.tr[i])) * ((sqrt(2 * const.pi) * exp(ss.A ^ 2 / 2)) / (2 * ss.A)) + (1 / (ss.d.tr[i] * ss.a * (1 - ss.a)))

  }

  ss.del.coa <- (1 - (1 / ss.d.coa)) * ((sqrt(2 * const.pi) * exp(ss.A ^ 2 / 2)) / (2 * ss.A)) + (1 / (ss.d.coa * ss.a * (1 - ss.a)))

  ttree <- extractTTree(ctree)

  record <- vector('list', mcmcIterations / thinning)

  trace.r <- array(dim = c(mcmcIterations, tr.dim))
  trace.p <- array(dim = c(mcmcIterations, tr.dim))
  trace.pi <- array(dim = c(mcmcIterations, pi.dim))
  trace.rho <- array(dim = c(mcmcIterations, rho.dim))

  trace.kappa <- numeric(mcmcIterations)
  trace.lambda <- numeric(mcmcIterations)

  grid <- seq(dateT, min(ttree$ttree[, 1]) - 0.5 * 1 - grid.delta, by = - grid.delta)

  ext.r <- init.r[tr.demes]
  ext.p <- init.p[tr.demes]
  ext.pi <- init.pi[pi.demes]

  ext.rho <- init.rho[rho.demes]

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

  pPTree <- log_lik_ptree_given_ctree(ctree, init.kappa, init.lambda, NA)


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


      #############################
      ## Make sure parms are correct
      #############################


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

      if (i %% thinning == 0) {

        tree_acc_rates["add"] <- tree_acc_count["add"] / tree_prop_count["add"]
        tree_acc_rates["remove"] <- tree_acc_count["remove"] / tree_prop_count["remove"]
        tree_acc_rates["move"] <- tree_acc_count["move"] / tree_prop_count["move"]

        record[[i / thinning]]$tree.prop.counts <- tree_prop_count
        record[[i / thinning]]$tree.acc.rates <- tree_acc_rates

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
                (dgamma(parms.prop.tr[["r"]], shape = r.shape, scale = r.scale, log = T) -
                   dgamma(parms.curr.tr[["r", k]], shape = r.shape, scale = r.scale, log = T)) +
                (dbeta(parms.prop.tr[["p"]], shape1 = p.shape1, shape2 = p.shape2, log = T) -
                   dbeta(parms.curr.tr[["p", k]], shape1 = p.shape1, shape2 = p.shape2, log = T))

            }

            if (k <= pi.dim) {

              ss.alpha.tr <- ss.alpha.tr +
                (dbeta(parms.prop.tr[["pi"]], shape1 = pi.shape1, shape2 = pi.shape2, log = T) -
                   dbeta(parms.curr.tr[["pi", k]], shape1 = pi.shape1, shape2 = pi.shape2, log = T))

            }

            if (k <= rho.dim) {

              ss.alpha.tr <- ss.alpha.tr +
                (dbeta(parms.prop.tr[["rho"]], shape1 = rho.shape1, shape2 = rho.shape2, log = T) -
                   dbeta(parms.curr.tr[["rho", k]], shape1 = rho.shape1, shape2 = rho.shape2, log = T))

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

            trace.r[i, k] <- parms.curr.tr[["r", k]]
            trace.p[i, k] <- parms.curr.tr[["p", k]]

          }

          if (k <= pi.dim) {

            trace.pi[i, k] <- parms.curr.tr[["pi", k]]

          }

          if (k <= rho.dim) {

            trace.rho[i, k] <- parms.curr.tr[["rho", k]]

          }

          if (i == 1) {

            mcmc.mu.tr[[k]] <- 0.5 * (parms.init.tr[, k] + parms.curr.tr[, k])

            mcmc.cov.tr[[k]][updt_idx, updt_idx] <- ((1 / (ss.v0.tr[k] + ss.d.tr[k] + 3)) * (parms.init.tr[, k] %*% t(parms.init.tr[, k]) + parms.curr.tr[, k] %*% t(parms.curr.tr[, k]) + (ss.v0.tr[k] + ss.d.tr[k] + 1) * mcmc.cov.tr[[k]] - 2 * mcmc.mu.tr[[k]] %*% t(mcmc.mu.tr[[k]])))[updt_idx, updt_idx]

          } else if (ss.f[i] == ss.f[i - 1]) {

            mcmc.mu.new.tr <- ((i - ss.f[i]) / (i - ss.f[i] + 1)) * mcmc.mu.tr[[k]] +
              (1 / (i - ss.f[i] + 1)) * parms.curr.tr[, k]

            mcmc.cov.tr[[k]][updt_idx, updt_idx] <- ((1 / (i - ss.f[i] + ss.v0.tr[k] + ss.d.tr[k] + 2)) * ((i - ss.f[i] + ss.v0.tr[k] + ss.d.tr[k] + 1) * mcmc.cov.tr[[k]] + parms.curr.tr[, k] %*% t(parms.curr.tr[, k]) + (i - ss.f[i]) * mcmc.mu.tr[[k]] %*% t(mcmc.mu.tr[[k]]) - (i - ss.f[i] + 1) * mcmc.mu.new.tr %*% t(mcmc.mu.new.tr)))[updt_idx, updt_idx]

            mcmc.mu.tr[[k]] <- mcmc.mu.new.tr

          } else {

            rem.el <- ss.f[i] - 1

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

            mcmc.mu.new.tr <- mcmc.mu.tr[[k]] + (1 / (i - ss.f[i] + 1)) * (parms.curr.tr[, k] - parms.rem.tr)

            mcmc.cov.tr[[k]][updt_idx, updt_idx] <- (mcmc.cov.tr[[k]] + (1 / (i - ss.f[i] + ss.v0.tr[k] + ss.d.tr[k] + 2)) * (parms.curr.tr[, k] %*% t(parms.curr.tr[, k]) - parms.rem.tr %*% t(parms.rem.tr) + (i - ss.f[i] + 1) * (mcmc.mu.tr[[k]] %*% t(mcmc.mu.tr[[k]]) - mcmc.mu.new.tr %*% t(mcmc.mu.new.tr))))[updt_idx, updt_idx]

            mcmc.mu.tr[[k]] <- mcmc.mu.new.tr

          }

          ss.lam.tr[k] <- max(c(ss.min, ss.lam.tr[k] * exp((ss.del.tr[k] / (ss.nstart + i)) * (min(c(1, exp(ss.alpha.tr))) - ss.a))))

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

      record[[i / thinning]]$off.r <- parms.curr.tr["r", ]
      record[[i / thinning]]$off.p <- parms.curr.tr["p", ]
      record[[i / thinning]]$pi <- parms.curr.tr["pi", ]
      record[[i / thinning]]$rho <- parms.curr.tr["rho", ]

      record[[i / thinning]]$w.shape <- w.shape
      record[[i / thinning]]$w.scale <- w.scale
      record[[i / thinning]]$ws.shape <- ws.shape
      record[[i / thinning]]$ws.scale <- ws.scale
      record[[i / thinning]]$mcmc.cov.tr <- mcmc.cov.tr
      record[[i / thinning]]$ss.lam.tr <- ss.lam.tr
      record[[i / thinning]]$mcmc.cov.coa <- mcmc.cov.coa
      record[[i / thinning]]$ss.lam.coa <- ss.lam.coa

      record[[i / thinning]]$source <- rec_ctree$ctree[rec_ctree$ctree[which(rec_ctree$ctree[, 4] == 0), 2], 4]
      if (record[[i / thinning]]$source <= length(rec_ctree$nam)) {

        record[[i / thinning]]$source <- rec_ctree$nam[record[[i / thinning]]$source]

      } else {

        record[[i / thinning]]$source <- 'Unsampled'

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

  if (update.kappa) {

    cnames[d] <- "lambda"

    d <- d + 1

  }

  fulltrace <- array(dim = c(mcmcIterations, mcmcdim), dimnames = list(NULL, cnames))

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

  class(record)<-'resTransPhylo'
  return(record)

}
