#' Calculate calculate log T for given host and deme
#'
#' @param ttree Transmission tree
#' @param obs Observation times
#' @param grid Discrete grid over which to evaluate functions
#' @param omega Exclusion probabilities
#' @param omega_bar Offspring exclusion probabilities
#' @param pit Observation probabilities
#' @param off.r Shape parameter for the number of offspring
#' @param off.p Probability parameter for the number of offspring
#' @param pi Probability of host being sampled
#' @param ws.shape Shape parameter of primary sampling time distribution
#' @param ws.scale Scale parameter of primary sampling time distribution
#' @param obs.start Start time of outbreak sampling
#' @param obs.end Stop time of outbreak sampling
#' @param grid.delta Discrete time step
#' @param hosts Host for which likelihood is calculated
#' @param deme Deme for which likelihood is calculated
dyn_T <- function(ttree, obs, grid, omega, omega_bar, pit, off.r, off.p, pi, ws.shape,
                          ws.scale, obs.start, obs.end, grid.delta = 1 / 365, host, deme) {

  inf_time <- ttree[host, 1]

  tidx <- 1 + floor((obs.end + 1e-10 - inf_time) / grid.delta)

  omega_int <- omega[tidx, deme] + ((inf_time - grid[tidx]) / (grid[tidx + 1] - grid[tidx])) * (omega[tidx + 1, deme] - omega[tidx, deme])
  omega_bar_int <- omega_bar[tidx, deme] + ((inf_time - grid[tidx]) / (grid[tidx + 1] - grid[tidx])) * (omega_bar[tidx + 1, deme] - omega_bar[tidx, deme])
  pit_int <- pit[tidx, deme] + ((inf_time - grid[tidx]) / (grid[tidx + 1] - grid[tidx])) * (pit[tidx + 1, deme] - pit[tidx, deme])

  inc_off <- length(which(ttree[, 3] == host))

  sum_lim <- qnbinom(1 - 1e-10, size = off.r[deme], prob = off.p[deme])

  alpha_sum <- sum(dnbinom(inc_off:sum_lim, size = off.r[deme], prob = off.p[deme]) *
                     choose(inc_off:sum_lim, inc_off) *
                     omega_bar_int ^ (0:(sum_lim - inc_off)))

  # Host is observed
  if (ttree[host, 2] > 0) {

    obs_time <- min(obs[which(obs[, 2] == host), 1])

    Tns <- log(pi[deme]) + dgamma(obs_time - inf_time, shape = ws.shape, scale = ws.scale, log = T) -
      log(1 - omega_int) + log(alpha_sum) + lfactorial(inc_off)

  # Host is unobserved
  } else {

    Tns <- log(1 - pit_int) - log(1 - omega_int) + log(alpha_sum) + lfactorial(inc_off)

  }

  return(Tns)

}


#' Calculate calculate log U for given hosts and demes
#'
#' @param ttree Transmission tree
#' @param grid Discrete grid over which to evaluate functions
#' @param omega Exclusion probabilities
#' @param w.shape Shape parameter of generation time distribution
#' @param w.scale Scale parameter of generation time distribution
#' @param obs.end Stop time of outbreak sampling
#' @param grid.delta Discrete time step
#' @param host1 Infector for which likelihood is calculated
#' @param deme1 Infector deme for which likelihood is calculated
#' @param host2 Infected for which likelihood is calculated
#' @param deme2 Infected deme for which likelihood is calculated
#' @param pm Transition probability matrix between demes
dyn_U <- function(ttree, grid, omega, w.shape, w.scale, obs.end, grid.delta = 1 / 365, host1, deme1, host2, deme2, pm) {

  inf_time1 <- ttree[host1, 1]
  inf_time2 <- ttree[host2, 1]

  tidx2 <- 1 + floor((obs.end + 1e-10 - inf_time2) / grid.delta)

  omega_int2 <- omega[tidx2, deme2] + ((inf_time2 - grid[tidx2]) / (grid[tidx2 + 1] - grid[tidx2])) * (omega[tidx2 + 1, deme2] - omega[tidx2, deme2])

  Uns <- log(1 - omega_int2) + dgamma(inf_time2 - inf_time1, shape = w.shape, scale = w.scale, log = T) +
    log(pm[deme1, deme2])

  return(Uns)

}



#' Calculate transmission tree log likelihood
#'
#' @param ttree Transmission tree
#' @param grid Discrete grid over which to evaluate functions
#' @param fn_list Precalculated discrete approximations of exclusion probabilities
#' @param off.r Shape parameter for the number of offspring
#' @param off.p Probability parameter for the number of offspring
#' @param pi Probability of host being sampled
#' @param w.shape Shape parameter of generation time distribution
#' @param w.scale Scale parameter of generation time distribution
#' @param ws.shape Shape parameter of primary sampling time distribution
#' @param ws.scale Scale parameter of primary sampling time distribution
#' @param obs.start Start time of outbreak sampling
#' @param obs.end Stop time of outbreak sampling
#' @param grid.delta Discrete time step
#' @param ndemes Number of demes
#' @param pm Transmission probability matrix between demes
#' @param demes.prior Prior probability for the deme of the root host
log_lik_ttree_multiparm <- function(ttree, grid, fn_list, off.r, off.p, pi, w.shape, w.scale, ws.shape,
                          ws.scale, obs.start, obs.end, grid.delta, ndemes, pm, demes.prior) {


  obs <- ttree$obs
  demes <- ttree$demes
  ttree <- ttree$ttree

  omega <- fn_list$omega
  omega_bar <- fn_list$omega_bar
  pit <- fn_list$pit
  gamma_prob <- fn_list$gamma_prob

  nhosts <- length(ttree[, 1])

  dyn_L <- array(-Inf, dim = c(nhosts, ndemes))
  lul <- numeric(ndemes)

  host_order <- order(ttree[, 1], decreasing = T)

  for (i in host_order) {

    # Leaf
    if (length(which(ttree[, 3] == i)) == 0) {

      if (!is.na(demes[i])) {

        if (demes[i] > 0) {

          dyn_L[i, demes[i]] <- dyn_T(ttree, obs, grid, omega, omega_bar, pit, off.r, off.p, pi, ws.shape,
                                      ws.scale, obs.start, obs.end, grid.delta, i, demes[i])

        } else {

          for (d in 1:ndemes) {

            dyn_L[i, d] <- dyn_T(ttree, obs, grid, omega, omega_bar, pit, off.r, off.p, pi, ws.shape,
                                        ws.scale, obs.start, obs.end, grid.delta, i, d)

          }

        }

      } else {

        for (d in 1:ndemes) {

          dyn_L[i, d] <- dyn_T(ttree, obs, grid, omega, omega_bar, pit, off.r, off.p, pi, ws.shape,
                               ws.scale, obs.start, obs.end, grid.delta, i, d)

        }

      }

      # Not a leaf
    } else {

      children <- which(ttree[, 3] == i)

      if (!is.na(demes[i])) {

        if (demes[i] > 0) {

          dyn_L[i, demes[i]] <- dyn_T(ttree, obs, grid, omega, omega_bar, pit, off.r, off.p, pi, ws.shape,
                                      ws.scale, obs.start, obs.end, grid.delta, i, demes[i])

          for (j in children) {

            for (d2 in 1:ndemes) {

              lul[d2] <- dyn_U(ttree, grid, omega, w.shape, w.scale, obs.end, grid.delta, i, demes[i], j, d2, pm) +
                dyn_L[j, d2]

            }

            m <- max(lul)

            dyn_L[i, demes[i]] <- dyn_L[i, demes[i]] +
              log(sum(exp(lul - m))) + m

          }

        } else {

          for (d in 1:ndemes) {

            dyn_L[i, d] <- dyn_T(ttree, obs, grid, omega, omega_bar, pit, off.r, off.p, pi, ws.shape,
                                 ws.scale, obs.start, obs.end, grid.delta, i, d)

            for (j in children) {

              for (d2 in 1:ndemes) {

                lul[d2] <- dyn_U(ttree, grid, omega, w.shape, w.scale, obs.end, grid.delta, i, d, j, d2, pm) +
                  dyn_L[j, d2]

              }

              m <- max(lul)

              dyn_L[i, d] <- dyn_L[i, d] +
                log(sum(exp(lul - m))) + m

            }

          }

        }

      } else {

        for (d in 1:ndemes) {

          dyn_L[i, d] <- dyn_T(ttree, obs, grid, omega, omega_bar, pit, off.r, off.p, pi, ws.shape,
                               ws.scale, obs.start, obs.end, grid.delta, i, d)

          for (j in children) {

            for (d2 in 1:ndemes) {

              lul[d2] <- dyn_U(ttree, grid, omega, w.shape, w.scale, obs.end, grid.delta, i, d, j, d2, pm) +
                dyn_L[j, d2]

            }

            m <- max(lul)

            dyn_L[i, d] <- dyn_L[i, d] +
              log(sum(exp(lul - m))) + m

          }

        }

      }

    }

  }

  root_host <- host_order[length(host_order)]

  lpl <- log(demes.prior) + dyn_L[root_host, ]

  rootm <- max(lpl)

  log_likelihood <- log(sum(exp(lpl - rootm))) + rootm

  return(list(loglik = log_likelihood, dyn_L = dyn_L))

}

#' Calculate transmission tree likelihood
#'
#' @param ttree Transmission tree
#' @param grid Discrete grid over which to evaluate functions
#' @param fn_list Precalculated discrete approximations of exclusion probabilities
#' @param off.r Shape parameter for the number of offspring
#' @param off.p Probability parameter for the number of offspring
#' @param pi Probability of host being sampled
#' @param w.shape Shape parameter of generation time distribution
#' @param w.scale Scale parameter of generation time distribution
#' @param ws.shape Shape parameter of primary sampling time distribution
#' @param ws.scale Scale parameter of primary sampling time distribution
#' @param obs.start Start time of outbreak sampling
#' @param obs.end Stop time of outbreak sampling
#' @param grid.delta Discrete time step
#' @param ndemes Number of demes
#' @param pm Transmission probability matrix between demes
#' @param demes.prior Prior probability for the deme of the root host
#' @param dyn_L Existing likelihood matrix from the dynamic programming algorithm
#' @param hosts Hosts over which likelihood is calculated
log_lik_ttree_multiparm_part <- function(ttree, grid, fn_list, off.r, off.p, pi, w.shape, w.scale, ws.shape,
                                    ws.scale, obs.start, obs.end, grid.delta, ndemes, pm, demes.prior,
                                    dyn_L, hosts) {


  obs <- ttree$obs
  demes <- ttree$demes
  ttree <- ttree$ttree

  omega <- fn_list$omega
  omega_bar <- fn_list$omega_bar
  pit <- fn_list$pit
  gamma_prob <- fn_list$gamma_prob

  nhosts <- length(ttree[, 1])

  lul <- numeric(ndemes)

  todo <- hosts

  while(length(todo) > 0) {

    todo <- todo[order(ttree[todo, 1], decreasing = T)]

    i <- todo[1]

    dyn_L[i, ] <- rep(-Inf, ndemes)

    # Leaf
    if (length(which(ttree[, 3] == i)) == 0) {

      if (!is.na(demes[i])) {

        if (demes[i] > 0) {

          dyn_L[i, demes[i]] <- dyn_T(ttree, obs, grid, omega, omega_bar, pit, off.r, off.p, pi, ws.shape,
                                        ws.scale, obs.start, obs.end, grid.delta, i, demes[i])

        } else {

          for (d in 1:ndemes) {

            dyn_L[i, d] <- dyn_T(ttree, obs, grid, omega, omega_bar, pit, off.r, off.p, pi, ws.shape,
                                   ws.scale, obs.start, obs.end, grid.delta, i, d)

          }

        }

      } else {

        for (d in 1:ndemes) {

          dyn_L[i, d] <- dyn_T(ttree, obs, grid, omega, omega_bar, pit, off.r, off.p, pi, ws.shape,
                                 ws.scale, obs.start, obs.end, grid.delta, i, d)

        }

      }

      # Not a leaf
    } else {

      children <- which(ttree[, 3] == i)

      if (!is.na(demes[i])) {

        if (demes[i] > 0) {

          dyn_L[i, demes[i]] <- dyn_T(ttree, obs, grid, omega, omega_bar, pit, off.r, off.p, pi, ws.shape,
                                        ws.scale, obs.start, obs.end, grid.delta, i, demes[i])

          for (j in children) {

            for (d2 in 1:ndemes) {

              lul[d2] <- dyn_U(ttree, grid, omega, w.shape, w.scale, obs.end, grid.delta, i, demes[i], j, d2, pm) +
                dyn_L[j, d2]

            }

            m <- max(lul)

            dyn_L[i, demes[i]] <- dyn_L[i, demes[i]] +
              log(sum(exp(lul - m))) + m

          }

        } else {

          for (d in 1:ndemes) {

            dyn_L[i, d] <- dyn_T(ttree, obs, grid, omega, omega_bar, pit, off.r, off.p, pi, ws.shape,
                                   ws.scale, obs.start, obs.end, grid.delta, i, d)

            for (j in children) {

              for (d2 in 1:ndemes) {

                lul[d2] <- dyn_U(ttree, grid, omega, w.shape, w.scale, obs.end, grid.delta, i, d, j, d2, pm) +
                  dyn_L[j, d2]

              }

              m <- max(lul)

              dyn_L[i, d] <- dyn_L[i, d] +
                log(sum(exp(lul - m))) + m

            }

          }

        }

      } else {

        for (d in 1:ndemes) {

          dyn_L[i, d] <- dyn_T(ttree, obs, grid, omega, omega_bar, pit, off.r, off.p, pi, ws.shape,
                                 ws.scale, obs.start, obs.end, grid.delta, i, d)

          for (j in children) {

            for (d2 in 1:ndemes) {

              lul[d2] <- dyn_U(ttree, grid, omega, w.shape, w.scale, obs.end, grid.delta, i, d, j, d2, pm) +
                dyn_L[j, d2]

            }

            m <- max(lul)

            dyn_L[i, d] <- dyn_L[i, d] +
              log(sum(exp(lul - m))) + m

          }

        }

      }


    }

    if (ttree[i, 3] > 0) {

      if (!(ttree[i, 3] %in% todo)) {

        todo <- c(todo, ttree[i, 3])

      }

    }

    todo <- todo[-1]

  }

  root_host <- which(ttree[, 3] == 0)

  lpl <- log(demes.prior) + dyn_L[root_host, ]

  rootm <- max(lpl)

  log_likelihood <- log(sum(exp(lpl - rootm))) + rootm

  return(list(loglik = log_likelihood, dyn_L = dyn_L))

}

