#' Calculate calculate T for given host and deme
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

    Tns <- ((pi[deme] * dgamma(obs_time - inf_time, shape = ws.shape, scale = ws.scale)) / (1 - omega_int)) *
      alpha_sum * factorial(inc_off)

  # Host is unobserved
  } else {

    Tns <- ((1 - pit_int) / (1 - omega_int)) * alpha_sum * factorial(inc_off)

  }

  return(Tns)

}


#' Calculate calculate U for given hosts and demes
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

  tidx1 <- 1 + floor((obs.end + 1e-10 - inf_time1) / grid.delta)
  tidx2 <- 1 + floor((obs.end + 1e-10 - inf_time2) / grid.delta)

  omega_int2 <- omega[tidx2, deme2] + ((inf_time2 - grid[tidx2]) / (grid[tidx2 + 1] - grid[tidx2])) * (omega[tidx2 + 1, deme2] - omega[tidx2, deme2])

  Uns <- (1 - omega_int2) * dgamma(inf_time2 - inf_time1, shape = w.shape, scale = w.scale) * pm[deme1, deme2]

  return(Uns)

}
