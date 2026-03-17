#' @name TransPhylo2
#' @title TransPhylo2
#'
#' @description TransPhylo2 is a R package extending TransPhylo.
#'
#' @author Jake Carson \email{jake.carson@warwick.ac.uk}
#'
#' @import stats
#' @import graphics
#' @import ape
#' @importFrom utils combn getFromNamespace
NULL


#' @useDynLib TransPhylo2, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL


.onUnload <- function (libpath) {
  library.dynam.unload("TransPhylo2", libpath)
}
NULL
