#' Build a consensus transmission tree from a MCMC output
#' @param record Output from inferTTree function
#' @param burnin Proportion of the MCMC output to be discarded as burnin
#' @return The consensus transmission tree
#' @export
consTTree = function(record,burnin=0.5)
{
  if (burnin>0) record=record[round(length(record)*burnin):length(record)]
  v=sapply(record, function(x) x$pTTree + x$pPTree)
  cons=record[[which.max(v)]]$ctree
  return(cons)
}
