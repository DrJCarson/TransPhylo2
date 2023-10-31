#' Converts an ape phylo object into a phylogenetic tree with host information
#'
#' @param tr phylo object
#' @param dateLastSample date of the last sample
#' @param host Optional host information, if not present the leaf names needs to contain this, eg 1.1 etc
#' @return phylogenetic tree
#'
#' @export
ptreeFromPhylo <- function(tr,dateLastSample,host) {
  if (!missing(host)) {
    r=rank(unique(host))
    names(r)=unique(host)
    h=r[as.character(host)]
    names(h)=names(host)
    tr$tip.label=sprintf('%d.%s',h[tr$tip.label],tr$tip.label)
  }
  ptree=ptreeFromPhylo(tr,dateLastSample = dateLastSample)
  ptree$host=as.numeric(unlist(strsplit(ptree$nam,'\\.'))[seq(1,length(ptree$nam)*2,2)]) #Extract host information from leaf names
  n=length(ptree$host) # Number of leaves
  w=(n+1):nrow(ptree$ptree) # Indexes of internal nodes
  o=order(ptree$host*1e10+ptree$ptree[1:n,1]) #Order needed for leaves: by host first and then from oldest to most recent
  o2=order(-ptree$ptree[w,1]) #Order needed for internal nodes: from most recent to oldest (root)
  o=c(o,n+o2)
  ot=o;ot[o]=1:length(o)
  ptree$host=ptree$host[o[1:n]] #Reorder host
  ptree$nam=ptree$nam[o[1:n]] #Reorder nam
  ptree$ptree=ptree$ptree[o,] #Reorder ptree
  ptree$ptree[w,2]=ot[ptree$ptree[w,2]]
  ptree$ptree[w,3]=ot[ptree$ptree[w,3]]
  class(ptree)<-'ptree'
  return(ptree)
}

