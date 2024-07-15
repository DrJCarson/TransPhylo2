#' Converts an ape phylo object into a phylogenetic tree with host information
#'
#' @param tr phylo object
#' @param dateLastSample date of the last sample
#' @param host Optional host information, if not present the leaf names needs to contain this, eg 1.1 etc
#' @param demes Optional location information, demes must be categorized as contiguous integers. NA can be used for missing data.
#' @return phylogenetic tree
#'
#' @export
ptreeFromPhylo <- function(tr, dateLastSample, host, demes) {

#  if (!missing(host)) {

#    r = rank(unique(host))

#    names(r) = unique(host)

#    h = r[as.character(host)]

##    names(h)=names(host)

#    tr$tip.label = sprintf('%d.%s', h[tr$tip.label], tr$tip.label)

#  }

  if (dateLastSample > 1900 && dateLastSample < 2100 && sum(tr$edge.length) < 1) {

    warning('Warning: input tree has small branch lengths. This needs to be a dated tree in the same unit as dateLastSample (eg years).\n')


  }

  n <- length(tr$tip.label)
  ed <- tr$edge
  le <- tr$edge.length

  tra <- c(1:n, (2*n - 1):(n + 1))
  ptree <- matrix(0, 2 * n - 1, 3)
  if (n == 1) {

    ptree[1, 1] = dateLastSample

  }

  for (i in 1:nrow(ed)) {

    father <- tra[ed[i, 1]]
    son <- tra[ed[i, 2]]
    if (ptree[father, 2] == 0) {

      ptree[father, 2] = son

    } else {

      ptree[father, 3] = son

    }

    ptree[son, 1]<-le[i]

  }

  todo <- 2 * n - 1
  while (length(todo) > 0) {

    t1 = todo[1]

    if (ptree[t1, 2] == 0) {

      todo = todo[-1]

      next

    }

    ptree[ptree[t1, 2], 1] <- ptree[ptree[t1, 2], 1] + ptree[t1, 1]
    ptree[ptree[t1, 3], 1] <- ptree[ptree[t1, 3], 1] + ptree[t1, 1]
    todo = c(todo[-1], ptree[t1, 2], ptree[t1, 3])

  }

  ptree[, 1] = ptree[, 1] - max(ptree[, 1]) + dateLastSample
  ptree[, 2:3] = ptree[, 3:2]

  l = list(ptree = ptree, nam = tr$tip.label)
  class(l) <- 'ptree'

  ptree <- l

#  ptree$host=as.numeric(unlist(strsplit(ptree$nam,'\\.'))[seq(1,length(ptree$nam)*2,2)]) #Extract host information from leaf names

  if (missing(host)) {

    ptree$host <- 1:length(tr$tip.label)

  } else {

    ptree$host <- host

  }

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

  if (!missing(demes)) {

    ptree$demes = demes[o[1:n]]

  }


  class(ptree)<-'ptree'
  return(ptree)

}

#' Converts a phylogenetic tree into an ape phylo object
#' @param ptree phylogenetic tree
#' @return phylo object
#' @examples
#' phyloFromPTree(extractPTree(simulateOutbreak()))
#' @export
phyloFromPTree <- function(ptree) {
  nam=ptree$nam
  ptree=ptree$ptree
  n<-ceiling(nrow(ptree)/2)
  if (n==1) return(ape::read.tree(text='(1);'))
  tr<-list()
  tr$Nnode<-n-1
  tr$tip.label<-nam
  tr$edge<-matrix(0,n*2-2,2)
  tr$edge.length<-rep(0,n*2-2)
  iedge<-1
  root<-which(ptree[,1]==min(ptree[,1]))
  tra<-c(1:n,root,setdiff((n+1):(2*n-1),root))
  tra2<-1:length(tra)
  tra[tra]<-tra2
  for (i in (n+1):(2*n-1)) {
    tr$edge[iedge,]<-c(tra[i],tra[ptree[i,3]])
    tr$edge.length[iedge]<-ptree[ptree[i,3],1]-ptree[i,1]
    iedge<-iedge+1
    tr$edge[iedge,]<-c(tra[i],tra[ptree[i,2]])
    tr$edge.length[iedge]<-ptree[ptree[i,2],1]-ptree[i,1]
    iedge<-iedge+1
  }
  tr$root.time=min(ptree[,1])
  class(tr)<-'phylo'
  tr=ape::reorder.phylo(tr, order = "cladewise")
  return(tr)
}
