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
