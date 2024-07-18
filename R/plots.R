#' Plot a transmission tree in an economic format
#' @param ttree Transmission tree
#' @param showLabels Boolean for whether or not to show the labels
#' @param maxTime Maximum value of time to show on x axis
#' @param cex Expansion factor
#' @export
plotTTree_summ = function(ttree, showLabels = T, maxTime = NA, cex = 1) {

  nam <- ttree$nam
  ttree <- ttree$ttree

  #Determine ys
  n <- nrow(ttree)

  ys <- rep(0, n)
  scale <- rep(1, n)

  todo <- c(which(ttree[, 3] == 0))
  while (length(todo) > 0) {

    f <- which(ttree[, 3] == todo[1])
    o <- rank(-ttree[f, 1], ties.method = "first")
    f[o] <- f

    for (i in f) {

      ys[i] <- ys[todo[1]] + scale[todo[1]] * which(f == i) / (length(f) + 1)
      scale[i] <- scale[todo[1]] / (length(f) + 1)
      todo <- c(todo, i)

    }

    todo <- todo[-1]

  }

  ys <- rank(ys)


  #Do the plot
  oldpar <- par('yaxt', 'bty')
  on.exit(par(oldpar))
  par(yaxt='n', bty='n')

  mi <- min(ttree[which(!is.na(ttree[, 1])), 1])
  ma <- max(ttree[which(!is.na(ttree[, 1])), 1])

  if (!is.na(maxTime)) {

    ma <- maxTime

  }

  plot(c(), c(), xlim = c(mi, ma), ylim=c(0, n + 1), xlab = '', ylab = '')

  for (i in 1:n) {

    if (ttree[i, 3] != 0) {

      arrows(ttree[ttree[i, 3], 1], ys[ttree[i, 3]], ttree[i, 1], ys[i], length = 0)

    }

    if (showLabels && ttree[i, 2] > 0) {

      text(ttree[i, 1], ys[i], i, pos = 4, cex = cex)

    }

  }

  for (i in 1:n) {

    points(ttree[i, 1], ys[i], pch = 21, bg = ifelse(ttree[i, 2] == 0, 'white', 'black'), cex = cex)

  }

  return(invisible(ttree))

}


#' Plot a transmission tree in a detailed format
#'
#' @param ttree Transmission tree
#' @param w.shape Shape parameter of the Gamma probability density function representing the generation time
#' @param w.scale Scale parameter of the Gamma probability density function representing the generation time
#' @param showLabels Whether or not to show the labels
#' @param maxTime Maximum value of time to show on x axis
#' @param cex Expansion factor
#' @export
plotTTree_det = function(ttree, w.shape, w.scale, showLabels = TRUE, maxTime = NA, cex=1) {

  nam <- ttree$nam
  obs <- ttree$obs
  ttree <- ttree$ttree

  n <- nrow(ttree)

  #Determine ys
  ys <- rep(0, n)
  scale <- rep(1, n)

  todo <- c(which(ttree[, 3] == 0))
  while (length(todo) > 0) {

    f <- which(ttree[, 3] == todo[1])
    o <- rank(-ttree[f, 1])
    f[o] <- f

    for (i in f) {

      ys[i] <- ys[todo[1]] + scale[todo[1]] * which(f == i) / (length(f) + 1)
      scale[i] <- scale[todo[1]] / (length(f) + 1)
      todo <- c(todo, i)

    }

    todo <- todo[-1]

  }

  ys <- rank(ys)

  oldpar <- par('yaxt', 'bty')
  on.exit(par(oldpar))
  par(yaxt='n', bty='n')

  mi <- min(ttree[, 1])

  if (is.na(maxTime)) {

    ma <- max(obs[, 1])

  }  else {

    ma <- maxTime

  }

  xstep <- (ma - mi) / 2000

  plot(c(), c(), xlim = c(mi - (ma - mi) * 0.05, ma + (ma - mi) * 0.05), ylim=c(0, n + 1), xlab = '', ylab = '')

  maxcol <- max(dgamma(seq(0, ma - mi, xstep), shape = w.shape, scale = w.scale))

  for (i in 1:n) {

    as <- seq(ttree[i, 1], ma, xstep)
    bs <- rep(ys[i], length(as))
    cs <- abs((maxcol - dgamma(as - ttree[i, 1], shape = w.shape, scale = w.scale)) / maxcol)
    cs <- grDevices::gray(cs)

    segments(as, bs, x1 = as + xstep, col = cs)

    host_obs <- which(obs[, 2] == i)

    if (length(host_obs) > 0) {

      obs_times <- obs[host_obs, 1]

      points(obs_times, rep(ys[i], length(obs_times)), col = 'red')

    }

    if (showLabels) {

      text(ma + (ma - mi) * 0.05, ys[i], i, cex = cex)

    }

    if (ttree[i, 3] == 0) {

      next

    }

    arrows(ttree[i, 1], ys[ttree[i, 3]], ttree[i, 1], ys[i], length = 0.1)

  }

  return(invisible(ttree))

}


#' Plot MCMC traces
#'
#' @param record Output from inferTTree function
#' @param burnin Proportion of the MCMC output to be discarded as burnin
#' @return Returns invisibly the first parameter
#' @export
plotTraces <- function(record, burnin = 0) {

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  par(mfrow = c(3, 2))

  record <- record[max(1, round(length(record) * burnin)):length(record)]


  rho_res <- sapply(record, function(x) x$rho)

  if (length(unique(rho_res)) > 1) {

    plot(rho_res, ylab = 'rho',
         xlab = 'MCMC iterations', type = 'l')

  } else {

    plot(sapply(record, function(x) x$pTTree + x$pPTree), ylab = 'Posterior probability',
         xlab = 'MCMC iterations', type = 'l')

  }

  plot(sapply(record, function(x) x$pi), ylab = 'Sampling proportion pi',
       xlab = 'MCMC iterations', type = 'l')

  plot(sapply(record, function(x) x$off.r), ylab = 'off.r',
       xlab = 'MCMC iterations', type = 'l')

  plot(sapply(record, function(x) x$off.p), ylab = 'off.p',
         xlab = 'MCMC iterations', type = 'l')

  plot(sapply(record, function(x) x$kappa), ylab = 'Within-host initial population kappa',
       xlab = 'MCMC iterations', type = 'l')

  plot(sapply(record, function(x) x$lambda), ylab = 'Within-host population growth rate lambda',
       xlab = 'MCMC iterations', type = 'l')

  return(invisible(record))

}



#' Plot both phylogenetic and transmission trees using colors on the phylogeny
#' @param tree Combined phylogenetic/transmission tree
#' @param showLabels Whether or not to show the labels
#' @param showStars Whether or not to show stars representing transmission events
#' @param cols Colors to use for hosts
#' @param maxTime Maximum time to show on the x axis
#' @param cex Expansion factor
#' @return Returns invisibly the first parameter
#' @examples
#' plotCTree(simulateOutbreak())
#' @export
plotCTree = function(tree,showLabels=TRUE,showStars=TRUE,cols=NA,maxTime=NA,cex=1)  {
  nam=tree$nam
  tree=tree$ctree
  nsam <- sum(tree[ ,2]+tree[ ,3] == 0)
  nh <- nrow(tree)-3*nsam+1
  ntot <- nsam+nh
  oldpar <- par('yaxt','bty','xpd')
  on.exit(par(oldpar))
  par(yaxt='n',bty='n',xpd=T)
  plot(0,0,type='l',xlim=c(min(tree[,1]),ifelse(is.na(maxTime),max(tree[,1]),maxTime)),ylim=c(0,nsam+1),xlab='',ylab='')
  host <- tree[ ,4]
  if (ntot>1) {
    if (is.na(cols[1])) grDevices::palette(grDevices::rainbow(min(1024,ntot)))#Need as many unique colors as there are hosts. If there are more than 1024 hosts, colors are recycled.
    else grDevices::palette(cols)
  }

  #Determine ys for leaves
  root<-which(host==0)
  ys <- matrix(0, nsam, 1)
  todo <- cbind(root,0,0.5,1)#Matrix of nodes to do,with associated starting x and y coordinates and scale
  while (nrow(todo) > 0)  {
    w <- todo[1,1]
    x <- todo[1,2]
    y <- todo[1,3]
    scale <- todo[1,4]
    if (tree[w,2] == 0 && tree[w,3] == 0)  {
      #Leaf node
      ys[w] <- y
    } else if (tree[w,3] == 0)  {
      #Transmission node
      todo <- rbind(todo,cbind(tree[w,2],tree[w,1],y,scale,deparse.level=0))
    } else {
      #Binary node
      todo <- rbind(todo,cbind(tree[w,2],tree[w,1],y + scale/2,scale/2,deparse.level=0),cbind(tree[w,3],tree[w,1],y-scale/2,scale/2,deparse.level=0))
    }
    todo <- rbind(todo[-1,])
  }
  ys<-rank(ys)

  #Determine ys for non-leaves
  for (i in ((nsam+1):nrow(tree))) {
    children <- c()
    todo <- i
    while (length(todo)>0) {
      children=c(children,todo[1])
      todo=c(todo[-1],setdiff(tree[todo[1],2:3],0))
    }
    ys[i] <- mean(ys[children[which(children<=nsam)]])
  }

  todo <- cbind(root,tree[root,1])
  while (nrow(todo) > 0)  {
    w <- todo[1,1]
    x <- todo[1,2]
    y <- ys[w]
    col=host[w]
    if (tree[w,2] == 0 && tree[w,3] == 0)  {
      #Leaf node
      lines(c(x,tree[w,1]),c(y,y),col=col,lwd=2)
      if (showLabels) text(tree[w,1],y,nam[w],cex=cex,pos=4)
    } else if (tree[w,3] == 0)  {
      #Transmission node
      lines(c(x,tree[w,1]),c(y,y),col=col,lwd=2)
      #points(tree[w,1],y,col = 'red',pch=8)
      todo <- rbind(todo,cbind(tree[w,2],tree[w,1]))
    } else {
      #Binary node
      lines(c(x,tree[w,1]),c(y,y),col=col,lwd=2)
      lines(c(tree[w,1],tree[w,1]),cbind(ys[tree[w,2]],ys[tree[w,3]]),col=col,lwd=2)
      todo <- rbind(todo,cbind(tree[w,2],tree[w,1]),cbind(tree[w,3],tree[w,1]))
    }
    todo <- rbind(todo[-1,])
  }

  todo <- cbind(root,tree[root,1])
  while (nrow(todo) > 0 && showStars)  {
    w <- todo[1,1]
    x <- todo[1,2]
    y <- ys[w]
    col=host[w]
    if (tree[w,2] == 0 && tree[w,3] == 0)  {
      #Leaf node
    } else if (tree[w,3] == 0)  {
      #Transmission node
      points(tree[w,1],y,col = 'red',pch=8)
      todo <- rbind(todo,cbind(tree[w,2],tree[w,1]))
    } else {
      #Binary node
      todo <- rbind(todo,cbind(tree[w,2],tree[w,1]),cbind(tree[w,3],tree[w,1]))
    }
    todo <- rbind(todo[-1,])
  }
  return(invisible(tree))
}



#' Plot both phylogenetic and transmission trees using colors on the phylogeny
#' @param tree Combined phylogenetic/transmission tree
#' @param showLabels Whether or not to show the labels
#' @param showStars Whether or not to show stars representing transmission events
#' @param cols Colors to use for hosts
#' @param maxTime Maximum time to show on the x axis
#' @param cex Expansion factor
#' @export
plotdemesCTree = function(tree,showLabels=TRUE,showStars=TRUE,cols=NA,maxTime=NA,cex=1)  {
  nam=tree$nam
  demes=tree$demes
  tree=tree$ctree
  nsam <- sum(tree[ ,2]+tree[ ,3] == 0)
  nh <- nrow(tree)-3*nsam+1
  ntot <- nsam+nh
  oldpar <- par('yaxt','bty','xpd')
  on.exit(par(oldpar))
  par(yaxt='n',bty='n',xpd=T)
  plot(0,0,type='l',xlim=c(min(tree[,1]),ifelse(is.na(maxTime),max(tree[,1]),maxTime)),ylim=c(0,nsam+1),xlab='',ylab='')
  host <- tree[ ,4]
  if (ntot>1) {
    if (is.na(cols[1])) grDevices::palette(grDevices::rainbow(min(1024,length(unique(demes)))))#Need as many unique colors as there are hosts. If there are more than 1024 hosts, colors are recycled.
    else grDevices::palette(cols)
  }

  #Determine ys for leaves
  root<-which(host==0)
  ys <- matrix(0, nsam, 1)
  todo <- cbind(root,0,0.5,1)#Matrix of nodes to do,with associated starting x and y coordinates and scale
  while (nrow(todo) > 0)  {
    w <- todo[1,1]
    x <- todo[1,2]
    y <- todo[1,3]
    scale <- todo[1,4]
    if (tree[w,2] == 0 && tree[w,3] == 0)  {
      #Leaf node
      ys[w] <- y
    } else if (tree[w,3] == 0)  {
      #Transmission node
      todo <- rbind(todo,cbind(tree[w,2],tree[w,1],y,scale,deparse.level=0))
    } else {
      #Binary node
      todo <- rbind(todo,cbind(tree[w,2],tree[w,1],y + scale/2,scale/2,deparse.level=0),cbind(tree[w,3],tree[w,1],y-scale/2,scale/2,deparse.level=0))
    }
    todo <- rbind(todo[-1,])
  }
  ys<-rank(ys)

  #Determine ys for non-leaves
  for (i in ((nsam+1):nrow(tree))) {
    children <- c()
    todo <- i
    while (length(todo)>0) {
      children=c(children,todo[1])
      todo=c(todo[-1],setdiff(tree[todo[1],2:3],0))
    }
    ys[i] <- mean(ys[children[which(children<=nsam)]])
  }

  todo <- cbind(root,tree[root,1])
  while (nrow(todo) > 0)  {
    w <- todo[1,1]
    x <- todo[1,2]
    y <- ys[w]


  if (length(demes[host[w]]) > 0) {

      if (is.na(demes[host[w]])) {

        col = 'grey80'

      } else {

        col=demes[host[w]]

      }

  } else {

    col = 'grey80'

  }

    if (tree[w,2] == 0 && tree[w,3] == 0)  {
      #Leaf node
      lines(c(x,tree[w,1]),c(y,y),col=col,lwd=2)
      if (showLabels) text(tree[w,1],y,nam[w],cex=cex,pos=4)
    } else if (tree[w,3] == 0)  {
      #Transmission node
      lines(c(x,tree[w,1]),c(y,y),col=col,lwd=2)
      #points(tree[w,1],y,col = 'red',pch=8)
      todo <- rbind(todo,cbind(tree[w,2],tree[w,1]))
    } else {
      #Binary node
      lines(c(x,tree[w,1]),c(y,y),col=col,lwd=2)
      lines(c(tree[w,1],tree[w,1]),cbind(ys[tree[w,2]],ys[tree[w,3]]),col=col,lwd=2)
      todo <- rbind(todo,cbind(tree[w,2],tree[w,1]),cbind(tree[w,3],tree[w,1]))
    }
    todo <- rbind(todo[-1,])
  }

  todo <- cbind(root,tree[root,1])
  while (nrow(todo) > 0 && showStars)  {
    w <- todo[1,1]
    x <- todo[1,2]
    y <- ys[w]


  if (length(demes[host[w]]) > 0) {

      if (is.na(demes[host[w]])) {

        col = 'grey80'

      } else {

        col=demes[host[w]]

      }

  } else {

    col = 'grey80'

  }

    if (tree[w,2] == 0 && tree[w,3] == 0)  {
      #Leaf node
    } else if (tree[w,3] == 0)  {
      #Transmission node
      points(tree[w,1],y,col = 'red',pch=8)
      todo <- rbind(todo,cbind(tree[w,2],tree[w,1]))
    } else {
      #Binary node
      todo <- rbind(todo,cbind(tree[w,2],tree[w,1]),cbind(tree[w,3],tree[w,1]))
    }
    todo <- rbind(todo[-1,])
  }
  return(invisible(tree))
}
