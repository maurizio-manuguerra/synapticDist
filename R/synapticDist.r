test <- function(dists=0:20){
  i0=4
  ii=(1:6)[-4]
  m=length(ii)
  fits.list=list()
  for (i in 1:m){
    cat("Processing neuron ", i, "of ", m, "neurons.\n")
    Ds=matrix(0,nrow=length(dists),ncol=m)
    Ds[,i] = dists
    Ds=split(Ds,row(Ds))
    fits.list[[i]]=do_fits(i0,ii,Ds)
  }
  names(fits.list) = 1:m
  return(fits.list)
}


#' @title do_fits function
#' for a range of distances Ds runs bulkem2 on matrices Nx(M+1), where M is the number of neurons inferring on the target.
#' the output is a list of the same length of Ds, where each element is the results of the analysis performed on the data sets obtained using Ds in data.gen
#' @param i0 the target neuron (index)
#' @param ii the inferring neurons (a vector of indexes)
#' @param Ds a list of M-vectors, giving the distances of the M neurons from the target
#' @param random.inits the number of random initializations to star the EM algorithm. Defauls to 3.
#' @param plot a logical variable. If TRUE (the default) a plot of the values of the log-likelihoods for each vector of Ds id created.
do_fits <- function(i0, ii, Ds, random.inits = 3, plot=T){
  t0 <- data.sim[which(data.sim[,2]==i0),1]
  ti <- lapply(ii, function(i)data.sim[which(data.sim[,2]==i),1])
  xlist <- lapply(Ds, function(d)data.gen(t0, ti, d))
  fits <- bulkem2(datasets=xlist, num.components=length(ii)+1, random.inits = random.inits, use.gpu = F, verbose=FALSE)
  if (plot) do_plots(fits, fits.names=paste(Ds))
  return(fits)
}

do_plots <- function(fits, fits_names=NULL){
  n <- length(fits)
  cols <- rainbow(n)
  if (is.null(fits_names)) fits_names=names(fits)
  if (is.null(fits_names)) fits_names=paste(1:n)
  fits_pch <- 1:n
  par(mar=c(5.1, 4.1, 4.1, 12.1))
  plot(sapply(fits,function(x)x$llik), col=cols, pch=fits_pch, ylab="log-likelihood", main=paste(length(ti)," neurons inferring on #",i0,sep=''))
  legend("topright", inset=c(-0.3*length(Ds),0), legend=fits_names, col=cols, pch=fits_pch, bg="white", xpd=TRUE)
  par(mar=c(5.1, 4.1, 4.1, 2.1))
}

test.IG.exp <- function(i1,i2, num.components=2, random.inits=1){
  #library(bulkem, lib.loc='/Library/Frameworks/R.framework/Versions/3.2/Resources/')
  #devtools::install("/Users/manuguerra/Documents/git/bulkem")
  library(statmod)
  source("/Users/manuguerra/Documents/git/bulkem/R/bulkem.R")
  source("/Users/manuguerra/Documents/git/bulkem/R/invgaussmixEM.R")
  Ds <- 0:20
  t1 <- data.sim[which(data.sim[,2]==i1),1]
  t2 <- data.sim[which(data.sim[,2]==i2),1]
  xlist=lapply(Ds, function(d)data.gen(t1, t2, d, num.components=num.components))

  #fits <- bulkem2(datasets=xlist)
  #print('--- fits ---')
  #print(fits)
  fits <- bulkem3(datasets=xlist, num.components=num.components, random.inits=random.inits, use.gpu=FALSE, verbose=TRUE)
  plot(Ds,sapply(fits,function(x)x$llik))
  lines(Ds,sapply(fits,function(x)x$llik))
  return(fits)
}

test.IG.gamma <- function(i1,i2, num.components=2, random.inits=1){
  #library(bulkem, lib.loc='/Library/Frameworks/R.framework/Versions/3.2/Resources/')
  #devtools::install("/Users/manuguerra/Documents/git/bulkem")
  library(statmod)
  source("/Users/manuguerra/Documents/git/bulkem/R/bulkem.R")
  source("/Users/manuguerra/Documents/git/bulkem/R/invgaussmixEM.R")
  Ds <- 0:20
  t1 <- data.sim[which(data.sim[,2]==i1),1]
  t2 <- data.sim[which(data.sim[,2]==i2),1]
  xlist=lapply(Ds, function(d)data.gen(t1, t2, d, num.components=num.components))

  #fits <- bulkem2(datasets=xlist)
  #print('--- fits ---')
  #print(fits)
  fits <- bulkem4(datasets=xlist, num.components=num.components, random.inits=random.inits, use.gpu=FALSE, verbose=TRUE)
  plot(Ds,sapply(fits,function(x)x$llik))
  lines(Ds,sapply(fits,function(x)x$llik))
  return(fits)
}


check.fit <- function(fit, mfrow=T){
  require(statmod)
  m=ncol(fit$x)
  if(mfrow) par(mfrow=c(m,1))
  #belongs_to <- apply(fit$member.prob,1,which.max)
  belongs_to <- apply(fit$member.prob, 1, function(x)sample(1:m, size=1, prob=x))
  for (j in 1:m){
    end.scale = as.integer(qinvgauss(.95,mean=fit$mu[j],shape=fit$lambda[j]))
    ii <- which(belongs_to == j)
    f2 <- fit$x[ii,j]
    hist(f2[f2<end.scale],end.scale, main=j)
    xx=0:end.scale
    fx=dinvgauss(xx,mean=fit$mu[j],shape=fit$lambda[j])
    scale <- length(ii)
    lines(xx,fx*scale, col="red")
  }
  if(mfrow) par(mfrow=c(1,1))
}

check.fit.IG.exp <- function(fit){
  require(statmod)
  m=length(fit$mu)
  par(mfrow=c(m,1))
  belongs_to <- apply(fit$member.prob,1,which.max)
  end.scale = c(
    as.integer(qinvgauss(.95,mean=fit$mu[1],shape=fit$lambda[1])),
    as.integer(qexp(.95,rate=1/fit$mu[2]))
  )
  fx= matrix(c(
    dinvgauss(xx,mean=fit$mu[1],shape=fit$lambda[1]),
    dexp(xx,rate=1/fit$mu[2])
    ), ncol=2, byrow=F)
  for (j in 1:m){
    ii <- which(belongs_to == j)
    scale <- length(ii)
    if (scale>0){
      f2 <- fit$x[ii,j]
      hist(f2[f2<end.scale[j]],end.scale[j], main=j)
      xx=0:end.scale[j]
      lines(xx,fx[,j]*scale, col="red")
    }
  }
  par(mfrow=c(1,1))
}


setup <- function(home.dir='.',code.dir='./R/', data.dir='./Data/', res.dir='./Results/'){
	library(Matrix)
	#library(gputools)
	library(glmnet)
	library(speedglm)
	library(statmod)
	library(mixtools)
	#Test the script fitting a subset of the generated neurons.
	#####
	setwd(code.dir)
	source('functions.r')
	#setwd(data.dir)
	setwd(home.dir)
	rr <<- list_rats_recs2(data.dir)
	nrecs <<- length(rr[[1]])
	code.dir <<- code.dir
	data.dir <<- data.dir
	res.dir <<- res.dir
	home.dir <<- home.dir
	tolerance <<- 1
}

list_rats_recs2 <- function(dir0='../../data'){
	rat=NULL
	rec=NULL
	l1=list.files(path=dir0)
	for (irat in 1:length(l1)){
		dir1=paste(dir0,l1[irat],sep='/')
		l2=list.files(path=dir1)
		rat=c(rat,rep(l1[irat],length(l2)))
		rec=c(rec,l2)
	}
	return(list(rat,rec))
}

which.machine <- function(){
	t1 <- try(system("uname -n", intern = TRUE))
	if (t1=="Maurizio-Manuguerras-Mac-mini.local"){
		m="mac"
	} else if (t1=="master185") {
		m="cuda"
	} else if (t1=="manuguerra-ubuntu"	) {
		m="server"
	} else if (t1=="hpcc.efs.mq.edu.au"	) {
		m="cluster"
	} else if (t1=="iMac-Home.local") {
		m="imac"
	} else {
		m=t1
	}
	return(m)
}

#' @title Returns the interspike intervals of m neurons firing (ti) on a specific neuron (t0)
#' Takes 1+m series of spike times,  and computes a matrix nX(m+1), where m is the number of neuron firing on neur 0, and n is the number of isi of neur 0.
#' d is an m-vector
#' Cond1: If more than a neur1 spike occur in a given isi of neur0, only the last spike is considered.
#' Cond2: If no spikes from neur1 occur in a given isi of neur0, t_12 = 0.
#' TODO: optimize the loop
data.gen <- function(t0, ti, d){
  m=length(ti)
  ti=lapply(1:m, function(i)ti[[i]]+d[i]) #shift/delay
  l0 <- length(t0)
  isi <- matrix(0, nrow=l0-1, ncol=m+1)
  isi[,1] <- t0[2:l0]-t0[1:(l0-1)]
  for (j in 1:m){
    t1 <- ti[[j]]
    ii <- findInterval(t1,t0)
  #%%%%%%%%%%%%%%%%%%%%%%
    extremes <- (ii==0 | ii == length(t0))
    t1 <- t1[!extremes]
    ii <- ii[!extremes]
  #%%%%%%%%%%%%%%%%%%%%%%
    dupl <- duplicated(ii, fromLast=T) #See description (cond1)
    t1 <- t1[!dupl]
    ii <- ii[!dupl]
  #%%%%%%%%%%%%%%%%%%%%%%
    l1 <- length(ii)
    isi[ii,j+1] <- t0[ii+1] - t1
  #%%%%%%%%%%%%%%%%%%%%%%
  }
  return(isi)
}




AIC.neuromixEM <- function(fits){
  sapply(fits, function(x) -2*x$llik + 2*(1+length(x$mu[x$mu!=0])+length(x$lambda[x$lambda!=0])))
}
BIC.neuromixEM <- function(fits){
  sapply(fits, function(x) -2*x$llik + log(nrow(x$x))*(1+length(x$mu[x$mu!=0])+length(x$lambda[x$lambda!=0])))
}
AICem.neuromixEM <- function(fits){
  sapply(fits, function(x) -2*x$Q)
}
AICem.ww.neuromixEM <- function(fits){ #http://www.smu.edu/~/media/Site/Dedman/Departments/Statistics/TechReports/TR303.ashx?la=en
  sapply(fits, function(x) {m=ncol(x$x); d=2; return(-2*x$llik + 2*((m-1)+d*m+(d*m*(m-1)/2)))})
}

compareICs <- function(fits.list, criterion="AIC"){
  if (criterion=='AIC')
    ics <- sapply(fits.list, AIC.neuromixEM)
  else if (criterion=='BIC')
    ics <- sapply(fits.list, BIC.neuromixEM)
  else if (criterion=='AICem')
    ics <- sapply(fits.list, AICem.neuromixEM)
  else if (criterion=='AICem.ww')
    ics <- sapply(fits.list, AICem.ww.neuromixEM)
  else
    stop("Criterion not implemented")
  n <- ncol(ics)
  cols <- rainbow(n)
  ics_names <- colnames(ics)
  if (any(is.null(ics_names))) ics_names <- paste(1:n)
  matplot(ics, col=cols, lty=1, t='l')
  legend("topleft", col=cols, ics_names, bg="white", lwd=1)
}
