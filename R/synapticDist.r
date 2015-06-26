test_extended <- function(ii, i0, d, random.inits = 3, plot=T){
  t0 <- data.sim[which(data.sim[,2]==i0),1]
  ti <- lapply(ii, function(i)data.sim[which(data.sim[,2]==i),1])
  isi=data.gen.extended(ti,t0,d)
  xlist=list(isi)
  fits <- bulkem2(datasets=xlist, num.components=length(ii)+1, random.inits = random.inits, use.gpu = F)
  #check.fit(fits[[1]], mfrow=F)
  return(fits[[1]])
}

#' @title test function
#' Runs do_fits for a range of num.components and put results in a list ready for compareICs
test <- function(i1=1, i2=4, max.num.components=2, random.inits = 3, plot=T){
  fits.list = lapply(1:max.num.components, function(i)do_fits(i1, i2, num.components = i, random.inits = random.inits, plot=plot))
  return(fits.list)
}


do_fits <- function(i1,i2, num.components=2, random.inits=1, plot=T){
  #library(bulkem, lib.loc='/Library/Frameworks/R.framework/Versions/3.2/Resources/')
  #devtools::install("/Users/manuguerra/Documents/git/bulkem")
  #library(statmod)
  #source("/Users/manuguerra/Documents/git/bulkem/R/bulkem.R")
  #source("/Users/manuguerra/Documents/git/bulkem/R/invgaussmixEM.R")
  Ds <- 0:20
  t1 <- data.sim[which(data.sim[,2]==i1),1]
  t2 <- data.sim[which(data.sim[,2]==i2),1]
  xlist=lapply(Ds, function(d)data.gen(t1, t2, d, num.components=num.components))

  #fits <- bulkem2(datasets=xlist)
  #print('--- fits ---')
  #print(fits)
  fits <- bulkem2(datasets=xlist, num.components=num.components, random.inits=random.inits, use.gpu=FALSE, verbose=TRUE)
  if (plot){
    plot(Ds,sapply(fits,function(x)x$llik), ylab="log-likelihood", main=paste(num.components,"components"))
    lines(Ds,sapply(fits,function(x)x$llik))
  }
  return(fits)
}

test2 <- function(i1,i2, random.inits=1, plot=T){
  #library(bulkem, lib.loc='/Library/Frameworks/R.framework/Versions/3.2/Resources/')
  #devtools::install("/Users/manuguerra/Documents/git/bulkem")
  library(statmod)
  source("/Users/manuguerra/Documents/git/bulkem/R/bulkem.R")
  source("/Users/manuguerra/Documents/git/bulkem/R/invgaussmixEM.R")
  Ds <- 0:20
  t1 <- data.sim[which(data.sim[,2]==i1),1]
  t2 <- data.sim[which(data.sim[,2]==i2),1]
  xlist=lapply(Ds, function(d)data.gen2(t1, t2, d))

  #fits <- bulkem2(datasets=xlist)
  #print('--- fits ---')
  #print(fits)
  fits <- bulkem2(datasets=xlist, num.components=num.components, random.inits=random.inits, use.gpu=FALSE, verbose=TRUE)
  if (plot){
    plot(Ds,sapply(fits,function(x)x$llik))
    lines(Ds,sapply(fits,function(x)x$llik))
  }
  return(fits)
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

#' @title Returns data ready to feed bulkem2
#' Takes two series of spike times, plus a temporal distance betweek neuron 1 and neuron 2, and computes t_22 and t_12.
#' Cond1: If more than a neur1 spike occur in a given isi of neur2, only the last spike is considered.
#' Cond2: If no spikes from neur1 occur in a give isi of neur2, t_12 = 0.
data.gen <- function(t1, t2, d, num.components=2){
  t1=t1+d #shift/delay
  l2 <- length(t2)
  isi22 <- t2[2:l2]-t2[1:(l2-1)]
  ii <- findInterval(t1,t2)
  #%%%%%%%%%%%%%%%%%%%%%%
  extremes <- (ii==0 | ii == length(t2))
  t1 <- t1[!extremes]
  ii <- ii[!extremes]
  #%%%%%%%%%%%%%%%%%%%%%%
  dupl <- duplicated(ii, fromLast=T) #See description (cond1)
  t1 <- t1[!dupl]
  ii <- ii[!dupl]
  #%%%%%%%%%%%%%%%%%%%%%%
  l1 <- length(ii)
  isi12 <- rep(0, l2-1) #See description (cond2)
  isi12[ii] <- t2[ii+1] - t1
  #%%%%%%%%%%%%%%%%%%%%%%
  return(matrix(c(isi22, rep(isi12, num.components-1)), ncol=num.components))
}

#' @title Returns data ready to feed bulkem2
#' Takes two series of spike times, plus a temporal distance betweek neuron 1 and neuron 2, and computes t_22 and t_12.
#' Cond1: If more than a neur1 spike occur in a give isi of neur2, all the spikes are considered.
#' Cond2: If no spikes from neur1 occur in a give isi of neur2, t_12 = 0.
data.gen2 <- function(t1, t2, d){
  t1=t1+d #shift/delay
  l2 <- length(t2)
  isi22 <- t2[2:l2]-t2[1:(l2-1)]
  ii <- findInterval(t1,t2)
  #%%%%%%%%%%%%%%%%%%%%%%
  extremes <- (ii==0 | ii == length(t2))
  t1 <- t1[!extremes]
  ii <- ii[!extremes]
  #%%%%%%%%%%%%%%%%%%%%%%
  isi_mat <- matrix(0, nrow=nrows, ncol=ncols+1)
  isi_mat[,1] <- isi22
  #FIXME: not efficient
  for (i in unique(ii)) {
    isi_seq <- t2[i+1] - t1[which(ii==i)]
    isi_mat[i, 2:(length(isi_seq)+1)] <- isi_seq
  }
  #%%%%%%%%%%%%%%%%%%%%%%
  return(isi_mat)
}


#' @title Returns the interspike intervals of m neurons firing (ti) on a specific neuron (t0)
#' Takes 1+m series of spike times,  and computes a matrix nX(m+1), where m is the number of neuron firing on neur 0, and n is the number of isi of neur 0.
#' d is an m-vector
#' Cond1: If more than a neur1 spike occur in a given isi of neur2, only the last spike is considered.
#' Cond2: If no spikes from neur1 occur in a give isi of neur2, t_12 = 0.
data.gen.extended <- function(ti, t0, d){
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
