test <- function(){
  #library(bulkem, lib.loc='/Library/Frameworks/R.framework/Versions/3.2/Resources/')
  #devtools::install("/Users/manuguerra/Documents/git/bulkem")
  library(statmod)
  source("/Users/manuguerra/Documents/git/bulkem/R/bulkem.R")
  source("/Users/manuguerra/Documents/git/bulkem/R/invgaussmixEM.R")
  i1 <- 1
  i2 <- 4
  Ds <- 1:10
  t1 <- data.sim[which(data.sim[,2]==i1),1]
  t2 <- data.sim[which(data.sim[,2]==i2),1]
  xlist=lapply(Ds, function(d)data.gen(t1, t2, d))

  fits <- bulkem2(datasets=xlist)
  print('--- fits ---')
  print(fits)
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
#' Cond1: If more than a neur1 spike occur in a give isi of neur2, only the last spike is considered.
#' Cond2: If no spikes from neur1 occur in a give isi of neur2, t_12 = 0.
data.gen = function(t1, t2, d){
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
  return(matrix(c(isi22, isi12), ncol=2))
}
