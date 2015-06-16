setup <- function(home.dir='~/Documents/git/synapticDist/',code.dir='~/Documents/git/synapticDist/R/', data.dir='~/Documents/git/synapticDist/Data/', res.dir='~/Documents/git/synapticDist/Results/'){
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

############################################################################
############################################################################
correlation <- function(){
	cor.res.dir = paste(res.dir,'cor/',sep='')
	dir.create(cor.res.dir,recursive=T)
	for (irec in recs){ #missing rec=22 and j from 162 to 166 - stop and restart from irec=47 - skip 47
		rat=rr[[1]][irec]
		rec=rr[[2]][irec]
		x1 = read_data(rat=rat,rec=rec,dir0=data.dir)
		print(paste(irec,rat,rec))
		for (j in 1:ncol(x1)){
	#		print(j)
			cat(j," .. ",fill=T)
			t1 = make.lagged.matrix(x1[,j],max.lag=100)
			try(system.time(c1 <- gpuCor(x1[1000:nrow(x1),],t1[1000:nrow(t1),])))
			if (exists('c1')) try(save(c1,file=paste(cor.res.dir,"cor.",irec,".",j,".R",sep='')))
			rm(c1)
		}
		#try(system.time(c0 <- cor(x1[1000:nrow(x1),],t1[1000:nrow(t1),])))
	}
}

correlation.spikes <- function(max.lag=100){
		cor.res.dir = paste(res.dir,'cor/',sep='')
		dir.create(cor.res.dir,recursive=T)
		x0=read.csv(paste(data.dir,"spikes.csv",sep=''))
		t=(x0[,1]-1)*1000+x0[,2]
		neur=x0[,3]
		x1=matrix(0,nrow=max(t),ncol=max(neur))
		for (i in 1:length(t)) x1[t[i],neur[i]] = 1
		####
		for (j in 1:ncol(x1)){
			#print(j)
			cat(j," .. ",fill=T)
			t1 = make.lagged.matrix(x1[,j],max.lag=max.lag)
			#try(system.time(c1 <- gpuCor(x1[1000:nrow(x1),],t1[1000:nrow(t1),])))
			try(system.time(c1 <- cor(x1[1000:nrow(x1),],t1[1000:nrow(t1),])))
			if (exists('c1')) {
				try(save(c1,file=paste(cor.res.dir,"cor.spikes.",j,".R",sep='')))
				rm(c1)
			}
		}
}


make.lagged.matrix <- function(x,max.lag=50){
	n=length(x)
	y=matrix(nrow=n,ncol=max.lag)
	for (i in 1:max.lag) y[,i]=c(rep(NA,i),x[1:(n-i)])
    return(y)
}



############################################################################
############################################################################
distances <- function(){
	cor.res.dir = paste(res.dir,'cor/',sep='')
	xdist=list()
	for (irec in recs){
		rat=rr[[1]][irec]
		rec=rr[[2]][irec]
		print(paste(irec, rat, rec))
		n=sum(read_cor_tal(rat,rec,data.dir))
		dd=matrix(nrow=n,ncol=n)
		for (j1 in 1:n){
			try(load(paste(cor.res.dir,"cor.",irec,".",j1,".R",sep='')))
			if (exists('c1')) {
				for (j2 in 1:n){
					#print(paste(j1,j2))
					try(dd[j1,j2] <- which.max(abs(c1$coefficients[j2,])))
				}
			}
		}
		xdist[[irec]]=dd
	}
	return(xdist)
}

distances.spikes <- function(n=200){
	cor.res.dir = paste(res.dir,'cor/',sep='')
	xdist=list()
	dd=matrix(nrow=n,ncol=n)
	for (j1 in 1:n){
		try(load(paste(cor.res.dir,"cor.spikes.",j1,".R",sep='')))
		if (exists('c1')) {
			for (j2 in 1:n){
				#print(paste(j1,j2))
				try(dd[j1,j2] <- which.max(abs(c1[j2,])))
			}
		}
	}
	xdist[[1]]=dd
	return(xdist)
}

read_cor_tal <- function(rat,rec,dir0){	
	dir=paste(rat,rec,sep='/')
	file=paste(dir0,dir,"/infoCT.txt",sep='')
	xtmp=scan(file,quiet = TRUE)
	return(xtmp)
}

############################################################################
############################################################################


glm.lasso.spikes <- function(xdist){
	glm.res.dir = paste(res.dir,'glm/',sep='')
	dir.create(glm.res.dir,recursive=T)
	####
	x0=read.csv(paste(data.dir,"spikes.csv",sep=''))
	t=(x0[,1]-1)*1000+x0[,2]
	neur=x0[,3]
	x1=matrix(0,nrow=max(t),ncol=max(neur))
	for (i in 1:length(t)) x1[t[i],neur[i]] = 1
	####
		n=ncol(x1)
		for (i in 1:n){
			cat(i,', ')
			#xstd=shift(x1,xdist[[irec]][,i])
			xstd=shift.and.sparse.and.integrate(x1,xdist[[1]][,i],ineur=i,lag=2)
			ystd=x1[1:nrow(xstd),i]
			###
			try(fit <- glmnet(xstd,ystd, family="poisson", alpha=0.99))
			try(cv.fit <- cv.glmnet(xstd,ystd, family="poisson", alpha=0.99))
			###
			#system.time(fit2 <- speedglm.wfit(ystd, xstd, family=poisson()))
			###
			#dfstd=as.data.frame(cbind(ystd,xstd))
			#fml = as.formula(paste("ystd ~", paste(paste("V",2:ncol(dfstd),sep=''),collapse='+')))
			#try(fit <- speedglm(fml, data=dfstd, family=poisson()))
			###
			#try(fit3 <- glm.fit(ystd, xstd, family=poisson()))
			###
			try(save(fit, file=paste(glm.res.dir,'fit.glmnet.rec_spikes.neur_',i,'.R',sep='')))
			try(save(cv.fit, file=paste(glm.res.dir,'cv.fit.glmnet.rec_spikes.neur_',i,'.R',sep='')))
			#try(save(fit, file=paste(glm.res.dir,'fit.speedglm.rec_',irec,'.neur_',i,'.R',sep='')))
		}		
}

glm.lasso <- function(xdist){
	glm.res.dir = paste(res.dir,'glm/',sep='')
	dir.create(glm.res.dir,recursive=T)
	for (irec in recs){
		rat=rr[[1]][irec]
		rec=rr[[2]][irec]
		x1 = read_data(rat=rat,rec=rec,data.dir)
		n=ncol(x1)
		print(paste(irec, rat, rec, ' - ',n,'neurons'))
		for (i in 1:n){
			cat(i,', ')
			#xstd=shift(x1,xdist[[irec]][,i])
			xstd=shift.and.sparse(x1,xdist[[irec]][,i])
			ystd=x1[1:nrow(xstd),i]
			###
			try(fit <- glmnet(xstd,ystd, family="poisson", alpha=0.99))
			try(cv.fit <- cv.glmnet(xstd,ystd, family="poisson", alpha=0.99))
			###
			#system.time(fit2 <- speedglm.wfit(ystd, xstd, family=poisson()))
			###
			#dfstd=as.data.frame(cbind(ystd,xstd))
			#fml = as.formula(paste("ystd ~", paste(paste("V",2:ncol(dfstd),sep=''),collapse='+')))
			#try(fit <- speedglm(fml, data=dfstd, family=poisson()))
			###
			#try(fit3 <- glm.fit(ystd, xstd, family=poisson()))
			###
			try(save(fit, file=paste(glm.res.dir,'fit.glmnet.rec_',irec,'.neur_',i,'.R',sep='')))
			try(save(cv.fit, file=paste(glm.res.dir,'cv.fit.glmnet.rec_',irec,'.neur_',i,'.R',sep='')))
			#try(save(fit, file=paste(glm.res.dir,'fit.speedglm.rec_',irec,'.neur_',i,'.R',sep='')))
		}
	}	
}



read_data <- function(rat,rec,dir0){
xtmp=read_cor_tal(rat,rec,dir0)
n=sum(xtmp)
dir=paste(rat,rec,sep='/')
 #First run: find m
m=0
for (i in 1:n){
	file=paste(dir0,dir,"/neuron",ifelse(trunc(i/100),"","0"),ifelse(trunc(i/10),"","0"),i,".txt",sep='')
	xtmp=scan(file,quiet = TRUE)
	m=max(m,xtmp)
}
mmax=m
 #Second run: create x
x=matrix(0,ncol=n,nrow=m)
for (i in 1:n){
	file=paste(dir0,dir,"/neuron",ifelse(trunc(i/100),"","0"),ifelse(trunc(i/10),"","0"),i,".txt",sep='')
	xtmp=scan(file,quiet = TRUE)
	x[xtmp,i]=1
}
return(x)
}

read_data_size <- function(rat,rec,dir0){
dir=paste(rat,rec,sep='/')
file=paste(dir0,dir,"/xall.txt",sep='')
xtmp=NULL
xtmp=try(scan(file,quiet = TRUE))
a=length(xtmp)
if (!(a>1)) a=0
return(a)
}



shift.and.sparse  <- function(x,d){ 
 #x: design matrix; d: vector of distances from i^th neuron
i=NULL
j=NULL
n=nrow(x)
for (k in 1:ncol(x)){
 #i0=max(1,(which(x[,k]!=0)-d[k]))
 i0=which(x[,k]!=0)+d[k]
	i0=i0[-c(1,which(i0>n))]
 i=c(i,i0)
 j=c(j,rep(k,length(i0)))
}
 # sort?
return(sparseMatrix(i=i,j=j))
}

shift.and.sparse.and.integrate  <- function(x,d,ineur,lag=2){   #Consider that post-synaptic neurons do not fire as soon as the voltage goes up (i.e. signals don't need to arrive exactly at the same time)
 #maybe need a different regression, on the differences of time. instead of (1, 0, 0, 1, 0)vs(1, 0, 0, 0, 1), need 3vs4....
 #x: design matrix; d: vector of distances from i^th neuron
y=x[,ineur]
len.y=length(y)

times=NULL
neurs=NULL
n=nrow(x)
for (k in 1:ncol(x)){
 #i0=max(1,(which(x[,k]!=0)-d[k]))
 i0=which(x[,k]!=0)+d[k]
	i0=i0[-c(1,which(i0>n))]
 times=c(times,i0)
 neurs=c(neurs,rep(k,length(i0)))
}
 #Integrate == move the 1 in x to the same position as in y if within lag.
t = times
for (jj in 1:length(t)){
	t0 = t[jj]
	l.limit = max((t0-lag),1)
	u.limit = min((t0+lag),len.y)
	if (y[t0] == 0 & sum(y[l.limit:u.limit])>0) {
		ifire = which(y[l.limit:u.limit]==1)[1]-(lag+1)
		times[jj]=t0+ifire
	}
}
######
 #print(paste(length(i),length(j)))
return(sparseMatrix(i=times,j=neurs))
}



shift  <- function(x,d){ 
 #x: design matrix; d: vector of distances from i^th neuron
n=nrow(x)
m=ncol(x)
x2=matrix(0,nrow=n,ncol=m)
for (k in 1:ncol(x)){
 #i0=max(1,(which(x[,k]!=0)-d[k]))
 i0=which(x[,k]!=0)+d[k]
	i0=i0[-c(1,which(i0>n))]
	x2[i0,k]=1
}
return(x2)
}



############################################################################
############################################################################


glm.lasso.res <- function(){
	glm.res.dir = paste(res.dir,'glm/',sep='')
	lasso.results=list()
	for (irec in recs){
		rat=rr[[1]][irec]
		rec=rr[[2]][irec]
		print(paste(irec, rat, rec))
		n=sum(read_cor_tal(rat,rec,data.dir))
		#pdf(paste(glm.res.dir,'plot.cv.fit_',irec,'.R',sep=''))
		Active.index.coef=list()
		for (i in 1:n){
			load(file=paste(glm.res.dir,'fit.glmnet.rec_',irec,'.neur_',i,'.R',sep=''))
			load(file=paste(glm.res.dir,'cv.fit.glmnet.rec_',irec,'.neur_',i,'.R',sep=''))
			#do something here
			#plot(cv.fit,main=i)
			#Coefficients <- coef(fit, s = cv.fit$lambda.min)
			Coefficients1 <- coef(fit, s = cv.fit$lambda.1se)
			#Active.Index <- which(Coefficients != 0)
			Active.Index1 <- which(Coefficients1 != 0)				
			## In fit (Coefficients1 and Active.Index1) there are n+1 possible indexes as the first is the intercept -> discard
			#Coefficients  <- Coefficients[2:length(Coefficients)]
			Coefficients1 <- Coefficients1[2:length(Coefficients1)]
			Active.Index1 <- Active.Index1[2:length(Active.Index1)]-1			
			#Active.Coefficients  <- Coefficients[Active.Index]
			Active.Coefficients1  <- Coefficients1[Active.Index1]
			#Active.Index
			#Active.Coefficients
			Active.index.coef[[i]]=data.frame(Active.Index1,Active.Coefficients1)
			####
			disconnected.neurs <- (1:n)[-Active.Index1]
			xdist[[irec]][i,disconnected.neurs] <<- NA
			####
		}
		#dev.off()
		lasso.results[[irec]] = Active.index.coef
	}	
	#return(lasso.results)
}

glm.lasso.spikes.res <- function(){
	glm.res.dir = paste(res.dir,'glm/',sep='')
	lasso.results=list()
		n=200
		#pdf(paste(glm.res.dir,'plot.cv.fit_',irec,'.R',sep=''))
		Active.index.coef=list()
		for (i in 1:n){
			load(file=paste(glm.res.dir,'fit.glmnet.rec_spikes.neur_',i,'.R',sep=''))
			load(file=paste(glm.res.dir,'cv.fit.glmnet.rec_spikes.neur_',i,'.R',sep=''))
			#do something here
			#plot(cv.fit,main=i)
			#Coefficients <- coef(fit, s = cv.fit$lambda.min)
			Coefficients1 <- coef(fit, s = cv.fit$lambda.1se)
			#Active.Index <- which(Coefficients != 0)
			Active.Index1 <- which(Coefficients1 != 0)		
				## In fit (Coefficients1 and Active.Index1) there are n+1 possible indexes as the first is the intercept -> discard
				#Coefficients  <- Coefficients[2:length(Coefficients)]
				Coefficients1 <- Coefficients1[2:length(Coefficients1)]
				Active.Index1 <- Active.Index1[2:length(Active.Index1)]-1
			#Active.Coefficients  <- Coefficients[Active.Index]
			Active.Coefficients1  <- Coefficients1[Active.Index1]
			#Active.Index
			#Active.Coefficients
			Active.index.coef[[i]]=data.frame(Active.Index1,Active.Coefficients1)
			####
			disconnected.neurs <- (1:n)[-Active.Index1]
			xdist[[1]][i,disconnected.neurs] <<- NA
			####
		}
		#dev.off()
		lasso.results[[1]] = Active.index.coef
	return(lasso.results)
}

############################################################################
############################################################################


search.bricks <- function(x,L){
	nneur=ncol(x)
	times.all=list()
	for (i in 1:nneur) times.all[[i]]=which(x[,i]==1)
	bricks=list()
	for (ineur in 1:nneur){
		cat(ineur,', ')
		#Neur i spike on neur ineur at these times:
		for (i in (1:nneur)[-ineur]) {
			d=xdist[[irec]][ineur,i]
			if (!is.na(d)){
				connection.times <- times.all[[ineur]][which(times.all[[ineur]] %in% (times.all[[i]]-d))]
				for (t in connection.times){
					bricks=c(bricks,list(matrix(c(i,ineur,(t-d),t),nrow=2,ncol=2,byrow=T)))
				}
			}
		}
	}
	t.dest <- sapply(bricks,function(x){return(x[2,2])})
	bricks.order <<- order(t.dest,decreasing=T)
	bricks <<- bricks
	#bricks.ordered=bricks[bricks.order]
	#return(list(bricks,bricks.order))
}

############################################################################
############################################################################


search.poly.groups.3 <- function(grp.min.len=2){
		bo=bricks.order
		groups=list()
		while (length(bo)>=grp.min.len){
			head=bo[1]
			bo=bo[-1]
			rel=find.relatives(list(head,bo))
			#print(rel[[1]])
			if (length(rel[[1]])>=grp.min.len) groups=c(groups,list(rel[[1]]))
			bo=rel[[2]]
		}
		return(groups)
}


find.relatives <- function(list.group.bo){ 
	group=list.group.bo[[1]]
	bo=list.group.bo[[2]]
	## find the time interval of group
	t.max <- max(sapply(group,function(x){return(bricks[[x]][2,2])}))+100
	t.min <- min(sapply(group,function(x){return(bricks[[x]][2,1])}))-100
	## search possible bricks from bo list (100ms far from the interval): candidates
	to.move=NULL
	for (i in 1:200){
		candidate=bricks[[bo[i]]]
		if (is.null(candidate)) break
		if (candidate[2,1]<t.min | candidate[2,2]>t.max) break
		test=compare(candidate,group)
		if (test) to.move=c(to.move,i)
	}
	if (length(to.move)>0){
		group=c(group,bo[to.move]) #copy successful candidates in group
		bo = bo[-to.move] #delete moved bricks from bo indexes
		list.group.bo <- find.relatives(list(group,bo))
	}
	return(list.group.bo)
}


eq.any.col <- function(a,b){
	return(identical(a[,1],b[,1]) | identical(a[,1],b[,2]) | identical(a[,2],b[,1]) | identical(a[,2],b[,2]))
}


compare <- function(a,group){
	for (ib in group){
		if (eq.any.col(a,bricks[[ib]])) {
			return(TRUE)
		}
	}
	return(FALSE)
}


find.relatives <- function(list.group.bo){ 
	group=list.group.bo[[1]]
	bo=list.group.bo[[2]]
	## find the time interval of group
	t.max <- max(sapply(group,function(x){return(bricks[[x]][2,2])}))+100
	t.min <- min(sapply(group,function(x){return(bricks[[x]][2,1])}))-100
	## search possible bricks from bo list (100ms far from the interval): candidates
	to.move=NULL
	for (i in 1:200){
		candidate=bricks[[bo[i]]]
		if (is.null(candidate)) break
		if (candidate[2,1]<t.min | candidate[2,2]>t.max) break
		test=compare(candidate,group)
		if (test) to.move=c(to.move,i)
	}
	if (length(to.move)>0){
		group=c(group,bo[to.move]) #copy successful candidates in group
		bo = bo[-to.move] #delete moved bricks from bo indexes
		list.group.bo <- find.relatives(list(group,bo))
	}
	return(list.group.bo)
}

############################################################################
############################################################################


simil.groups.by.similarity <- function(bricks, groups,irec=irec,grp.min.len=6,min.simil=0.2,do.plot=F){
	n <- max(sapply(bricks,function(x){return(x[1,])})) 
	group.matrix=list()
	lg=length(groups)
	long.enough.groups=NULL
	for (i in 1:lg){
		if (length(groups[[i]])>=grp.min.len) long.enough.groups=c(long.enough.groups,i)
	}
	if (length(long.enough.groups)<=1) return()	
	lg=length(long.enough.groups)
	cat(paste(lg,"long enough groups to indagate in recording",irec),"\n")
	cat("----------------------------------------------------","\n")
	simil = data.frame()
	for (i in 1:(lg-1)){
		#cat(i," ")
		for (j in (i+1):lg){
			ig=groups[[long.enough.groups[i]]]
			jg=groups[[long.enough.groups[j]]]
			len.grp.1=length(ig)
			len.grp.2=length(jg)
			simil0=0
			for (i2 in ig){
				for (j2 in jg){
					bi=bricks[[i2]]
					bj=bricks[[j2]]
					bi[2,]=bi[2,]-bi[2,1]
					bj[2,]=bj[2,]-bj[2,1]
					if (identical(bi,bj)) simil0=simil0+1
				}
			}
			simil0=simil0/(min(len.grp.1,len.grp.2))
			#print(c(i,j,simil0))
			if (simil0 > min.simil) simil = rbind(simil, c(long.enough.groups[i],long.enough.groups[j],simil0))
		}
	}
	if (nrow(simil)==0) return("No simil groups!")
	simil=simil[order(simil[,3],decreasing=T),]
	###
	cat("\n")
	cat('# By similarity:',"\n")
	cat(paste(nrow(simil),'pairs of selected groups.'),"\n")
	print(summary(simil[,3]))
	cat('##########',"\n")
	###
	if (do.plot){
 	  if (nrow(simil)>0){
		pdf(paste("simil.grp.by.simil.rec",irec,".pdf",sep=''))
		par(mfrow=c(2,1))
		cat("Plotting similarities..","\n")
		for (ii in 1:nrow(simil)){
			i=simil[ii,1]
			j=simil[ii,2]
			simil0=simil[ii,3]
			###
			#print(paste("simil:",i,j))
			i1=long.enough.groups[i]
			j1=long.enough.groups[j]
			#1: i
			ig1=groups[[i1]]
			br1=bricks[ig1]
			t1.min <- min(sapply(br1,function(x){return(x[2,1])}))
			t1.max <- max(sapply(br1,function(x){return(x[2,2])}))
			t1.int <- t1.max-t1.min
			#2: j
			ig2=groups[[j1]]
			br2=bricks[ig2]
			t2.min <- min(sapply(br2,function(x){return(x[2,1])}))
			t2.max <- max(sapply(br2,function(x){return(x[2,2])}))
			t2.int <- t2.max-t2.min	
			#both
			t.int=max(t1.int,t2.int)	
			#print(paste(br1[2,],t1.min,br1[1,]))
			n.min <- 0 
			n.max <- max(sapply(bricks,function(x){return(x[1,])})) 
			plot(x=br1[[1]][2,]-t1.min, y=br1[[1]][1,], t='b', xlim=c(0,t.int), ylim=c(n.min,n.max), pch=19, xlab="Time [ms]", ylab="Neuron #",main=paste(i1,'- Similarities=',simil0))
			if (length(br1)>0){
				for (jj in 2:length(br1)){
					lines(x=(br1[[jj]][2,]-t1.min), y=br1[[jj]][1,])
					points(x=(br1[[jj]][2,]-t1.min), y=br1[[jj]][1,],pch=19)
				}
			}
			plot(x=br2[[1]][2,]-t2.min, y=br2[[1]][1,], t='b', xlim=c(0,t.int), ylim=c(n.min,n.max), pch=19, xlab="Time [ms]", ylab="Neuron #",main=paste(j1,'- Similarities=',simil0))
			if (length(br2)>0){
				for (jj in 2:length(br2)){
					lines(x=(br2[[jj]][2,]-t2.min), y=br2[[jj]][1,])
					points(x=(br2[[jj]][2,]-t2.min), y=br2[[jj]][1,],pch=19)
				}
			}
		}
		par(mfrow=c(1,1))
		dev.off()
	  }
	}
	return(simil)
}



plot.groups <- function(groups,irec=NULL,grp.min.len=6){
	pdf(paste("groups.",irec,".pdf",sep=''))
	lg=length(groups)
	to.remove=NULL
	for (i in 1:lg){
		if (length(groups[[i]])<grp.min.len) to.remove=c(to.remove,i)
	}
	if (length(to.remove)>0) groups=groups[-to.remove]
	n.min <- 0 
	n.max <- max(sapply(bricks,function(x){return(x[1,])})) 
	lg=length(groups)
	if (lg>0){
	  for (i in 1:lg){
		gr=groups[[i]]
		br=bricks[gr]
		t.min <- min(sapply(br,function(x){return(x[2,1])}))
		t.max <- max(sapply(br,function(x){return(x[2,2])}))
		t.int <- t.max-t.min
		#n.min <- min(sapply(br,function(x){return(x[1,])}))
		#n.max <- max(sapply(br,function(x){return(x[1,])}))
		plot(x=(br[[1]][2,]-t.min), y=br[[1]][1,], t='b', xlim=c(0,t.int), ylim=c(n.min,n.max),pch=19,xlab="Time [ms]",ylab="Neuron #")
		for (j in 2:length(br)){
			lines(x=(br[[j]][2,]-t.min), y=br[[j]][1,])
			points(x=(br[[j]][2,]-t.min), y=br[[j]][1,],pch=19)
		}
		#dev.off()
	  }
	}
	dev.off()
}

###count how many times each pair of neurons is connected
plot.number.connections <- function(){
	num.connect = matrix(0,nrow=nneur,ncol=nneur)
	for (i in 1:length(bricks.ordered)){
		i1=bricks.ordered[[i]][1,1]
		i2=bricks.ordered[[i]][1,2]
		num.connect[i1,i2]=num.connect[i1,i2]+1
	} 
	#and plot..
	filled.contour(num.connect, col=colorpanel(20, "white", "grey10"), nlevels=15)
	dev.off()
	###
} 
############################################################
############################################################
prev_spikes_distr <- function(irec,max.dist=100){
	rat=rr[[1]][irec]
	rec=rr[[2]][irec]
	x1 = read_data(rat=rat,rec=rec,dir0=data.dir)
	times=list()
	Nneur=ncol(x1)
	for (j in 1:Nneur){
		times[[j]]=which(x1[,j]==1)
	}
	neurons=1 #1:Nneur
	prev.t=list()
	for (i in neurons){
		pt=matrix(nrow=length(times[[i]]),ncol=Nneur)
		for (j in 1:Nneur){
			pt[,j] = previous.delta(times[[j]],times[[i]],ceiling=max.dist)
		}
		prev.t[[i]] = pt
	}
	return(prev.t)
}

multi.plot <- function(x){
	attach(mtcars)
	n=ncol(x[[1]])
	mm=max(x[[1]],na.rm=T)
	breaks=0:mm
	h=matrix(nrow=mm,ncol=n)
	for (j in 1:n) {
		y=hist(prev.t[[1]][,j],breaks=breaks,plot=F)$counts
		h[,j]=y
	}
	ymax=max(h[,2:n])
	par(mfrow=c(2,2))
	q=ceiling(n/4)
	to.plot=(1:q)*4
	ll=length(to.plot)
	for (i in 1:ll) {
		for (j in to.plot[i]:(to.plot[i]+3)){
 #			if (j<=n) hist(prev.t[[1]][,j],breaks=breaks, main=j, ylim=c(0,mm))
			if (j<=n) plot(h[,j],t='h', main=j, ylim=c(0,50))
		}
		readline("Press enter")
	}
}


previous <- function(t1,t2){  
 #Given two ordered list of spike times t1 and t2, returns an array (same length as t2) with the times in t1 which are before the corresponding (i.e. same index) times in t2. If there aren't any or the time has already been used, NA.
 #Basically, find which spikes in t1 caused spikes in t2
	pre=NULL
	for (i in t2){
		inds=which(t1<i)
		if (length(inds)>0){
			pre=c(pre,t1[inds[length(inds)]]) #length faster than max(inds)
			t1=t1[-inds]
		} else {
			pre=c(pre,NA)
		}
	}
	return(pre)
}

previous.delta <- function(t1,t2,ceiling=1000){  
 #Given two ordered list of spike times t1 and t2, returns an array (same length as t2) with the times in t1 which are before the corresponding (i.e. same index) times in t2. If there aren't any or the time has already been used, NA.
 #Basically, find which spikes in t1 caused spikes in t2
	pre=NULL
	for (i in t2){
		inds=which(t1<=i)
		if (length(inds)>0){
			delta=i-t1[inds[length(inds)]]
			if (delta>ceiling) delta=NA
			pre=c(pre,delta) #length faster than max(inds)
			t1=t1[-inds]
		} else {
			pre=c(pre,NA)
		}
	}
	return(pre)
}

following <- function(t1,t2){  
 #Given two ordered list of spike times t1 and t2, returns an array (same length as t1) with the times in t2 which are after the corresponding (i.e. same index) times in t1. If there aren't any or the time has already been used, NA.
 #Basically, find which spikes in t2 are caused by spikes in t1
	pre=NULL
	for (i in t1){
		inds=which(t2<i)
		if (length(inds)>0) t2=t2[-inds]
		if (length(t2)>0) pre=c(pre,t2[1]) 		
	}
	return(pre)
}

find.PGs.original.idea <- function(x, D, max.dist=100, Nneur=NULL){ 
	# x is a matrix [Nspikes, 2]. 1st col: neur no. 2nd col: time
	# D is the matrix of distances
	# max.dist is the maximum distance between neurons that is cosidered sensible
	if (is.null(Nneur)) Nneur <- max(x[,1])
	Nspikes <- i <- nrow(x) #last spike
	ineur=x[i,1]
	it=x[i,2]
	min.t=it-max.dist
	x.redox=x[which(x[,2]>min.t),]
	Nspikes.redox=nrow(x.redox)
	if (Nspikes.redox>1){
		for (ii in 1:(Nspikes.redox-1)){
			if (are.communicating(x.redox[ii],x.redox[Nspikes.redox],D)) update.PG()
		}
	}
}

 ### Based on http://stackoverflow.com/questions/21499742/fast-minimum-distance-interval-between-elements-of-2-logical-vectors-take-2 #2014.02.03
 # x2 and x1 are the times of two neurons (with or without considering their distance) and are sorted.
 # the functions outputs, for each element of x2, the minimum distance with a smaller element of x1.
 # "i" is as long as x2. Each element is the index of the element in x1 which is closer to (and smaller than) the correspond element in x2

distances.backward <- function(x1, x2,lag=0) {
	if (lag!=0) x1 <- x1+lag
	i <- findInterval(x2, x1, all.inside = TRUE) 
	distance=x2-x1[i]
	return(matrix(c(i,(1:length(i)),distance),byrow=F,ncol=3))
}

### 2014.02.17
 # For each spike in x1, computes the time till the next spike in x2, possibly with a lag to take into account the distance between the two neurons
 # Info on neurons which spike at the same time is lost. We indagate causation and only consider future spikes.
distances.forward <- function(x1, x2, lag=0) {
	if (lag!=0) x1 <- x1+lag
  	i <- findInterval(x1, x2, all.inside = FALSE)+1
	distance=x2[i]-x1
	return(matrix(c((1:length(i)),i,distance),byrow=F,ncol=3))
}

### 2014.02.10
 # Takes the registration number and returns the data in format [nspikes,2] (1: times in ms; 2: neuron number)
read.data <-function(irec){
	rat=rr[[1]][irec]
	rec=rr[[2]][irec]
	x <- read_data(rat=rat,rec=rec,dir0=data.dir)
	a<-which( x==1, arr.ind=T )
	return(a)
}

### 2014.02.03
find.links <- function(a=NULL,D=NULL,irec=17){
	if (is.null(a)) a=read.data(irec)
	links = NULL
	n=max(a[,2])
	neurs=1:n
	aa=list()
	for (i in 1:n) aa[[i]]=a[which(a[,2]==i),1]
	if (is.null(D)) D=matrix(ceiling(runif(n*n,0,50)),nrow=n,ncol=n) #for testing purposes
	for (i1 in neurs){
		for (i2 in neurs[-i1]){
			if (!is.na(D[i1,i2])){
				synaptic.dist <- D[i1,i2]
				print(paste("Analysing",i1,"and",i2))
				n1 <- aa[[i1]]
				n2 <- aa[[i2]]
				dd <- distances(n1,n2) 
				dd <- matrix(dd[which(dd[,3]<=(synaptic.dist + tolerance) & dd[,3]>=(synaptic.dist + tolerance)),], ncol=3) 
				if (nrow(dd)>0){
				  for (j in 1:nrow(dd)) {
					link <- c(i1,n1[dd[j,1]],i2,n2[dd[j,2]])
					#print(new.record)
					links <- rbind(links,link)
				  }
				}
			}
		}
	}
	return(links)
}



### 2014.02.04
find.PGs <- function(links, pg.min.length=3){
	ol <<- links[order(links[,4]),]
	search.inds=1:nrow(ol)
	PGs=list()
	repeat{
		if (length(search.inds)<2) break
		pg.ind <- tail(search.inds,1) #nrow(ol)
		level<<-1
		#print(paste("Evaluating", pg.ind))
		pg.ind <- connect.links(pg.ind,exclude=pg.ind)
		#print(level)
		level<<-1
		if (length(pg.ind)>=pg.min.length) {
			PGs=c(PGs,list(pg.ind))
			print(pg.ind)
			#break
		}
		search.inds <- setdiff(search.inds, pg.ind)		
	}
	return(PGs)
}

### 2014.02.04
connect.links <- function(pg.ind,exclude){
	level <<- level+1
	# This is to protect from stack/memory problems
	#Better, before running R, run: "ulimit -s 16384" or more.
	if (level>700) return(pg.ind)
	for (ipg in pg.ind){
		before <- which(ol[,3]==ol[ipg,1] & ol[,4]==ol[ipg,2])
		before <- setdiff(before, exclude)
		after  <- which(ol[,1]==ol[ipg,3] & ol[,2]==ol[ipg,4])
		after <- setdiff(after, exclude)
		if (length(before)>0) {
			pg.ind = union(pg.ind, before)
			pg.before <- connect.links(before,exclude=pg.ind)
			pg.ind = union(pg.ind, pg.before)
		}
		if (length(after)> 0) {
			pg.ind = union(pg.ind, after)
			pg.after  <- connect.links(pg.ind,exclude=pg.ind)
			pg.ind = union(pg.ind, pg.after)
		}
	}
	return(pg.ind)
}

###useless as the delta is decided in the distance matrix D # 2014.02.03
PG.2.deltaPG <- function(pg){
	return(cbind(pg[,1],pg[,3],pg[,4]-pg[,2]))
}

### 2014.02.03
plot.links <- function(pg,win=1000){
	opg=pg[order(pg[,2]),]
	Nspikes=nrow(opg)
	Nneurs=max(opg[,1])
	for (i1 in 1:floor(Nspikes/win)){
		tmax=win*i1
		plot(opg[1,c(2,4)],opg[1,c(1,3)], ylim=c(0,Nneurs), xlim=c(tmax-win, tmax), xlab="ms", ylab="Neuron #", t='l')
		ind=which(opg[,4]<=tmax & opg[,4]>(tmax-win))
		for (i in ind) lines(opg[i,c(2,4)],opg[i,c(1,3)])
		readline("Press enter")
	}
}

### 2014.02.04
plot.PGs <- function(pg){
	Nlinks=length(pg)
	links=ol[pg,]
	tlim=range(c(links[,2],links[,4]))
	nlim=range(c(links[,1],links[,3]))
	plot(c(links[,2],links[,4]),c(links[,1],links[,3]), ylim=nlim, xlim=tlim, xlab="ms", ylab="Neuron #",t='p',pch=19)
 #	plot(links[1,c(2,4)],links[1,c(1,3)], ylim=nlim, xlim=tlim, xlab="ms", ylab="Neuron #", t='l')
	for (i in 1:Nlinks) lines(links[i,c(2,4)],links[i,c(1,3)])
}




##likelihood for an oriented 1-->2 couple of neuronal point processes #2014.02.05
loglik <- function(t1,t2,d,bump=0.1){ 
	n1=length(t1)
	n2=length(n2)
	lam0 = n2/1200
	lam1 = lam0 + bump
	t1check <- t1 + d
	fire <- sum(t1check %in% t2)
	nofire <- n-fire
	pfire=exp(-lam1*d)-exp(-lam1*(d+1))
	loglik <- fire*log(pfire) + nofire*log(1-pfire) 
	return(loglik)
}

###2014.02.05
find.D <- function(a=NULL,irec=17,bseq=NULL,dseq=NULL){
	#bseq: sequence of possible bumbs, i.e. increase in rate of Poisson process when the distance is right.
	if (is.null(bseq)) bseq=c(0e+00, 1e-04, 1e-03, 1e-02, 1e-01, 3e-01, 5e-01, 7e-01, 1e+00, 2e+00, 3e+00, 4e+00, 5e+00, 6e+00, 7e+00, 8e+00, 9e+00, 1e+01)
	if (is.null(dseq)) dseq=1:20
	if (is.null(a)) a=read.data(irec)
	n=max(a[,2])
	neurs=1:n
	aa=list()
	D=matrix(NA,nrow=n,ncol=n)
	for (i in 1:n) aa[[i]]=a[which(a[,2]==i),1]
	for (i1 in neurs){
		for (i2 in neurs[-i1]){
			print(paste("Analysing",i1,"and",i2))
			n1 <- aa[[i1]]
			n2 <- aa[[i2]]
			maxlik=maxbump=NULL
			for (d in dseq){
				llik=NULL
				for (bump in bseq) llik=c(llik,loglik(n1,n2,d,bump=bump))
				maxlik <- c(maxlik,max(llik))
				maxbump = c(maxbump,bseq[which.max(llik)])
			}
			no.ind=which(maxbump==0)
			try(D[i1,i2] <- dseq[which.max(maxlik[-no.ind])])
		}
	}
	return(D)
}




 ### correlation between point processes in the format [nspikes, 2]: 1:times, 2:neuron number #2014.02.11
t.cor <- function(a){
	n.neurs <- max(a[,2])
	t.max <- max(a[,1])
	rho <- matrix(nrow=n.neurs, ncol=n.neurs)
	times <- sapply(1:n.neurs, function(x,y)a[which(y[,2]==x),1],a)
	nsp <- sapply(times, length)
 #	aij <- sapply(times, function(t)sapply(times,function(x)sum(x %in% t)))
 #	rho <- (t.max*as.real(aij))/(sqrt(nsp*as.real((t.max-nsp))) %*% sqrt(nsp*as.real((t.max-nsp))))
	for (i in 1:n.neurs){
		aij <- sapply(times,function(x)sum(x %in% times[[i]]))
		rho[i,]=(t.max*as.real(aij))/sqrt(t.max*as.real(nsp[i])-as.real(nsp[i])^2)/sqrt(t.max*as.real(nsp)-as.real(nsp)^2)
	}
	#	for (j in 1:n.neurs){
	#		aij <- a[which(a[,2]==i | a[,2]==j),1]
	#		rho[i,j] <- (t.max*sum(duplicated(aij)))/sqrt(t.max*nsp[i]-nsp[i]^2)/sqrt(t.max*nsp[j]-nsp[j]^2)
	#		print(paste(i,j,rho[i,j]))
	#	}
	#}
	return(rho)
}

###2014.02.11
pre.distr <- function(sx,ineur,max.lag=20){
	n.neurs=max(sx[,2])
	inds=which(sx[,2]==ineur)
	times=sx[inds,1]
	counts=matrix(0,nrow=max.lag,ncol=n.neurs)
	for (it in times){
		ii=as.matrix(sx[which(sx[,1]<it & sx[,1]>=(it-max.lag) & sx[,2]!=ineur),])
		ii[,1]=it-ii[,1]
		#ii=unique(ii)
		counts[ii]=counts[ii]+1
	}
	return(counts)
}

###2014.02.11
pre.distr2 <- function(sx,ineur,max.lag=30){
	n.neurs=max(sx[,2])
	tt=sapply(1:n.neurs,function(i)sx[which(sx[,2]==i),1])
	dd=sapply(tt,function(t)distances(t,tt[[ineur]])[,3])
	#dd=dd[all(dd>0 & dd<=max.lag),]
	return(dd)
}

###2014.05.14
pre.distr2.ij <- function(sx,i,j,dist=0,max.lag=30){
	ti=sx[which(sx[,2]==i),1]
	tj=sx[which(sx[,2]==j),1]
	dd=distances.backward(ti,tj,dist)[,3]
	#dd=dd[all(dd>0 & dd<=max.lag),]
	return(dd)
}

###2014.02.17
post.distr <- function(sx,ineur,max.lag=30){
	n.neurs=max(sx[,2])
	inds=which(sx[,2]==ineur)
	times=sx[inds,1]
	counts=matrix(0,nrow=max.lag,ncol=n.neurs)
	for (it in times){
		print(paste(it,tail(times,1)))
		ii=as.matrix(sx[which(sx[,1]>it & sx[,1]<=(it+max.lag) & sx[,2]!=ineur),])
		ii[,1]=ii[,1]-it
		#ii=unique(ii)
		counts[ii]=counts[ii]+1
	}
	return(counts)
}

###2014.02.17
post.distr2 <- function(sx,ineur,max.lag=30){
	n.neurs=max(sx[,2])
	tt=sapply(1:n.neurs,function(i)sx[which(sx[,2]==i),1])
	dd=sapply(tt,function(t)distances.forward(tt[[ineur]],t)[,3])
	#dd=dd[all(dd>0 & dd<=max.lag),]
	return(dd)
}

###2014.05.13
post.distr2.ij <- function(sx,i,j,dist=0,max.lag=30){
	ti=sx[which(sx[,2]==i),1]
	tj=sx[which(sx[,2]==j),1]
	dd=distances.forward(ti,tj,dist)[,3]
	#dd=dd[all(dd>0 & dd<=max.lag),]
	return(dd)
}

###2014.02.18
loglik1 <- function(pars,x,d=0){ #shape a; scale s.
	#2 gamma distr, the second with d delay
	#pars: mean and std
	s1=pars[2]^2/pars[1]
	a1=pars[1]/s1
	#r2=pars[3]
	s2=pars[4]^2/pars[3]
	a2=pars[3]/s2
 #	d=pars[5]
	alpha=pars[5]
	if (alpha<=0 | alpha>=1 | d<0 | a2<1) return(10000000)
	f1=dgamma(x,shape=a1,scale=s1)
	f2=dgamma(x-d,shape=a2,scale=s2)
	return(-sum(log(alpha*f1+(1-alpha)*f2)))
}

###2014.10.17
loglik1.tmp <- function(pars, d=0, data){ #shape a; scale s.
	#2 gamma distr, the second with d delay
	#pars: mean and std
	x=data
	s1=pars[2]^2/pars[1]
	a1=pars[1]/s1
	#r2=pars[3]
	s2=pars[4]^2/pars[3]
	a2=pars[3]/s2
 #	d=pars[5]
	alpha=pars[5]
	if (alpha<=0 | alpha>=1 | d<0 | a2<1) return(10000000)
	f1=dgamma(x,shape=a1,scale=s1)
	f2=dgamma(x-d,shape=a2,scale=s2)
	return(-sum(log(alpha*f1+(1-alpha)*f2)))
}

###2014.11.03
## Conditional likelihood of isi12 given isi21 with IG model
loglik2 <- function(pars, d=0, data){ 
	#mean m; lambda l; dispersion=1/lambda
	#2 IG distr, the second with d delay
	#pars: mean and std
	#data is a two cols array (1st col: isi12; 2nd col:isi2)
	x=data
	m1=pars[1]
	l1=pars[1]^3/pars[2]^2
	m2=pars[3]
	l2=pars[3]^3/pars[4]^2
	alpha=pars[5]
	if (alpha<=0 | alpha>=1 | d<0) return(10000000)
	f1num=dinvgauss(x[,2],mean=m1,dispersion=1/l1)
	p1den=1-pinvgauss(x[,2]-x[,1],mean=m1,dispersion=1/l1)
	f2=dinvgauss(x[,1]-d,mean=m2,dispersion=1/l2)
	return(-sum(log(alpha*f1num/p1den+(1-alpha)*f2)))
}

###2014.11.03
## Conditional likelihood of isi12 given isi21 with IG model - distr. of isi2 known
loglik2.bis <- function(pars, pars0, shifted_data){ 
	#mean m; lambda l; dispersion=1/lambda
	#2 IG distr, the second with d delay
	#pars: mean and std
	#data is a two cols array (1st col: isi12; 2nd col:isi2)
	x=shifted_data
	m1=pars0[1]
	l1=pars0[1]^3/pars0[2]^2
	m2=pars[1]
	l2=pars[1]^3/pars[2]^2
	alpha=pars[3]
	if (alpha<0 | alpha>1 | pars[2]<=0) return(10000000)
	f1num=dinvgauss(x[,2],mean=m1,dispersion=1/l1)
	p1den=1-pinvgauss(x[,2]-x[,1],mean=m1,dispersion=1/l1)
	f2=dinvgauss(x[,1],mean=m2,dispersion=1/l2)
	pen=0
	if (pars[2]<2) pen=10/pars[2]
	loglik=-sum(log(alpha*f1num/p1den+(1-alpha)*f2))+pen
	return(loglik)
}



###2014.02.18
loglik1.bis <- function(pars1,pars0,x,d){ #shape a; scale s.
	#2 gamma distr, the second with d delay
	#pars: mean and std
	s1=pars0[2]^2/pars0[1]
	a1=pars0[1]/s1
 #	r2=pars[1]
	s2=pars1[2]^2/pars1[1]
	a2=pars1[1]/s2
 #	d=pars[5]
	alpha=pars1[3]
	#if (alpha<=0 | alpha>=1 | d<=0 | r2<=0) return(10000000)
	if (alpha<=0 | alpha>=1 | a2<1) return(10000000)
	f1=dgamma(x,shape=a1,scale=s1)
	f2=dgamma(x-d,shape=a2,scale=s2)
	return(-sum(log(alpha*f1+(1-alpha)*f2)))
}

###2014.02.19
take.off.f2 <- function(fit0,fit1,x,d){
	pars0=fit0$par
	pars1=fit1$par
	s1=pars0[2]^2/pars0[1]
	a1=pars0[1]/s1
	#	r2=pars[1]
	s2=pars1[2]^2/pars1[1]
	a2=pars1[1]/s2
	alpha=pars1[3]
	f1=alpha*dgamma(x,shape=a1,scale=s1)
	f2=(1-alpha)*dgamma(x-d,shape=a2,scale=s2)
	return(x[f1>f2])
}

###2014.02.19
loglik0 <- function(pars,x){ #shape a; scale s.
	#2 gamma distr, the second with d delay
	#pars: mean and std
	s1=pars[2]^2/pars[1]
	a1=pars[1]/s1
	logf1=dgamma(x,shape=a1,scale=s1,log=TRUE)
	return(-sum(logf1))
}

###2014.02.19
 # Not adjusting x for distance between neurons: the next spike was not really the next spike, especially with high rates of the second neuron
model.comparison.old <- function(pars,x,d,plot=T){
	aic1=NULL
	pars0=pars[1:2]
	pars1=pars
	fit0=optim(pars0,loglik0,x=x,lower=c(50,50),upper=c(2000,2000),method="L-BFGS-B")
	fit=list()
	fit=c(fit,list(fit0))
	for (id in d){
		fit1=optim(pars1,loglik1,x=x,d=id,lower=c(50,50,0.01,0.01,0.01),upper=c(2000,2000,20,20,1),method="L-BFGS-B")
		fit=c(fit,list(fit1))
		aic0=2*fit0$value+2*length(fit0$par)
		aic1=c(aic1,2*fit1$value+2*(length(fit0$par)+length(fit1$par)))
		cat("Distance = ",id,"\n")
		cat("AIC_0:",aic0,"\n")
		cat("AIC_1:",tail(aic1,1),"\n")
		cat("-----------------------------------------------------------","\n")
		if (plot){
			nn=100
			par(mfrow=c(2,1))
			hist(x[x<nn],breaks=(0:(nn+1)),main="Data")
			if (tail(aic1,1)<aic0){
				xx=rdensity(fit1$par,id)
				hist(xx[xx<nn],breaks=(0:(nn+1)),main=paste("H1 with distance d=",id))
			} else {
				pars=fit1$par
				pars[length(pars)]=1
				xx=rdensity(pars,id)
				hist(xx[xx<nn],breaks=(0:(nn+1)),main=paste("H0 (H1 with distance d=",id,")"))
			}
			par(mfrow=c(1,1))
		}
	}
	best.model=which.min(aic1)
	plot(d,aic1,main=paste("Best model with d=",d[best.model]),ylim=c(min(c(aic0,aic1)),max(c(aic0,aic1))))
	lines(d,aic1,col='red')
	lines(d[c(1,length(d))],rep(aic0,2))
	return(fit)
}

###2014.02.19
 # Not adjusting x for distance between neurons: the next spike was not really the next spike, especially with high rates of the second neuron
 # Takes the 1st gamma process constant when fitting the second and fits iteratively
model.comparison.refit.old <- function(pars,x,d, plot=T){
	aic1=NULL
	pars0=pars[1:2]
	pars1=pars[3:5]
	fit0=fit0.tmp=optim(pars0,loglik0,x=x,lower=c(50,50),upper=c(2000,2000),method="L-BFGS-B")
	fit=list()
	fit=c(fit,list(fit0))
	for (id in d){
		fit1=optim(pars1,loglik1.bis,pars0=fit0.tmp$par,x=x,d=id,lower=c(0.01,0.01,0.01),upper=c(20,20,1),method="L-BFGS-B")
		for (i in 1:10){
			xredox=take.off.f2(fit0.tmp,fit1,x,id)
			fit0.tmp=optim(fit0.tmp$par,loglik0,x=xredox,lower=c(50,50),upper=c(2000,2000),method="L-BFGS-B")
			pars1.prev=fit1$par
			fit1=optim(pars1,loglik1.bis,pars0=fit0.tmp$par,x=x,d=id,lower=c(0.01,0.01,0.01),upper=c(20,20,1),method="L-BFGS-B")
			if (sum(abs(fit1$par-pars1.prev))<1E-5) break()
		}
		fit=c(fit,list(fit1))
		aic0=2*fit0$value+2*length(fit0$par)
		aic1=c(aic1,2*fit1$value+2*(length(fit0$par)+length(fit1$par)))
		cat("Distance = ",id,"\n")
		cat("AIC_0:",aic0,"\n")
		cat("AIC_1:",tail(aic1,1),"\n")
		cat("-----------------------------------------------------------","\n")	
		if (plot){
			nn=100
			par(mfrow=c(2,1))
			hist(x[x<nn],breaks=(0:(nn+1)),main="Data")
			if (tail(aic1,1)<aic0){
				xx=rdensity(c(fit0$par,fit1$par),id)
				hist(xx[xx<nn],breaks=(0:(nn+1)),main=paste("H1 with distance d=",id))
			} else {
				pars=c(fit0$par,fit1$par)
				pars[length(pars)]=1
				xx=rdensity(pars,id)
				hist(xx[xx<nn],breaks=(0:(nn+1)),main=paste("H0 (H1 with distance d=",id,")"))
			}
			par(mfrow=c(1,1))
		}
	}
	best.model=which.min(aic1)
	plot(d,aic1,main=paste("Best model with d=",d[best.model]),ylim=c(min(c(aic0,aic1)),max(c(aic0,aic1))))
	lines(d,aic1,col='red')
	lines(d[c(1,length(d))],rep(aic0,2))
	return(fit)
}

###2014.05.13
 # distances already considered when making x
model.comparison <- function(pars,x.in,d,plot=T){
	aic1=NULL
	pars0=pars[1:2]
	pars1=pars
	#Should be x0, not x[[1]]
	fit0=optim(pars0,loglik0,x=x.in[[1]],lower=c(50,50),upper=c(2000,2000),method="L-BFGS-B")
	fit=list()
	fit=c(fit,list(fit0))
	for (id in d){
		## x.in is a list
		x=x.in[[id]]
		fit1=optim(pars1,loglik1,x=x,d=0,lower=c(50,50,0.01,0.01,0.01),upper=c(2000,2000,20,20,1),method="L-BFGS-B")
		fit=c(fit,list(fit1))
		aic0=2*fit0$value+2*length(fit0$par)
		aic1new=2*fit1$value+2*(length(fit0$par)+length(fit1$par))
		aic1=c(aic1,aic1new)
		cat("Distance = ",id,"\n")
		cat("AIC_0:",aic0,"\n")
		cat("AIC_1:",tail(aic1,1),"\n")
		cat("loglik:",fit0$value,"\n")
		cat("loglik:",fit1$value,"\n")
		cat("-----------------------------------------------------------","\n")
		if (plot){
			nn=100
			par(mfrow=c(2,1))
			hist(x[x<nn],breaks=(0:(nn+1)),main="Data")
			if (tail(aic1,1)<aic0){
				xx=rdensity(fit1$par,d=0)
				hist(xx[xx<nn],breaks=(0:(nn+1)),main=paste("H1 with distance d=",id))
			} else {
				pars=fit1$par
				pars[length(pars)]=1
				xx=rdensity(pars,d=0)
				hist(xx[xx<nn],breaks=(0:(nn+1)),main=paste("H0 (H1 with distance d=",id,")"))
			}
			par(mfrow=c(1,1))
		}
	}
	best.model=which.min(aic1)
	plot(d,aic1,main=paste("Best model with d=",d[best.model]),ylim=c(min(c(aic0,aic1)),max(c(aic0,aic1))))
	lines(d,aic1,col='red')
	lines(d[c(1,length(d))],rep(aic0,2))
	return(fit)
}

###2014.05.13
 # distances already considered when making x
 # Takes the 1st gamma process constant when fitting the second and fits iteratively
model.comparison.refit <- function(pars,x.in,d, plot=T){
	aic1=NULL
	pars0=pars[1:2]
	pars1=pars[3:5]
	fit0=fit0.tmp=optim(pars0,loglik0,x=x.in[[1]],lower=c(50,50),upper=c(2000,2000),method="L-BFGS-B")
	fit=list()
	fit=c(fit,list(fit0))
	for (id in d){
		## x.in is a list
		## we took care of the distance already, so now d=0 in loglik1.bis
		x=x.in[[id]]
		fit1=optim(pars1,loglik1.bis,pars0=fit0.tmp$par,x=x,d=0,lower=c(0.01,0.01,0.01),upper=c(20,20,1),method="L-BFGS-B")
		for (i in 1:10){
			xredox=take.off.f2(fit0.tmp,fit1,x,0)
			fit0.tmp=optim(fit0.tmp$par,loglik0,x=xredox,lower=c(50,50),upper=c(2000,2000),method="L-BFGS-B")
			pars1.prev=fit1$par
			fit1=optim(pars1,loglik1.bis,pars0=fit0.tmp$par,x=x,d=0,lower=c(0.01,0.01,0.01),upper=c(20,20,1),method="L-BFGS-B")
			if (sum(abs(fit1$par-pars1.prev))<1E-5) break()
		}
		fit=c(fit,list(fit1))
		aic0=2*fit0$value+2*length(fit0$par)
		aic1=c(aic1,2*fit1$value+2*(length(fit0$par)+length(fit1$par)))
		cat("Distance = ",id,"\n")
		cat("AIC_0:",aic0,"\n")
		cat("AIC_1:",tail(aic1,1),"\n")
		cat("-----------------------------------------------------------","\n")	
		if (plot){
			nn=100
			par(mfrow=c(2,1))
			hist(x[x<nn],breaks=(0:(nn+1)),main="Data")
			if (tail(aic1,1)<aic0){
				xx=rdensity(c(fit0$par,fit1$par),d=0)
				hist(xx[xx<nn],breaks=(0:(nn+1)),main=paste("H1 with distance d=",id))
			} else {
				pars=c(fit0$par,fit1$par)
				pars[length(pars)]=1
				xx=rdensity(pars,d=0)
				hist(xx[xx<nn],breaks=(0:(nn+1)),main=paste("H0 (H1 with distance d=",id,")"))
			}
			par(mfrow=c(1,1))
		}
	}
	best.model=which.min(aic1)
	plot(d,aic1,main=paste("Best model with d=",d[best.model]),ylim=c(min(c(aic0,aic1)),max(c(aic0,aic1))))
	lines(d,aic1,col='red')
	lines(d[c(1,length(d))],rep(aic0,2))
	return(fit)
}

#2014.06.24
# distances already considered when making x
model.comparison.over.num.obs <- function(pars,x0,x.in,d,plot=T){
	aic1=NULL
	pars0=pars[1:2]
	pars1=pars
	#Should be x0, not x[[1]]
	lx=length(x0)
	fit0=optim(pars0,loglik0,x=x0,lower=c(50,50),upper=c(2000,2000),method="L-BFGS-B")
	aic0=2/lx*fit0$value+2*length(fit0$par) 
	fit=list()
	fit=c(fit,list(fit0))
	for (id in d){
		## x.in is a list
		x=x.in[[id]]
		lx=length(x) ###
		fit1=optim(pars1,loglik1,x=x,d=0,lower=c(50,50,0.01,0.01,0.01),upper=c(2000,2000,20,20,1),method="L-BFGS-B")
		fit=c(fit,list(fit1))
		#aic0=2*fit0$value+2*length(fit0$par)
		###
		#aic1new=2*fit1$value+2*(length(fit0$par)+length(fit1$par))
		aic1new=2/lx*fit1$value+2*length(fit1$par) ###
		aic1=c(aic1,aic1new)
		cat("Distance = ",id,"\n")
		cat("AIC_0:",aic0,"\n")
		cat("AIC_1:",tail(aic1,1),"\n")
		cat("-----------------------------------------------------------","\n")
		if (plot){
			nn=100
			par(mfrow=c(2,1))
			hist(x[x<nn],breaks=(0:(nn+1)),main="Data")
			if (tail(aic1,1)<aic0){
				xx=rdensity(fit1$par,d=0)
				hist(xx[xx<nn],breaks=(0:(nn+1)),main=paste("H1 with distance d=",id))
			} else {
				pars=fit1$par
				pars[length(pars)]=1
				xx=rdensity(pars,d=0)
				hist(xx[xx<nn],breaks=(0:(nn+1)),main=paste("H0 (H1 with distance d=",id,")"))
			}
			par(mfrow=c(1,1))
		}
	}
	best.model=which.min(aic1)
	plot(d,aic1,main=paste("Best model with d=",d[best.model]),ylim=c(min(c(aic0,aic1)),max(c(aic0,aic1))))
	lines(d,aic1,col='red')
	lines(d[c(1,length(d))],rep(aic0,2))
	return(fit)
}

###2014.06.24
 # distances already considered when making x
 # Takes the 1st gamma process constant when fitting the second and fits iteratively
model.comparison.refit.over.num.obs <- function(pars,x0,x.in,d, plot=T){
	aic1=NULL
	pars0=pars[1:2]
	pars1=pars[3:5]
	lx=length(x0)
	fit0=fit0.tmp=optim(pars0,loglik0,x=x0,lower=c(50,50),upper=c(2000,2000),method="L-BFGS-B")
	aic0=2/lx*fit0$value+2*length(fit0$par)
	fit=list()
	fit=c(fit,list(fit0))
	for (id in d){
		## x.in is a list
		## we took care of the distance already, so now d=0 in loglik1.bis
		x=x.in[[id]]
		lx=length(x)
		fit1=optim(pars1,loglik1.bis,pars0=fit0.tmp$par,x=x,d=0,lower=c(0.01,0.01,0.01),upper=c(20,20,1),method="L-BFGS-B")
		for (i in 1:10){
			xredox=take.off.f2(fit0.tmp,fit1,x,0)
			fit0.tmp=optim(fit0.tmp$par,loglik0,x=xredox,lower=c(50,50),upper=c(2000,2000),method="L-BFGS-B")
			pars1.prev=fit1$par
			fit1=optim(pars1,loglik1.bis,pars0=fit0.tmp$par,x=x,d=0,lower=c(0.01,0.01,0.01),upper=c(20,20,1),method="L-BFGS-B")
			if (sum(abs(fit1$par-pars1.prev))<1E-5) break()
		}
		fit=c(fit,list(fit1))
		#aic0=2*fit0$value+2*length(fit0$par)
		#aic1=c(aic1,2*fit1$value+2*(length(fit0$par)+length(fit1$par)))
		aic1new=2/lx*fit1$value+2*length(fit1$par) ###
		aic1=c(aic1,aic1new)
		cat("Distance = ",id,"\n")
		cat("AIC_0:",aic0,"\n")
		cat("AIC_1:",tail(aic1,1),"\n")
		cat("-----------------------------------------------------------","\n")	
		if (plot){
			nn=100
			par(mfrow=c(2,1))
			hist(x[x<nn],breaks=(0:(nn+1)),main="Data")
			if (tail(aic1,1)<aic0){
				xx=rdensity(c(fit0$par,fit1$par),d=0)
				hist(xx[xx<nn],breaks=(0:(nn+1)),main=paste("H1 with distance d=",id))
			} else {
				pars=c(fit0$par,fit1$par)
				pars[length(pars)]=1
				xx=rdensity(pars,d=0)
				hist(xx[xx<nn],breaks=(0:(nn+1)),main=paste("H0 (H1 with distance d=",id,")"))
			}
			par(mfrow=c(1,1))
		}
	}
	best.model=which.min(aic1)
	plot(d,aic1,main=paste("Best model with d=",d[best.model]),ylim=c(min(c(aic0,aic1)),max(c(aic0,aic1))))
	lines(d,aic1,col='red')
	lines(d[c(1,length(d))],rep(aic0,2))
	return(fit)
}


## EM fit of mixture of distributions

### EM fit of mixture of gammas with only one distribution lagged - d given #2014.02.20

my.gammamixEM <- function (x, d=0, lambda = NULL, alpha = NULL, beta = NULL, k = 2, epsilon = 1e-08, maxit = 1000, maxrestarts = 20, verb = FALSE,method="Nelder-Mead") 
{
    x <- as.vector(x)
    tmp <- gammamix.init(x = x, lambda = lambda, alpha = alpha, beta = beta, k = k)
    lambda <- tmp$lambda
    alpha <- tmp$alpha
    beta <- tmp$beta
    theta <- c(alpha, beta)
    k <- tmp$k
    iter <- 0
    mr <- 0
    diff <- epsilon + 1
    n <- length(x)
    dens <- NULL
    dens <- function(lambda, theta, k) {
		k=2 #works only with k=2
        temp <- NULL
        alpha = theta[1:k]
        beta = theta[(k + 1):(2 * k)]
        #for (j in 1:k) { #works only with k=2
		#	if (j==1) {
				j=1
            	temp = cbind(temp, dgamma(x, shape = alpha[j], scale = beta[j]))
		#	} else if (j==2) {
				j=2
				temp = cbind(temp, dgamma(x-d, shape = alpha[j], scale = beta[j]))
		#	}
        #}
        temp = t(lambda * t(temp))
        temp
    }
    old.obs.ll <- sum(log(apply(dens(lambda, theta, k), 1, sum)))
    ll <- old.obs.ll
    par.chain <- c(lambda,theta)
    gamma.ll <- function(theta, z, lambda, k) {
    	inds=which(z[,1]!=0 & z[,2]!=0)
    	return(-sum(z[inds,] * log(dens(lambda, theta, k)[inds,])))
    }
	while(diff>epsilon && iter<maxit){
        dens1 = dens(lambda, theta, k)
        z = dens1/apply(dens1, 1, sum)
        lambda.hat = apply(z, 2, mean)
        #out = try(suppressWarnings(nlm(gamma.ll, p = theta, lambda = lambda.hat, k = k, z = z)), silent = TRUE)
        out = try(suppressWarnings(optim(theta, gamma.ll, lambda = lambda.hat, k = k, z = z, method=method)), silent = TRUE)
        #print(str(out))
        if (class(out) == "try-error") {
#            cat("Note: Choosing new starting values.", "\n")
#            if (mr == maxrestarts) stop(paste("Try different number of components?", "\n"))
#            mr <- mr + 1
#            tmp <- gammamix.init(x = x, k = k)
#            lambda <- tmp$lambda
#            alpha <- tmp$alpha
#            beta <- tmp$beta
#            theta <- c(alpha, beta)
#            k <- tmp$k
#            iter <- 0
#            diff <- epsilon + 1
#            old.obs.ll <- sum(log(apply(dens(lambda, theta, k), 1, sum)))
#            ll <- old.obs.ll
			theta.hat=rep(0,4)
			alpha.hat = theta.hat[1:k]
            beta.hat = theta.hat[(k + 1):(2 * k)]
            new.obs.ll <- old.obs.ll
            diff <- 0
            ll <- c(ll, old.obs.ll)
            lambda = lambda.hat
            theta = theta.hat
            par.chain=rbind(par.chain,c(lambda,theta))
            alpha = alpha.hat
            beta = beta.hat
            iter = iter + 1
        }
        else {
            #theta.hat = out$estimate
            theta.hat=out$par
            alpha.hat = theta.hat[1:k]
            beta.hat = theta.hat[(k + 1):(2 * k)]
            new.obs.ll <- sum(log(apply(dens(lambda.hat, theta.hat, k), 1, sum)))
            diff <- new.obs.ll - old.obs.ll
            old.obs.ll <- new.obs.ll
            ll <- c(ll, old.obs.ll)
            lambda = lambda.hat
            theta = theta.hat
            par.chain=rbind(par.chain,c(lambda,theta))
            alpha = alpha.hat
            beta = beta.hat
            iter = iter + 1
            if (verb) {
                cat("iteration =", iter, " log-lik diff =", diff, " log-lik =", new.obs.ll, "\n")
            }
        }
    }
    if (iter == maxit) {
        cat("WARNING! NOT CONVERGENT!", "\n")
    }
    cat("number of iterations=", iter, "\n")
    theta = rbind(alpha, beta)
    rownames(theta) = c("alpha", "beta")
    colnames(theta) = c(paste("comp", ".", 1:k, sep = ""))
    colnames(par.chain) = c("lambda1","lambda2","alpha1","alpha2","beta1","beta2")
    par.chain=par.chain[,c(1,3,5,2,4,6)]
    a = list(x = x, lambda = lambda, gamma.pars = theta, loglik = new.obs.ll, posterior = z, all.loglik = ll, n.iters=iter, par.chain=par.chain, z=z, ft = "gammamixEM")
    class(a) = "mixEM"
    a
}

my.log.gammamixEM <- function (x, d=0, lambda = NULL, alpha = NULL, beta = NULL, k = 2, epsilon = 1e-08, maxit = 1000, maxrestarts = 20, verb = FALSE) 
{
    x <- as.vector(x)
    tmp <- gammamix.init(x = x, lambda = lambda, alpha = alpha, beta = beta, k = k)
    lambda <- tmp$lambda
    alpha <- tmp$alpha
    beta <- tmp$beta
    theta <- c(alpha, beta)
    k <- tmp$k
    iter <- 0
    mr <- 0
    diff <- epsilon + 1
    n <- length(x)
    dens <- NULL
    dens <- function(lambda, log.theta, k) {
    	theta=exp(log.theta)
		k=2 #works only with k=2
        temp <- NULL
        alpha = theta[1:k]
        beta = theta[(k + 1):(2 * k)]
        #for (j in 1:k) { #works only with k=2
		#	if (j==1) {
				j=1
            	temp = cbind(temp, dgamma(x, shape = alpha[j], scale = beta[j]))
		#	} else if (j==2) {
				j=2
				temp = cbind(temp, dgamma(x-d, shape = alpha[j], scale = beta[j]))
		#	}
        #}
        temp = t(lambda * t(temp))
        temp
    }
    old.obs.ll <- sum(log(apply(dens(lambda, log(theta), k), 1, sum)))
    ll <- old.obs.ll
    par.chain <- c(lambda,theta)
    log.gamma.ll <- function(log.theta, z, lambda, k) -sum(z * log(dens(lambda, log.theta, k)))
	while(diff>epsilon && iter<maxit){
        dens1 = dens(lambda, log(theta), k)
        z = dens1/apply(dens1, 1, sum)
        lambda.hat = apply(z, 2, mean)
        log.theta=log(theta)
        out = try(suppressWarnings(nlm(log.gamma.ll, p = log.theta, lambda = lambda.hat, k = k, z = z)), silent = TRUE)
        if (class(out) == "try-error") {
            cat("Note: Choosing new starting values.", "\n")
            if (mr == maxrestarts) stop(paste("Try different number of components?", "\n"))
            mr <- mr + 1
            tmp <- gammamix.init(x = x, k = k)
            lambda <- tmp$lambda
            alpha <- tmp$alpha
            beta <- tmp$beta
            theta <- c(alpha, beta)
            k <- tmp$k
            iter <- 0
            diff <- epsilon + 1
            old.obs.ll <- sum(log(apply(dens(lambda, log(theta), k), 1, sum)))
            ll <- old.obs.ll
        }
        else {
            theta.hat = exp(out$estimate)
            alpha.hat = theta.hat[1:k]
            beta.hat = theta.hat[(k + 1):(2 * k)]
            new.obs.ll <- sum(log(apply(dens(lambda.hat, log(theta.hat), k), 1, sum)))
            diff <- new.obs.ll - old.obs.ll
            old.obs.ll <- new.obs.ll
            ll <- c(ll, old.obs.ll)
            lambda = lambda.hat
            theta = theta.hat
            par.chain=rbind(par.chain,c(lambda,theta))
            alpha = alpha.hat
            beta = beta.hat
            iter = iter + 1
            if (verb) {
                cat("iteration =", iter, " log-lik diff =", diff, " log-lik =", new.obs.ll, "\n")
            }
        }
    }
    if (iter == maxit) {
        cat("WARNING! NOT CONVERGENT!", "\n")
    }
    cat("number of iterations=", iter, "\n")
    theta = rbind(alpha, beta)
    rownames(theta) = c("alpha", "beta")
    colnames(theta) = c(paste("comp", ".", 1:k, sep = ""))
    colnames(par.chain) = c("lambda1","lambda2","alpha1","alpha2","beta1","beta2")
    par.chain=par.chain[,c(1,3,5,2,4,6)]
    a = list(x = x, lambda = lambda, gamma.pars = theta, loglik = new.obs.ll, posterior = z, all.loglik = ll, n.iters=iter, par.chain=par.chain, ft = "gammamixEM")
    class(a) = "mixEM"
    a
}

### EM fit of mixture of gammas with d given #2014.06.24

my.gammamixEM.dist <- function (x, d, lambda = NULL, alpha = NULL, beta = NULL, k = 2, epsilon = 1e-08, maxit = 1000, maxrestarts = 20, verb = FALSE) 
{
    x <- as.vector(x)
    i.mix = which(x>max(d))
    i.no_mix = which(x<=max(d))
	d=d
    tmp <- gammamix.init(x = x, lambda = lambda, alpha = alpha, beta = beta, k = k)
    lambda <- tmp$lambda
    alpha <- tmp$alpha
    beta <- tmp$beta
    theta <- c(alpha, beta)
    k <- tmp$k
    iter <- 0
    mr <- 0
    diff <- epsilon + 1
    n <- length(x)
    dens <- NULL
    dens <- function(lambda, theta, k) {
		k=2 #works only with k=2
        temp <- NULL
        alpha = theta[1:k]
        beta = theta[(k + 1):(2 * k)]
        #for (j in 1:k) { #works only with k=2
		#	if (j==1) {
				j=1
            	temp = cbind(temp, dgamma(x-d[j], shape = alpha[j], scale = beta[j]))
		#	} else if (j==2) {
				j=2
				temp = cbind(temp, dgamma(x-d[j], shape = alpha[j], scale = beta[j]))
		#	}
        #}
        lambda.mix=cbind(rep(lambda[1],length(x)),rep(lambda[2],length(x)))
        lambda.mix[i.no_mix,]=cbind(rep(1,length(i.no_mix)),rep(0,length(i.no_mix)))
        temp=lambda.mix*temp
        #temp = t(lambda * t(temp))
        temp
    }
    old.obs.ll <- sum(log(apply(dens(lambda, theta, k), 1, sum)))
    ll <- old.obs.ll
    gamma.ll <- function(theta, z, lambda, k) -sum(z * log(dens(lambda, theta, k)))
    check=F
    check=try((diff>epsilon && iter<maxit),silent=T)
    count=0
	while(check){
		count=count+1
		print(count)
        dens1 = dens(lambda, theta, k)
        z = dens1/apply(dens1, 1, sum)
        #z[i.no_mix,]=cbind(rep(1,length(i.no_mix)),rep(0,length(i.no_mix)))
        lambda.hat = apply(z[i.mix,], 2, mean)
        out = try(suppressWarnings(nlm(gamma.ll, p = theta, lambda = lambda.hat, k = k, z = z)), silent = TRUE)
        if (class(out) == "try-error") {
            cat("Note: Choosing new starting values.", "\n")
            if (mr == maxrestarts) stop(paste("Try different number of components?", "\n"))
            mr <- mr + 1
            tmp <- gammamix.init(x = x, k = k)
            lambda <- tmp$lambda
            alpha <- tmp$alpha
            beta <- tmp$beta
            theta <- c(alpha, beta)
            k <- tmp$k
            iter <- 0
            diff <- epsilon + 1
            old.obs.ll <- sum(log(apply(dens(lambda, theta, k), 1, sum)))
            ll <- old.obs.ll
        }
        else {
            theta.hat = out$estimate
            alpha.hat = theta.hat[1:k]
            beta.hat = theta.hat[(k + 1):(2 * k)]
            new.obs.ll <- sum(log(apply(dens(lambda.hat, theta.hat, k), 1, sum)))
            #print(dens(lambda.hat, theta.hat, k))
            #print(theta.hat)
            #print(lambda.hat)
            diff <- new.obs.ll - old.obs.ll
            old.obs.ll <- new.obs.ll
            ll <- c(ll, old.obs.ll)
            lambda = lambda.hat
            theta = theta.hat
            alpha = alpha.hat
            beta = beta.hat
            iter = iter + 1
            if (verb) {
                cat("iteration =", iter, " log-lik diff =", diff, " log-lik =", new.obs.ll, "\n")
            }
        }
        check=F
        check=try((diff>epsilon && iter<maxit),silent=T)
    }
    if (iter == maxit) {
        cat("WARNING! NOT CONVERGENT!", "\n")
    }
    cat("number of iterations=", iter, "\n")
    theta = rbind(alpha, beta)
    rownames(theta) = c("alpha", "beta")
    colnames(theta) = c(paste("comp", ".", 1:k, sep = ""))
    a = list(x = x, lambda = lambda, gamma.pars = theta, loglik = new.obs.ll, posterior = z, all.loglik = ll, ft = "gammamixEM")
    class(a) = "mixEM"
    a
}
### EM fit of mixture of gammas with d given and iterative fitting method #2014.07.15

my.gammamixEM.refit.dist <- function (x, d, lambda = NULL, alpha = NULL, beta = NULL, k = 2, epsilon = 1e-08, maxit = 1000, maxrestarts = 20, verb = FALSE, method="Nelder-Mead") 
{
    x <- as.vector(x)
    i.mix = which(x>max(d))
    i.no_mix = which(x<=max(d))
	d=d
    tmp <- gammamix.init(x = x, lambda = lambda, alpha = alpha, beta = beta, k = k)
    lambda <- tmp$lambda
    alpha <- tmp$alpha
    beta <- tmp$beta
    theta <- c(alpha, beta)
    k <- tmp$k
    iter <- 0
    mr <- 0
    diff <- epsilon + 1
    n <- length(x)
    dens <- NULL
    dens <- function(theta, lambda, k) {
		k=2 #works only with k=2
        temp <- NULL
        alpha = theta[1:2]
        beta = theta[3:4]
				j=1
            	temp = cbind(temp, dgamma(x-d[j], shape = alpha[j], scale = beta[j]))
				j=2
				temp = cbind(temp, dgamma(x-d[j], shape = alpha[j], scale = beta[j]))
        #lambda.mix=cbind(rep(lambda[1],length(x)),rep(lambda[2],length(x)))
        #lambda.mix[i.no_mix,]=cbind(rep(1,length(i.no_mix)),rep(0,length(i.no_mix)))
        #temp=lambda.mix*temp
        temp = t(lambda * t(temp))
        temp
    }
    dens1 <- function(theta1, x1) {
        alpha = theta1[1]
        beta = theta1[2]
		j=1
        temp = dgamma(x1-d[j], shape = alpha, scale = beta)
        #temp = lambda1 * temp
        temp = -sum(log(temp))
        temp
    }
    dens2 <- function(theta2, x2) {
        temp <- NULL
        alpha = theta2[1]
        beta = theta2[2]
				j=2
            	temp = cbind(temp, dgamma(x2-d[j], shape = alpha, scale = beta))
        #temp = lambda2 * temp
        temp = -sum(log(temp))
        temp
    }
    old.obs.ll <- sum(log(apply(dens(theta, lambda, k), 1, sum)))
    ll <- old.obs.ll
    gamma.ll <- function(theta, z, lambda, k) {
    	ind0=which(z[,1]!=0 & z[,2]!=0)
    	ind1=which(z[,1]==0 & z[,2]!=0)
    	ind2=which(z[,1]!=0 & z[,2]==0)
    	-sum(z[ind0,] * log(dens(theta, lambda, k)[ind0,])) -sum(z[ind1,]*log(dens(theta, lambda, k)[ind1,2])) -sum(z[ind2,]*log(dens(theta, lambda, k)[ind2,1]))
    }
    gamma.ll1 <- function(theta1, z1, x1, lambda1) {
    	if (any(theta1<=0)) {
    		ll1=10E6
    	} else {
	    	de=dens1(theta1, lambda1, x1)
	    	ind=which(z1!=0 & de!=0)
	    	if (length(ind)==0) {
	    		ll1=0
	    	} else {
	    		ll1=-sum(z1[ind] * log(de[ind]))
	    	}
	    }
    	return(ll1)
    }
    gamma.ll2 <- function(theta2, z2, x2, lambda2) {
    	if (any(theta2<=0)) {
    		ll2=10E6
    	} else {
	    	de=dens2(theta2, lambda2, x2)
	    	ind=which(z2!=0 & de!=0)
	    	if (length(ind)==0) {
	    		ll2=0
	    	} else {
	    		ll2=-sum(z2[ind] * log(de[ind]))
	    	}
	    }
    	return(ll2)
    }
	#tmpfun <- function(t,l){-sum(log(dens1(l,t,x)))}
	#out=try(suppressWarnings(nlm(tmpfun,p=theta[c(1,3)],0.5,x)), silent = TRUE)
	#print(str(out))
	#theta[c(1,3)]=out$estimate
	while(diff>epsilon && iter<maxit){
		#E-step
        dens0 = dens(theta, lambda, k)
        z <<- dens0/apply(dens0, 1, sum)
        #theta.mat=matrix(0,nrow=100,ncol=4)
        theta.mat=NULL
        lik00=NULL
      #for (i in 1:100){
      	#rnd.chk=(z[,1]<z[,2])
        rnd.chk=(runif(length(x))<z[,1]) #rnd version
        sub1.x.ind=which(rnd.chk)
        sub2.x.ind=which(!rnd.chk)
        x1=x[sub1.x.ind]
        x2=x[sub2.x.ind]
        z1=z[sub1.x.ind,1]
        z2=z[sub2.x.ind,2]
        #M-step 1
        lambda.hat <<- apply(z[i.mix,], 2, mean)
        lambda1 = lambda.hat[1]
        lambda2 = lambda.hat[2]
      #  i=0
      ### Fitting ###  
      #while (if (i>5) sum(abs((theta.mat[i,]-theta.mat[i-1,])/theta.mat[i,]))>0.1 else TRUE){
      #	i=i+1
#      for (i in 1:100){
        #out1 = nlm(gamma.ll1, p = theta[c(1,3)], x1=x1,lambda1 = lambda2,  z1 = z[,1])
        #theta1 = out1$estimate
        #out1 = optim(theta[c(1,3)], gamma.ll1, z1 = z1, x1=x1, lambda1 = lambda1)
        out1 = optim(theta[c(1,3)], dens1, x1=x1)
        theta1 = out1$par
        #----------
        #out2 = nlm(gamma.ll2, p = theta[c(2,4)], x2=x2,lambda2 = lambda2,  z2 = z[,2])
        #theta2 = out2$estimate
        #out2 = optim(theta[c(2,4)], gamma.ll2, z2 = z2, x2=x2, lambda2 = lambda2)
        out2 = optim(theta[c(2,4)], dens2, x2=x2)
        theta2= out2$par
        #----------
        #theta.hat=rep(0,4)
        #theta.hat[c(1,3)]=theta1
        #theta.hat[c(2,4)]=theta2
        theta = c(theta1[1],theta2[1],theta1[2],theta2[2])
        #out = try(suppressWarnings(nlm(gamma.ll, p = theta, lambda = lambda.hat, k = k, z = z)), silent = TRUE)
        #out = try(suppressWarnings(optim(theta, gamma.ll, lambda = lambda.hat, k = k, z = z, method=method)), silent = TRUE)
        #theta.mat[i,]=theta.last
       # theta.mat=rbind(theta.mat,theta)
        #ll00=-sum(log(dens(theta, lambda.hat, k)))
       # lik00=c(lik00, ll00)
       #}
        #theta.hat=apply(theta.mat,2,mean)
       # print(summary(dens(theta.hat, lambda.hat, k)))
       # print(-sum(log(dens(theta, lambda.hat, k))))
       # print(paste(i,which.max(lik00)))
        #print(lik00)
        #theta.hat=theta.mat[which.max(lik00),]
		#out=try(suppressWarnings(nlm(gamma.ll, p = theta.hat, z=z, lambda = lambda.hat, k = k)), silent = TRUE)
		### End of Fitting ###
        if (class(out1) == "try-error" | class(out2) == "try-error") {
            #cat("Note: Choosing new starting values.", "\n")
            #if (mr == maxrestarts) stop(paste("Try different number of components?", "\n"))
            #mr <- mr + 1
            #tmp <- gammamix.init(x = x, k = k)
            #lambda <- tmp$lambda
            #alpha <- tmp$alpha
            #beta <- tmp$beta
            #theta <- c(alpha, beta)
            #k <- tmp$k
            #iter <- 0
            #diff <- epsilon + 1
            #old.obs.ll <- sum(log(apply(dens(theta, lambda, k), 1, sum)))
            #ll <- old.obs.ll
            iter=maxit
        }
        else {
            #theta.hat = out$estimate
            #theta.hat = out$par
            theta.hat=theta
            alpha.hat = theta.hat[1:2]
            beta.hat = theta.hat[3:4]
            new.obs.ll <- sum(log(apply(dens(theta.hat, lambda.hat, k), 1, sum)))
            #print(new.obs.ll)
            if (!is.finite(new.obs.ll)) new.obs.ll=10E6
            diff <- new.obs.ll - old.obs.ll
            #print(paste(count,diff))
            #print(theta.hat)
            old.obs.ll <- new.obs.ll
            ll <- c(ll, old.obs.ll)
            lambda = lambda.hat
            theta <<- theta.hat
            alpha = alpha.hat
            beta = beta.hat
            iter = iter + 1
            if (verb) {
                cat("iteration =", iter, " log-lik diff =", diff, " log-lik =", new.obs.ll, "\n")
            }
        }
    }
    if (iter == maxit) {
        cat("WARNING! NOT CONVERGENT!", "\n")
    }
    cat("number of iterations=", iter, "\n")
    theta = rbind(alpha, beta)
    rownames(theta) = c("alpha", "beta")
    colnames(theta) = c(paste("comp", ".", 1:k, sep = ""))
    a = list(x = x, lambda = lambda, gamma.pars = theta, loglik = new.obs.ll, posterior = z, all.loglik = ll, n.iters= iter, z=z, ft = "gammamixEM")
    class(a) = "mixEM"
    a
}




my.optim <- function(likfun,theta,x,lambda,z, fun='optim'){
	if (fun=='nlm'){
		out = nlm(likfun, p = theta, x=x, lambda = lambda,  z = z)
    } else if (fun=='optim'){
        out = optim(theta, likfun, z = z, x=x, lambda = lambda)
    } else if (fun=='constrOptim'){
    	out = NA
	}
	return(out)
}


### Summary of EM fit of mixture of gammas with only one distribution lagged - d given #2014.02.24
summary.myEM <- function(fit,d,max.x=NULL,y.lim=1.3){
	x=fit$x
	if (is.null(max.x)) max.x=max(x)
	f12 <- function(x1,x2,fit) fit$lambda[1]*dgamma(x1, shape= fit$gamma.pars[1,1], scale= fit$gamma.pars[2,1]) + fit$lambda[2]*dgamma(x2, shape= fit$gamma.pars[1,2], scale= fit$gamma.pars[2,2])
	f1 <- function(x,fit) fit$lambda[1]*dgamma(x, shape= fit$gamma.pars[1,1], scale= fit$gamma.pars[2,1])
	f2 <- function(x,fit) fit$lambda[2]*dgamma(x, shape= fit$gamma.pars[1,2], scale= fit$gamma.pars[2,2])
	xx=seq(1,max.x)
	xxlag=xx-d
	if (0 %in% xxlag) xxlag[xxlag==0]=0.0001
	freqs12=f12(xxlag,xx,fit)
	freqs1=f1(xxlag,fit)		
	freqs2=f2(xx,fit)
	par(mfrow=c(2,1))	
	plot(xx, freqs12, ylim=c(0,y.lim*max(freqs12)),"l")
	lines(xx,freqs1,col='red')
	lines(xx,freqs2,col='red')
	data=x[x<max.x]
	hist(data, freq=F, breaks=round(length(data)/5),xlim=c(0,max.x))
	par(mfrow=c(1,1))	
}

### Summary of EM fit of mixture of gammas with d given #2014.02.24

summary.myEM.dist <- function(fit,d,max.x=NULL,y.lim=1.3){
	x=fit$x
	if (is.null(max.x)) max.x=max(x)
	f1 <- function(x,fit) fit$lambda[1]*dgamma(x, shape= fit$gamma.pars[1,1], scale= fit$gamma.pars[2,1])
	f2 <- function(x,fit) fit$lambda[2]*dgamma(x, shape= fit$gamma.pars[1,2], scale= fit$gamma.pars[2,2])
	xx=seq(1,max.x)
	x1lag=xx-d[1]
	x2lag=xx-d[2]
	if (0 %in% x1lag) x1lag[x1lag==0]=0.0001 #to solve singularity when dgamma(x=0,alpha,beta)=Inf
	if (0 %in% x2lag) x2lag[x2lag==0]=0.0001
	freqs1=f1(x1lag,fit)		
	freqs2=f2(x2lag,fit)
	freqs12=freqs1+freqs2
	par(mfrow=c(2,1))	
	plot(xx, freqs12, ylim=c(0,y.lim*max(freqs12)),"l")
	lines(xx,freqs1,col='red')
	lines(xx,freqs2,col='red')
	data=x[x<max.x]
	hist(data, freq=F, breaks=round(length(data)/5),xlim=c(0,max.x))
	par(mfrow=c(1,1))	
}

 

## Random numbers from mixture of gammas #2014.02.18
rdensity <- function(pars,d,n=10000){
		a1=pars[1]^2/pars[2]^2
		s1=pars[1]/a1
		a2=pars[3]^2/pars[4]^2
		s2=pars[3]/a2
		alpha=pars[5]
	x1=rgamma(round(n*alpha),shape=a1,scale=s1)
	x2=rgamma(round(n*(1-alpha)),shape=a2,scale=s2)+d
	return(sort(c(x1,x2)))
}

## Random numbers from mixture of gammas  #2014.02.26

my.rgamma <- function(n,lambda,mean,sd){
 n1=round(lambda[1]*n)
 n2=n-n1
 shape1=mean[1]^2/sd[1]^2
 scale1=mean[1]/shape1
 shape2=mean[2]^2/sd[2]^2
 scale2=mean[2]/shape2
 cat(shape1,scale1,shape2,scale2,"\n")
 return(sort(c(rgamma(n1,shape=shape1,scale=scale1),rgamma(n2,shape=shape2,scale=scale2))))
}


## Random numbers from lagged mixture of gammas or inverse gaussians #2014.06.23

rmixture <- function(distr='gamma',n=1000,pi=c(0.95,0.05),means=c(200,8),sds=c(1000,2),d=c(0,0), plot.lim=NULL){	
	n1=round(n*pi[1])
	n2=n-n1
	data.latent <<- 1:n1
	if (distr=='gamma'){
		scales=sds^2/means
		shapes=means^2/sds^2
		data.simple=c(d[1]+rgamma(n1,shape=shapes[1],scale=scales[1]),d[2]+rgamma(n2,shape=shapes[2],scale=scales[2]))
	} else if (distr=='invgauss'){
		lambdas=means^3/sds^2
		data.simple=c(d[1]+rinvgauss(n1,mean=means[1],shape=lambdas[1]), d[2]+rinvgauss(n2,mean=means[2],shape=lambdas[2]))
	}
	par(mfrow=c(2,1))
	x=seq(0,max(data.simple),by=0.1)
	x[x==0]=0.001
	if (is.null(plot.lim)) plot.lim=range(x)
	subdata=data.simple[data.simple <= plot.lim[2] & data.simple >= plot.lim[1]]
	nn=length(subdata)
	hist(subdata,round(nn/5),xlim=plot.lim,xlab='ISI')
	if (distr=='gamma'){
		x1=dgamma(x-d[1],shape=shapes[1],scale=scales[1])
		x2=dgamma(x-d[2],shape=shapes[2],scale=scales[2])
	} else if (distr=='invgauss'){
		x1=dinvgauss(x-d[1],mean=means[1],shape=lambdas[1])
		x2=dinvgauss(x-d[2],mean=means[2],shape=lambdas[2])
	}
	if (is.null(plot.lim)) plot.lim=range(x)
	plot(x,pi[1]*x1+pi[2]*x2,t='l',ylab='Density',xlab='ISI',xlim=plot.lim)
	par(mfrow=c(1,1))
	return(data.simple)
}

## Fano factor #2014.05.13

fano.factor <- function(x,wins, plot=T){ #x are absolute times, i.e. data[which(data[,2]==ineur2),1]
	fano=NULL
	for (w in wins){
		print(paste("Computing fano factor for a",w,"ms window"))
		br = seq(0,max(x)+w,by=w) #???
		counts = hist(x,breaks=br,plot=F)$counts
		ave=mean(counts)
		stdev=sd(counts)
		print(paste("Mean:",ave))
		print(paste("SDev:",stdev))
		fano = c(fano,(stdev^2)/ave)
	}
	if (plot) {
		plot(wins,fano,t='l',log='x')
		points(wins,fano)
	}
	return(fano)
}


## Gamma and inverse gaussian models and parameter translations

models <- function(imodel){
	ind=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
	means1=c(200,600,800,1000,1200,20,800,800,800,20,20,20,20,20)
	means2=c(1,2,3,4,5,80,2,2,2,60,60,60,60,60)
	means=cbind(means1,means2)
	sds1=c(120,500,500,500,500,5,1000,1000,1000,5,5,5,5,5)
	sds2=c(0.8,1,1,1,1,10,1,1,1,10,10,10,10,10)
	sds=cbind(sds1,sds2)
	alpha1=c(.5,.6,.7,.8,.9,.5,.5,.97,.90,.5,.9,.8,.95,.5)
	alpha2=1-alpha1
	alphas=cbind(alpha1,alpha2)
	return(list(means=means[imodel,],sds=sds[imodel,],alphas=alphas[imodel,]))
}

models.trans.pars <- function(pars,distr='gamma'){
	if (distr=='gamma'){
		means=pars[1,]*pars[2,]
		sds=sqrt(pars[1,])*pars[2,]
	} else if (distr=='invgauss'){
		means=pars[1,]
		sds=sqrt(means^3/pars[2,])
	}
	return(list(means=means,sds=sds)) #,alphas=fit$lambda))
}

models.trans.pars.inv <- function(msa,distr='gamma'){
	means=msa$means
	sds=msa$sds
	if (distr=='gamma'){
		scales=sds^2/means
		shapes=means^2/sds^2
		pars=rbind(shapes,scales)
	} else if (distr=='invgauss'){
		lambdas=means^3/sds^2
		pars=rbind(means,lambdas)
	}
	colnames(pars) = c(paste("comp", ".", 1:2, sep = ""))
	return(pars)
}

ll <- function(x, d, lambda, alpha, beta){ #for gamma distr
    i.no_mix = which(x<=max(d))
    temp <- NULL
	temp = cbind(temp, dgamma(x-d[1], shape = alpha[1], scale = beta[1]))
	temp = cbind(temp, dgamma(x-d[2], shape = alpha[2], scale = beta[2]))
    lambda.mix=cbind(rep(lambda[1],length(x)),rep(lambda[2],length(x)))
 #    lambda.mix[i.no_mix,]=cbind(rep(1,length(i.no_mix)),rep(0,length(i.no_mix)))
    temp=lambda.mix*temp
    temp=apply(temp,1,sum)
    temp=sum(log(temp))
    return(temp)
}

ll.map <- function(n=100, x, d, lambda, alpha, beta, fix=1, focus.on.peak=F){ #fix is the mixuture component to fix
	change=3-fix
	#n1=round(n/2)
	#n2=n-n1
	if (focus.on.peak){
		xpts=c(alpha[change]*0.8,alpha[change]*1.2)
		ypts=c(beta[change]*0.8,beta[change]*1.2)
	} else {
		xpts=c(min(alpha)*0.1,max(alpha)*10)
		ypts=c(min(beta)*0.01,max(beta)*1)
	}
	#
	n1=max(round(n*(alpha[change]-xpts[1])/(xpts[2]-xpts[1])),10)
	n2=n-n1
	xgrid=exp(seq(log(xpts[1]),log(xpts[2]),length.out=n))
	#
	n1=max(round(n*(beta[change]-ypts[1])/(ypts[2]-ypts[1])),10)
	n2=n-n1
	ygrid=exp(seq(log(ypts[1]),log(ypts[2]),length.out=n))
	#
	llmap=NULL
	loglik=matrix(nrow=n*n,ncol=3)
	i=0
	for (alpha0 in xgrid){
		for (beta0 in ygrid){
			i=i+1
			alpha1=alpha
			beta1=beta
			alpha1[change]=alpha0
			beta1[change]=beta0
			ll.comp=ll(x,d,lambda,alpha1,beta1)
			cat(i, "/", n*n," - ", lambda[1], alpha0, beta0, ll.comp, "\n")
			loglik[i,]=c(alpha1[change],beta1[change],ll.comp)
		}
	}
	return(loglik)
}

ll.map.cross <- function(n=100, x, d, lambda, alpha, beta, fix=1, focus.on.peak=F){ #fix is the mixuture component to fix
	change=3-fix
	#n1=round(n/2)
	#n2=n-n1
	if (focus.on.peak){
		xpts=c(alpha[change]*0.8,alpha[change]*1.2)
		ypts=c(beta[change]*0.8,beta[change]*1.2)
	} else {
		xpts=c(min(alpha)*0.5,max(alpha)*10)
		ypts=c(min(beta)*0.01,max(beta)*1)
	}
	#
	n1=max(round(n*(alpha[change]-xpts[1])/(xpts[2]-xpts[1])),10)
	n2=n-n1
	xgrid1= alpha[change]-(alpha[change]-xpts[1])*seq(alpha[change]-xpts[1],0,length.out=n1)^2/(alpha[change]-xpts[1])^2
	xgrid2= alpha[change]+(xpts[2]-alpha[change])*seq(0,xpts[2]-alpha[change],length.out=n2)^2/(xpts[2]-alpha[change])^2
	xgrid=unique(c(xgrid1,xgrid2))
	#
	n1=max(round(n*(beta[change]-ypts[1])/(ypts[2]-ypts[1])),10)
	n2=n-n1
	ygrid1= beta[change]-(beta[change]-ypts[1])*seq(beta[change]-ypts[1],0,length.out=n1)^2/(beta[change]-ypts[1])^2
	ygrid2= beta[change]+(ypts[2]-beta[change])*seq(1,ypts[2]-beta[change],length.out=n2)^2/(ypts[2]-beta[change])^2
	ygrid=unique(c(ygrid1,ygrid2))	
	#
	llmap=NULL
	loglik=matrix(nrow=n*n,ncol=3)
	i=0
	for (alpha01 in xgrid){
		for (alpha02 in xgrid){
			i=i+1
			alpha1=c(alpha01,alpha02)
			beta1=beta
			ll.comp=ll(x,d,lambda,alpha1,beta1)
			cat(i, "/", n*n," - ", ll.comp, "\n")
			loglik[i,]=c(alpha01,alpha02,ll.comp)
		}
	}
	return(loglik)
}

### Function to plot nice 3d graphs

plot.nice.3d<-function(ll){
	scale.x<-function(x){
		rr=range(x)[2]-range(x)[1]
		return(((x-min(x))/rr))
	}
	ind=which(!is.na(ll[,3]) & !is.infinite(ll[,3]))
	x=unique(log(ll[ind,1]))
	x=scale.x(x)
	z=unique(log(ll[ind,2])) ##weird parametrization for rgl: z is called y
	z=scale.x(z)
	y=ll[ind,3]
	y=scale.x(y)
	y=matrix(y,length(x),length(z),byrow=T)
	require(rgl,quietly=T)
	ylim <- range(y)
	ylen <- ylim[2] - ylim[1] + 1
	jet.colors <- 
      colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                      "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
	colorzjet <- jet.colors(100)
	open3d()
	rgl.surface(x,z,y, color=colorzjet[ findInterval(y, seq(min(y), max(y), length=100))])
	axes3d()
}

