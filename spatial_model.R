##### Bayesian Spatial model for mapping two traits independently #####

library(mvtnorm)
library(pscl)
library(matrixStats)
library(boot)
library(fields)
#library(MASS)

################ a few functions for the GP regression #################
myknn <- function(dist,m){
	i=which(dist==0)
	if(i==1) {im<-0}
      	else{
		if(m>=(i-1)) im<-1:(i-1)
	      else 	im<-order(dist[1:i-1])[1:m]
	}
	return(im)
	}

sigen<-function(i,imvec,n){
	si<-0
	if(i<n){
		for(j in (i+1):n){
			if(0<sum(unlist(imvec[j])==i)) si<-c(si,j)
			}
		if(1<length(si)) {si<-si[-1]}
		}	
	return(si)
	}

ufgen<-function(i,distmat,imvec,phi,m,n){
	loc<-c(imvec[[i]],i)
	l<-length(loc)-1
	cov<-exp(-phi*distmat[loc,loc])
	u<-c(rep(0,m),1)
	if(0<loc[1]){	
		u[1:l]<-chol2inv(chol(cov[1:l,1:l]))%*%cov[1:l,l+1]
		#if(i<m+1) u<-c(u,rep(0,m+1-i))
		u[m+1]<-max(1e-12,1-t(cov[1:l,l+1])%*%u[1:l])
		}
	return(u)
	}

cijgen<-function(j,i,ufmat,imvec){
	loc<-which(imvec[[j]]==i)
	return(ufmat[loc,j])
	}

tijgen<-function(j,i,ufmat,imvec,W){
	loc1<-imvec[[j]]
	loc2<-which(loc1==i)
	tij<-W[j]-t(W[loc1])%*%ufmat[

1:length(loc1),j]+ufmat[loc2,j]*W[i]
	return(tij)
	}

ewgen<-function(i,ufmat,imvec,W){
	ew<-W[i]
	loc<-imvec[[i]]
	if(i>1) 	ew<-ew-t(W[loc])%*%ufmat[1:length(loc),i]
	return(ew)
	}

outknn <- function(dist,m){
	order(dist)[1:m]
	}

nngpkrig=function(ulocout,ulocs,s,wsmallmat,phivec,sigsqvec,locs,locsmap,m,cart){
	pos=which(ulocs==ulocout)
	if(length(pos) > 0) wout=wsmallmat[,pos] else{
		sout=as.vector(locs[locsmap[[ulocout]][1],])
		if(cart==0) sout=latlontocart(sout)
		distvec=as.vector(rdist(t(sout),s))
		im=outknn(distvec,m)
		wt=sapply(phivec,function(phi) as.vector(solve(exp(-phi*as.matrix(dist(s[im,]))))%*%exp(-phi*distvec[im])))
		mu<-sapply(1:length(phivec), function(j) as.vector(wsmallmat[j,im]%*%wt[,j]))
		v<-sigsqvec*(1-sapply(1:length(phivec),function(j) as.vector(t(exp(-phivec[j]*distvec[im])%*%wt[,j]))))
		wout=mu+sqrt(v)*rnorm(length(phivec))
	}
	wout
}		

## function to calculate the logit Jacobian in log-scale 
logJlogit<-function(x,a,b){log(x-a)+log(b-x)}
####################################################################

#### Main function for GP regression ####
#### Arguments: a list dt containing the traits (y1,y2), predictor matrices (X1,X2) and other accessory information ####
#### N is the number of MCMC iterations, N1 is number of burn-in iterations for the MCMC #### 
#### Inits is the set of initial values of the parameters ####
model_ind_spatial=function(dt,N,N1,Nout,inits,const,m,cart,spat){
	t0=Sys.time()
	#set.seed(1)

	locin=dt$locin
	n=length(locin)
	X=as.matrix(dt$X[locin,])
	y=dt$y[locin]
	p=ncol(X)

	##setting initial values
	w=rep(0,n)
	tausq=inits$tausq
	sigsq=inits$sigsq
	phi=inits$phi
	aphi=const$aphi
	bphi=const$bphi
	dphi=bphi-aphi
	ltune=const$ltune

	betamat=matrix(0,N,p)
	wmat=matrix(0,N,n)
	tausqvec=sigsqvec=phivec=flagvec=rep(0,N)

	## precalculating some expressions to save computation time
	# v=ginv(t(X)%*%X) ## Moore-Penrose pseudo-inverse
	v=solve(t(X)%*%X+1e-5*diag(p))
	cholv=chol(v)
	mu=v%*%t(X)
	#mat=rmvnorm(N,rep(0,p),v,method="svd") # svd is the only error free method...
	mat=as.matrix(rmvnorm(N,rep(0,p),diag(p))%*%cholv)

	if(spat==1){
	## mapping the spatial effects 
	ulocs=unique(dt$locsmapinv[locin])
	nsmall=length(ulocs)
	wsmallmat=matrix(0,N,nsmall)
	wsmall=rep(0,nsmall)
	Xsmall=t(sapply(ulocs,function(x)  colMeans(dt$X[as.vector(intersect(dt$locsmap[[x]],locin)),,drop=F])))
	if(p==1) Xsmall=t(Xsmall)
	#Xsmall=t(sapply(ulocs,function(x) dt$X[dt$locsmap[[x]][1],]))
	ysmall=sapply(ulocs,function(x) mean(dt$y[intersect(dt$locsmap[[x]],locin)]))
	count=sapply(ulocs,function(x) length(intersect(dt$locsmap[[x]],locin)))
	s=t(sapply(ulocs,function(x) as.vector(dt$locs[dt$locsmap[[x]][1],])))
	if(cart==0) s=t(apply(s,1,latlontocart))
	#return(s)
	distmat=as.matrix(dist(s))  ## 6400 miles is the radius of the earth but we are not multiplying
	row.names(distmat)=colnames(distmat)=NULL

	imvec<- apply(distmat,1,myknn,m)
	sivec<-sapply(1:nsmall,sigen,imvec,nsmall)
	ufmat<-sapply(1:nsmall,ufgen,distmat,imvec,phi,m,nsmall)
	
	#wsmallmat=as.matrix(rmvnorm(N,rep(0,nsmall),covwsmall))
	invmap=sapply(locin,function(x) which(ulocs==dt$locsmapinv[x]))

	}

	if(spat==1) dev.new()
	## MCMC (Gibbs) sampler 	
	for(i in 1:N){
		## updating the coefficient vectors beta1 and beta2
		beta=as.vector(mu%*%(y-w)+sqrt(tausq)*mat[i,])
		betamat[i,]=beta
		res1=y-as.vector(X%*%beta)

		## updating the variance components tausq1 and tausq2
		res=res1-w
		tausq=rigamma(1,2+n/2,0.1+0.5*t(res)%*%res)
		tausqvec[i]=tausq
		
		if(spat==1){

		## updating the spatial component w
		for(k in 1:nsmall){
			nk=count[k]
			a<-(nk/tausq)+1/(ufmat[m+1,k]*sigsq)
			loc<-imvec[[k]]
			b<-nk*(ysmall[k]-Xsmall[k,,drop=F]%*%beta)/tausq
			if(loc[1]>0) b<-b+t(wsmall[loc])%*%ufmat[1:length(loc),k]/(ufmat[m+1,k]*sigsq)
			if(0<(sivec[[k]][1])){
				cijvec<-sapply(sivec[[k]],cijgen,k,ufmat,imvec)
				tijvec<-sapply(sivec[[k]],tijgen,k,ufmat,imvec,wsmall)
				a<-a+sum(cijvec*cijvec/ufmat[m+1,sivec[[k]]])/sigsq
				b<-b+sum(cijvec*tijvec/ufmat[m+1,sivec[[k]]])/sigsq
				}
			wsmall[k]<-rnorm(1,b/a,1/sqrt(a))
			}
		#wsmall=rep(0,nsmall)
		#wmean=sapply(ulocs,function(x) mean(res1[intersect(dt$locsmap[[x]],locin)]))
		#res1small=ysmall-as.vector(Xsmall%*%beta)
		#res2small=ysmall-as.vector(Xsmall1%*%beta)
		#covwsmall=solve(diag(count)/tausq+solve(exp(-phi*distmat))/sigsq)
		#wsmall=as.vector(rmvnorm(1,as.vector(covwsmall%*%res1small*count/tausq),covwsmall))
		#wsmall2=covwsmall%*%res2small*count/tausq+wsmallmat[i,]
		#print(sum((wsmall-wsmall2)^2))
		wsmallmat[i,]=wsmall
		w=wsmall[invmap]
		#w=rep(0,n)
		wmat[i,]=w

		## updating spatial variance
		ew<-sapply(1:nsmall,ewgen,ufmat,imvec,wsmall)
		ew<-sum(ew*ew/ufmat[m+1,])
		sigsq<-rigamma(1,2+nsmall/2,0.2+0.5*ew)
		#sigsq=0.01
		sigsqvec[i]=sigsq

		## updating the spatial range
		phinew<-aphi+dphi*inv.logit(rnorm(1,logit((phi-aphi)/dphi),ltune))
		ufmatnew<-sapply(1:nsmall,ufgen,distmat,imvec,phinew,m,nsmall)
		ewnew<-sapply(1:nsmall,ewgen,ufmatnew,imvec,wsmall)
		ewnew<-sum(ewnew*ewnew/ufmatnew[m+1,])
		lacceptance<-0.5*(sum(log(ufmat[m+1,])-log(ufmatnew[m+1,]))+(ew-ewnew)/sigsq)+logJlogit(phinew,aphi,bphi)-logJlogit(phi,aphi,bphi)
		#print(lacceptance)	
		flag<-0
		un<-runif(1,0,1)
		if(log(un)<= lacceptance){
			flag=1
			phi=phinew
			ufmat<-ufmatnew
			}
		flagvec[i]=flag
		#phi=1
		phivec[i]=phi
		}else{
			sigsqvec[i]=0
			phivec[i]=0
			wmat[i,]=rep(0,n)
			flagvec[i]=0
			}		

		## printing the progress of the MCMC and checking acceptance rate is reasonable (roughly between (0.2,0.4))
		if((i %% 500)==0) {
			
			if(spat==1) plot(cumsum(flagvec[1:i])/(1:i))
			#print(all(dt$w[locin]==w))
			#print(head(ufmat))
			print(i)
			#print(round(colMeans(betamat[1:i,]),3))
		
			rate=mean(flagvec[1:i])
			if(rate<0.3) ltune=ltune/2
			if(rate > 0.5) ltune=ltune*2
			}			

		}

	## storing all the posterior samples
	#post=list(beta=betamat,tausq=tausqvec,sigsq=sigsqvec,phi=phivec,w=wmat)
	#if(spat==1) post$wsmall=wsmallmat

	## removing pre burn-in samples
	betamat=betamat[N1:N,,drop=F]
	wmat=wmat[N1:N,]
	if(spat==1)	wsmallmat=wsmallmat[N1:N,]
	tausqvec=tausqvec[N1:N]
	sigsqvec=sigsqvec[N1:N]
	phivec=phivec[N1:N]
	deltamat=betamat+wmat%*%X%*%solve(t(X)%*%X+1e-5*diag(p))

	## setting trait and predictor names for displaying the results table
	if(is.null(dt$names)) names=1:p else names=dt$names
	cnames=c(names,"tausq","sigsq","phi")
	cnameslow=sapply(cnames,function(x) paste(x,"_low",sep=""))
	cnamesup=sapply(cnames,function(x) paste(x,"_up",sep=""))
	cnamesfull=as.vector(rbind(cnames,cnameslow,cnamesup))

	index=3*(sapply(names,function(x) which(cnames == x))-1)+1
	indexfull=as.vector(rbind(index,index+rep(1,p),index+rep(2,p)))
	varindex=length(cnamesfull)-(8:0)

	## table containing parameter estimates and 95% confidence intervals
	tabest=tabestadj=matrix(NA,1,length(cnamesfull))
	qvec=c(0.5,0.025,0.975)
	tabest[1,indexfull]=as.vector(apply(betamat,2,qntlgen))
	tabest[1,varindex]=c(quantile(tausqvec,qvec),quantile(sigsqvec,qvec),quantile(phivec,qvec))
	colnames(tabest)=cnamesfull
	if(is.null(dt$traitnames)) rname=1 else rname=dt$traitnames
	row.names(tabest)=rname

	tabestadj[1,indexfull]=as.vector(apply(deltamat,2,qntlgen))
	tabestadj[1,varindex]=c(quantile(tausqvec,qvec),quantile(sigsqvec,qvec),quantile(phivec,qvec))
	colnames(tabestadj)=cnamesfull
	row.names(tabestadj)=rname

	## computing model evaluation metrics: gpd, Rsq, DIC, RMSPE
	N2=N-N1+1

	## generating posterior samples for the traits
	fitmat=X%*%t(betamat)+t(wmat)+t(sqrt(tausqvec)*matrix(rnorm(n*N2),N2,n))
	fit=apply(fitmat,1,mean)
	fitsd=apply(fitmat,1,sd)
	
	## gpd score
	gpd=gpdgen_ind_nonspatial(y,fitmat)

	## pseudo Rsq
	Rsq=cor(y,fit)^2
  
	## dic
	dic=dicgen_ind_spatial(y,X,betamat,wmat,tausqvec)

	## MSPE
	locout=dt$locout
	nout=length(locout)
	Xout=dt$X[locout,,drop=F]
	yout=dt$y[locout]

	## generating posterior predictive samples
	outfitmat=Xout%*%t(betamat)+t(sqrt(tausqvec)*matrix(rnorm(nout*N2),N2,nout))

	if(spat==1){
	## out of sample spatial kriging
	ulocsout=unique(dt$locsmapinv[locout])
	nsmallout=length(ulocsout)
	wsmalloutmat=sapply(ulocsout,nngpkrig,ulocs,s,wsmallmat,phivec,sigsqvec,dt$locs,dt$locsmap,m,cart)
	woutmat=wsmalloutmat[,sapply(locout,function(x) which(ulocsout==dt$locsmapinv[x]))]
	outfitmat=outfitmat+t(woutmat)
	}
	subsample=sample(N2,Nout,replace=F)
	outfitmat=outfitmat[,subsample]
	outfitmat=exp(outfitmat) ### transforming to original exponential scale
	outfit=apply(outfitmat,1,mean)
	outfitsd=apply(outfitmat,1,sd)

	## coverage probability of out-of-sample samples
	cp=mean(sapply(1:nout,function(j) (exp(yout[j]) > quantile(outfitmat[j,],0.025))&(exp(yout[j]) < quantile(outfitmat[j,],0.975))))

	## mspe for one trait hold out and both traits hold out 
	rmspe=sqrt(mean((exp(yout)-outfit)^2))
	rmspelog=sqrt(mean((yout-log(outfit))^2))

	## table containing model evaluation metrics
	tabeval=matrix(c(Rsq,dic,gpd,rmspe,rmspelog,cp),1,6)
	rownames(tabeval)=c(rname)
	colnames(tabeval)=c("Rsq","DIC","GPD","RMSPE","RMSPEL","CP")

	t1=Sys.time()
	## output containing everything
	
	betamat=betamat[subsample,,drop=F]
	tausqvec=tausqvec[subsample]
	sigsqvec=sigsqvec[subsample]
	phivec=phivec[subsample]
	if(spat==1) wsmallmat=wsmallmat[subsample,]
	wmat=wmat[subsample,]
	deltamat=deltamat[subsample,]	

	postout=list(betamat=betamat,tausqvec=tausqvec,sigsqvec=sigsqvec,phivec=phivec,ltune=ltune,
		deltamat=deltamat,fit=fit,fitsd=fitsd,outfit=outfit,fitmat=fitmat,outfitmat=outfitmat,
		outfitsd=outfitsd,tabest=tabest,tabeval=tabeval,tabestadj=tabestadj,time=t1-t0)
	if(spat==1) postout=c(postout,list(woutmat=woutmat,wsmallmat=wsmallmat,s=s,wmat=wmat))

	postout
	}

## latitude longitude to cartesian co-ordinates
latlontocart=function(s){
	s=s/180*pi
	c(cos(s[1])*cos(s[2]),sin(s[1])*cos(s[2]),sin(s[2]))
	}


## function to generate 2.5%, 50% and 97.5% quantiles
qntlgen<-function(v){
	qnvec=c(0.5,0.025,0.975)
	quantile(v,probs=qnvec)
	}

## function to calculate gpd score
gpdgen_ind_nonspatial=function(y,fitmat){
	m=rowMeans(fitmat)
	v=rowVars(fitmat)
	mean((y-m)^2)+mean(v)
	}	

## function to calculate dic score
dicgen_ind_spatial=function(y,X,betamat,wmat,tausqvec){
	n=length(y)
	ymeanmat=X%*%t(betamat)+t(wmat)
	avgdev=n*mean(log(tausqvec))+mean(colSums((ymeanmat-y)^2)/tausqvec)
	beta=colMeans(betamat)
	w=colMeans(wmat)
	ymean=as.vector(X%*%beta)+w
	tausq=mean(tausqvec)
	devavg=n*log(tausq)+sum((y-ymean)^2)/tausq
	2*avgdev-devavg
	}
