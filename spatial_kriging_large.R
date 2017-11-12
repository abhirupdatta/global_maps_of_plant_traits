#library(Matrix)
#library(plotrix)
#library(corrplot)
library(fields)
library(parallel)

### CHANGES--- NOW GIVES MCMC BASED KRIGING RESULTS USING 5000 SAMPLES

### function to create the processed file for predictor data at  
### new locations from the raw file pred_locs
### it removes NAs and scales the PFT weights
process_preds=function(predfname,intercept_names){
	#### loading the raw data file for predictors at new locations
	pred_locs=read.csv(predfname)

	#### removing NAs ####
	filtered_preds=pred_locs[complete.cases(pred_locs),]

	### scaling the weights so that they sum up to one, this section needs to be modified if we are using 15 PFTs
	filtered_preds$totweight=rowSums(filtered_preds[,intercept_names])
	filtered_preds=filtered_preds[which(filtered_preds$totweight >0),]
	#filtered_preds[,c("pftN" ,  "pftB" ,  "pftS" ,  "pftG")]=filtered_preds[,c("pftN" ,  "pftB" ,  "pftS" ,  "pftG")]/filtered_preds$totweight

	#row.names(filtered_preds)=NULL
	
	### saving the processed file, so that the processing step can be avoided in future
	write.csv(filtered_preds,"filtered_preds.csv",row.names=F)
	rm(list="pred_locs")
	filtered_preds

    }

#### kriging at all new locations
#### Xout= explanatory variables at new locations
#### locsout= new locations, subpost = fit file containing parameter estimates for the data subgrp
#### Nout=number of samples, m=number of neighbors, 
#### cart=indicator if data is in cartersian coordinates (always 0 for us)
#### spat=spatial model or not
#### weights = matrix of abundance for each pft at each new location
#### seed = random seed
#### range of quantiles from which the samples will be stored, currently it is (25%,75%)
masskrig=function(Xout,locsout,subpost,Nout,m,cart,spat,weights,seed,qrange){
    
    set.seed(seed)
    Nin=length(subpost$tausq)
    #Nout=Nin/(qrange[2]-qrange[1])
    #subsample=sample(Nin,Nout,replace=F)
    subsample=1:Nin
    betamat=subpost$betamat[subsample,]
    tausqvec=subpost$tausqvec[subsample]
    nout=nrow(Xout)
    
    ## generating posterior predictive samples
    print(dim(Xout))
    print(dim(betamat))
    
    locsoutnew=locsout
    if(cart==0) locsoutnew=t(apply(locsout,1,latlontocart))
    locsin=subpost$s
    wsmallmat=subpost$wsmallmat[subsample,]
    sigsqvec=subpost$sigsqvec[subsample]
    phivec=subpost$phivec[subsample]
    
    #no_cores <- detectCores() - 1
    #cl <- makeCluster(no_cores)
    #clusterExport(cl, varlist=c("Xout","betamat","tausqvec","nout","Nin","Nout","weights","locsoutnew","cart","spat",
    #  "locsin","wsmallmat","sigsqvec","phivec","m","nngpkrig2"), envir=environment())
    
    #outfitmat=Xout%*%t(betamat)+t(sqrt(tausqvec)*matrix(rnorm(nout*Nin),Nin,nout))
    #outfitmat[which(weights==0),]=0                                       
    #rowMeans(exp(outfitmat))
    #v=mu=w=rep(0,nout)
    #if(spat==1){
    ## out of sample spatial kriging
    
    youtmat=t(sapply(1:nout,nngpkrig2,Xout,betamat,tausqvec,weights,locsoutnew,locsin,wsmallmat,sigsqvec,phivec,Nin,Nout,m))
    fitmat=cbind(locsout,weights,youtmat)
    colnames(fitmat)=c("lon","lat","weights","mean","sd","logmean","logsd","wmean","wsd","freqwmean","freqwsd",
                       sapply(1:Nout, function(x) paste("sample",x,sep="")))
    fitmat
    
}


## latitude longitude to cartesian co-ordinates
latlontocart=function(s){
  s=s/180*pi
  c(cos(s[1])*cos(s[2]),sin(s[1])*cos(s[2]),sin(s[2]))
}

nngpkrig2=function(k,Xout,betamat,tausqvec,weights,locsout,locsin,wsmallmat,sigsqvec,phivec,Nin,Nout,m){
    Nin=as.numeric(Nin)
    ymean=ysd=lymean=lysd=wmean=wsd=flwmean=flwsd=0
    yout=rep(0,Nout)
    wout=mu=v=rep(0,Nin)
    #wout=mu=v=rep(0,Nout)
    if(weights[k]>0){
        outfit=Xout[k,]%*%t(betamat)+sqrt(tausqvec)*rnorm(Nin)
        locout=locsout[k,]
        distvec=as.vector(rdist(t(locout),locsin))
        im=order(distvec)[1:m]
        #print(im)
        if(spat==1){
            wtmat=sapply(phivec,function(phi) as.vector(solve(exp(-phi*as.matrix(dist(locsin[im,]))))%*%exp(-phi*distvec[im])))
            mu<-sapply(1:Nout, function(j) sum(wsmallmat[j,im]*wtmat[,j]))
            v<-sigsqvec*(1-sapply(1:length(phivec),function(j) sum(exp(-phivec[j]*distvec[im])*wtmat[,j]) ))
            v=sqrt(pmax(0,v))
            wout=mu+v*rnorm(length(phivec))
        }
        ly=outfit+wout
        lymean=mean(ly)
        lysd=sd(ly)
        yout=exp(ly)
        ymean=mean(yout)
        ysd=sd(yout)
        wmean=mean(wout)
        wsd=sd(wout)
        flwmean=mean(mu)
        flwsd=mean(v)
        set.seed(1)
        yout=yout[sample(Nin,Nout,replace=FALSE)]
        #if((k%%100)==0) print(paste0("kriging loc: ",k))
    }
    c(ymean,ysd,lymean,lysd,wmean,wsd,flwmean,flwsd,yout)
}
