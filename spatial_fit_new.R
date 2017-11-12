library(Matrix)
library(plotrix)
library(corrplot)

## Processes the rawdata, adds explanatory variables, scales and centers them,
## adds intercept, filters missing data and redundant columns
## (parameters are described in the fnction call in spatial_main.R)
#process=function(traitnames,envvar,rawdata,glob_preds){
#	pred_nms <- names(glob_preds)
#	d <- dim(rawdata)
#	soil_clim <- matrix(nrow=d[1],ncol=10)	## matrix to store all predictors 
#	for (i in seq(1,10)){
#	  soil_clim[,i] <- as.vector(glob_preds[[i+4]][loc])
#	}
#	soil_clim=cbind(1,soil_clim)	## adding intercept
#	colnames(soil_clim) <- c("Intercept",pred_nms[5:14])
#	rawdata <- cbind(rawdata,soil_clim,loc)	## append predictor variables to rawdata
#
#	allpreds <- c("Intercept",envvar) ## adds intercept to the list of predictors	
#	varlist=c("Observation","Lat","Lon",traitnames,allpreds,pftcol,"loc")
#
#	### subsetting the data where the trait and predictor variables are available
#	traitdata=subset(rawdata[,varlist],!is.na(rawdata[,traitnames]))  ## filters missing traits
#	traitdata$ObsNum=1:nrow(traitdata)
#	
#	logtraitnames=sapply(traitnames,function(x) paste("log",x,sep=""))
#	traitdata[,logtraitnames]=log(traitdata[,traitnames])
#	
#	meanvec=sdvec=vector()
#	p=length(envvar)
#	if(p > 0){
#		traitdata=subset(traitdata,apply(!is.na(as.matrix(traitdata[,envvar])),1,all)) ## filters missing predictors
#		meanvec=apply(as.matrix(traitdata[,envvar]),2,mean)
#		sdvec=apply(as.matrix(traitdata[,envvar]),2,sd)
#		traitdata[,envvar]=t((t(as.matrix(traitdata[,envvar]))-meanvec)/sdvec) ##Z-scoring to make scales similar
#		}
#
#	list(traitdata=traitdata,meanvec=meanvec,sdvec=sdvec)
#	}

## format try30 data for spatial model
process = function(traitnames,envvar,allpreds,try30_pred,meanvec,sdvec){
  varlist <- c("LAT_site","LON_site",traitnames,allpreds,"dists","PFTs")
  traitdata <- subset(try30_pred[,varlist],!is.na(try30_pred[,traitnames]))
  #traitdata$Intercept <- rep(1,length(traitdata$Observation))
  traitdata$ObsNum <- 1:nrow(traitdata)
  logtraitnames=sapply(traitnames,function(x) paste("log",x,sep=""))
  traitdata[,logtraitnames]=log(traitdata[,traitnames])
  p=length(envvar)
  if(p > 0){
    traitdata=subset(traitdata,apply(!is.na(as.matrix(traitdata[,envvar])),1,all)) ## filters missing predictors
    #meanvec=apply(as.matrix(traitdata[,envvar]),2,mean)
    #sdvec=apply(as.matrix(traitdata[,envvar]),2,sd)
    traitdata[,envvar]=t((t(as.matrix(traitdata[,envvar]))-meanvec)/sdvec) ##Z-scoring to make scales similar
  }
  
  l=list(traitdata=traitdata)
  if(p>0) l=c(l,list(meanvec=meanvec,sdvec=sdvec))
  l
  }

## Processes each subdata (corresponding to pft grps)
## Runs the spatial model (parameters are described in the fnction call in spatial_main.R)
model_spatial_sub=function(subd,allpreds,
  logtraitnames,N,N1,Nout,inits,const,m,cart,spat,seed,uniqlocout){

	locs2=with(data=subd, aggregate(LAT_site,by=list(subd$dists), FUN=length))
	locs2=locs2[order(locs2$x),]	
	
	locsmap=sapply(1:nrow(locs2),function(i) which(subd$dists==locs2[i,1]))
	locsmapinv=sapply(1:nrow(subd),function(i) which(locs2[,1]==subd[i,"dists"]))

	##locout=sample(1:nrow(subd),nrow(subd)/10)
	locout=which(subd$dist %in% uniqlocout)
	print(nrow(subd))
	print(length(locout))
	locin=setdiff(1:nrow(subd),locout)
  
	#nPFT=length(intercept_names)
	#intercept_cols=intercept_names[which(colSums(subd[,intercept_names])>0)]
	#print(intercept_cols)
	#allpreds=c(intercept_cols,envvar)
	X=as.matrix(subd[,allpreds])
	
	y=as.vector(subd[,logtraitnames])

	sublist=list(y=y,X=X,locin=locin,locout=locout,traitnames=logtraitnames,names=allpreds,
		locsmap=locsmap,locsmapinv=locsmapinv,locs=as.matrix(subd[,c("LON_site","LAT_site")]))

	subpost=model_ind_spatial(sublist,N,N1,Nout,inits,const,m,cart,spat)
	subpost=c(subpost,list(seed=seed,locin=locin,locout=locout,locs2=locs2,locsmap=locsmap,locsmapinv=locsmapinv,subpreds=allpreds))
	subpost

	}