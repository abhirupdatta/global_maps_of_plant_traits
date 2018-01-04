library(plyr)
library(lattice)
library(latticeExtra)
library(maps)
library(matrixStats)

main=function(trait,maindatafile,outpct,model,level,spat,inits,const,N,N1,cart,m,Nout){
      
  print("Data processing")
  if(!dir.exists(trait)) dir.create(trait)
  set.seed(seed)
  
  trydt=read.csv(maindatafile)
  dt=trydt[,setdiff(colnames(trydt),
    c("X","ObservationID","PFTs","correct.names","dists"))]
  
  species=trydt$correct.names
  location=dists=trydt$dists
  PFTs=trydt$PFTs
  
  #### this is the species-location-PFT level aggregated dataset
  if(!file.exists(paste0("agg_",maindatafile))){
    agg_try30=aggregate(x=dt, FUN=mymean, by=list(species=species,dists=dists,PFTs=PFTs))
    agg_try30$superPFT=sapply(as.character(agg_try30$PFTs),superpftgen)
    row.names(agg_try30)=NULL
    write.csv(agg_try30,paste0("agg_",maindatafile),row.names=FALSE)
    rm(list=c("dt"))
  }else{
    agg_try30=read.csv(paste0("agg_",maindatafile))
  }
  
  agg_try30=subset(agg_try30,!is.na(get(trait)))
  agg_try30$global=1  
  
  #### holdout locations ####
  outlocs=stratifiedholdoutlocations(agg_try30,trait,outpct)$outlocs
  
  #### categorical models fits ####
  print("Categorical model")
  model="cat"
  tr_all <- data.frame(sla = trydt$sla,
    lnm = trydt$lnm,
    lpm = trydt$lpm,
    loc = trydt$dists, 
    lon = trydt$LON_site, lat = trydt$LAT_site, 
    spec = trydt$correct.names, PFT=trydt$PFTs, 
    myT=trydt$myT, rad=trydt$rad, wet=trydt$wet, 
    MI=trydt$MI, pH=trydt$pH, cly=trydt$cly)
  tr_all <- tr_all[!is.na(tr_all$spec),]
  tr_all$superPFT=sapply(as.character(tr_all$PFT),superpftgen)
  
  # omitted locations
  loc_in <- setdiff(unique(tr_all$loc),outlocs)
  pl <- which(tr_all$loc %in% loc_in)
  tr_trait <- tr_all[pl,]
  tr_trait$all="global"
  
  grp=ifelse(level=="global","all",ifelse(level=="superpft","superPFT","PFT"))
  pft_trait_om <- ddply(tr_trait,.(grp=get(grp),spec),summarise,
    sla=mean(sla,na.rm=T),lnm=mean(lnm,na.rm=T),
    lpm=mean(lpm,na.rm=T) )
  
  cat <- vector("list",length=uniqlength(pft_trait_om$grp))
  if(level=="global"){names(cat)="global"}else{ if(level=="superpft"){
    names(cat)=superPFTs}else{names(cat)=PFTnames}}
  
  for(i in 1:length(names(cat))){
    pl <- which(pft_trait_om$grp==names(cat)[i])
    tmp <- na.exclude(pft_trait_om[pl,trait])
    cat[[i]]$vals <- tmp[1:length(tmp)]
    cat[[i]]$mean <- mean(tmp)
  }
  catmeans=sapply(cat,function(x) x$mean)
  catmeans1=catmeans
  catlow=sapply(cat,function(x) quantile(x$vals,0.025))
  cathigh=sapply(cat,function(x) quantile(x$vals,0.975))
  outdata=subset(agg_try30,dists %in% outlocs)[,c("PFTs","superPFT","global",trait)]
  grp=ifelse(level=="global","global",ifelse(level=="superpft","superPFT","PFTs"))
  # re-group categories to preserve outdata ordering
  if (level != "global") {
    catmeans = catmeans[levels(outdata[,grp])]
    catlow = catlow[paste0(levels(outdata[,grp]),'.2.5%')]
    cathigh = cathigh[paste0(levels(outdata[,grp]),'.97.5%')]
  }
  outdata$mean=sapply(outdata[,grp],function(i) catmeans[i])
  outdata$low=sapply(outdata[,grp],function(i) catlow[i])
  outdata$high=sapply(outdata[,grp],function(i) cathigh[i])
  
  #### adding Rsq numbers for cat model with the in-sample data ####
  indata=subset(agg_try30,!(dists %in% outlocs))[,c("PFTs","superPFT","global",trait)]
  indata$mean=sapply(indata[,grp],function(i) catmeans[i])
  
  evalgen(outdata,indata,trait,model,spat,level)
  
  mapmat=filtered_preds[,c("lon","lat")]
  colname1=paste(model,spat,level,trait,sep="_")
  
  filtered_preds$global=1
  filtered_preds[,PFTnames]=filtered_preds[,PFTnames]/100
  filtered_preds[,superPFTs]=sapply(1:4, function(i) rowSums(filtered_preds[,pftgrps[[i]]]))
  if(level=="global"){ intercepts="global"} else{ if(level=="pft"){
    intercepts=pftgrps} else intercepts=superPFTs}
  
  mapweight=as.matrix(cbind(0,filtered_preds[,unlist(intercepts)]))
  mapmat[,colname1]=as.vector(mapweight%*%c(0,as.vector(catmeans1)))/rowSums(mapweight)
  
  rm(list=c("trydt"))
  
  ######### spatial model ##########
  print("Spatial model")
  model="spat"
  colname1=paste(model,spat,level,trait,sep="_")
  agg_try30[,PFTnames]=t(sapply(agg_try30$PFTs,pfttobin))
  
  logtrait=sapply(trait,function(x) paste("log",x,sep="")) ## Column name for log(plant trait)
  agg_try30[,logtrait]=log(agg_try30[,trait])
  # if(level=="global"){ intercepts="global"} else{ if(level=="pft"){
  #   intercepts=PFTnames} else intercepts=superPFTs}
  # allpreds=c(intercepts,envvar) ## adds PFT specific intercept to the list of predictors
  
  agg_try30[,superPFTs]=t(apply(agg_try30[,PFTnames],1,function(x) {
    sapply(1:4,function(i) sum(x[pftgrps[[i]]]))}))
  
  ## Current column which decides the PFT groupings    
  source(filefit)
  source(filemodel)
  
  meanvec=as.vector(colMeans(filtered_preds[,envvar]))
  sdvec=as.vector(apply(filtered_preds[,envvar],2,sd))
  
  agg_try30[,envvar]=t((t(as.matrix(agg_try30[,envvar]))-meanvec)/sdvec)
  
  grp=ifelse(level=="global","global","superPFT")
  
  ## creating the data tables for pft subgroups (4 superpfts as of now, can be 15 pfts as well)
  subdata=split(agg_try30,agg_try30[,grp])
  ngrp=length(subdata)
  
  ## loading holdout locations
  uniqlocout=outlocs
  
  ## running the model on each subdata
  ## the plots printed for each run is the acceptance rate for the metropolis sampler
  ## the overall acceptance rate should be between 0.2 and 0.5
  ## the number printed is the MCMC iteration number for each run
  
  allpreds=lapply(intercepts,c,envvar)
  if(level=="global"){perm=1}else{perm=sapply(1:4,
    function(i) which(unique(subdata[[i]][,grp])==superPFTs))}
  
  dev.new()
  subpost1=outfitmat=infitmat=list()
  print("Running the spatial model. MCMC iterations will be displayed")
  for(i in 1:ngrp){
    subpost1[[i]]=model_spatial_sub(subdata[[i]],allpreds[[perm[i]]],
      logtrait,N,N1,N-N1+1,inits,const,m,cart,spat,seed,uniqlocout)
  }
  
  for(i in 1:ngrp){
    outqntl=t(apply(subpost1[[i]]$outfitmat,1,quantile,c(0.025,0.975)))
    outfitmat[[i]]=cbind(subpost1[[i]]$outfit,outqntl)
    colnames(outfitmat[[i]])=c("mean","low","high")
    outfitmat[[i]]=cbind(subdata[[i]][subpost1[[i]]$locout,c("dists","PFTs",
      "superPFT","global",trait,logtrait)],outfitmat[[i]])
    
    infitmat[[i]]=cbind(subdata[[i]][subpost1[[i]]$locin,c("dists","PFTs",
      "superPFT","global",trait,logtrait)],exp(subpost1[[i]]$fit))
    colnames(infitmat[[i]])[ncol(infitmat[[i]])]="mean"
    
    }
  outfit=Reduce('rbind',outfitmat)
  infit=Reduce('rbind',infitmat)
  evalgen(outfit,infit,trait,model,spat,level)
  
  print(paste0("Model evluation metrics are stored in eval_large.csv in ",trait," folder"))
  #print("Saving barplots of RMSPE and Coverage probability")
  
  #mybarplot(trait) # this needs to come after all the models have been run
  
  source(filekrig)
  
  filtered_preds[,envvar]=t((t(filtered_preds[,envvar])-meanvec)/sdvec)
  Xoutlist=lapply(intercepts,function(cols) {m1=as.matrix(filtered_preds[,cols]);
    m1=m1/rowSums(m1);cbind(m1,as.matrix(filtered_preds[,envvar]))})
  
  #qrange=c(0.25,0.75)	## Upper and lower quantiles of the samples, we are only going to store the samples within this range
  qrange=c(0,1)	
  locsout=as.matrix(filtered_preds[,c("lon","lat")])  ### longitude then latitude to be compatible with the model estimation
  
  ## matrix of PFT abundance weights
  if(grp=="global"){
    weightmat=matrix(rowSums(filtered_preds[,superPFTs]),ncol=1)
    grpnum=1
      }else{
      weightmat=filtered_preds[,superPFTs]
      grpnum=match(superPFTs,sapply(1:4, function(i) unique(subdata[[i]]$superPFT)))
      }
  totweights=rowSums(weightmat)
  ### filenames to store the krigged distributions for each PFT and each location
  if(grp=="global") grplabels=grp else grplabels=superPFTs
  filenames=sapply(grplabels,function(x,spat) paste(trait,"/pred_",level,"_",spat,"_",x,"_large.csv",sep=""),spat)
  mapnames=sapply(grplabels,function(x,spat) paste(trait,"/pred_",level,"_",spat,"_",x,"_large.pdf",sep=""),spat)
  pred_tables=list()
  
  ### kriging ###
  #pred_tables=list()
  sampcols=paste0('sample',1:Nout)
  print("Kriging at 53900 locations for creating maps")
  colname2=paste(model,spat,level,trait,sep="_")
  for (i in 1:ngrp){
    print(i)
    pred_table=masskrig(Xoutlist[[i]],locsout,subpost1[[grpnum[i]]],Nout,m,cart,spat,weightmat[,i],seed,qrange)
    #pred_table[which(pred_table[,"sd"]==0),c("mean",sampcols)]=0
    pred_tables[[i]]=as.data.frame(pred_table)
    write.csv(pred_table,filenames[[i]],row.names=F,quote=F)
    map_plotter(pred_tables[[i]],paste0(trait,"/",grplabels[i],"_mean_large.pdf"),paste0("Spatial model: ",grplabels[i]," mean"),
                "mean",at="range",palette=terrain.colors,weights=weightmat[,i]/totweights)
    map_plotter(pred_tables[[i]],paste0(trait,"/",colname1,"_",grplabels[i],"_sd_large.pdf"),paste0("Spatial model: ",grplabels[i]," sd"),
                "sd",at="range",palette=terrain.colors,weights=weightmat[,i]/totweights)
    if(spat==1) map_plotter(pred_tables[[i]],paste0(trait,"/",colname1,"_",grplabels[i],"_wmean_large.pdf"),paste0("Spatial model: ",grplabels[i]," wmean"),
                "wmean",at="posrange",palette=terrain.colors)
    if(spat==1) map_plotter(pred_tables[[i]],paste0(trait,"/",colname1,"_",grplabels[i],"_wsd_large.pdf"),paste0("Spatial model: ",grplabels[i]," wsd"),
                "wsd",at="posrange",palette=terrain.colors)
    map_plotter(pred_tables[[i]],paste0(trait,"/",colname1,"_",grplabels[i],"_logmean_large.pdf"),paste0("Spatial model: ",grplabels[i]," logmean"),
                "logmean",at="posrange",palette=terrain.colors)
    map_plotter(pred_tables[[i]],paste0(trait,"/",colname1,"_",grplabels[i],"_logsd_large.pdf"),paste0("Spatial model: ",grplabels[i]," logsd"),
                "logsd",at="posrange",palette=terrain.colors)
    if(spat==1) map_plotter(pred_tables[[i]],paste0(trait,"/",colname1,"_",grplabels[i],"_freq_logmean_large.pdf"),paste0("Spatial model: ",grplabels[i]," freq_logmean"),
                "freqwmean",at="posrange",palette=terrain.colors)
    if(spat==1) map_plotter(pred_tables[[i]],paste0(trait,"/",colname1,"_",grplabels[i],"_freq_logsd_large.pdf"),paste0("Spatial model: ",grplabels[i]," freq_logsd"),
                "freqwsd",at="posrange",palette=terrain.colors)
    }
  
  weightedfits=lapply(1:ngrp,function(i) pred_tables[[i]][,
    paste0("sample",1:Nout)]*weightmat[,i])
  avgfits=cbind(locsout,Reduce('+',weightedfits)/totweights)
  avgfits=as.data.frame(avgfits)
  #avgfits$mean=rowMeans(avgfits[,sampcols])
  avgfits$mean=Reduce('+',lapply(1:ngrp,function(i) pred_tables[[i]][,"mean"]*weightmat[,i]))/totweights
  avgfits$sd=sqrt(Reduce('+',lapply(1:ngrp,function(i) pred_tables[[i]][,"sd"]^2*weightmat[,i]^2))/(totweights^2))
  write.csv(avgfits,paste0(trait,"/",colname2,"_large.csv"))
  
  mapmat[,colname2]=avgfits$mean
  mapmat[,paste0(colname2,"_sd")]=avgfits$sd
  
  #traitrange=seq(quantile(agg_try30[,trait],c(0.05,0.95))[1],
  #  quantile(agg_try30[,trait],c(0.5,0.95))[2],length=21)
  
  qntl=quantile(c(mapmat[,colname1],mapmat[,colname2]),c(0.025,0.975))
  traitrange=seq(qntl[1],qntl[2],length=21)
  
  print(paste0("Maps are created and stored in ",trait," folder"))
  map_plotter(mapmat,paste0(trait,"/",colname1,"_large.pdf"),paste0("Categorical model ",trait),
    colname1,at=traitrange,palette=terrain.colors)
  map_plotter(mapmat,paste0(trait,"/",colname2,"_large.pdf"),paste0("Spatial model ",trait),
    colname2,at=traitrange,palette=terrain.colors)
  map_plotter(mapmat,paste0(trait,"/",colname2,"_sd_large.pdf"),paste0("s.d. of Spatial model ",trait),
    paste0(colname2,"_sd"),at="qntl",palette=terrain.colors)
  
  ## name of the R workspace 
  mydate=gsub(":","_",gsub(" ","_",date()))
  fitresults=paste(trait,"/large_fitresults_",model,"_",level,
    "_",spat,"_",paste(envvar,collapse="_"),"_",mydate,".Rdata",sep="")
  
  print(paste0("Saving fits in ",trait," folder"))
  ## storing all important objects related to model fitting
  save(subpost1,outlocs,trait,level,seed,outpct,N,N1,Nout,inits,const,m,outdata,indata,subdata,
    outfit,infit,pred_tables,pred_table,mapmat,nPFT,PFTnames,pftgrps,superPFTs,intercepts,allpreds,
    cat,catmeans,weightmat,Xoutlist,perm,grpnum,meanvec,sdvec,logtrait,envvar,spat,cart,file=fitresults)  
  
}

####### function to save barplots of RMSPE and Cov. Prob. for all models
mybarplot=function(trait){
  tab=as.matrix(read.csv(paste0(trait,"/eval_large.csv"),row.names=1) )
  row.names(tab)=sapply(row.names(tab), function(x) unlist(strsplit(x,"_"))[2])  
  
  filename=paste0(trait,"/rmspe_14pft.png")
  png(filename,height=1024,width=1920)
  barplot(t(tab[20:38,paste0(c("cat_1_pft_","spat_1_pft_"),trait)]),
    beside=TRUE,args.legend = list(x="topleft",bty='n',cex=2),cex.names=1.2,
    legend=c("Cat.","Spat."),main=paste0(trait,": RMSPE of models using 14 PFTs"),cex.main=2)
  dev.off()
  
  filename=paste0(trait,"/RMSPE_4superpft.png")
  png(filename,height=1024,width=1920)
  barplot(t(tab[20:38,paste0(c("cat_1_superpft_","spat_1_superpft_"),trait)]),
    beside=TRUE,args.legend = list(x="topleft",bty='n',cex=2),cex.names=1.2,
    legend=c("Cat.","Spat."),main=paste0(trait,": RMSPE of models using 4 super PFTs"),cex.main=2)
  dev.off()
  
  filename=paste0(trait,"/RMSPE_global.png")
  png(filename,height=1024,width=1920)
  barplot(t(tab[20:38,paste0(c("cat_1_global_","spat_1_global_"),trait)]),
    beside=TRUE,args.legend = list(x="topleft",bty='n',cex=2),cex.names=1.2,
    legend=c("Cat.","Spat."),main=paste0(trait,": RMSPE of global models using no PFT information"),cex.main=2)
  dev.off()
  
  filename=paste0(trait,"/cov_prob_14pft.png")
  png(filename,height=1024,width=1920)
  barplot(t(tab[39:57,paste0(c("cat_1_pft_","spat_1_pft_"),trait)]),
    beside=TRUE,args.legend = list(x="topleft",bty='n',cex=2),cex.names=1.2,
    legend=c("Cat.","Spat."),main=paste0(trait,": Coverage prob. of models using 14 PFTs"),cex.main=2)
  abline(h=0.95)
  dev.off()
  
  filename=paste0(trait,"/cov_prob_4superpft.png")
  png(filename,height=1024,width=1920)
  barplot(t(tab[39:57,paste0(c("cat_1_superpft_","spat_1_superpft_"),trait)]),
    beside=TRUE,args.legend = list(x="topleft",bty='n',cex=2),cex.names=1.2,
    legend=c("Cat.","Spat."),main=paste0(trait,": Coverage prob. of models using 4 super PFTs"),cex.main=2)
  abline(h=0.95)
  
  dev.off()
  filename=paste0(trait,"/cov_prob_global.png")
  png(filename,height=1024,width=1920)
  barplot(t(tab[39:57,paste0(c("cat_1_global_","spat_1_global_"),trait)]),
    beside=TRUE,args.legend = list(x="topleft",bty='n',cex=2),cex.names=1.2,
    legend=c("Cat.","Spat."),main=paste0(trait,": Coverage prob. of global models using no PFT information"),cex.main=2)
  abline(h=0.95)
  dev.off()
  
  
}

mymean=function(x) { if(all(is.na(x))) y=NA else y=mean(x,na.rm=TRUE);y}

myequal=function(x,y){ if(is.na(x)){if(is.na(y)) z=1 else z=0} else{
  if(is.na(y)) z=0 else z=1*(abs(x-y)<1e-5)};z}

uniqlength=function(x) length(unique(x))

#### choosing holdout locations stratified by latitude ####
stratifiedholdoutlocations=function(agg_try30,trait,outpct){
  #agg_try30_sub=subset(agg_try30,!is.na(get(trait)))
  agg_tropic=subset(agg_try30,(abs(LAT_site) < 23.5))
  agg_temp=subset(agg_try30,(abs(LAT_site) > 23.5) & (abs(LAT_site) < 66.5))
  agg_boreal=subset(agg_try30,(abs(LAT_site) > 66.5))
  
  locs_tropic=holdoutlocations(agg_tropic,trait,outpct)
  locs_temp=holdoutlocations(agg_temp,trait,outpct)
  locs_boreal=holdoutlocations(agg_boreal,trait,outpct)
  
  outlocs=c(locs_tropic$outlocs,locs_temp$outlocs,locs_boreal$outlocs)
  list(outlocs=outlocs,tropic=locs_tropic,temp=locs_temp,boreal=locs_boreal)
}

#### choosing holdout locations ####
holdoutlocations=function(tab,trait,outpct){
  tab_sub=subset(tab,!is.na(get(trait)))
  numobs=with(data=tab_sub,aggregate(LAT_site, FUN=length, by=list(dists=dists)))
  numobs=numobs[order(numobs$x,decreasing = TRUE),]
  numout=round(nrow(numobs)*outpct)
  num=nrow(numobs)
  numspecloc=sum(numobs$x)
  outlocs=sample(numobs$dists,numout)
  numspeclocout=sum(numobs$x[which(numobs$dists %in% outlocs)])
  ratio=round(numspeclocout/numspecloc,2)
  print(ratio)
  list(outlocs=outlocs,num=nrow(numobs),numout=numout,numspecloc=numspecloc,
    numspeclocout=numspeclocout,outpctspecloc=ratio)
}

#### superpfts
superpftgen=function(pft) {num=as.numeric(unlist(strsplit(pft,"T"))[2]);
c(rep("Needleleafs",3),rep("Broadleafs",5),rep("Shrubs",3),rep("Grasses",3))[num]}

#### creates binary vectors based on PFT status
pfttobin=function(pftx){
  y=rep(0,nPFT);
  y[as.numeric(unlist(strsplit(as.character(pftx),"T"))[2])]=1;
  y
}

#### Rsq (calculated at the logarithmic scale which is used for the regression) ####
rsqgen=function(tab,trait){
  cor(log(tab$mean),log(tab[,trait]))^2
}

#### rmspe ####
rmspegen=function(tab,trait){
  sqrt(mean(((tab[,trait])-(tab$mean))^2))
}

#### cp ####
cpgen=function(tab,trait){
  mean((tab[,trait]> tab$low)&(tab[,trait]< tab$high))
}

#### function which calculates evaluation numbers ####
evalgen=function(tab,intab,trait,model,spat,level){
  print(paste0("Running model evaluation for ",model," model"))
  grps=c("global","superPFT","PFTs")
  rmspe=sapply(grps,function(g) sapply(split(tab,tab[,g]),rmspegen,trait))
  rmspe$superPFT=rmspe$superPFT[superPFTs]
  rmspe$PFTs=rmspe$PFTs[PFTnames]
  
  cp=sapply(grps,function(g) sapply(split(tab,tab[,g]),cpgen,trait))
  cp$superPFT=cp$superPFT[superPFTs]
  cp$PFTs=cp$PFTs[PFTnames]

  print("Rsq numbers may be NA for the categorical model which will give warnings")
  rsq=sapply(grps,function(g) sapply(split(intab,intab[,g]),rsqgen,trait))
  rsq$superPFT=rsq$superPFT[superPFTs]
  rsq$PFTs=rsq$PFTs[PFTnames]
  
  m=unlist(rsq)
  names(m)=NULL
  #m
  #rm(list=setdiff(ls(),c("agg_try30","grps","rsqgen")))
  
  colname=paste(model,spat,level,trait,sep="_")
  path=paste0(trait,"/eval_large.csv")
  if(!file.exists(path)){
    evaltab=matrix(c(unlist(rsq),unlist(rmspe),unlist(cp)),ncol=1)
    row.names(evaltab)=unlist(sapply(c("rsq","rmspe","cp"),function(x) paste0(x,"_",
      c("global",superPFTs,PFTnames))))
    colnames(evaltab)=colname
  }else{
    evaltab=read.csv(path,row.names = 1)
    evaltab[,colname]=c(unlist(rsq),unlist(rmspe),unlist(cp))
  }
  write.csv(evaltab,path)
}

map_plotter=function(tab,figname,title,colname,at,palette=terrain.colors,weights=1){
  #tab=read.csv(filename)
  
  locs <- matrix(nrow=nrow(tab),ncol=2)
  locs[,1] <- as.vector(tab$lon)
  locs[,2] <- as.vector(tab$lat)
  vals=tab[,colname]*weights
  
  worldmap <- map('world',plot=F)
  world.df <- data.frame(lon=worldmap$x,lat=worldmap$y)
  if (length(at)==1) { 
  if(at=="range") { at=seq(min(vals),max(vals),length=21) }
  else if(at=="evenrange") { at=seq(min(even(vals)),max(even(vals)),length=21) }
  else if(at=="qntl") { at=unique(quantile(vals,p=(1:19)/20)) }
  else if(at=="qr") { at=seq(round(100*quantile(vals,p=0.1)),round(100*quantile(vals,p=0.9)),length=128)/100 }
  else if(at=="posrange") { at=seq(max(min(vals),min(vals[which(vals>0)]),na.rm=TRUE),max(vals),length=21) }
  }
  print(paste0("Range of values: (",round(min(vals),1),",",round(max(vals),1),")"))
  print(paste0(round(100*mean((vals<max(at))&(vals>min(at)))),
    " % of the model fits are within the empirical range"))
  vals=pmin(pmax(vals,min(at)),max(at))
  pdf(figname)
  #at=seq(0,round(max(vals))+1,1)
  #at=quantile(vals,p=(0:10)/10)
  #at=seq(min(vals),max(vals),length=21)
  print(levelplot(vals~locs[,1]+locs[,2],
    xlab = "Longitude", ylab="Latitude",at=at, 
    main = list(title,cex=1.9),col.regions=rev(palette(128)),
    add=T) + 
      xyplot(lat~lon,world.df,
        type='l',lty=1,lwd=1,col="black"))
  dev.off()
  
}