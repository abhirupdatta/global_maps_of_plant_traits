##### caller file ####
##### specifies all initial values and settings ####
##### runs the model ####

envvar=c("myT","rad","MI","pH","cly")
outpct=0.15 ### percentage of locations left out

## model fitting parameters
inits=list(tausq=2,sigsq=1,phi=100)  ## initial values for the spatial model
const=list(aphi=50,bphi=2000,ltune=1)  ## list of constants used by spatial model

#seed=sample(.Random.seed,1) ## random seed, this will be stored to reproduce the results in the future
#set.seed(seed)

N=10000 ## total number of MCMC iterations
N1=5001 ## number of MCMC iterations used for final results
Nout=100 ## number of samples used to create the trait distribution
m=5 ## number of neighbors used in Nearest Neighbor Gaussian Process
cart=0 ## If the Locations are in cartesian coordinates (always zero for our data)
spat=1 ## Spatial (spat=1) or non-spatial (spat=0) (just regression) model

maindatafile="try30_pred6.csv" ## need to rename SLA,LNM,LPM to lowercase
filtered_preds=read.csv("filtered_preds.csv")
filemodel="spatial_model.R" ## File containing the spatial model
filefit="spatial_fit_new.R" ## Name of file used to fit the model to the data 
filekrig="spatial_kriging_large.R" ## Name of file which uses the fitted model to predict at new locations 
seed=1
filemain="allmodels_large.R"

superPFTs=c("Needleleafs","Broadleafs","Shrubs","Grasses")
PFTnames=paste0("PFT",1:14)
nPFT=length(PFTnames)

pftgrps=list()
pftgrps[[1]]=c("PFT1","PFT2","PFT3")
pftgrps[[2]]=c("PFT4","PFT5","PFT6","PFT7","PFT8")
pftgrps[[3]]=c("PFT9","PFT10","PFT11")
pftgrps[[4]]=c("PFT12","PFT13","PFT14")

# note the modest change in map_plotter
settings=expand.grid(c("sla","lnm","lpm"),c("pft","superpft","global"),c(1,0))
#i=as.numeric(Sys.getenv('SGE_TASK_ID'))

i=1 ### run this for i =1, 2, 3, ..., 18 for all the runs
trait=as.character(settings[i,1])
level=as.character(settings[i,2])
spat=as.numeric(settings[i,3])

source(filemain)
main(trait,maindatafile,outpct,model,level,spat,inits,const,N,N1,cart,m,Nout)
mybarplot(trait)
# 
# main("sla",maindatafile,outpct,model,"pft",spat,inits,const,N,N1,cart,m,Nout)
# main("sla",maindatafile,outpct,model,"superpft",spat,inits,const,N,N1,cart,m,Nout)
# main("sla",maindatafile,outpct,model,"global",spat,inits,const,N,N1,cart,m,Nout)
# mybarplot("sla")
#  
# main("lnm",maindatafile,outpct,model,"pft",spat,inits,const,N,N1,cart,m,Nout)
# main("lnm",maindatafile,outpct,model,"superpft",spat,inits,const,N,N1,cart,m,Nout)
# main("lnm",maindatafile,outpct,model,"global",spat,inits,const,N,N1,cart,m,Nout)
# mybarplot("lnm")
# 
# main("lpm",maindatafile,outpct,model,"pft",spat,inits,const,N,N1,cart,m,Nout)
# main("lpm",maindatafile,outpct,model,"superpft",spat,inits,const,N,N1,cart,m,Nout)
# main("lpm",maindatafile,outpct,model,"global",spat,inits,const,N,N1,cart,m,Nout)
# mybarplot("lpm")