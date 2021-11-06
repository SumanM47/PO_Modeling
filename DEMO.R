rm(list=ls())

## compile the C codices
system("R CMD SHLIB ./auxiliary/aux_c.c -lm")

## load the compiled functions
##UNIX or Mac
dyn.load(paste(getwd(),"/auxiliary/aux_c.so",sep=""))
##Windows
# dyn.load(paste(getwd(),"/auxiliary/aux_c.dll",sep=""))

## load the R functions
source(paste(getwd(),"/auxiliary/aux_r.R",sep=""))


## load libararies
library(stringr)
library(sp)
library(spatstat)

## generate fake data

data.ppp <- multiprocppp(win=owin(),
                         lambdapar=150,
                         noff=2,mu0vec=c(3.5,2.75),hvec=c(0.03,0.02),
                         nxtra=1,lambdaxtra=95)

## Run it on the fake data

CPUtime <- system.time(M_out <- NS_MCMC_new_fixp(data.ppp,c("parent"),c("offspring1","offspring2"),bdtype=2,jitter=FALSE,B=100,hclimp=0.025,iters=100000,burn=80000,thin=2,store_res=FALSE,outfile=NULL))


