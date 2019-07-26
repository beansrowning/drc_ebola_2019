require(pomp)
require(plyr)
require(reshape2)
require(magrittr)
options(stringsAsFactors=FALSE)

require(foreach)
require(doMPI)
require(iterators)

source("ebola.R")
noexport <- c("ebolaModel")

cl <- startMPIcluster()
registerDoMPI(cl)

bake <- function (file, expr) {
  if (file.exists(file)) {
    readRDS(file)
    } else {
      val <- eval(expr)
      saveRDS(val,file=file)
      val
      }
  }

## trajectory matching simulation study

tic <- Sys.time()

po <- ebolaModel(country="Guinea",timestep=0.1)
params <- parmat(coef(po),3)
params["k",] <- c(0,0.2,0.5)
paramnames <- names(coef(po))

nsims <- 500

bake(file="sims.rds",{
  simulate(po,params=params,nsim=nsims,seed=208335746L,
           as.data.frame=TRUE,obs=TRUE) %>%
    rename(c(time="week")) %>%
    mutate(k=params["k",((as.integer(sim)-1)%%ncol(params))+1])
  }) -> simdat

pompUnload(po)

bake(file="tm-sim-profiles-R0.rds",{
  foreach(simul=1:nsims,
          .combine=rbind,.inorder=FALSE,
          .noexport=noexport,
          .options.mpi=list(chunkSize=10,seed=1598260027L,info=TRUE)) %:%
    foreach(type=c("raw","cum"),.combine=rbind,.inorder=FALSE) %dopar%
    {
      dat <- subset(simdat,sim==simul,select=c(week,cases,deaths))
      tm <- ebolaModel(country="Guinea",data=dat,type=as.character(type))
      st <- params[,(simul-1)%%3+1]
      true.k <- unname(st["k"])
      true.R0 <- unname(st["R0"])
      true.rho <- unname(st["rho"])
      st["k"] <- st["k"]+1e-6
      tm <- traj.match(tm,start=st,est=c("R0","k","rho"),transform=TRUE)
      if (coef(tm,"rho")==1) coef(tm,"rho") <- 0.999
      if (coef(tm,"rho")==0) coef(tm,"rho") <- 0.001
      tm <- traj.match(tm,est=c("R0","k","rho"),method='subplex',transform=TRUE)
      pompUnload(tm)
      data.frame(sim=simul,type=as.character(type),
                 true.k=true.k,true.R0=true.R0,true.rho=true.rho,
                 as.list(coef(tm)),loglik=logLik(tm),
                 conv=tm$convergence)
      } -> fits

  foreach (fit=iter(fits,by="row"),
           .noexport=noexport,
           .combine=rbind,.inorder=FALSE) %:%
    foreach (r0=seq(from=0.7,to=3,length=200),
             .combine=rbind,.inorder=FALSE,
             .options.mpi=list(chunkSize=200,seed=1598260027L,info=TRUE)) %dopar%
    {
      dat <- subset(simdat,sim==fit$sim,select=c(week,cases,deaths))
      tm <- ebolaModel(country="Guinea",data=dat,type=as.character(fit$type))
      coef(tm) <- unlist(fit[paramnames])
      coef(tm,"R0") <- r0
      if (coef(tm,"rho")==1) coef(tm,"rho") <- 0.999
      if (coef(tm,"rho")==0) coef(tm,"rho") <- 0.001
      tm <- traj.match(tm,est=c("k","rho"),transform=TRUE)
      if (coef(tm,"rho")==1) coef(tm,"rho") <- 0.999
      if (coef(tm,"rho")==0) coef(tm,"rho") <- 0.001
      tm <- traj.match(tm,est=c("k","rho"),transform=TRUE,method='subplex')
      pompUnload(tm)
      data.frame(sim=fit$sim,type=fit$type,
                 true.k=fit$true.k,true.R0=fit$true.R0,true.rho=fit$true.rho,
                 as.list(coef(tm)),loglik=logLik(tm),
                 conv=tm$convergence)
      }
  }) -> profiles

bake(file="tm-sim-fits.rds",{

  ddply(profiles,~type+sim+true.k,subset,loglik==max(loglik)) -> starts

  foreach(fit=iter(starts,by="row"),.combine=rbind,.inorder=FALSE,
          .noexport=noexport,
          .options.mpi=list(chunkSize=30,seed=1598260027L,info=TRUE)) %dopar%
    {
      dat <- subset(simdat,sim==fit$sim,select=c(week,cases,deaths))
      tm <- ebolaModel(country="Guinea",data=dat,type=as.character(fit$type))
      coef(tm) <- unlist(fit[paramnames])
      if (coef(tm,"rho")==1) coef(tm,"rho") <- 0.999
      if (coef(tm,"rho")==0) coef(tm,"rho") <- 0.001
      tm <- traj.match(tm,est=c("k","rho"),transform=TRUE)
      if (coef(tm,"rho")==1) coef(tm,"rho") <- 0.999
      if (coef(tm,"rho")==0) coef(tm,"rho") <- 0.001
      tm <- traj.match(tm,est=c("k","rho"),transform=TRUE,method='subplex')
      pompUnload(tm)
      data.frame(sim=fit$sim,type=fit$type,
                 true.k=fit$true.k,true.R0=fit$true.R0,true.rho=fit$true.rho,
                 as.list(coef(tm)),loglik=logLik(tm),
                 conv=tm$convergence)
      }
  }) -> fits

toc <- Sys.time()
print(toc-tic)

## trajectory matching with least squares simulation study

tic <- Sys.time()

bake(file="ls-sim-profiles-R0.rds",{

  foreach(simul=1:nsims,
          .combine=rbind,.inorder=FALSE,
          .noexport=noexport,
          .options.mpi=list(chunkSize=10,seed=1598260027L,info=TRUE)) %:%
    foreach(type=c("raw","cum"),.combine=rbind,.inorder=FALSE) %dopar%
    {
      dat <- subset(simdat,sim==simul,select=c(week,cases,deaths))
      tm <- ebolaModel(country="Guinea",data=dat,type=as.character(type),
                       least.sq=TRUE)
      st <- params[,(simul-1)%%3+1]
      true.k <- unname(st["k"])
      true.R0 <- unname(st["R0"])
      true.rho <- unname(st["rho"])
      st["k"] <- 10
      tm <- traj.match(tm,start=st,est=c("R0","k","rho"),transform=TRUE)
      if (coef(tm,"rho")==1) coef(tm,"rho") <- 0.999
      if (coef(tm,"rho")==0) coef(tm,"rho") <- 0.001
      tm <- traj.match(tm,est=c("R0","k","rho"),method='subplex',transform=TRUE)
      pompUnload(tm)
      data.frame(sim=simul,type=as.character(type),
                 true.k=true.k,true.R0=true.R0,true.rho=true.rho,
                 as.list(coef(tm)),loglik=logLik(tm),
                 conv=tm$convergence)
      } -> fits

  foreach (fit=iter(fits,by="row"),
           .noexport=noexport,
           .combine=rbind,.inorder=FALSE) %:%
    foreach (r0=seq(from=0.7,to=3,length=200),
             .combine=rbind,.inorder=FALSE,
             .options.mpi=list(chunkSize=200,seed=1598260027L,info=TRUE)) %dopar%
    {
      dat <- subset(simdat,sim==fit$sim,select=c(week,cases,deaths))
      tm <- ebolaModel(country="Guinea",data=dat,type=as.character(fit$type),
                       least.sq=TRUE)
      coef(tm) <- unlist(fit[paramnames])
      coef(tm,"R0") <- r0
      if (coef(tm,"rho")==1) coef(tm,"rho") <- 0.999
      if (coef(tm,"rho")==0) coef(tm,"rho") <- 0.001
      tm <- traj.match(tm,est=c("k","rho"),transform=TRUE)
      if (coef(tm,"rho")==1) coef(tm,"rho") <- 0.999
      if (coef(tm,"rho")==0) coef(tm,"rho") <- 0.001
      tm <- traj.match(tm,est=c("k","rho"),transform=TRUE,method='subplex')
      pompUnload(tm)
      data.frame(sim=fit$sim,type=fit$type,
                 true.k=fit$true.k,true.R0=fit$true.R0,true.rho=fit$true.rho,
                 as.list(coef(tm)),loglik=logLik(tm),
                 conv=tm$convergence)
      }
  }) -> profiles

bake(file="ls-sim-fits.rds",{

  ddply(profiles,~type+sim+true.k,subset,loglik==max(loglik)) -> starts

  foreach(fit=iter(starts,by="row"),.combine=rbind,.inorder=FALSE,
          .noexport=noexport,
          .options.mpi=list(chunkSize=30,seed=1598260027L,info=TRUE)) %dopar%
    {
      dat <- subset(simdat,sim==fit$sim,select=c(week,cases,deaths))
      tm <- ebolaModel(country="Guinea",data=dat,type=as.character(fit$type),
                       least.sq=TRUE)
      coef(tm) <- unlist(fit[paramnames])
      if (coef(tm,"rho")==1) coef(tm,"rho") <- 0.999
      if (coef(tm,"rho")==0) coef(tm,"rho") <- 0.001
      tm <- traj.match(tm,est=c("k","rho"),transform=TRUE)
      if (coef(tm,"rho")==1) coef(tm,"rho") <- 0.999
      if (coef(tm,"rho")==0) coef(tm,"rho") <- 0.001
      tm <- traj.match(tm,est=c("k","rho"),transform=TRUE,method='subplex')
      pompUnload(tm)
      data.frame(sim=fit$sim,type=fit$type,
                 true.k=fit$true.k,true.R0=fit$true.R0,true.rho=fit$true.rho,
                 as.list(coef(tm)),loglik=logLik(tm),
                 conv=tm$convergence)
      }
  }) -> fits

toc <- Sys.time()
print(toc-tic)

closeCluster(cl)
mpi.quit()
