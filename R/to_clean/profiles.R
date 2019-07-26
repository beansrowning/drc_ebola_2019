require(pomp)
require(plyr)
require(reshape2)
require(magrittr)
options(stringsAsFactors=FALSE)

require(foreach)
require(doMPI)
require(iterators)

source("ebola.R")

foreach (country=c("SierraLeone","Liberia","Guinea","WestAfrica"),
         .inorder=TRUE,.combine=c) %:%
  foreach (type=c("raw","cum"),.inorder=TRUE,.combine=c) %do%
  {
    ebolaModel(country=country,type=type,na.rm=TRUE)
    } -> models
dim(models) <- c(2,4)
dimnames(models) <- list(c("raw","cum"),
                         c("SierraLeone","Liberia","Guinea","WestAfrica"))

noexport <- c("models")

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

## trajectory matching: R0 profile

bake(file="tm-fits-R0.rds",{

  starts <- profileDesign(R0=seq(1,3,length=200),
                          upper=c(k=1),
                          lower=c(k=1e-8),
                          nprof=40)

  foreach (start=iter(starts,by='row'),
           .combine=rbind,.inorder=FALSE,
           .noexport=noexport,
           .options.mpi=list(chunkSize=100,seed=2016138277L,info=TRUE)
           ) %:%
    foreach (type=c("raw","cum"),.combine=rbind,.inorder=FALSE) %:%
    foreach (country=c("SierraLeone","Liberia","Guinea","WestAfrica"),
             .combine=rbind,.inorder=FALSE) %dopar%
    {
      tm <- models[type,country][[1]]
      tic <- Sys.time()
      coef(tm,names(start)) <- unname(unlist(start))
      coef(tm,"rho") <- 0.2
      tm <- traj.match(tm,est=c("k","E_0","I_0"),transform=TRUE)
      if (coef(tm,"k")==0) coef(tm,"k") <- 1e-8
      if (coef(tm,"E_0")==0) coef(tm,"E_0") <- 1e-12
      if (coef(tm,"I_0")==0) coef(tm,"I_0") <- 1e-12
      tm <- traj.match(tm,method='subplex',control=list(maxit=1e5))
      toc <- Sys.time()
      etime <- toc-tic
      units(etime) <- "hours"
      data.frame(country=country,type=type,as.list(coef(tm)),
                 loglik=logLik(tm),conv=tm$convergence,
                 etime=as.numeric(etime))
      } %>% mutate(sum=S_0+E_0+I_0+R_0,
                   S_0=round(N*S_0/sum),
                   E_0=round(N*E_0/sum),
                   I_0=round(N*I_0/sum),
                   R_0=round(N*R_0/sum)) %>%
    subset(conv %in% c(0,1),select=-sum) %>%
    unique()

  }) -> profR0

## trajectory matching: k profile

bake(file="tm-fits-k.rds",{

  starts <- profileDesign(k=seq(0,1,length=100),
                          upper=c(R0=1),
                          lower=c(R0=3),
                          nprof=40)

  foreach (start=iter(starts,by='row'),
           .combine=rbind,.inorder=FALSE,
           .noexport=noexport,
           .options.mpi=list(chunkSize=100,seed=2016138277L,info=TRUE)
           ) %:%
    foreach (type=c("raw","cum"),.combine=rbind,.inorder=FALSE) %:%
    foreach (country=c("SierraLeone","Liberia","Guinea","WestAfrica"),
             .combine=rbind,.inorder=FALSE) %dopar%
    {
      tm <- models[type,country][[1]]
      tic <- Sys.time()
      coef(tm,names(start)) <- unname(unlist(start))
      coef(tm,"rho") <- 0.2
      tm <- traj.match(tm,est=c("R0","E_0","I_0"),transform=TRUE)
      if (coef(tm,"E_0")==0) coef(tm,"E_0") <- 1e-12
      if (coef(tm,"I_0")==0) coef(tm,"I_0") <- 1e-12
      tm <- traj.match(tm,method='subplex',control=list(maxit=1e5))
      toc <- Sys.time()
      etime <- toc-tic
      units(etime) <- "hours"
      data.frame(country=country,type=type,as.list(coef(tm)),
                 loglik=logLik(tm),conv=tm$convergence,
                 etime=as.numeric(etime))
      } %>% mutate(sum=S_0+E_0+I_0+R_0,
                   S_0=round(N*S_0/sum),
                   E_0=round(N*E_0/sum),
                   I_0=round(N*I_0/sum),
                   R_0=round(N*R_0/sum)) %>%
    subset(conv %in% c(0,1),select=-sum) %>%
    unique()

  }) -> profk

## All trajectory matching computations
ldply(list(R0=profR0,k=profk),.id='profile') -> profTM

## Iterated filtering, R0 profile

bake(file="if-fits-R0_a.rds",{

  profTM %>% subset(profile=="R0") %>%
    ddply(~country+type,subset,
          is.finite(loglik)&loglik>max(loglik)-20) %>%
    ddply(~country+type+R0,subset,
          loglik==max(loglik),
          select=-c(loglik,etime,conv,profile)) -> pars

  foreach (start=iter(pars,by='row'),
           .combine=rbind,.inorder=FALSE,
           .options.mpi=list(chunkSize=1,seed=1264624821L),
           .noexport=noexport) %dopar%
    {
      tic <- Sys.time()

      country <- as.character(start$country)
      type <- as.character(start$type)
      st <- unlist(subset(start,select=-c(country,type)))

      po <- models[type,country][[1]]
      coef(po,names(st)) <- st
      if (coef(po,"E_0")==0) coef(po,"E_0") <- 1e-5
      if (coef(po,"I_0")==0) coef(po,"I_0") <- 1e-5

      mf <- mif(po, Nmif=10,
                rw.sd = c(k=0.02,E_0=1,I_0=1),
                ivps = c("E_0","I_0"),
                Np = 2000,
                var.factor = 2,
                method = "mif2",
                cooling.type = "hyperbolic",
                cooling.fraction = 0.5,
                transform = TRUE,
                verbose = FALSE)
      mf <- continue(mf, Nmif = 50, cooling.fraction = 0.1)

      ## Runs 10 particle filters to assess Monte Carlo error in likelihood
      pf <- replicate(10,pfilter(mf,Np=5000,max.fail=Inf))
      ll <- sapply(pf,logLik)
      ll <- logmeanexp(ll, se = TRUE)
      nfail <- sapply(pf,getElement,"nfail")
      toc <- Sys.time()
      etime <- toc-tic
      units(etime) <- "hours"

      data.frame(country=country,type=type,as.list(coef(mf)),
                 loglik = ll[1],
                 loglik.se = ll[2],
                 nfail.min = min(nfail),
                 nfail.max = max(nfail),
                 etime = as.numeric(etime))
      }
  }) -> profR0

## Filter once more on maxima

bake(file="if-fits-R0.rds",{

  profR0 %>% subset(is.finite(loglik)&nfail.max==0) %>%
    ddply(~country+type+R0,subset,rank(-loglik)<=5) %>%
    subset(select=-c(loglik,loglik.se,nfail.max,nfail.min,etime)) -> pars

  foreach (start=iter(pars,by='row'),
           .combine=rbind,.inorder=FALSE,
           .options.mpi=list(chunkSize=1,seed=1264624821L),
           .noexport=noexport) %dopar%
    {
      tic <- Sys.time()

      country <- as.character(start$country)
      type <- as.character(start$type)
      st <- unlist(subset(start,select=-c(country,type)))

      po <- models[type,country][[1]]
      coef(po,names(st)) <- unname(st)

      ## Runs 10 particle filters to assess Monte Carlo error in likelihood
      pf <- try(replicate(10,pfilter(po,Np=5000,max.fail=Inf)))

      toc <- Sys.time()
      etime <- toc-tic
      units(etime) <- "hours"

      ll <- sapply(pf,logLik)
      ll <- logmeanexp(ll, se = TRUE)
      nfail <- sapply(pf,getElement,"nfail")

      data.frame(country=country,type=type,as.list(coef(po)),
                 loglik = ll[1],
                 loglik.se = ll[2],
                 nfail.min = min(nfail),
                 nfail.max = max(nfail),
                 etime = as.numeric(etime))
      }
  }) -> profR0

## Iterated filtering, k profile

bake(file="if-fits-k_a.rds",{

  profTM %>% subset(profile=="k") %>%
    ddply(~country+type,subset,
          is.finite(loglik)&loglik>max(loglik)-20) %>%
    ddply(~country+type+k,subset,
          loglik==max(loglik),
          select=-c(loglik,etime,conv,profile)) -> pars

  foreach (start=iter(pars,by='row'),
           .combine=rbind,.inorder=FALSE,
           .options.mpi=list(chunkSize=1,seed=1264624821L),
           .noexport=noexport) %dopar%
    {
      tic <- Sys.time()

      country <- as.character(start$country)
      type <- as.character(start$type)
      st <- unlist(subset(start,select=-c(country,type)))

      po <- models[type,country][[1]]
      coef(po,names(st)) <- st
      if (coef(po,"E_0")==0) coef(po,"E_0") <- 1e-5
      if (coef(po,"I_0")==0) coef(po,"I_0") <- 1e-5

      mf <- mif(po, Nmif=10,
                rw.sd = c(R0=0.02,E_0=1,I_0=1),
                ivps = c("E_0","I_0"),
                Np = 2000,
                var.factor = 2,
                method = "mif2",
                cooling.type = "hyperbolic",
                cooling.fraction = 0.5,
                transform = TRUE,
                verbose = FALSE)
      mf <- continue(mf, Nmif = 50, cooling.fraction = 0.1)

      ## Runs 10 particle filters to assess Monte Carlo error in likelihood
      pf <- replicate(10,pfilter(mf,Np=5000,max.fail=Inf))
      ll <- sapply(pf,logLik)
      ll <- logmeanexp(ll, se = TRUE)
      nfail <- sapply(pf,getElement,"nfail")
      toc <- Sys.time()
      etime <- toc-tic
      units(etime) <- "hours"

      data.frame(country=country,type=type,as.list(coef(mf)),
                 loglik = ll[1],
                 loglik.se = ll[2],
                 nfail.min = min(nfail),
                 nfail.max = max(nfail),
                 etime = as.numeric(etime))
      }
  }) -> profk

## Filter once more on maxima

bake(file="if-fits-k.rds",{

  profk %>% subset(is.finite(loglik)&nfail.max==0) %>%
    ddply(~country+type+R0,subset,rank(-loglik)<=5) %>%
    subset(select=-c(loglik,loglik.se,nfail.max,nfail.min,etime)) -> pars

  foreach (start=iter(pars,by='row'),
           .combine=rbind,.inorder=FALSE,
           .options.mpi=list(chunkSize=1,seed=1264624821L),
           .noexport=noexport) %dopar%
    {
      tic <- Sys.time()

      country <- as.character(start$country)
      type <- as.character(start$type)
      st <- unlist(subset(start,select=-c(country,type)))

      po <- models[type,country][[1]]
      coef(po,names(st)) <- unname(st)

      ## Runs 10 particle filters to assess Monte Carlo error in likelihood
      pf <- try(replicate(10,pfilter(po,Np=5000,max.fail=Inf)))

      toc <- Sys.time()
      etime <- toc-tic
      units(etime) <- "hours"

      ll <- sapply(pf,logLik)
      ll <- logmeanexp(ll, se = TRUE)
      nfail <- sapply(pf,getElement,"nfail")

      data.frame(country=country,type=type,as.list(coef(po)),
                 loglik = ll[1],
                 loglik.se = ll[2],
                 nfail.min = min(nfail),
                 nfail.max = max(nfail),
                 etime = as.numeric(etime))
      }
  }) -> profk

ldply(list(R0=profR0,k=profk),.id='profile') -> profIF

ldply(list(det=profTM,stoch=profIF),.id='model') %>%
  saveRDS(file='profiles.rds')

closeCluster(cl)
mpi.quit()
