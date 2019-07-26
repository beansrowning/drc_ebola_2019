require(pomp)
require(plyr)
require(reshape2)
require(magrittr)
options(stringsAsFactors=FALSE)

set.seed(988077383L)

require(foreach)
require(doMPI)
require(iterators)

source("ebola.R")

horizon <- 13

foreach (country=c("SierraLeone"),.inorder=TRUE,.combine=c) %:%
  foreach (type=c("raw","cum"),.inorder=TRUE,.combine=c) %do%
  {
    M1 <- ebolaModel(country=country,type=type,
                     timestep=0.01,nstageE=3,na.rm=TRUE)
    M2 <- ebolaModel(country=country,type="raw",
                     timestep=0.01,nstageE=3,na.rm=TRUE)
    time(M2) <- seq(from=1,to=max(time(M1))+horizon,by=1)
    M3 <- ebolaModel(country=country,type="raw",
                     timestep=0.01,nstageE=3,na.rm=TRUE)
    time(M3) <- seq(from=max(time(M1))+1,to=max(time(M1))+horizon,by=1)
    timezero(M3) <- max(time(M1))
    list(M1,M2,M3)
    } -> models
dim(models) <- c(3,2,1)
dimnames(models) <- list(c("fit","det.forecast","stoch.forecast"),
                         c("raw","cum"),c("SierraLeone"))

noexport <- c("models")

## Weighted quantile function
wquant <- function (x, weights, probs = c(0.025,0.5,0.975)) {
  idx <- order(x)
  x <- x[idx]
  weights <- weights[idx]
  w <- cumsum(weights)/sum(weights)
  rval <- approx(w,x,probs,rule=1)
  rval$y
  }

starts <- c(Guinea="2014-01-05",Liberia="2014-06-01",SierraLeone="2014-06-08")

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

readRDS("profiles.rds") %>%
  ddply(~country+type+model,subset,loglik>max(loglik)-6,
        select=-c(conv,etime,loglik.se,nfail.min,nfail.max,profile)) -> mles

mles %>% melt(id=c("country","type","model"),variable.name='parameter') %>%
  ddply(~country+type+model+parameter,summarize,
        min=min(value),max=max(value)) %>%
  subset(parameter!="loglik") %>%
  melt(measure=c("min","max")) %>%
  acast(country~type~model~parameter~variable) -> ranges

mles %>% ddply(~country+type+model,subset,loglik==max(loglik),select=-loglik) %>%
  mutate(k=round(k,4),rho=round(rho,4),R0=round(R0,4),E_0=3*round(E_0/3)) %>%
  unique() %>%
  arrange(country,type,model) -> mles

### DETERMINISTIC MODELS

bake(file="forecasts_det.rds",{
  foreach (country=c("SierraLeone"),
           .inorder=TRUE,.combine=rbind) %:%
    foreach (type=c("raw","cum"),nsamp=c(1000,3000),
             .inorder=TRUE,.combine=rbind) %do%
    {

      params <- sobolDesign(lower=ranges[country,type,'det',,'min'],
                            upper=ranges[country,type,'det',,'max'],
                            nseq=nsamp)

      foreach(p=iter(params,by='row'),
              .inorder=FALSE,
              .combine=rbind,
              .noexport=noexport,
              .options.multicore=list(set.seed=TRUE),
              .options.mpi=list(chunkSize=10,seed=1568335316L,info=TRUE)
              ) %dopar%
        {
          M1 <- models["fit",type,country][[1]]
          M2 <- models["det.forecast",type,country][[1]]
          ll <- logLik(traj.match(M1,start=unlist(p)))
          x <- trajectory(M2,params=unlist(p))
          p <- parmat(unlist(p),20)
          rmeasure(M2,x=x,times=time(M2),params=p) %>%
            melt() %>%
            mutate(time=time(M2)[time],
                   period=ifelse(time<=max(time(M1)),"calibration","projection"),
                   loglik=ll)
        } %>%
        subset(variable=="cases",select=-variable) %>%
        mutate(weight=exp(loglik-mean(loglik))) %>%
        arrange(time,rep) -> sims

      ess <- with(subset(sims,time==max(time)),weight/sum(weight))
      ess <- 1/sum(ess^2)
      cat("ESS det",country,type,"=",ess,"\n")

      sims %>%
        ddply(~time+period,summarize,prob=c(0.025,0.5,0.975),
              quantile=wquant(value,weights=weight,probs=prob)) %>%
        mutate(prob=mapvalues(prob,from=c(0.025,0.5,0.975),
                              to=c("lower","median","upper"))) %>%
        dcast(period+time~prob,value.var='quantile') %>%
        mutate(country=country,type=type)
      }
  }) -> fc_tm

### STOCHASTIC MODEL

bake(file="forecasts_stoch.rds",{
  foreach (country=c("SierraLeone"),
           .inorder=TRUE,.combine=rbind) %:%
    foreach (type=c("raw","cum"),nsamp=c(200,200),
             .inorder=TRUE,.combine=rbind) %do%
    {

      params <- sobolDesign(lower=ranges[country,type,'stoch',,'min'],
                            upper=ranges[country,type,'stoch',,'max'],
                            nseq=nsamp)

      foreach(p=iter(params,by='row'),
              .inorder=FALSE,
              .combine=rbind,
              .noexport=noexport,
              .options.multicore=list(set.seed=TRUE),
              .options.mpi=list(chunkSize=1,seed=1568335316L,info=TRUE)
              ) %dopar%
        {
          M1 <- models["fit",type,country][[1]]
          M2 <- models["stoch.forecast",type,country][[1]]
          pf <- pfilter(M1,params=unlist(p),Np=2000,save.states=TRUE)
          pf$saved.states %>% tail(1) %>% melt() %>%
            acast(variable~rep,value.var='value') %>%
            apply(2,function (x) {
              setNames(c(x["S"],sum(x[c("E1","E2","E3")]),x["I"],x["R"]),
                       c("S_0","E_0","I_0","R_0"))}) -> x
          pp <- parmat(unlist(p),ncol(x))
          pp[rownames(x),] <- x
          simulate(M2,params=pp,obs=TRUE) %>%
            melt() %>%
            mutate(time=time(M2)[time],
                   period=ifelse(time<=max(time(M1)),"calibration","projection"),
                   loglik=logLik(pf))
        } %>% subset(variable=="cases",select=-variable) %>%
        mutate(weight=exp(loglik-mean(loglik))) %>%
        arrange(time,rep) -> sims

      ess <- with(subset(sims,time==max(time)),weight/sum(weight))
      ess <- 1/sum(ess^2)
      cat("ESS stoch",country,type,"=",ess,"\n")

      sims %>% ddply(~time+period,summarize,prob=c(0.025,0.5,0.975),
                     quantile=wquant(value,weights=weight,probs=prob)) %>%
        mutate(prob=mapvalues(prob,from=c(0.025,0.5,0.975),
                              to=c("lower","median","upper"))) %>%
        dcast(period+time~prob,value.var='quantile') %>%
        mutate(country=country,type=type)
      }
  }) -> fc_if

ldply(list(stoch=fc_if,det=fc_tm),.id='model') %>%
  ddply(~country,mutate,
        model=factor(model,levels=c("stoch","det")),
        date=as.Date(starts[unique(as.character(country))])+7*(time-1)) %>%
  saveRDS(file='forecasts.rds')

closeCluster(cl)
mpi.quit()
