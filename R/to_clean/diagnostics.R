library(pomp)
library(plyr)
library(reshape2)
library(magrittr)
options(stringsAsFactors=FALSE)

require(foreach)
require(doMC)
require(iterators)

source("ebola.R")

readRDS("profiles.rds") %>%
  ddply(~country+type+model,subset,loglik==max(loglik)) %>%
  subset(type=="raw"&model=="stoch") -> mles

time1 <- c(Guinea="2014-01-05",Liberia="2014-06-01",
           SierraLeone="2014-06-08",WestAfrica="2014-01-05")

registerDoMC(4)

foreach (mle=iter(mles,by='row'),.combine=rbind) %dopar%
  {
    country=as.character(mle$country)
    type=as.character(mle$type)
    M <- ebolaModel(country=country,type=type,
                    na.rm=TRUE,nstage=3,timestep=0.01)
    p <- unlist(subset(mle,select=-c(country,type,model,profile,
                                     loglik,loglik.se,
                                     nfail.min,nfail.max,conv,etime)))
    coef(M,names(p)) <- unname(unlist(p))

    t0 <- as.Date(time1[country])

    simulate(M,nsim=10,as.data.frame=TRUE,obs=TRUE,include.data=TRUE,seed=2186L) %>%
      mutate(date=t0+time*7,country=country,type=type)
    } %>% saveRDS(file="diagnostics-sim.rds")

foreach (mle=iter(mles,by='row'),.combine=rbind) %dopar%
  {
    country=as.character(mle$country)
    type=as.character(mle$type)
    M <- ebolaModel(country=country,type=type,
                    na.rm=TRUE,nstage=3,timestep=0.01)
    p <- unlist(subset(mle,select=-c(country,type,model,profile,
                                     loglik,loglik.se,
                                     nfail.min,nfail.max,conv,etime)))
    coef(M,names(p)) <- unname(unlist(p))

    probe(M,probes=list(probe.acf(var="cases",lags=1,type="correlation")),
          nsim=500,seed=1878812716L) %>%
      as.data.frame() -> pb
    pb %>% mutate(sim=rownames(pb),
                  data=ifelse(sim=="data","data","simulation"),
                  type=type,
                  country=country)
    } %>% saveRDS(file="diagnostics-probes.rds")

## Additional diagnostics

## Run probes for each country
## Custom probe: exponential growth rate
probe.trend <- function (y) {
  cases <- y["cases",]
  df <- data.frame(week=seq_along(cases),cases=cases)
  fit <- lm(log1p(cases)~week,data=df)
  unname(coef(fit)[2])
  }

foreach (mle=iter(mles,by='row'),.combine=rbind) %dopar%
  {
    country=as.character(mle$country)
    type=as.character(mle$type)
    M <- ebolaModel(country=country,type=type,
                    na.rm=TRUE,nstage=3,timestep=0.01)
    p <- unlist(subset(mle,select=-c(country,type,model,profile,
                                     loglik,loglik.se,
                                     nfail.min,nfail.max,conv,etime)))
    coef(M,names(p)) <- unname(unlist(p))

    ## remove an exponential trend, give residuals on the log scale
    dm <- model.matrix(lm(log1p(cases)~time,data=as.data.frame(M)))
    rm <- diag(nrow(dm))-dm%*%solve(crossprod(dm))%*%t(dm)
    detrend <- function (x) log1p(x)%*%rm

    probe(M,probes=list(
      probe.acf(var="cases",lags=1,type="correlation"),
      sd=probe.sd(var="cases",transform=log1p),
      probe.quantile(var="cases",prob=c(0.9)),
      d=probe.acf(var="cases",lags=c(1,2,3),type="correlation",
                  transform=detrend),
      trend=probe.trend),
      nsim=2000,seed=2186L
      ) %>% as.data.frame() -> pb
    pb %>% mutate(sim=rownames(pb),
                  kind=ifelse(sim=="data","data","simulation"),
                  type=type,
                  country=country)
    } %>%
  melt(id=c("country","type","kind","sim"),variable.name="probe") %>%
  arrange(country,type,probe,kind,sim) %>%
  saveRDS(file="diagnostics-addl-probes.rds")
