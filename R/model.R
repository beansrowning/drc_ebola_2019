library(pomp)
library(dplyr)
library(readr)
library(tidyr)

drc_ebola_data <- read.csv("data/health_zone_counts.csv", stringsAsFactors = FALSE) %>%
  mutate(
    week = epiweek(report_date),
    year = epiyear(report_date)
  ) %>%
  group_by(country, week, year) %>%
  mutate(
    date = min(report_date)
  ) %>%
  group_by(country, week, date, year) %>%
  summarize(
    cases = sum(total_cases_change, na.rm = TRUE)
  ) %>%
  ungroup()

drc_ebola_data_2018 <- drc_ebola_data %>%
  filter(year == 2018) %>%
  mutate(
    week = week - min(week) + 1
  )

drc_ebola_data <- drc_ebola_data %>%
  filter(year == 2019) %>%
  mutate(
    week = week + nrow(drc_ebola_data_2018)
  ) %>%
  bind_rows(drc_ebola_data_2018) %>%
  select(-year) %>%
  arrange(week) %>%
  mutate(country = "DRC") %>%
  mutate(deaths = NA)


# Population from 2018 World Bank estimates (https://data.worldbank.org/indicator/SP.POP.TOTL?locations=CD)
populations <- c(DRC = 84068091)

## Parameter transformations

paruntrans <- Csnippet('
  double *IC = &S_0;
  double *TIC = &TS_0;
  TR0 = log(R0);
  Trho = logit(rho);
  Tk = log(k);
  to_log_barycentric(TIC,IC,4);
')

partrans <- Csnippet('
  double *IC = &S_0;
  double *TIC = &TS_0;
  TR0 = exp(R0);
  Trho = expit(rho);
  Tk = exp(k);
  from_log_barycentric(TIC,IC,4);
')

##  Measurement model: hierarchical model for cases
## p(C_t | H_t): Negative binomial with mean rho*H_t and variance rho*H_t*(1+k*rho*H_t)
dObs <- Csnippet('
  double f;
  if (k > 0.0)
    f = dnbinom_mu(nearbyint(cases),1.0/k,rho*N_EI,1);
  else
    f = dpois(nearbyint(cases),rho*N_EI,1);
  lik = (give_log) ? f : exp(f);
')

rObs <- Csnippet('
  if (k > 0) {
    cases = rnbinom_mu(1.0/k,rho*N_EI);
    deaths = rnbinom_mu(1.0/k,rho*cfr*N_IR);
  } else {
    cases = rpois(rho*N_EI);
    deaths = rpois(rho*cfr*N_IR);
  }')

### measurement model for ordinary least-squares
dObsLS <- Csnippet('
  double f;
  f = dnorm(cases,rho*N_EI,k,1);
  lik = (give_log) ? f : exp(f);
')

rObsLS <- Csnippet('
  cases = rnorm(rho*N_EI,k);
  deaths = NA_REAL;
')

## Process model simulator
rSim <- Csnippet('
  double lambda, beta;
  double *E = &E1;
  beta = R0 * gamma; // Transmission rate
  lambda = beta * I / N; // Force of infection
  int i;

  // Transitions
  // From class S
  double transS = rbinom(S, 1.0 - exp(- lambda * dt)); // No of infections
  // From class E
  double transE[nstageE]; // No of transitions between classes E
  for(i = 0; i < nstageE; i++){
    transE[i] = rbinom(E[i], 1.0 - exp(- nstageE * alpha * dt));
  }
  // From class I
  double transI = rbinom(I, 1.0 - exp(- gamma * dt)); // No of transitions I->R

  // Balance the equations
  S -= transS;
  E[0] += transS - transE[0];
  for(i=1; i < nstageE; i++) {
    E[i] += transE[i-1] - transE[i];
  }
  I += transE[nstageE - 1] - transI;
  R += transI;
  N_EI += transE[nstageE - 1]; // No of transitions from E to I
  N_IR += transI; // No of transitions from I to R
')

## Deterministic skeleton (an ODE), used in trajectory matching
skel <- Csnippet('
  double lambda, beta;
  double *E = &E1;
  double *DE = &DE1;
  beta = R0 * gamma; // Transmission rate
  lambda = beta * I / N; // Force of infection
  int i;

  // Balance the equations
  DS = - lambda * S;
  DE[0] = lambda * S - nstageE * alpha * E[0];
  for (i=1; i < nstageE; i++)
    DE[i] = nstageE * alpha * (E[i-1]-E[i]);
  DI = nstageE * alpha * E[nstageE-1] - gamma * I;
  DR = gamma * I;
  DN_EI = nstageE * alpha * E[nstageE-1];
  DN_IR = gamma * I;
')

ebolaModel <- function(
  country = c("DRC"),
  data = NULL,
  timestep = 0.01, nstageE = 3L,
  type = c("raw","cum"), na.rm = FALSE, least.sq = FALSE) {

  type <- match.arg(type)
  ctry <- match.arg(country)
  pop <- unname(populations[ctry])

  ## Incubation period is supposed to be Gamma distributed with shape parameter 3
  ## and mean 11.4 days.  The discrete-time formula is used to calculate the
  ## corresponding alpha (cf He et al., Interface 2010).
  ## Case-fatality ratio is fixed at 0.7 (cf WHO Ebola response team, NEJM 2014)
  incubation_period <- 11.4 / 7
  infectious_period <- 7 / 7
  index_case <- 10 / pop
  dt <- timestep
  nstageE <- as.integer(nstageE)

  globs <- paste0("static int nstageE = ", nstageE, ";");

  theta <- c(
    N = pop,
    R0 = 1.4,
    alpha = -1 / (nstageE * dt) * log(1 - nstageE * dt / incubation_period),
    gamma = -log(1 - dt / infectious_period) / dt,
    rho = 0.2, 
    cfr = 0.7,
    k = 0,
    S_0 = 1 - index_case,
    E_0 = index_case / 2 - 5e-9,
    I_0 = index_case / 2 - 5e-9,
    R_0 = 1e-8
  )

  if (is.null(data)) {
    if (ctry == "WestAfrica") {
      dat <- dat %>%
        group_by(week) %>%
        summarize(
          cases = sum(cases, na.rm = TRUE),
          deaths = sum(deaths, na.rm = TRUE)
        )
      } else {
        dat <- dat %>%
          filter(county == ctry) %>%
          select(-county)
        }
    } else {
      dat <- data
      }

  if (na.rm) {
    dat <- dat %>%
      filter(!is.na(cases)) %>%
      mutate(week = week - min(week) + 1)
    }

  if (type=="cum") {
    dat <- dat %>%
      mutate(
        cases = cumsum(cases),
        deaths = cumsum(deaths)
      )
    }

  ## Create the pomp object
  pomp(
    data = dat[c("week", "cases", "deaths")],
    times = "week",
    t0 = 0,
    params = theta,
    globals = globs,
    obsnames = c("cases", "deaths"),
    statenames = c("S", "E1", "I", "R", "N_EI", "N_IR"),
    zeronames = if (type=="raw") c("N_EI", "N_IR") else character(0),
    paramnames=c("N","R0","alpha","gamma","rho","k","cfr",
                 "S_0","E_0","I_0","R_0"),
    nstageE=nstageE,
    dmeasure=if (least.sq) dObsLS else dObs,
    rmeasure=if (least.sq) rObsLS else rObs,
    rprocess=discrete.time.sim(step.fun=rSim,delta.t=timestep),
    skeleton=skel,
    skeleton.type="vectorfield",
    parameter.transform=partrans,
    parameter.inv.transform=paruntrans,
    initializer = function (params, t0, nstageE, ...) {
      all.state.names <- c("S", paste0("E",1:nstageE),"I","R","N_EI","N_IR")
      comp.names <- c("S",paste0("E",1:nstageE),"I","R")
      x0 <- setNames(numeric(length(all.state.names)),all.state.names)
      frac <- c(params["S_0"],rep(params["E_0"]/nstageE,nstageE),params["I_0"],params["R_0"])
      x0[comp.names] <- round(params["N"]*frac/sum(frac))
      x0
      }
    ) -> po
  }