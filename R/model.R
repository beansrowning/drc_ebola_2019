# ==================================================================== #
# SEIR Ebola model (from King, et al. 2015)                            #
# Sean Browning                                                        #
# ==================================================================== #
library(pomp)
library(dplyr)
library(lubridate)

# === Define folder locations =============================================
tmp_dir <- normalizePath("tmp", winslash = "/", mustWork = FALSE)
data_dir <- normalizePath("data", winslash = "/", mustWork = TRUE)

# Make sure we have a tmp_dir
# and that we don't track any files
if (!dir.exists(tmp_dir)) {
  dir.create(tmp_dir)
  cat("*.*", file = file.path(tmp_dir, ".gitignore"))
}
# === Load Ebola counts ===================================================
# NOTE:
# - This pulls from a live CSV that is updated by update_ebola_counts.R
# - incident counts are back-computed from cumulative counts, so totals will not
#   match up.
# - A better fit might be achieved with taking cumulative counts as-is, or
#   refining calc to use only confirmed cases, etc.
drc_ebola_data <- read.csv(file.path(data_dir, "drc_model_counts.csv"), stringsAsFactors = FALSE) %>%
  select(
    country, week, date,
    cases = new_cases, deaths = new_deaths
  ) %>%
  slice(-1) %>%
  mutate(
    date = as_date(date),
    week = week - 1
  ) %>%
  as_tibble()

# Population from 2018 World Bank estimates (https://data.worldbank.org/indicator/SP.POP.TOTL?locations=CD)
population <- c(DRC = 84068091)

# ==== Measurement model ==================================================
# Hierarchical model for cases
# p(C_t | H_t): Negative binomial with mean rho*H_t and variance rho*H_t*(1+k*rho*H_t)
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
                       type = c("raw", "cum"), na.rm = FALSE, least.sq = FALSE) {
  type <- match.arg(type)
  ctry <- match.arg(country)
  pop <- unname(population[ctry])

  ## Incubation period is supposed to be Gamma distributed with shape parameter 3
  ## and mean 11.4 days.  The discrete-time formula is used to calculate the
  ## corresponding alpha (cf He et al., Interface 2010).
  ## Case-fatality ratio is fixed at 0.7 (cf WHO Ebola response team, NEJM 2014)
  incubation_period <- 11.4 / 7
  infectious_period <- 7 / 7
  index_case <- 10 / pop
  dt <- timestep
  nstageE <- as.integer(nstageE)

  globs <- paste0("static int nstageE = ", nstageE, ";")

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

  if (na.rm) {
    data <- data %>%
      filter(!is.na(cases)) %>%
      mutate(week = week - min(week) + 1)
  }

  if (type == "cum") {
    data <- data %>%
      mutate(
        cases = cumsum(cases),
        deaths = cumsum(deaths)
      )
  }

  ## Create the pomp object
  pomp(
    data = data[c("week", "cases", "deaths")],
    times = "week",
    t0 = 0,
    params = theta,
    globals = globs,
    cdir = tmp_dir,
    obsnames = c("cases", "deaths"),
    statenames = c("S", "E1", "I", "R", "N_EI", "N_IR"),
    accumvars = if (type == "raw") c("N_EI", "N_IR") else character(0),
    paramnames = c(
      "N", "R0", "alpha", "gamma", "rho", "k", "cfr",
      "S_0", "E_0", "I_0", "R_0"
    ),
    nstageE = nstageE,
    dmeasure = if (least.sq) dObsLS else dObs,
    rmeasure = if (least.sq) rObsLS else rObs,
    rprocess = discrete_time(step.fun = rSim, delta.t = timestep),
    skeleton = vectorfield(skel),
    partrans = parameter_trans(
      log = c("R0", "k"),
      logit = "rho",
      barycentric = c("S_0", "E_0", "I_0", "R_0")
    ),
    rinit = function(S_0, E_0, I_0, R_0, N, t0, nstageE, ...) {
      all.state.names <- c("S", paste0("E", 1:nstageE), "I", "R", "N_EI", "N_IR")
      comp.names <- c("S", paste0("E", 1:nstageE), "I", "R")
      x0 <- setNames(numeric(length(all.state.names)), all.state.names)
      frac <- c(S_0, rep(E_0 / nstageE, nstageE), I_0, R_0)
      x0[comp.names] <- round(N * frac / sum(frac))
      x0
    }
  )
}