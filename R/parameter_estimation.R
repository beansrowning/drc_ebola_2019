# ================================================================================= #
# Parameter Estimation (adapted from King, et al. 2015)                             #
# (For desktop/SMP node, PSOCK/FORK parallelism)                                    #
# Sean Browning                                                                     #
# ================================================================================= #
library(pomp)
library(subplex)
library(dplyr)
library(tidyr)
library(magrittr)
library(foreach)
library(doParallel)
library(iterators)

# Load Ebola data and models
source("R/model.R")

# Run Preliminary model
models <- c(
  raw = ebolaModel(country = "DRC", data = drc_ebola_data, type = "raw", na.rm = TRUE),
  cum = ebolaModel(country = "DRC", data = drc_ebola_data, type = "cum", na.rm = TRUE)
)

# Start cluster
# BUG: POMP models do not seem to be thread-safe.
# Concurrent loading of shared libraries throws error
clust <- makeCluster(4, ifelse(.Platform$OS.type == "windows", "PSOCK", "FORK"))
registerDoParallel(clust)

# Define helper function to save results to file
bake <- function (file, expr) {
  if (file.exists(file)) {
    readRDS(file)
    } else {
      val <- eval(expr)
      saveRDS(val, file = file)
      val
      }
  }

# === Trajectory matching the R0 profile ========================================================

profR0 <- bake(file = "output/R0_fits.rds", {

  starts <- profileDesign(
    R0 = seq(1, 3, length = 200),
    upper = c(k = 1),
    lower = c(k = 1e-8),
    nprof = 40
  )

  out <- foreach(
    start = iter(starts, by = "row"),
    .combine = rbind,
    .inorder = FALSE
    ) %:%
    foreach(
      type = c("raw", "cum"),
      .combine = rbind,
      .inorder = FALSE,
      .packages = c("pomp", "subplex")
      ) %dopar% {

      tm <- models[[type]]

      tic <- Sys.time()
      
      coef(tm, names(start)) <- unname(unlist(start))

      coef(tm, "rho") <- 0.2

      tm <- traj_objfun(
        tm,
        est = c("k", "E_0", "I_0")
      )

      if (coef(tm, "k") == 0) {
        coef(tm, "k") <- 1e-8
      } 
      if (coef(tm, "E_0") == 0) {
        coef(tm, "E_0") <- 1e-12
      } 
      if (coef(tm, "I_0") == 0) {
        coef(tm,"I_0") <- 1e-12
      }

      tm_optim <- subplex(
        fn = tm,
        par = coef(tm, c("k", "E_0", "I_0")),
        control = list(maxit = 1e5)
      )

      toc <- Sys.time()

      # Compute system time
      etime <- toc - tic
      
      units(etime) <- "hours"
      
      data.frame(
        country = "DRC",
        type = type,
        as.list(coef(tm)),
        loglik = logLik(tm),
        conv = tm_optim$convergence,
        etime = as.numeric(etime)
      )
      } 
      
  # Final transformation
  out <- out %>% mutate(
    sum = S_0 + E_0 + I_0 + R_0,
    S_0 = round(N * S_0 / sum),
    E_0 = round(N * E_0 / sum),
    I_0 = round(N * I_0 / sum),
    R_0 = round(N * R_0 / sum)
  ) %>%
  filter(conv %in% 0:1) %>%
  select(-sum) %>%
  distinct()

  out
})

# === Trajectory Matching the k profile ===================================================
profk <- bake(file = "output/k_fits.rds", {

  starts <- profileDesign(
    k = seq(0, 1, length = 100),
    upper = c(R0 = 1),
    lower = c(R0 = 3),
    nprof = 40
  )

  out <- foreach(
    start = iter(starts, by = "row"),
    .combine = rbind,
    .inorder = FALSE
    ) %:%
    foreach(
      type = c("raw", "cum"),
      .combine = rbind,
      .inorder = FALSE) %dopar% {

      tm <- models[[type]][[1]]

      tic <- Sys.time()
      
      coef(tm, names(start)) <- unname(unlist(start))

      coef(tm, "rho") <- 0.2

      tm <- traj_objfun(
        tm,
        est = c("R0", "E_0", "I_0"),
        transform = TRUE
      )

      if (coef(tm, "E_0") == 0) {
        coef(tm, "E_0") <- 1e-12
      } 
      if (coef(tm, "I_0") == 0) {
        coef(tm,"I_0") <- 1e-12
      }

      tm_optim <- subplex(
        fn = tm,
        par = coef(tm, c("R0", "E_0", "I_0")),
        control = list(maxit = 1e5)
      )

      toc <- Sys.time()

      # Compute system time
      etime <- toc - tic
      
      units(etime) <- "hours"
      
      data.frame(
        country = "DRC",
        type = type,
        as.list(coef(tm)),
        loglik = logLik(tm),
        conv = tm_optim$convergence,
        etime = as.numeric(etime)
      )
      } 
      
  # Final transformation
  out <- out %>% mutate(
    sum = S_0 + E_0 + I_0 + R_0,
    S_0 = round(N * S_0 / sum),
    E_0 = round(N * E_0 / sum),
    I_0 = round(N * I_0 / sum),
    R_0 = round(N * R_0 / sum)
  ) %>%
  filter(conv %in% 0:1) %>%
  select(-sum) %>%
  distinct()

  out
})


# Combine all trajectory matching
profTM <- profR0 %>%
  mutate(profile = "RO") %>%
  bind_rows(
    mutate(profk, profile = "k")
  )

# === Iterated Filtering on the R0 profile ================================================= #
profR0_if <- bake(file = "output/R0_fits_if.rds", {
  pars <- profTM %>%
    filter(profile == "R0") %>%
    group_by(country, type) %>%
    filter(is.finite(loglik) & loglik > max(loglik) - 20) %>%
    group_by(R0, add = TRUE) %>%
    filter(loglik == max(loglik)) %>%
    ungroup() %>%
    select(-c(loglik, etime, conv, profile))


  foreach (
    start = iter(pars, by = "row"),
    .combine = rbind,
    .inorder = FALSE
  ) %dopar% {
      tic <- Sys.time()

      type <- as.character(start$type)

      st <- start %>%
        select(-country, -type) %>%
        unlist()

      po <- models[[type]]

      coef(po, names(st)) <- st

      if (coef(po, "E_0") == 0) {
        coef(po, "E_0") <- 1e-5
      }
      if (coef(po, "I_0") == 0) {
        coef(po, "I_0") <- 1e-5
      }

      mf <- mif2(
        po, 
        Nmif = 10,
        rw.sd = c(k = 0.02, E_0 = 1,I_0 = 1),
        Np = 2000,
        cooling.type = "hyperbolic",
        cooling.fraction.50 = 0.5,
        verbose = FALSE
      )

      mf <- continue(
        mf,
        Nmif = 50,
        cooling.fraction.50 = 0.1
      )

      ## Runs 10 particle filters to assess Monte Carlo error in likelihood
      pf <- replicate(
        10,
        pfilter(mf, Np = 5000, max.fail = Inf)
      
      )
      ll <- sapply(pf, logLik)
      ll <- logmeanexp(ll, se = TRUE)
      nfail <- sapply(pf, getElement, "nfail")

      toc <- Sys.time()

      etime <- toc - tic

      units(etime) <- "hours"

      data.frame(
        type = type,
        as.list(coef(mf)),
        loglik = ll[1],
        loglik.se = ll[2],
        nfail.min = min(nfail),
        nfail.max = max(nfail),
        etime = as.numeric(etime)
      )
      }
  })

## Filter once more on maxima

profR0_if <- bake(file = "output/R0_fits_if_maxima.rds",{

  pars <- profR0_if %>%
    filter(is.finite(loglik) & nfail.max == 0) %>%
    group_by(type, R0) %>%
    filter(rank(-loglik) <= 5) %>%
    ungroup() %>%
    select(-c(loglik, loglik.se, nfail.max, nfail.min, etime))

  foreach(
    start = iter(pars, by = "row"),
    .combine = rbind,
    .inorder = FALSE) %dopar% {
      tic <- Sys.time()

      type <- as.character(start$type)

      st <- start %>%
        select(-country, -type) %>%
        unlist()

      po <- models[[country]]

      coef(po, names(st)) <- unname(st)

      ## Runs 10 particle filters to assess Monte Carlo error in likelihood
      pf <- try(replicate(10, pfilter(po, Np = 5000, max.fail = Inf)))

      toc <- Sys.time()
      etime <- toc - tic
      units(etime) <- "hours"

      ll <- sapply(pf, logLik)
      ll <- logmeanexp(ll, se = TRUE)
      nfail <- sapply(pf, getElement, "nfail")

      data.frame(
        country = country,
        type = type,
        as.list(coef(po)),
        loglik = ll[1],
        loglik.se = ll[2],
        nfail.min = min(nfail),
        nfail.max = max(nfail),
        etime = as.numeric(etime)
      )
      }
  })

## Iterated filtering, k profile

profk_if <- bake(file = "output/k_fits_if.rds", {

  pars <- profTM %>%
    filter(profile == "k") %>%
    group_by(country, type) %>%
    filter(is.finite(loglik) & loglik > max(loglik) - 20) %>%
    group_by(k, add = TRUE) %>%
    filter(loglik == max(loglik)) %>%
    ungroup() %>%
    select(-c(loglik, etime, conv, profile))

  foreach (
    start = iter(pars, by = "row"),
    .combine = rbind,
    .inorder = FALSE
  ) %dopar% {
      tic <- Sys.time()

      type <- as.character(start$type)

      st <- start %>%
        select(-county, -type) %>%
        unlist()

      po <- models[[type]]

      coef(po, names(st)) <- st

      if (coef(po, "E_0") == 0) {
        coef(po, "E_0") <- 1e-5
      }
      if (coef(po, "I_0") == 0) {
        coef(po, "I_0") <- 1e-5
      }

      mf <- mif2(
        po,
        Nmif = 10,
        rw.sd = c(R0 = 0.02, E_0 = 1, I_0 = 1),
        ivps = c("E_0", "I_0"),
        Np = 2000,
        cooling.type = "hyperbolic",
        cooling.fraction.50 = 0.5,
        verbose = FALSE
      )

      mf <- continue(
        mf,
        Nmif = 50,
        cooling.fraction.50 = 0.1
      )

      ## Runs 10 particle filters to assess Monte Carlo error in likelihood
      pf <- replicate(10, pfilter(mf, Np = 5000, max.fail = Inf))
      ll <- sapply(pf, logLik)
      ll <- logmeanexp(ll, se = TRUE)
      nfail <- sapply(pf, getElement, "nfail")
      toc <- Sys.time()
      etime <- toc - tic
      units(etime) <- "hours"

      data.frame(
        country = country,
        type = type,
        as.list(coef(mf)),
        loglik = ll[1],
        loglik.se = ll[2],
        nfail.min = min(nfail),
        nfail.max = max(nfail),
        etime = as.numeric(etime)
      )
      }
  })

## Filter once more on maxima

profk_if <- bake(file = "output/k_fits_if_maxima.rds", {

  pars <- profk_if %>%
    filter(is.finite(loglik) & nfail.max == 0) %>%
    group_by(type, R0) %>%
    filter(rank(-loglik) <= 5) %>%
    ungroup() %>%
    select(-c(loglik, loglik.se, nfail.max, nfail.min, etime))

  foreach(
    start = iter(pars, by = "row"),
    .combine = rbind,
    .inorder = FALSE
  ) %dopar% {
      tic <- Sys.time()

      type <- as.character(start$type)

      st <- start %>%
        select(-country, -type) %>%
        unlist()

      po <- models[[country]]
      coef(po,names(st)) <- unname(st)

      ## Runs 10 particle filters to assess Monte Carlo error in likelihood
      pf <- try(replicate(10, pfilter(po, Np = 5000, max.fail = Inf)))

      toc <- Sys.time()
      etime <- toc - tic
      units(etime) <- "hours"

      ll <- sapply(pf, logLik)
      ll <- logmeanexp(ll, se = TRUE)
      nfail <- sapply(pf, getElement, "nfail")

      data.frame(
        country = country,
        type = type,
        as.list(coef(po)),
        loglik = ll[1],
        loglik.se = ll[2],
        nfail.min = min(nfail),
        nfail.max = max(nfail),
        etime = as.numeric(etime)
      )
      }
  })

# Combine iterative filtering runs
profIF <- profR0_if %>%
  mutate(profile = "R0") %>%
  bind_rows(mutate(profk_if, profile = "k"))

# Combine both models and save
all_models <- profIF %>%
  mutate(model = "Stochastic") %>%
  bind_rows(
    profTM %>%
      mutate(model = "Deterministic")
  )

saveRDS("output/model_profiles.RDS")
closeCluster(clust)