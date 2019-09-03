# ================================================================================ #
# Ebola SEIR Simulation (From King et al., 2015)                                   #
# NOTE: SMP/PSOCK parallelism                                                      #
# Sean Browning                                                                    #
# ================================================================================ #

library(pomp)
library(dplyr)
library(doParallel)
library(foreach)
library(iterators)

# Load model and data
source("R/model.R")
source("R/util.R")

out_dir <- normalizePath("output/model_fits", winslash = "/", mustWork = TRUE)

# === Set up cluster ==============================================
clust <- make_smp_cluster()
registerDoParallel(clust)

# === Trajectory matching simulation study =========================================

tic <- Sys.time()

po <- ebolaModel(country = "DRC", data = drc_ebola_data, timestep = 0.1)
params <- parmat(coef(po), 3)
params["k", ] <- c(0, 0.2, 0.5)
paramnames <- names(coef(po))

nsims <- 500

sim_dat <- bake(file = file.path(out_dir, "sims.rds"), {
  simulate(
    po,
    params = params,
    nsim = nsims,
    seed = 208335746L,
    format = "data.frame"
  ) %>% 
  rename(time = week) %>%
  mutate(
    sim = as.integer(gsub("^(\\d*)_(\\d*)$", "\\2", .id)),
    param_set = as.integer(gsub("^(\\d*)_(\\d*)$", "\\1", .id)),
    k = params["k", param_set]
  )
})

pompUnload(po)

profiles <- bake(file = file.path(out_dir, "sim_profiles_R0_deter.rds"), {
  fits <- foreach(
    dat = isplit(sim_dat, sim_dat[[".id"]]),
    .combine = rbind,
    .inorder = FALSE
  ) %:%
    foreach(
      type = c("raw", "cum"),
      .combine = rbind,
      .inorder = FALSE,
      .packages = c("pomp", "subplex", "dplyr")
    ) %dopar% {

      tm <- ebolaModel(
        country = "DRC",
        data = select(dat$value, week, cases, deaths),
        type = as.character(type)
      )

      st <- params[, unique(dat$value[["param_set"]])]
      true.k <- unname(st["k"])
      true.R0 <- unname(st["R0"])
      true.rho <- unname(st["rho"])
      st["k"] <- st["k"] + 1e-6

      tm <- traj_objfun(
        tm,
        start = st,
        est = c("R0", "k", "rho")
      )

      if (coef(tm, "rho") == 1) coef(tm, "rho") <- 0.999
      if (coef(tm, "rho") == 0) coef(tm, "rho") <- 0.001

      tm_optim <- subplex(
        fn = tm,
        par = coef(tm, c("R0", "k", "rho"))
      )

      pompUnload(tm)

      data.frame(
        sim = dat$key[[1]],
        param_set = unique(dat$value[["param_set"]]),
        type = as.character(type),
        true.k = true.k,
        true.R0 = true.R0,
        true.rho = true.rho,
        as.list(coef(tm)),
        loglik = logLik(tm),
        conv = tm_optim$convergence
      )
    }

  foreach(
    fit = iter(fits, by = "row"),
    .combine = rbind,
    .inorder = FALSE
  ) %:%
    foreach(
      r0 = seq(from = 0.7, to = 3, length = 200),
      .combine = rbind,
      .inorder = FALSE
    ) %dopar% {
      dat <- sim_dat %>%
        filter(sim == fit$sim, param_set = fit$param_set) %>%
        select(week, cases, deaths)

      tm <- ebolaModel(
        country = "DRC",
        data = dat,
        type = as.character(fit$type)
      )

      coef(tm) <- unlist(fit[paramnames])
      coef(tm, "R0") <- r0
      if (coef(tm, "rho") == 1) coef(tm, "rho") <- 0.999
      if (coef(tm, "rho") == 0) coef(tm, "rho") <- 0.001

      tm <- traj_objfun(
        tm,
        est = c("k", "rho")
      )

      if (coef(tm, "rho") == 1) coef(tm, "rho") <- 0.999
      if (coef(tm, "rho") == 0) coef(tm, "rho") <- 0.001
      
      tm_optim <- subplex(
        fn = tm,
        par = coef(tm, c("k", "rho"))
      )
      
      # Re-run parameter values through
      # stateful function to make sure values are
      # accurate
      tm(tm_optim$par)

      pompUnload(tm)

      data.frame(
        sim = fit$sim,
        param_set = fit$param_set,
        type = fit$type,
        true.k = fit$true.k,
        true.R0 = fit$true.R0,
        true.rho = fit$true.rho,
        as.list(coef(tm)),
        loglik = logLik(tm),
        conv = tm_optim$convergence
      )
    }
})

fits <- bake(file = file.path(out_dir, "tm-sim-fits.rds"), {
  starts <- profiles %>%
    group_by(type, sim, true.k) %>%
    filter(loglik == max(loglik)) %>%
    ungroup()

  foreach(
    fit = iter(starts, by = "row"),
    .combine = rbind,
    .inorder = FALSE,
  ) %dopar% {
    dat <- sim_dat %>%
      filter(sim == fit$sim, param_set == fit$param_set) %>%
      select(week, cases, deaths)

    tm <- ebolaModel(
      country = "DRC",
      data = dat,
      type = as.character(fit$type)
    )

    coef(tm) <- unlist(fit[paramnames])
    if (coef(tm, "rho") == 1) coef(tm, "rho") <- 0.999
    if (coef(tm, "rho") == 0) coef(tm, "rho") <- 0.001


    tm <- traj_objfun(
      tm,
      est = c("k", "rho")
    )


    if (coef(tm, "rho") == 1) coef(tm, "rho") <- 0.999
    if (coef(tm, "rho") == 0) coef(tm, "rho") <- 0.001

    tm_optim <- subplex(
      fn = tm,
      par = coef(tm, est = c("k", "rho"))
    )

    # Re-run parameter values through
    # stateful function to make sure values are
    # accurate
    tm(tm_optim$par)

    pompUnload(tm)

    data.frame(
      sim = fit$sim,
      param_set = fit$param_set,
      type = fit$type,
      true.k = fit$true.k,
      true.R0 = fit$true.R0,
      true.rho = fit$true.rho,
      as.list(coef(tm)),
      loglik = logLik(tm),
      conv = tm_optim$convergence
    )
  }
})

toc <- Sys.time()
print(toc - tic)

# === Trajectory matching with least squares simulation ================================

tic <- Sys.time()

profiles_ls <- bake(file = "ls-sim-profiles-R0.rds", {
  fits <- foreach(
    dat = isplit(sim_dat, sim_dat[[".id"]]),
    .combine = rbind,
    .inorder = FALSE
  ) %:%
    foreach(
      type = c("raw", "cum"),
      .combine = rbind,
      .inorder = FALSE
    ) %dopar% {

      tm <- ebolaModel(
        country = "DRC",
        data = select(dat$value, week, cases, deaths),
        type = as.character(type),
        least.sq = TRUE
      )

      st <- params[, unique(dat$value[["param_set"]])]
      true.k <- unname(st["k"])
      true.R0 <- unname(st["R0"])
      true.rho <- unname(st["rho"])
      st["k"] <- 10

      tm <- traj_objfun(
        tm,
        start = st,
        est = c("R0", "k", "rho")
      )
      
      if (coef(tm, "rho") == 1) coef(tm, "rho") <- 0.999
      if (coef(tm, "rho") == 0) coef(tm, "rho") <- 0.001
      
      tm_optim <- subplex(
        fn = tm,
        par = coef(tm, c("R0", "k", "rho"))
      )

      # Re-run parameter values through
      # stateful function to make sure values are
      # accurate
      tm(tm_optim$par)

      pompUnload(tm)

      data.frame(
        sim = dat$key[[1]],
        param_set = unique(dat$value[["param_set"]]),
        type = as.character(type),
        true.k = true.k,
        true.R0 = true.R0,
        true.rho = true.rho,
        as.list(coef(tm)),
        loglik = logLik(tm),
        conv = tm_optim$convergence
      )
    }

  foreach(
    fit = iter(fits, by = "row"),
    .combine = rbind,
    .inorder = FALSE
  ) %:%
    foreach(
      r0 = seq(from = 0.7, to = 3, length = 200),
      .combine = rbind,
      .inorder = FALSE,
    ) %dopar% {
      dat <- sim_dat %>%
        filter(sim == fit$sim) %>%
        select(week, cases, deaths)

      tm <- ebolaModel(
        country = "DRC",
        data = dat,
        type = as.character(fit$type),
        least.sq = TRUE
      )

      coef(tm) <- unlist(fit[paramnames])
      coef(tm, "R0") <- r0

      if (coef(tm, "rho") == 1) coef(tm, "rho") <- 0.999
      if (coef(tm, "rho") == 0) coef(tm, "rho") <- 0.001

      tm <- traj_objfun(
        tm,
        est = c("k", "rho")
      )

      if (coef(tm, "rho") == 1) coef(tm, "rho") <- 0.999
      if (coef(tm, "rho") == 0) coef(tm, "rho") <- 0.001
      
      tm_optim <- subplex(
        fn = tm,
        par = coef(tm, c("k", "rho")),
        control = list(maxit = 1e5)
      )

      # Re-run parameter values through
      # stateful function to make sure values are
      # accurate
      tm(tm_optim$par)

      pompUnload(tm)

      data.frame(
        sim = fit$sim,
        type = fit$type,
        true.k = fit$true.k,
        true.R0 = fit$true.R0,
        true.rho = fit$true.rho,
        as.list(coef(tm)),
        loglik = logLik(tm),
        conv = tm_optim$convergence
      )
    }
})

fits_ls <- bake(file = "ls-sim-fits.rds", {
  starts <- profiles %>%
    group_by(type, sim, true.k) %>%
    filter(loglik == max(loglik))
    ungroup()

  foreach(
    fit = iter(starts, by = "row"),
    .combine = rbind,
    .inorder = FALSE,
  ) %dopar% {
    dat <- sim_dat %>%
      filter(sim == fit$sim, param_set = fit$param_set) %>%
      select(week, cases, deaths)

    tm <- ebolaModel(
      country = "DRC",
      data = dat,
      type = as.character(fit$type),
      least.sq = TRUE
    )
    if (coef(tm, "rho") == 1) coef(tm, "rho") <- 0.999
    if (coef(tm, "rho") == 0) coef(tm, "rho") <- 0.001

    tm <- traj_objfun(
      tm,
      est = c("k", "rho")
    )

    if (coef(tm, "rho") == 1) coef(tm, "rho") <- 0.999
    if (coef(tm, "rho") == 0) coef(tm, "rho") <- 0.001
    
    tm_optim <- subplex(
      fn = tm,
      par = coef(tm, c("k", "rho")),
      control = list(maxit = 1e5)
    )

    # Re-run parameter values through
    # stateful function to make sure values are
    # accurate
    tm(tm_optim$par)

    pompUnload(tm)

    data.frame(
      sim = fit$sim,
      type = fit$type,
      param_set = fit$param_set,
      true.k = fit$true.k,
      true.R0 = fit$true.R0,
      true.rho = fit$true.rho,
      as.list(coef(tm)),
      loglik = logLik(tm),
      conv = tm_optim$convergence
    )
  }
})

toc <- Sys.time()
print(toc - tic)
