# ========================================================================== #
# Helperfuncs                                                                #
# Sean Browning                                                              #
# ========================================================================== #

library(parallel)

# === Simple helper to return a FORK/PSOCK cluster ========================
make_smp_cluster <- function() {
  # NOTES: On HPC runs, I set N_CORES env var
  # prior to running
  n_cores <- ifelse(
    Sys.getenv("N_CORES") == character(1),
    parallel::detectCores() - 1,
    as.numeric(Sys.getenv("N_CORES")) - 1
  )

  message(sprintf("Making cluster with %d cores", n_cores))

  clust <- parallel::makeCluster(n_cores, ifelse(.Platform$OS.type == "windows", "PSOCK", "FORK"))

  return(clust)
}

# === Helper function to save results ====================================
bake <- function(file, expr) {
  if (file.exists(file)) {
    readRDS(file)
  } else {
    val <- eval(expr)
    saveRDS(val, file = file)
    val
  }
}
