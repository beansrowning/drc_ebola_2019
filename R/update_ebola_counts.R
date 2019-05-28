# =========================================================================== #
# Update Ebola case counts from HDX data source                               #
# Sean Browning                                                               #
# =========================================================================== #

# === Lib =================================================================== #
library(rhdx) # devtools::install_github("dickoa/rhdx")
library(readr)

# Set dirs
data_folder <- normalizePath("data", winslash = "/", mustWork = TRUE)
tmp_folder <- tempdir()

# point data source at prod data
set_rhdx_config(hdx_site = "prod")
rhdx_cache$cache_path_set(normalizePath(".cache", winslash = "/"))

# === Pull Ebola case and death count dataset ============================== #
# https://data.humdata.org/dataset/ebola-cases-and-deaths-drc-north-kivu
dataset <- pull_dataset("383945db-762c-46e2-bec6-07adc41fbd16")

health_zone_figures <- dataset %>%
  get_resource(2) %>%
  read_resource(folder = tmp_folder) %>%
  write_csv(file.path(data_folder,"health_zone_counts.csv"))

# === Download shapefiles ================================================== #
# https://data.humdata.org/dataset/democratic-republic-of-congo-health-boundaries
shape_dataset <- pull_dataset("9690f4f5-d9e5-469a-b9dd-ae59b629a853")

districts <- shape_dataset %>%
  get_resource(1)

health_zones <- shape_dataset %>%
  get_resource(3)

districts$download(folder = file.path(data_folder, "shape_file"), filename = "drc_districts.geojson")
health_zones$download(folder = file.path(data_folder, "shape_file"), filename = "drc_health_zones.geojson")


