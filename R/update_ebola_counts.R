# =========================================================================== #
# Update Ebola case counts from HDX data source                               #
# Sean Browning                                                               #
# =========================================================================== #

# === Lib =================================================================== #
library(rhdx) # devtools::install_github("dickoa/rhdx")
library(readr)
library(dplyr)
library(lubridate)

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
  write_csv(file.path(data_folder, "health_zone_counts.csv"))

# Cumulative cases and deaths country-wide for model
drc_ebola_data <- health_zone_figures %>%
  mutate(date = as_date(report_date)) %>%
  group_by(country, date) %>%
  summarize(
    cases = sum(total_cases, na.rm = TRUE),
    deaths = sum(total_deaths, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  arrange(date) %>%
  mutate(
    epi_week = epiweek(date),
    year = epiyear(date)
  ) %>%
  group_by(epi_week, year) %>%
  mutate(sort_date = min(date)) %>%
  ungroup() %>%
  arrange(sort_date) %>%
  group_by(sort_date) %>%
  mutate(week = group_indices()) %>%
  ungroup() %>%
  select(country, week, date, cases, deaths) %>%
  arrange(date) %>%
  group_by(week) %>%
  filter(date == max(date)) %>%
  ungroup() %>%
  arrange(date) %>%
  mutate(
    new_cases = cases - lag(cases, default = NA),
    new_deaths = deaths - lag(deaths, default = NA)
  ) %>%
  write_csv(file.path(data_folder, "drc_model_counts.csv"))


# === Download shapefiles ================================================== #
# https://data.humdata.org/dataset/democratic-republic-of-congo-health-boundaries
shape_dataset <- pull_dataset("9690f4f5-d9e5-469a-b9dd-ae59b629a853")

districts <- shape_dataset %>%
  get_resource(1)

health_zones <- shape_dataset %>%
  get_resource(3)

districts$download(folder = file.path(data_folder, "shape_file"), filename = "drc_districts.geojson")
health_zones$download(folder = file.path(data_folder, "shape_file"), filename = "drc_health_zones.geojson")

# === Download denomonator data ============================================ #
# https://data.humdata.org/dataset/dr-congo-health-0
regional_denom <- pull_dataset("70ad1012-189c-4b49-b9f0-b97a71ec4c7b")

all_denom <- regional_denom %>%
  get_resource(1) %>%
  read_resource() %>%
  as_tibble() %>%
  select(
    province = PROVINCE,
    health_zone = Nom,
    pop = Population
  ) %>%
  write_csv(file.path(data_folder, "health_zone_population.csv"))