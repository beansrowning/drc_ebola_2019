# =================================================================================== #
# Basic visualization of DRC Ebola outbreak                                           #
# MPI-based cluster render (also PSOCK for desktop usage)                             #
# Sean Browning                                                                       #
# =================================================================================== #

# === Lib =========================================================================== #
library(dplyr)
library(tidyr)
library(lubridate)
library(janitor)
library(ggplot2)
library(sf)
library(Cairo)
library(patchwork)
library(ggthemes)
library(ggmap)
library(gifski)
library(parallel)
# library(Rmpi)
# library(doMPI)
library(doParallel)
library(iterators)
library(foreach)

# === Make Cluster ================================================================
# PSOCK (for home use)
clust <- makeCluster(detectCores(), "PSOCK")
registerDoParallel(clust)

# MPI (For high performance cluster)
# clust <- startMPIcluster()
# registerDoMPI(clust)

# === Read in data and shapefile ==================================================
# BUG: read_csv was being terrible at guessing types, so I reverted to base
# NOTE: A few health_zone names need to be re-mapped manually to join properly
ebola_data <- read.csv("data/health_zone_counts.csv", stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  mutate_at(vars(matches("date")), ymd) %>%
  mutate(
    health_zone =  case_when(
      health_zone == "Mangurujipa" ~ "Manguredjipa",
      health_zone == "Rwampara (Bunia)" ~ "Rwampara",
      health_zone == "Kayina" ~ "Kayna",
      health_zone == "N/A" ~ NA_character_,
      TRUE ~ health_zone
    )
  )

shp <- st_read(dsn = "data/shape_file/drc_districts.geojson", layer = "drc_districts") %>%
  clean_names() %>%
  filter(adm0_viz_name == "Democratic Republic of the Congo") %>%
  mutate(
    adm2_vis_name = tools::toTitleCase(tolower(adm2_name))
  )

# Bounding box for entire map to pull basemap
bbox_shp <- shp %>%
  st_bbox() %>%
  unname()

# Bounding box for the inset map
bbox_inset <- c(xmin = 27.5, ymin = -0.9, xmax = 31, ymax = 2.6)

# Testing just plotting the features in the inset plot
sf_inset <- st_as_sfc(st_bbox(bbox_inset, crs = 4326))

shp <- shp %>%
  st_intersection(sf_inset)

# Color scale for choropleth
col_scale <- c("#00a650", "#318a4a", "#5f6f40", "#8e5236", "#be382d", "#ee1c25")

# Unicode caption
caption_unicode <- "Source: Minist\U00E8re de la Sant\U00E9, DRC"

# Extra tailing frames for the latest date
n_extra <- 20

# Time delay per frame
n_delay <- 0.1

# === Summarize and combine data =================================================

adm2_names <- shp %>%
  as.data.frame() %>%
  pull(adm2_vis_name)

case_tmp <- ebola_data %>%
  select(report_date, adm2_vis_name = health_zone, total_cases) %>%
  mutate(adm2_vis_name = factor(adm2_vis_name, levels = adm2_names)) %>%
  complete(report_date, crossing(adm2_vis_name)) %>%
  group_by(adm2_vis_name) %>%
  arrange(report_date) %>%
  fill(total_cases, .direction = "down") %>%
  ungroup()

case_count <- case_tmp %>%
  group_by(report_date) %>%
  summarize(total_cases = sum(total_cases, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(new_cases = total_cases - lag(total_cases, default = NA))


report_dates <- ebola_data %>%
  arrange(report_date) %>%
  distinct(report_date) %>%
  pull()

shp_all <- shp %>%
  left_join(
    case_tmp,
    by = "adm2_vis_name"
  ) %>%
  mutate(
    total_cases = case_when(
      total_cases == 0 ~ "0",
      total_cases > 0 & total_cases <= 10 ~ "[1,10]",
      total_cases > 10 & total_cases <= 50 ~ "(10,50}",
      total_cases > 50 & total_cases <= 100 ~ "(50,100]",
      total_cases > 100 & total_cases <= 200 ~ "(100,200]",
      total_cases > 200 ~ "200+"
    ),
    total_cases = factor(total_cases, levels = c("0", "[1,10]", "(10,50}", "(50,100]", "(100,200]", "200+"), ordered = TRUE)
  )

# Pull base maps
# Stamen maps, colored only in inset map looks better IMO
main_base_map <- get_stamenmap(bbox_shp, zoom = 6, maptype = "toner-lite")
inset_base_map <- get_stamenmap(bbox = unname(bbox_inset), zoom = 9, maptype = "terrain-background")

# === Create zipped iter of our data =================================================== #
# MPI does not share memory, so each worker will only receive the chunk of data
# it needs to create its frame
data_split <- split(shp_all, shp_all$report_date)
data_split <- c(data_split, rep(data_split[length(data_split)], n_extra))

report_dates <- c(report_dates, rep(report_dates[length(report_dates)], n_extra))

data_iter <- iter(data_split)
date_iter <- iter(report_dates)

message(sprintf("%s - Starting PNG render via MPI", Sys.time()))

# === Iterate over every report date and render a separate PNG ========================= #
# NOTE: Could be made much faster with foreach / parallelized code
foreach(
  data = data_iter,
  rpt_date = date_iter,
  i = seq_along(report_dates),
  .packages = c(
    "ggplot2", "sf", "Cairo", "ggforce", "dplyr", "patchwork", "ggmap", "ggthemes"
  ),
  .noexport = c(
    "shp_all", "adm2_names", "ebola_data", "shp_names",
    "report_dates", "case_tmp"
  ),
  .export = c(
    "case_count", "main_base_map", "inset_base_map",
    "col_scale", "caption_unicode", "bbox_inset"
  ),
  .inorder = FALSE,
  .combine = function(...) {return(NULL)}
) %dopar% {

  case_plot <- case_count %>%
    ggplot(aes(report_date, new_cases)) +
      geom_bar(aes(fill = report_date == rpt_date), stat = "identity", show.legend = FALSE) +
      scale_fill_manual(
        values = c(`TRUE` = "Red", `FALSE` = "Grey50")
      ) +
      scale_x_date(
        date_breaks = "1 month",
        date_labels = "%b"
      ) +
      labs(
        x = "Date Reported",
        y = "New Cases (n)"
      )
  
  map_main <- main_base_map %>%
    ggmap() +
    geom_sf(
      data = data,
      aes(fill = total_cases),
      inherit.aes = FALSE,
      show.legend = FALSE
      ) +
      geom_rect(
        aes(xmin = bbox_inset[1], xmax = bbox_inset[3], ymin = bbox_inset[2], ymax = bbox_inset[4]),
        fill = "transparent",
        size = 1,
        color = "black"
      ) +
      theme_map() +
      theme(
        panel.grid.major = element_line("transparent")
      ) +
      scale_fill_manual(values = col_scale, na.value = col_scale[1])

  map_inset <- inset_base_map %>%
    ggmap() +
    geom_sf(
      data = data,
      aes(fill = total_cases),
      alpha = 0.7,
      inherit.aes = FALSE
    ) +
      coord_sf(xlim = bbox_inset[c(1,3)], ylim = bbox_inset[c(2,4)], expand = FALSE) +
      theme_map() +
      theme(
        panel.grid.major = element_line("transparent"),
        legend.position = "right"
      ) +
      scale_fill_manual(
        values = col_scale,
        breaks = levels(pull(data, total_cases)),
        na.value = col_scale[1],
        drop = FALSE
      ) +
      labs(
        fill = "Total Cases"
      )

    plot_all <- (map_main + map_inset + plot_layout(ncol = 2, widths = c(2, 1))) / case_plot + 
      plot_annotation(
        title = "Total Ebola Cases in Democratic Republic of Congo by District",
        subtitle = sprintf(
          "%s, n = %s  (confirmed + probable)",
          as.character(rpt_date, format = "%d %B %Y"),
          case_count %>%
            filter(report_date == rpt_date) %>%
            pull(total_cases) %>%
            format(big.mark = ",")
        ),
        caption = enc2utf8(caption_unicode)
        ) +
      plot_layout(nrow = 2, heights = c(3, 1))
  
  CairoPNG(sprintf("output/graphic_frame_%03d.png", i), width = 1280, height = 720)

  print(plot_all)

  dev.off()
}

# === Stop Cluster =================================================================
# MPI
# closeCluster(clust)
# PSOCK
stopCluster(clust)
message(sprintf("%s - Cluster closed, starting GIF render", Sys.time()))

# === Make the GIF =================================================================
gifski(
  list.files("output", pattern = "^.*\\.png$", full.names = TRUE),
  gif_file = "output/viz.gif",
  width = 1280,
  height = 720,
  delay = n_delay
)

message(sprintf("%s - GIF render complete, shutting down", Sys.time()))

# === Terminate MPI =======================================================================
# mpi.quit()