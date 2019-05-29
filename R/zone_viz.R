# =================================================================================== #
# Basic visualization of DRC Ebola outbreak                                           #
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
library(ggforce)
library(gifski)
library(progress)

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
  )

i <- 1

pb <- progress_bar$new(
  format = ":elapsedfull [:bar] :current/:total (:percent) :eta",
  total = length(report_dates)
)

# === Iterate over every report date and render a separate PNG ========================= #
# NOTE: Could be made much faster with foreach / parallelized code
for (rpt_date in report_dates) {
  case_plot <- case_count %>%
    ggplot(aes(report_date, new_cases)) +
      geom_bar(aes(fill = report_date == rpt_date), stat = "identity", show.legend = FALSE) +
      scale_fill_manual(
        values = c(`TRUE` = "Red", `FALSE` = "Grey50")
      ) +
      labs(
        x = "Date Reported",
        y = "New Cases (n)",
        caption = "Source: Ministère de la Santé, DRC"
      )
  
  faceted_map <- shp_all %>%
    filter(report_date == rpt_date) %>%
    ggplot() +
      geom_sf(aes(fill = total_cases)) +
      scale_fill_distiller(palette = "OrRd", limits = c(0, 600)) +
      labs(
        title = "Total Ebola Cases in Democratic Republic of Congo by District",
        subtitle = sprintf(
          "%s, n = %d",
          as.character(as_date(rpt_date), format = "%d %B %Y"),
          case_count %>%
            filter(report_date == rpt_date) %>%
            pull(total_cases)
        ),
        fill = "Total Cases"
      ) +
      facet_zoom(xy = center_lon > 27.5 & center_lat > -2)
  
  CairoPNG(sprintf("output/graphic_frame_%03d.png", i), width = 1680, height = 1000)

  print(faceted_map / case_plot)

  dev.off()

  i <- i + 1
  pb$tick()
}

# Make a GIF
gifski(list.files("output", pattern = "^.*\\.png$", full.names = TRUE), gif_file = "output/viz.gif", width = 1680, height = 1000, delay = 0.1)