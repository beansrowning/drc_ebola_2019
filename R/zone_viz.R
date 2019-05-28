# =================================================================================== #
# Basic visualization of DRC Ebola outbreak                                           #
# Sean Browning                                                                       #
# =================================================================================== #

# === Lib =========================================================================== #
library(dplyr)
library(lubridate)
library(janitor)
library(ggplot2)
library(sf)
library(Cairo)
library(patchwork)
library(ggforce)
library(gifski)

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
case_count <- ebola_data %>%
  group_by(report_date) %>%
  summarize(total_cases = sum(confirmed_cases, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(total_cases = total_cases - lag(total_cases, default = 0L))


report_dates <- ebola_data %>%
  arrange(report_date) %>%
  distinct(report_date) %>%
  pull()

shp_all <- shp %>%
  left_join(
    ebola_data %>%
      select(
        report_date,
        adm2_vis_name = health_zone,
        total_cases
      ),
    by = "adm2_vis_name"
  )

i <- 1

# === Iterate over every report date and render a separate PNG ========================= #
# NOTE: Could be made much faster with foreach / parallelized code
for (rpt_date in report_dates) {
  case_plot <- case_count %>%
    ggplot(aes(report_date, total_cases)) +
      geom_bar(aes(fill = report_date == rpt_date), stat = "identity", show.legend = FALSE) +
      scale_fill_manual(
        values = c(`TRUE` = "Red", `FALSE` = "Grey50")
      ) +
      labs(
        x = "Date Reported",
        y = "New Cases (n)",
        caption = "Source: Ministère de la Santé, DRC"
      )
  
  shp_current <- shp %>%
    left_join(
      ebola_data %>%
      group_by(health_zone) %>%
      filter(
        report_date <= rpt_date,
      ) %>%
      arrange(desc(report_date)) %>%
      slice(1) %>%
      ungroup() %>%
      select(
        report_date,
        adm2_vis_name = health_zone,
        total_cases
      ),
    by = "adm2_vis_name"
    )

  faceted_map <- shp_current %>%
    ggplot() +
      geom_sf(aes(fill = total_cases)) +
      scale_fill_distiller(palette = "OrRd", limits = c(0, 600)) +
      labs(
        title = "Total Ebola Cases in Democratic Republic of Congo by Health Zone",
        subtitle = as.character(as_date(rpt_date), format = "%d %B %Y"),
        fill = "Total Cases"
      ) +
      facet_zoom(xy = center_lon > 27.5 & center_lat > -2)
  
  CairoPNG(sprintf("output/graphic_frame_%03d.png", i), width = 1680, height = 1000)

  print(faceted_map / case_plot)

  dev.off()

  i <- i + 1
}

# Make a GIF
gifski(list.files("output", full.names = TRUE), gif_file = "output/viz.gif", width = 1680, height = 1000, delay = 0.2)