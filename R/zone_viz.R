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
library(gifski)
library(progress)
library(ggthemes)
library(ggmap)
library(patchwork)

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
  
# Bounding box for entire map
# To pull basemap
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
    )
  )

# Pull base maps
main_base_map <- get_stamenmap(bbox_shp, zoom = 6, maptype = "toner-lite")
inset_base_map <- get_stamenmap(bbox = c(27.5, -0.9, 31, 2.6), zoom = 9, maptype = "terrain-background")

i <- 1

pb <- progress_bar$new(
  format = ":elapsedfull [:bar] :current/:total (:percent) :eta",
  total = length(report_dates)
)

caption_unicode <- "Source: Minist\U00E8re de la Sant\U00E9, DRC"

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
        y = "New Cases (n)"
      )

  map_main <- main_base_map %>%
    ggmap() +
    geom_sf(
      data = shp_all %>% filter(report_date == rpt_date),
      aes(fill = total_cases),
      inherit.aes = FALSE,
      show.legend = FALSE
      ) +
      geom_rect(
        aes(xmin = 27.5, xmax = 31, ymin = -0.9, ymax = 2.6),
        fill = "transparent",
        size = 1,
        color = "black"
      ) +
      theme_map() +
      theme(
        panel.grid.major = element_line("transparent")
      ) +
      scale_fill_manual(values = col_scale)

  map_inset <- inset_base_map %>%
    ggmap() +
    geom_sf(
      data = shp_all %>% filter(report_date == rpt_date),
      aes(fill = total_cases),
      alpha = 0.7,
      inherit.aes = FALSE
    ) +
      coord_sf(xlim = c(27.5, 31), ylim = c(-0.9, 2.6), expand = FALSE) +
      theme_map() +
      theme(
        panel.grid.major = element_line("transparent"),
        legend.position = "right"
      ) +
      scale_fill_manual(
        values = col_scale,
        breaks = Filter(function(x) !is.na(x), unique(pull(shp_all, total_cases))),
        na.value = col_scale[1]
      ) +
      labs(
        fill = "Total Cases"
      )

    plot_all <- (map_main + map_inset + plot_layout(ncol = 2, widths = c(2, 1))) / case_plot + 
      plot_annotation(
        title = "Total Ebola Cases in Democratic Republic of Congo by District",
        subtitle = sprintf(
            "%s, n = %d",
            as.character(as_date(rpt_date), format = "%d %B %Y"),
            case_count %>%
              filter(report_date == rpt_date) %>%
              pull(total_cases)
          ),
        caption = enc2utf8(caption_unicode)
        ) +
      plot_layout(nrow = 2, heights = c(3, 1))

  CairoPNG(sprintf("output/graphic_frame_%03d.png", i), width = 1280, height = 720)

  print(plot_all)

  dev.off()

  i <- i + 1
  pb$tick()
}

# Make a GIF
gifski(list.files("output", pattern = "^.*\\.png$", full.names = TRUE), gif_file = "output/viz.gif", width = 1280, height = 720, delay = 0.1)