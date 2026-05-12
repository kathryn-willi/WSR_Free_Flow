library(tidyverse)
library(sf)
library(mapview)






watersheds_inside <- list.files("data/nwis_ws_backup/inside", full.names = TRUE) %>%
  map(~readRDS(.)) %>%
  bind_rows() %>%
  rowid_to_column("index") %>%
  st_transform(4326)

watersheds_outside <- st_read("data/nwis_ws_backup/outside/") %>% # also found at: "data/archive/shapefiles/watersheds.shp"
  filter(!site_no %in% watersheds_inside$site_no &
           site_no %in% str_remove(list.files("data/nwis_data/outside"), "\\.csv$")) %>%
  left_join(., read_csv("data/nwis_data/all_outside_gages.csv") %>% select(site_no, comid), by = "site_no") %>%
  select(site_no, comid) %>%
  st_transform(4326)

watersheds_wsr_raw <- st_read("data/wsr_watershed_raw.shp") 
watersheds_wsr <- list.files("data/archive/wsr_ws/", full.names = TRUE) %>%
  map_dfr(~readRDS(.))
mapview(watersheds_wsr)
