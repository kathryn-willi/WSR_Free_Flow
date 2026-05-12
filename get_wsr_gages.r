library(tidyverse)
library(sf)
library(dataRetrieval)
library(nhdplusTools)
library(mapview)

listNWIS <- function(aoi = NULL){
  
  # Grab NWIS by an area of interest:
  
  gage_sites <- vector("list", length = nrow(aoi))
  
  for (i in 1:nrow(aoi)){
    
    b_box <- sf::st_bbox(aoi[i,]) %>%
      as.vector() %>%
      round(., digits = 1) %>% 
      paste(collapse = ",")
    
    # Try to get data, skip if nothing returned or error occurs
    result <- tryCatch({
      dataRetrieval::whatNWISdata(bbox = b_box, parameterCd = "00060") %>%
        filter(stat_cd == "00003",
               site_tp_cd == "ST")
    }, error = function(e) {
      message("Error for bounding box ", i, ": ", e$message)
      return(NULL)
    })
    
    # Only add to list if result has rows
    if (!is.null(result) && nrow(result) > 0) {
      gage_sites[[i]] <- result
    } else {
      message("No data found for bounding box ", i)
      gage_sites[[i]] <- NULL
    }
    
  }
  
  # Remove NULL elements using tidyverse - keep only non-null elements
  gage_sites <- gage_sites %>%
    purrr::keep(~ !is.null(.x))
  
  # Check if any data was found
  if (length(gage_sites) == 0) {
    warning("No NWIS data found for any of the bounding boxes")
    return(NULL)
  }
  
  gage_sites <- dplyr::bind_rows(gage_sites) %>%
    sf::st_as_sf(coords = c('dec_long_va', 'dec_lat_va'), crs = 4269)
  
  if(st_crs(aoi) != st_crs(gage_sites)){
    
    gage_sites <- st_transform(gage_sites, st_crs(aoi))
    
  }
  
  inventory <- gage_sites %>%
    .[aoi,] %>%
    dplyr::select(c(site_no,
                    site_name = station_nm,
                    data_type_cd,
                    site_type_cd = site_tp_cd,
                    n_obs = count_nu,
                    begin_date,
                    end_date,
                    code = parm_cd)) %>%
    sf::st_join(., dplyr::select(aoi, WSR, state))
  
  return(inventory)
  
}

mainstem <- listNWIS(aoi = river_corridors) %>%
  mutate(location = "MAINSTEM") 

watershed <- listNWIS(aoi = final_watersheds %>% rename(WSR_RIVER_ = WSR, STATE = state)) %>%
  mutate(location = "WATERSHED") 

# Need >85% data availability for at least 20 years between water years 2000-2024
all_gages <- bind_rows(mainstem, watershed) %>%
  filter(begin_date <= "1999-10-01" & end_date >= "2025-09-30")

size <- readNWISsite(all_gages$site_no) %>%
  mutate(drain_area_km = drain_area_va * 2.58999) %>%
  filter(drain_area_km <= 1500)

good_gages <- all_gages %>%
  filter(site_no %in% size$site_no) %>%
  distinct(site_no)

for(i in 1:nrow(good_gages)){
  
  data <- readNWISdv(siteNumber = good_gages[i,], 
                     parameterCd = "00060",
                     startDate = "1999-10-01",
                     endDate = "2025-09-30")
  
  enough <- data %>%
    filter(!is.na(X_00060_00003),
           grepl("A|P", X_00060_00003_cd, ignore.case = FALSE))
  
  if(nrow(enough) >= 9497 * 0.85){
  write_csv(enough, paste0("data/nwis_data/inside/", good_gages[i,], ".csv"))
  }
  
}

gooder_gages <- all_gages %>% 
  filter(site_no %in% str_remove(list.files("data/nwis_data/inside/"), "\\.csv$")) %>%
  distinct(site_no, .keep_all = TRUE)


for(i in 1:nrow(gooder_gages)){
  comid <- get_nhdplus(AOI = gooder_gages[i,]) %>%
    pull(comid)
  try(getXYWatersheds(sf = gooder_gages[i,]) %>%
        mutate(site_no = gooder_gages[i,]$site_no,
               comid = comid) %>%
        saveRDS(paste0("data/nwis_ws_backup/inside/", gooder_gages[i,]$site_no, ".RDS")))
  Sys.sleep(3)
  
}

