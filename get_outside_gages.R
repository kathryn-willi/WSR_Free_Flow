library(tidyverse)
library(sf)
library(dataRetrieval)
library(nhdplusTools)
library(mapview)

nwis_watersheds <- st_read("data/archive/shapefiles/watersheds.shp")

nwis_sites <- nwis_watersheds %>% pull(site_no)

grabNWIS <- function(site_no = NULL){
  
  # Grab NWIS by an area of interest:
  
  gage_sites <- vector("list", length = length(site_no))
  
  for (i in 1:length(site_no)){
    
    # Try to get data, skip if nothing returned or error occurs
    result <- tryCatch({
      dataRetrieval::whatNWISdata(site = site_no[i], parameterCd = "00060") %>%
        dplyr::filter(stat_cd == "00003",
                      site_tp_cd == "ST") %>%
        dplyr::mutate(dplyr::across(everything(), as.character))
    }, error = function(e) {
      message("Error for site_no ", i, ": ", e$message)
      return(NULL)
    })
    
    # Only add to list if result has rows
    if (!is.null(result) && nrow(result) > 0) {
      gage_sites[[i]] <- result
    } else {
      message("No data found for site_no ", i)
      gage_sites[[i]] <- NULL
    }
    
  }
  
  # Remove NULL elements using tidyverse - keep only non-null elements
  gage_sites <- gage_sites %>%
    purrr::keep(~ !is.null(.x))
  
  # Check if any data was found
  if (length(gage_sites) == 0) {
    warning("No NWIS data found for any of the site_no's")
    return(NULL)
  }
  
  gage_sites <- dplyr::bind_rows(gage_sites) %>%
    sf::st_as_sf(coords = c('dec_long_va', 'dec_lat_va'), crs = 4269) %>%
    mutate(comid = 1)
  
  for(i in 1:nrow(gage_sites)){
    gage_sites$comid[i] <- discover_nhdplus_id(gage_sites[i,])
    
  }
  
  inventory <- gage_sites %>%
    dplyr::select(c(site_no,
                    site_name = station_nm,
                    comid,
                    data_type_cd,
                    site_type_cd = site_tp_cd,
                    n_obs = count_nu,
                    begin_date,
                    end_date,
                    code = parm_cd))
  
  return(inventory)
  
}

all_gages <- grabNWIS(site_no = nwis_sites) %>%
  mutate(location = "OUTSIDE") %>%
  # Need >85% data availability for at least 20 years between water years 2000-2024
  filter(begin_date <= "1999-10-01" & end_date >= "2025-09-30")

write_csv(all_gages, "data/nwis_data/all_outside_gages.csv")

size <- readNWISsite(all_gages$site_no) %>%
  mutate(drain_area_km = drain_area_va * 2.58999) %>%
  filter(drain_area_km <= 1500)

good_gages <- all_gages %>%
  filter(site_no %in% size$site_no) %>%
  distinct(site_no)

 for(i in 1:nrow(good_gages)) {
  
  data <- readNWISdv(siteNumber = good_gages[i,], 
                     parameterCd = "00060",
                     startDate = "1999-10-01",
                     endDate = "2025-09-30")
  
  # Dynamically find the first matching column names
  value_col <- grep("00060_00003$", names(data), value = TRUE)[1]
  code_col  <- grep("00060_00003_cd$", names(data), value = TRUE)[1]
  
  # Skip if missing expected columns
  if (is.na(value_col) || is.na(code_col)) {
    message("Skipping site ", good_gages[i,], " — expected columns not found.")
    next
  }
  
  enough <- data %>%
    filter(!is.na(.data[[value_col]]),
           grepl("A|P", .data[[code_col]], ignore.case = FALSE))
  
  if(nrow(enough) >= 9497 * 0.85){
    write_csv(enough, paste0("data/nwis_data/outside/", good_gages[i,], ".csv"))
  }
  
  print(i)
}



gooder_gages <- all_gages %>% 
  filter(site_no %in% str_remove(list.files("data/nwis_data/outside"), "\\.csv$")) %>%
  distinct(site_no, .keep_all = TRUE)

archive_wsr <- st_read("data/archive/shapefiles/watersheds.shp")
