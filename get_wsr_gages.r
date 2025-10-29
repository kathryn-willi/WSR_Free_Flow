library(tidyverse)
library(sf)
library(dataRetrieval)
library(nhdplusTools)
library(mapview)

rivers <- st_read("data/S_USA.WildScenicRiver_LN/S_USA.WildScenicRiver_LN.shp") %>%
  filter(!STATE %in% c("Puerto Rico", "Alaska")) %>%
  filter(!WSR_RIVER_ %in% c("Green Wild and Scenic River", "Missouri Wild and Scenic River", "Rio Grande Wild and Scenic River", "Snake Wild and Scenic River"))

# mapview(rivers) + rivers_buf

get_river_endpoints <- function(sf) {
  
  # explode all multilinestring to linestring
  lines_df <- st_cast(sf, "LINESTRING")
  
  # convert line tips to coordinates
  endpoints_list <- map_dfr(1:nrow(lines_df), function(i) {
    coords <- st_coordinates(lines_df[i, ])
    
    # start point (first coordinate)
    start_pt <- st_point(c(coords[1, "X"], coords[1, "Y"]))
    
    # end point (last coordinate)  
    end_pt <- st_point(c(coords[nrow(coords), "X"], coords[nrow(coords), "Y"]))
    
    data.frame(original_feature_id = st_drop_geometry(lines_df[i, "WSR_RIVER_"]),
               STATE = st_drop_geometry(lines_df[i, "STATE"]),
               point_type = c("start", "end"),
               geometry = st_sfc(start_pt, end_pt, crs = st_crs(sf)))
  })
  
  # make the coordinates into points
  endpoints <- st_as_sf(endpoints_list)
  
  return(endpoints)
}

all_tips <- get_river_endpoints(rivers) %>%
  rowid_to_column()

getXYWatersheds <- function(sf = NULL, coordinates = NULL, crs = NULL){
  
  if(is.null(sf)){
    
    # Create a data frame with a column named 'geometry'
    df <- tibble::tibble(long = coordinates[1],
                         lat = coordinates[2])
    
    aoi_raw <- sf::st_as_sf(df, coords = c("long", "lat"), crs = crs) 
    
  }
  
  if(is.null(coordinates)){
    
    aoi_raw <- sf
  }
  
  if(st_crs(aoi_raw)$epsg != 4326){
    
    aoi <- aoi_raw %>% st_transform(crs = 4326)
    
  } else {
    
    aoi <- aoi_raw
  }
  
  # Get the NHDPlus flowline near the site and create a 30m buffer around it
  # This is considered a "danger zone" to avoid placing points on the flowline
  # if they aren't in fact on the flowline
  flowline_danger_zone <- get_nhdplus(AOI = aoi, realization = "flowline") %>%
    st_buffer(30)
  
  # is point on a small tributary? (ie too small for NHD flowlines)
  if(is.na(as.numeric(st_intersects(aoi %>% st_buffer(100), flowline_danger_zone)))){
    
    
    # If the site is known to be a small watershed but is still very near "main" flowline...
    if(!is.na(as.numeric(st_intersects(aoi %>% st_buffer(35), flowline_danger_zone)))){
      
      # create a 35m buffer around the site location
      site_buffer <- aoi %>%
        st_buffer(35) 
      
      # get points from the buffer boundary while avoiding any along the flowline
      boundary_points <- site_buffer %>%
        st_boundary() %>%      
        st_cast("POINT") %>%    
        mutate(point_id = row_number()) %>%   
        st_difference(., flowline_danger_zone) %>%  # remove points that fall within the flowline buffer
        # reduce the number of points by taking only every 50th point
        filter(point_id %% 50 == 0) %>%
        bind_rows(aoi)  # add the original site point to the collection
      
      splits <- vector("list")
      
      # loop through each boundary point to get split catchments
      for(i in 1:nrow(boundary_points)){
        Sys.sleep(2)
        # Get the split catchment for each point
        # We're looking for areas without a catchmentID (NA values)
        splits[[i]] <- get_split_catchment(point = st_as_sfc(boundary_points[i,])) %>%
          filter(is.na(catchmentID)) %>%
          mutate(area = as.numeric(st_area(.)))  # Calculate the area
      }
      
      # Combine all splits and select the largest one
      splits <- bind_rows(splits) %>%
        filter(as.numeric(area) == max(as.numeric(area))) %>%  # Keep only the largest area
        select(-catchmentID,-id) 
      
      return(splits)
      
    } else {
      Sys.sleep(2)
      # If site is on small trib but far away from flowline
      mini_ws <- get_split_catchment(point = st_as_sfc(aoi)) %>%
        filter(is.na(catchmentID)) %>%
        mutate(area = as.numeric(st_area(.))) %>%
        select(-catchmentID,-id)
      
      if(st_crs(mini_ws) != st_crs(aoi_raw)){
        
        mini_ws <- mini_ws %>% st_transform(crs = st_crs(aoi_raw)$epsg)
        
      }
      
      return(mini_ws)
      
    }
    
  } else {
    
    # Site is along flowline...
    flowline <- get_nhdplus(AOI = aoi, realization = "flowline", t_srs = 4326)
    
    # "Snap" our site to the nearest NHD flowline feature
    nearest_points <- st_nearest_points(aoi, flowline)
    
    snapped_points_sf <- st_cast(nearest_points, "POINT")[2,]
    
    trace <- get_raindrop_trace(snapped_points_sf, direction = "down")
    Sys.sleep(2)
    raindrop <- sf::st_sfc(sf::st_point(trace$intersection_point[[1]][1:2]),
                           crs = 4326)
    # Clip/split our catchment to only include the portion of the
    # catchment upstream of our site and grab whole watershed 
    nhd_catch <- get_split_catchment(raindrop, upstream = T)[2,] %>%
      st_make_valid() %>%
      mutate(area = as.numeric(st_area(.))) %>%
      select(area)
    
    if(st_crs(nhd_catch) != st_crs(aoi_raw)){
      
      nhd_catch <- nhd_catch %>% st_transform(crs = st_crs(aoi_raw)$epsg)
      
    }
    
    return(nhd_catch)
    
  }
  
}  


ws <- vector("list", length = nrow(all_tips))

for(i in 1:nrow(all_tips)){
  
  try(ws[[i]] <- getXYWatersheds(sf = all_tips[i,]) %>%
        mutate(WSR = all_tips[i,]$WSR_RIVER_,
               state = all_tips[i,]$STATE) %>%
        saveRDS(paste0("data/wsr_ws_backup//", all_tips[i,]$rowid, ".RDS")))
  Sys.sleep(3)
  
}


# which ones didn't work?

missing <- all_tips %>% filter(!rowid %in% as.numeric(str_remove(list.files("data/wsr_ws_backup//"), "\\.RDS$")))

for(i in 1:nrow(missing)){
  
  try(ws[[i]] <- getXYWatersheds(sf = missing[i,]) %>%
        mutate(WSR = missing[i,]$WSR_RIVER_,
               state = missing[i,]$STATE) %>%
        saveRDS(paste0("data/ws_backup/", missing[i,]$rowid, ".RDS")))
  Sys.sleep(3)
  
}






test1 <- list.files("data/wsr_ws_backup//", full.names = TRUE) %>% #[1:500] %>%
  map(~readRDS(.)) %>%
  bind_rows() %>%
  filter(!WSR %in% c("White Clay Wild and Scenic River", "Great Egg Harbor Wild and Scenic River", 
                     "Flathead Wild and Scenic River", "Smith Wild and Scenic River")) %>%
  st_make_valid() %>%
  group_by(WSR, state) %>%
  summarize() %>%
  nngeo::st_remove_holes()

test2 <- list.files("data/wsr_ws_backup//", full.names = TRUE) %>%
  map(~readRDS(.)) %>%
  bind_rows() %>%
  filter(WSR %in% c("White Clay Wild and Scenic River", "Great Egg Harbor Wild and Scenic River", 
                    "Flathead Wild and Scenic River", "Smith Wild and Scenic River")) %>%
  filter(!st_is_empty(.)) %>%  
  st_make_valid() %>%
  st_cast("POLYGON") %>%
  group_by(WSR, state) %>%
  summarize(geometry = st_union(geometry)) %>%
  .[c(2,3,4,6),] %>%
  #filter(st_geometry_type(.) == "POLYGON")# c("POLYGON", "MULTIPOLYGON")) 
  #st_cast("MULTIPOLYGON")
  nngeo::st_remove_holes()

final_watersheds <- bind_rows(test1, test2) %>%
  ungroup() %>%
  st_make_valid() %>%
  group_by(WSR, state) %>%
  summarize() %>%
  nngeo::st_remove_holes()

river_corridors <- rivers %>%
  st_transform(crs = 32616) %>% 
  st_buffer(402.336, endCapStyle = "FLAT") %>% # quarter mile around lines
  rowid_to_column() %>%
  st_transform(4269)

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
    sf::st_join(., dplyr::select(aoi, WSR_RIVER_ = WSR, state))
  
  return(inventory)
  
}

mainstem <- listNWIS(aoi = river_corridors) %>%
  mutate(location = "MAINSTEM") %>%
  rename(WSR = WSR_RIVER_)

watershed <- listNWIS(aoi = final_watersheds) %>%
  mutate(location = "WATERSHED") %>%
  rename(STATE = state)

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
  
  data <- readNWISdv(siteNumber = good_gages[i,]$site_no, 
                     parameterCd = "00060",
                     startDate = "1999-10-01",
                     endDate = "2025-09-30")
  
  enough <- data %>%
    filter(!is.na(X_00060_00003),
           grepl("A|P", X_00060_00003_cd, ignore.case = FALSE))
  
  if(nrow(enough) >= 9497 * 0.85){
  write_csv(enough, paste0("data/nwis_data/inside/", good_gages[i,]$site_no, ".csv"))
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

