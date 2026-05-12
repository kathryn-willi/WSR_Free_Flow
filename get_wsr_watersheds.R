library(tidyverse)
library(sf)
library(dataRetrieval)
library(nhdplusTools)
library(mapview)

rivers <- st_read("data/S_USA.WildScenicRiver_LN/S_USA.WildScenicRiver_LN.shp")
  # readRDS("data/wsr_official_lines.RDS") %>%
  # # requires different data
  # filter(!grepl("Alaska|Puerto Rico", WSR_RIVER_, ignore.case = TRUE))

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
               point_type = c("start", "end"),
               geometry = st_sfc(start_pt, end_pt, crs = st_crs(sf)))
  })
  
  # make the coordinates into points
  endpoints <- st_as_sf(endpoints_list)
  
  return(endpoints)
}

all_tips <- get_river_endpoints(rivers) %>%
  rowid_to_column()

ws <- vector("list", length = nrow(all_tips))

getXYWatersheds <- function(sf = NULL, coordinates = NULL, crs = NULL, snap = FALSE){
  
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
  
  if(snap == FALSE){
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
          
          # Get the split catchment for each point
          # We're looking for areas without a catchmentID (NA values)
          splits[[i]] <- get_split_catchment(point = st_as_sfc(boundary_points[i,])) %>%
            filter(is.na(catchmentID)) %>%
            st_make_valid() %>%
            mutate(area = as.numeric(st_area(.)))  # Calculate the area
        }
        
        # Combine all splits and select the largest one
        splits <- bind_rows(splits) %>%
          filter(as.numeric(area) == max(as.numeric(area))) %>%  # Keep only the largest area
          select(-catchmentID,-id) 
        
        return(splits)
        
      } else {
        
        # If site is on small trib but far away from flowline
        mini_ws <- get_split_catchment(point = st_as_sfc(aoi)) %>%
          filter(is.na(catchmentID)) %>%
          st_make_valid() %>%
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
      
    }
    
  }
  
  if(snap == TRUE){
    
    # Site is along flowline...
    flowline <- get_nhdplus(AOI = aoi, realization = "flowline", t_srs = 4326)
    
    # "Snap" our site to the nearest NHD flowline feature
    nearest_points <- st_nearest_points(aoi, flowline)
    
    snapped_points_sf <- st_cast(nearest_points, "POINT")[2,]
    
    trace <- get_raindrop_trace(snapped_points_sf, direction = "down")
    
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
    
    
  }
  
  
  return(nhd_catch)
  
}

for(i in 1:nrow(all_tips)){
  
  try(ws[[i]] <- getXYWatersheds(sf = all_tips[i,], snap = TRUE) %>%
        mutate(WSR_RIVER_ = all_tips[i,]$WSR_RIVER_) %>%
        saveRDS(paste0("data/wsr_ws_backup/", all_tips[i,]$rowid, ".RDS")))
  Sys.sleep(3)
  
}

# which ones didn't work?
 
missing <- all_tips %>% filter(!rowid %in% as.numeric(str_remove(list.files("data/wsr_ws_backup//"), "\\.RDS$")))

for(i in 1:nrow(missing)){
  
  try(ws[[i]] <- getXYWatersheds(sf = missing[i,], snap = TRUE) %>%
        mutate(WSR_RIVER_ = missing[i,]$WSR_RIVER_) %>%
        saveRDS(paste0("data/wsr_ws_backup/", missing[i,]$rowid, ".RDS")))
  Sys.sleep(3)
  
}

black <- all_tips %>%
  filter(WSR_RIVER_ == "Black, Michigan",
         point_type == "start") 

get



test1 <- list.files("data/wsr_ws_backup/", full.names = TRUE) %>% #[1:500] %>%
  map(~readRDS(.)) %>%
  bind_rows() %>%
  filter(!WSR %in% c("White Clay Wild and Scenic River", "Great Egg Harbor Wild and Scenic River", 
                     "Flathead Wild and Scenic River", "Smith Wild and Scenic River")) %>%
  st_make_valid() %>%
  group_by(WSR, state) %>%
  summarize() %>%
  nngeo::st_remove_holes()

test2 <- list.files("data/wsr_ws_backup/", full.names = TRUE) %>%
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
  dplyr::mutate(WSR = WSR_RIVER_, state = STATE) %>%
  sf::st_transform(crs = 32616) %>% 
  st_buffer(402.336, endCapStyle = "FLAT") %>% # quarter mile around lines
  rowid_to_column() %>%
  st_transform(4269) %>%
  select(WSR, state)



