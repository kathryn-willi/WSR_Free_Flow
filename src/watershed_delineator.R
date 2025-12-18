# function to delineate each gage's watershed:
watershed_delineator <- function(site_list){
  
  # Filter to get a specific site from the dataset and transform to WGS84 coordinate system
  site <- all_sites %>%
    filter(index == site_list) %>%
    st_transform(4326)
  
  # these are watersheds that are contained within a catchment (very small). require
  # a different workflow:
  if (site$combo %in% c("SPACRECO", "10172200", "08329700", "SCHFLMCO")){
    
    # Get the NHDPlus flowline near the site and create a 30m buffer around it
    # This is considered a "danger zone" to avoid placing points on the flowline
    flowline_danger_zone <- get_nhdplus(AOI = site, realization = "flowline") %>%
      st_buffer(30)
    
    # Create a 35m buffer around the site location
    site_buffer <- site %>%
      st_buffer(35) 
    
    # Extract points from the buffer boundary while avoiding the flowline
    boundary_points <- site_buffer %>%
      st_boundary() %>%      # Extract the boundary (a LINESTRING)
      st_cast("POINT") %>%   # Convert the LINESTRING to individual POINTs
      mutate(point_id = row_number()) %>%  # Add an ID for each point
      st_difference(., flowline_danger_zone) %>%  # Remove points that fall within the flowline buffer
      # Reduce the number of points by taking only every 50th point
      filter(point_id %% 50 == 0) %>%
      bind_rows(site)  # Add the original site point to the collection
    
    # Initialize an empty list to store catchment splits
    splits <- vector("list")
    
    # Loop through each boundary point to get split catchments
    for(i in 1:nrow(boundary_points)){
      
      # Get the split catchment for each point
      # We're looking for areas without a catchmentID (NA values)
      splits[[i]] <- get_split_catchment(point = st_as_sfc(boundary_points[i,])) %>%
        filter(is.na(catchmentID)) %>%
        mutate(area = as.numeric(st_area(.)))  # Calculate the area
    }
    
    # Combine all splits and select the largest one
    splits <- bind_rows(splits) %>%
      filter(as.numeric(area) == max(as.numeric(area))) %>%  # Keep only the largest area
      select(-catchmentID,-id) %>%  # Remove unnecessary columns
      mutate(ws_flag = "SMALL") %>%  # Add a flag indicating this is a small watershed
      # Add site information and initialize area metrics as NA
      dplyr::mutate(index = site$index,
                    combo = site$combo,
                    comid = site$comid,
                    area_original_m2 = NA,
                    area_updated_m2 = NA,
                    area_dif_m2 = NA)
    # back it up:
    saveRDS(splits, paste0("data/raw_watersheds/", site$combo, ".RDS"))
    
    message(paste0(site$station_nm, " delineated!"))
    
    return(splits)
    
  } 
    
    # "Normal" approach for delineating watersheds:
    
    raindrop <- get_raindrop_trace(point = st_as_sfc(site), direction = "down")
    
    flowline <- get_nhdplus(comid = site$comid)
    
    if(raindrop[1,3]$comid == site$comid) {
      
      nearest_point <- st_as_sf(data.frame(x = raindrop$intersection_point[[1]][1],
                                           y = raindrop$intersection_point[[1]][2]),
                                coords = c("x", "y"),
                                crs = st_crs(raindrop)) %>%
        st_as_sfc()
      
      better_termination <- get_split_catchment(nearest_point, upstream = F) %>%
        filter(is.na(catchmentID)) %>%
        mutate(featureid = site$comid)
      
     # if the raindrop trace says that the gage is actually along a different
     # flowline...
      
    } else {
      
      message("One of the funky ones")
      
      nearest_points <- st_nearest_points(site, flowline, pairwise = T) %>%
        st_cast(., "POINT") %>%
        st_as_sf() %>%
        st_make_valid()  %>%
        mutate(distance = as.numeric(st_distance(., site))) %>%
        distinct() %>%
        arrange(distance) 
      
      nearest_point_1 <- nearest_points %>%
        .[1,] %>%
        st_as_sfc()
      
      nearest_point_2 <- nearest_points %>%
        .[2,] %>%
        st_as_sfc()
      
      better_termination_1 <- get_split_catchment(nearest_point_1, upstream = F) %>%
        filter(is.na(catchmentID)) 
      
      better_termination <- get_split_catchment(nearest_point_2, upstream = F) %>%
        filter(is.na(catchmentID)) %>%
        bind_rows(better_termination_1) %>%
        summarize() %>%
        mutate(featureid = as.numeric(flowline$comid))
      
    }
      
      # Weird ones that require even more work:
      if(site$combo == "VALBAYCO"){
        better_termination <- better_termination %>%
          mutate(featureid = 17034207)
      } else if(site$combo == "LVNGEOCO"){
        better_termination <- better_termination %>%
          mutate(featureid = 2885290)
    }
    
    # once termination point is figured out... start delineating!
    
    # use get_UT to list all comids that are upstream of our gage using the comid the
    # gage falls on:
    upstream <- nhdplusTools::get_UT(nhd, better_termination$featureid)
    
    # These ones are weird braided channels that the NHD doesn't
    # work super well with:
    if(site$combo == "NORLASCO"){
      more_upstream <- get_UT(nhd, "17876457")
      upstream <- c(upstream, more_upstream)  
    } else if(site$combo == "NORCONCO"){
      more_upstream <- get_UT(nhd, "17875505")
      upstream <- c(upstream, more_upstream)  
    } else if(site$combo == "TOMGUNCO"){
      more_upstream <- get_UT(nhd, "9772879")
      upstream <- c(upstream, more_upstream)  
    }
    
    
    # grab all the catchments associated with the upstream comids:
    nhd_catch <- nhdplusTools::get_nhdplus(comid = upstream,
                                           realization = 'catchment') 
    # another weird one where the NHD is broken and the catchment comid and nhd comid
    # aren't the same:
    if(site$combo == "06799445"){
      nhd_catch <- nhdplusTools::get_nhdplus(AOI = site,
                                             realization = 'catchment') %>%
        mutate(featureid = flowline$comid) %>%
        bind_rows(nhd_catch)
    }
    
    site_catch <- nhd_catch %>%
      filter(featureid == better_termination$featureid) %>%
      sf::st_area()
    
    split_catch <- better_termination %>%
      sf::st_area()
    
    # Weirdos:
    if(site$combo %in%c("WFKMOUCO","NORCONCO")){
      split_catch <- nhd_catch %>%
        filter(featureid == better_termination$featureid) %>%
        sf::st_area()
      better_termination <- nhd_catch %>%
        filter(featureid == better_termination$featureid)
    }
    
    dif <- as.numeric(site_catch) - as.numeric(split_catch)
    
    if(raindrop[1,3]$comid != site$comid && dif < 0) {
      split_catch <- 0 
      dif <- as.numeric(site_catch) - as.numeric(split_catch)
    }
    
    if(raindrop[1,3]$comid == site$comid && dif < 0) {
      split_catch <- site_catch 
      dif <- as.numeric(site_catch) - as.numeric(split_catch)
    }
    
    nhd_watershed <- nhd_catch %>%
      # remove dupes (precautionary step, not likely necessary)
      dplyr::distinct(featureid, .keep_all=TRUE) %>%
      dplyr::filter(featureid != better_termination$featureid) %>%
      dplyr::bind_rows(., better_termination) %>%
      st_make_valid() %>%
      # "dissolve" all the catchments into a single polygon
      dplyr::summarize() %>%
      # remove weird hole by-products that exist if the catchment boundaries don't
      # line up perfectly:
      nngeo::st_remove_holes() %>%
      # tack on the site name and comid to the watershed
      dplyr::mutate(index = site$index,
                    combo = site$combo,
                    comid = site$comid,
                    area_original_m2 = as.numeric(site_catch),
                    area_updated_m2 = as.numeric(split_catch),
                    area_dif_m2 = dif)
    
    # mapview(nhd_watershed) + site + flowline
    
    # back it up:
    saveRDS(nhd_watershed, paste0("data/raw_watersheds/", site$index, ".RDS"))
    
    message(paste0(site$station_nm, " delineated!"))
    
    return(nhd_watershed)
  }