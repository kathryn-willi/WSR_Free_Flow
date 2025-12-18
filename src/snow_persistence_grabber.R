snow_persistence_grabber <- function(watersheds){
  
  # sp_preview <- rast("data/snow_persistence_hammond/MOD10A2_SCI_2020.tif")
  
  # tm_shape(sp_preview) +
  #   tm_raster(palette = "viridis", title = "Snow Persistence") +
  #   tm_layout(frame = FALSE)
  
  # load all the .tif files into a list
  tif_files <- list.files("data/snow_persistence_hammond/", pattern = "\\.tif$", full.names = TRUE)
  
  # stack the snow persistence .tif files into a single raster stack
  raster_stack <- terra::rast(tif_files)
  
  # convert the shapefile to a 'terra' object (if necessary)
  if(crs(watersheds) != crs(raster_stack)){
    polygon <- st_transform(watersheds, crs(raster_stack))  # Align CRS
  } else {
  
  polygon <- watersheds
  
  }
  
  # convert the polygons to 'terra' vector format
  polygon_terra <- vect(polygon)
  
  # mask the raster stack to the watershed polygons
  masked_stack <- mask(raster_stack, polygon_terra)
  
  # extract mean SP across each watershed. weights = TRUE means get the area-weighted average
  mean_sp <- extract(masked_stack, polygon_terra, fun = mean, weights = TRUE)
  
  # calculate mean annual snow persistence for each cell across all years
  masked_mean_annual_sp <- app(masked_stack, fun = mean, na.rm = TRUE)
  
  # create a binary raster where 1 = SP > 60%, 0 = SP <= 60%
  sp_60_mask <- masked_mean_annual_sp > 60
  
  # calculate the fraction of each watershed with SP > 60%
  area60_fraction <- zonal(sp_60_mask, polygon_terra, fun = "mean", weights = TRUE, na.rm = TRUE) %>%
    rename(area60 = mean)
  
  # convert the results to a data frame listing each gage's SP
  watershed_sp <- as_tibble(mean_sp) %>%
    # bind_cols(st_drop_geometry(watersheds)) %>%
    # select(-ID) %>%
    pivot_longer(cols = c(contains("MOD"))) %>%
    group_by(ID) %>% 
    summarize(mean_sp_2001_2020 = mean(value)) %>%
    select(-ID) %>%
    bind_cols(watersheds, .) %>%
    # add area60 to the results
    bind_cols(area60_fraction)
  
  return(watershed_sp)
  
}