transbasin_finder <- function(site_list){
  
  site <- raw_watersheds %>%
    filter(index == site_list)
  
  flowlines_unnatural <- NA
  
  flowlines_unnatural <- readRDS(paste0("data/raw_modifications/", site$index, ".RDS"))
  
  if(!is.data.frame(flowlines_unnatural)){
    
    return(site %>% mutate(transbasin_diversion = NA) %>% select(index, transbasin_diversion) %>% st_drop_geometry())
    
  }
  
  # For linestring transformation step to work, need the watershed to be a polygon object:
  if (st_geometry_type(site) != "POLYGON") {
    # If not, cast to a Polygon... which will "explode" it into multiple. 
    # This is a rare thing... I think...
    site <- st_cast(site, "POLYGON") 
  }
  
  polyline <-  site %>% st_cast("LINESTRING")
  
  crossovers <- flowlines_unnatural %>%
    .[polyline,] %>%
    nrow()
  
  site <- site %>%
    mutate(transbasin_diversion = ifelse(crossovers > 0, "transbasin diversion", NA)) %>%
    select(index, transbasin_diversion) %>%
    st_drop_geometry()
  
  return(site)
  
}