# data is a sf data frame I guess
fetch_stream_mods <- function(site_list){
  
  data <- raw_watersheds %>%
    filter(index == site_list)
  
  nhd_plus_hr_url <- "https://hydro.nationalmap.gov/arcgis/rest/services/NHDPlus_HR/MapServer"
  
  # open the nhd_hr - which contains a bunch of layers
  nhd_hr <- arcgislayers::arc_open(nhd_plus_hr_url)
  
  # arcgislayers::list_items(nhd_hr)
  
  # 0: NHDPlusGage 
  # 1: NHDPlusSink
  # 2: NHDPoint
  # 3: NetworkNHDFlowline, this is where you want to start
  # 4: NonNetworkNHDFlowline
  # 5: FlowDirection
  # 6: NHDPlusWall
  # 7: NHDLine
  # 8: NHDArea
  # 9: NHDWaterbody, NOT THIS ONE
  # 10: NHDPlusCatchment, this is like a small watershed that feeds into the big ones
  # 11: NHDPlusBoundaryUnit
  # 12: WBDHU12
  
  nhd_hr_flowlines <- arcgislayers::get_layer(nhd_hr, 3)
  
  
  if (st_geometry_type(data) != "POLYGON") {
    # If not, cast to a Polygon
    data <- st_cast(data, "POLYGON")
  }
  
  # use bbox to return associated flowlines
  geospatial_aoi <- data %>% 
    # st_buffer(0.0000001) %>% 
    # Convert sf object to sfc object (required for downloading from the map server)
    st_as_sfc(.)
  
  nhd_flowlines <- NULL
  
  nhd_flowlines <- vector("list", length = length(geospatial_aoi)) #data frame of each water shed?
  
  # fields <- arcgislayers::list_fields(nhd_hr_flowlines)
  
  # TODO: update the query for only artificial stuff
  # query <- "Ftype IN (428, 420, 336, 468, 558)"
  
  # Here is where the intersection happens, which is really just a filter based on 
  # the shape of the watersheds. The flow lines seem to be cut off at the edge of the polygons,
  # because it is a filter instead of an intersection.
  for(i in 1:length(geospatial_aoi)){
    tryCatch({
      nhd_flowlines[[i]] <- arcgislayers::arc_select(nhd_hr_flowlines,
                                                     # where = query,
                                                     filter_geom = geospatial_aoi[i],
                                                     crs = st_crs(geospatial_aoi[i])) %>% 
        st_make_valid()},
      error = function(e){
        cat("Index ", i, " from input data failed.")
      }) 
  }
  
  nhd_flowlines <- nhd_flowlines %>% 
    keep(~!is.null(.))
  
  try(nhd_flowlines <- nhd_flowlines %>% 
        dplyr::bind_rows() %>%
        dplyr::distinct() %>%
        mutate(#natural = ifelse(ftype == 460, T, F),
          flowline_type = case_when(ftype == 460 ~ "natural",
                                    ftype == 558 ~ "artificial path",
                                    ftype == 468 ~ "drainageway",
                                    ftype == 336 ~ "canal ditch",
                                    ftype == 566 ~ "coastline",
                                    ftype == 334 ~ "connector",
                                    ftype == 428 ~ "pipeline",
                                    ftype == 420 ~ "underground conduit",
                                    .default = "natural")) %>%
        filter(!flowline_type %in% c("natural", "coastline", "connector", "artificial path")),
      silent = TRUE)
  
  if (is.data.frame(nhd_flowlines) && nrow(nhd_flowlines) > 0) {
    saveRDS(nhd_flowlines, paste0("data/raw_modifications/", data$index, ".RDS"))
  } else {
    message("Skipping saveRDS: nhd_flowlines is either not a data frame or is empty.")
  }
  
  
  # transform polygon to polyline and find those intersections instead
  return(nhd_flowlines)
  
}