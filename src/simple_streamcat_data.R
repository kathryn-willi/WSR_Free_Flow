simple_streamcat_data <- function(df){

streamcat_vars <- StreamCatTools::sc_get_data(metric = paste(selected_vars, collapse = ","),
                                              aoi = 'watershed',
                                              showAreaSqKm = TRUE,
                                              comid = df$comid) %>%
  # remove variables we don't particularly care about that get returned:
  select(-contains("areasqkm"))

streamcat_list_simple <- df %>%
  select(index, comid) %>%
  left_join(., streamcat_vars, by = "comid")

  return(streamcat_list_simple)
}
