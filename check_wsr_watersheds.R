# Load previous results if they exist
if (file.exists("data/watershed_review_results.RDS")) {
  results <- readRDS("data/watershed_review_results.RDS")
  cat("Loaded previous results. Continuing from where you left off.\n")
}

library(mapview)
library(dplyr)
library(sf)

# Load data
rivers <- readRDS("data/wsr_official_lines.RDS") %>%
  filter(!grepl("Alaska|Puerto Rico", standard_name, ignore.case = TRUE))

watersheds <- list.files("data/wsr_ws_backup", full.names = TRUE) %>%
  map_dfr(~readRDS(.)) %>%
  group_by(WSR_RIVER_) %>%
  summarize() %>%
  nngeo::st_remove_holes() #%>%
#rename(standard_name = WSR)

# Initialize results dataframe
results <- data.frame(
  standard_name = watersheds$standard_name,
  status = NA_character_,
  stringsAsFactors = FALSE
)

# Function to review watersheds
review_watersheds <- function(watersheds, rivers, results) {
  n <- nrow(watersheds)
  i <- 1
  
  while (i <= n) {
    # Get current watershed
    current_ws <- watersheds[i, ]
    current_name <- current_ws$standard_name
    
    # Get matching rivers
    current_rivers <- rivers %>% 
      filter(standard_name == current_name)
    
    # Create and display map
    cat("\n========================================\n")
    cat(sprintf("Watershed %d of %d: %s\n", i, n, current_name))
    cat("========================================\n")
    
    # Display map
    m <- mapview(current_ws, layer.name = "Watershed") + 
      mapview(current_rivers, color = "blue", layer.name = "Rivers")
    print(m)
    
    # Get user input
    response <- ""
    while (!response %in% c("y", "Y", "n", "N", "q", "Q", "b", "B")) {
      response <- readline(prompt = "Good watershed? (Y/N, B=back, Q=quit): ")
    }
    
    # Process response
    if (response %in% c("q", "Q")) {
      cat("\nQuitting review. Progress saved.\n")
      break
    } else if (response %in% c("b", "B")) {
      if (i > 1) {
        i <- i - 1
        cat("\nGoing back to previous watershed.\n")
      } else {
        cat("\nAlready at first watershed.\n")
      }
    } else {
      results$status[i] <- ifelse(response %in% c("y", "Y"), "Good", "Bad")
      cat(sprintf("Marked as: %s\n", results$status[i]))
      i <- i + 1
    }
  }
  
  return(results)
}

# Run the review
results <- review_watersheds(watersheds, rivers, results)

# Save results
saveRDS(results, "data/watershed_review_results.RDS")
write.csv(results, "data/watershed_review_results.csv", row.names = FALSE)

# Display summary
cat("\n========================================\n")
cat("REVIEW SUMMARY\n")
cat("========================================\n")
cat(sprintf("Total watersheds: %d\n", nrow(results)))
cat(sprintf("Good: %d\n", sum(results$status == "Good", na.rm = TRUE)))
cat(sprintf("Bad: %d\n", sum(results$status == "Bad", na.rm = TRUE)))
cat(sprintf("Not reviewed: %d\n", sum(is.na(results$status))))
cat("\nResults saved to:\n")
cat("  - data/watershed_review_results.RDS\n")
cat("  - data/watershed_review_results.csv\n")