# =============================================================================
# FUNC_worms_lifestyle.R
# =============================================================================
#
# Query WoRMS functional groups for species and add a `lifestyle` column.
# Results are cached to CSV so repeated runs don't re-query the API.
#
# FUNCTIONS:
#
#   query_worms_functional_groups(species_list, cache_file, delay)
#     Query WoRMS API for functional group traits for a vector of species.
#     Caches results to CSV. Subsequent calls skip already-queried species.
#
#   add_lifestyle_from_worms(site_level, worms_traits)
#     Join WoRMS functional groups to site_level data and map to
#     standardised lifestyle categories.
#
# USAGE:
#   source("Scripts/Final/FUNC_worms_lifestyle.R")
#
#   # Query all species (runs once, ~30 min for 1000 species)
#   all_species <- unique(na.omit(combined_all$site_level$Species))
#   worms_traits <- query_worms_functional_groups(
#     all_species, 
#     cache_file = "Processed_data/worms_functional_groups.csv"
#   )
#
#   # Add lifestyle column
#   combined_all$site_level <- add_lifestyle_from_worms(
#     combined_all$site_level, worms_traits
#   )
#
# =============================================================================

library(worrms)
library(dplyr)
library(tidyr)
library(purrr)


#' Query WoRMS functional group attributes for a list of species
#'
#' For each species, resolves the AphiaID via wm_name2id(), then retrieves
#' all attributes via wm_attr_data() and extracts the "Functional group"
#' trait. Results are cached to a CSV file. On subsequent runs, only new
#' (uncached) species are queried.
#'
#' @param species_list Character vector of species names
#' @param cache_file Path to CSV cache file. If it exists, cached results
#'   are loaded and only missing species are queried. Set NULL to disable.
#' @param delay Seconds between API calls (default 0.3 to respect rate limits)
#' @return Data frame with columns: Species, AphiaID, functional_group
#'
#' @examples
#' # First run (queries API):
#' traits <- query_worms_functional_groups(
#'   c("Clupea harengus", "Crassostrea gigas", "Mnemiopsis leidyi"),
#'   cache_file = "Processed_data/worms_functional_groups.csv"
#' )
#'
#' # Second run (loads from cache, only queries new species):
#' traits <- query_worms_functional_groups(
#'   c("Clupea harengus", "Crassostrea gigas", "Mytilus edulis"),
#'   cache_file = "Processed_data/worms_functional_groups.csv"
#' )
query_worms_functional_groups <- function(
    species_list,
    cache_file = "Processed_data/worms_functional_groups.csv",
    delay = 0.3
) {
  
  species_list <- unique(na.omit(species_list))
  
  # Load cache if it exists
  cached <- NULL
  if (!is.null(cache_file) && file.exists(cache_file)) {
    cached <- read.csv(cache_file, stringsAsFactors = FALSE)
    cat("Loaded", nrow(cached), "cached results from", cache_file, "\n")
    
    # Only query species not yet cached
    to_query <- setdiff(species_list, cached$Species)
    cat("Already cached:", length(species_list) - length(to_query), "\n")
    cat("New species to query:", length(to_query), "\n")
  } else {
    to_query <- species_list
    cat("No cache found. Querying all", length(to_query), "species\n")
  }
  
  if (length(to_query) == 0) {
    cat("All species already cached.\n")
    return(cached %>% filter(Species %in% species_list))
  }
  
  cat("Estimated time:", round(length(to_query) * delay * 2 / 60, 1), "minutes\n\n")
  
  results <- list()
  
  for (i in seq_along(to_query)) {
    sp <- to_query[i]
    cat(i, "/", length(to_query), ":", sp, "...")
    
    # Resolve AphiaID
    aphia_id <- tryCatch(
      wm_name2id(sp),
      error = function(e) NA_integer_
    )
    
    if (is.na(aphia_id)) {
      cat(" not found in WoRMS\n")
      results[[i]] <- tibble(
        Species = sp, AphiaID = NA_integer_, functional_group = NA_character_
      )
      Sys.sleep(delay)
      next
    }
    
    Sys.sleep(delay)
    
    # Get attributes
    attrs <- tryCatch(
      wm_attr_data(id = aphia_id),
      error = function(e) NULL
    )
    
    if (is.null(attrs) || nrow(attrs) == 0) {
      cat(" AphiaID", aphia_id, "- no attributes\n")
      results[[i]] <- tibble(
        Species = sp, AphiaID = aphia_id, functional_group = NA_character_
      )
      Sys.sleep(delay)
      next
    }
    
    # Extract functional group
    fg <- attrs %>%
      filter(grepl("unctional", measurementType, ignore.case = TRUE)) %>%
      pull(measurementValue)
    
    if (length(fg) == 0) {
      cat(" AphiaID", aphia_id, "- no functional group\n")
      fg_str <- NA_character_
    } else {
      fg_str <- paste(unique(fg), collapse = "; ")
      cat(" ", fg_str, "\n")
    }
    
    results[[i]] <- tibble(
      Species = sp, AphiaID = aphia_id, functional_group = fg_str
    )
    
    Sys.sleep(delay)
  }
  
  new_results <- bind_rows(results)
  
  # Combine with cache and save
  if (!is.null(cached)) {
    all_results <- bind_rows(cached, new_results)
  } else {
    all_results <- new_results
  }
  
  if (!is.null(cache_file)) {
    dir.create(dirname(cache_file), showWarnings = FALSE, recursive = TRUE)
    write.csv(all_results, cache_file, row.names = FALSE)
    cat("\nSaved", nrow(all_results), "results to", cache_file, "\n")
  }
  
  all_results %>% filter(Species %in% species_list)
}


#' Map WoRMS functional groups to standardised lifestyle categories
#'
#' Takes the raw functional_group strings from WoRMS and maps them to:
#'   Fish, Sessile invertebrate, Mobile invertebrate, Phytoplankton,
#'   Zooplankton, Macroalgae, Microeukaryote, Other
#'
#' @param worms_traits Output from query_worms_functional_groups()
#' @param custom_mapping Named character vector mapping WoRMS functional
#'   group strings to lifestyle categories. Default handles common groups.
#'   Override or extend as needed.
#' @return Input data frame with added `lifestyle` column
#'
#' @examples
#' # Default mapping:
#' traits <- classify_worms_lifestyle(worms_traits)
#'
#' # Custom mapping — reclassify jellyfish:
#' traits <- classify_worms_lifestyle(worms_traits, custom_mapping = c(
#'   "benthos" = "Sessile invertebrate",
#'   "nekton" = "Fish",
#'   "zooplankton" = "Zooplankton",
#'   "phytoplankton" = "Phytoplankton",
#'   "macroalgae" = "Macroalgae",
#'   "jellyfish" = "Zooplankton"
#' ))
classify_worms_lifestyle <- function(
    worms_traits,
    custom_mapping = c(
      "benthos"       = "Sessile invertebrate",
      "nekton"        = "Fish",
      "zooplankton"   = "Zooplankton",
      "phytoplankton" = "Phytoplankton",
      "macroalgae"    = "Macroalgae",
      "algae"         = "Macroalgae",
      "plants"        = "Other",
      "birds"         = "Other",
      "mammals"       = "Other",
      "reptiles"      = "Other",
      "fish"          = "Fish"
    )
) {
  
  worms_traits %>%
    mutate(
      # Take first functional group if multiple
      fg_primary = str_extract(functional_group, "^[^;]+"),
      fg_primary = str_trim(tolower(fg_primary)),
      lifestyle = custom_mapping[fg_primary],
      lifestyle = ifelse(is.na(lifestyle) & !is.na(functional_group), 
                         "Other", lifestyle)
    ) %>%
    select(-fg_primary)
}


#' Add WoRMS-derived lifestyle to site_level data
#'
#' Joins the classified WoRMS traits to the site_level data frame.
#'
#' @param site_level Data frame from combined_all$site_level
#' @param worms_traits Output from classify_worms_lifestyle()
#' @return site_level with `lifestyle` and `functional_group` columns added
add_lifestyle_from_worms <- function(site_level, worms_traits) {
  
  site_level %>%
    left_join(
      worms_traits %>% select(Species, functional_group, lifestyle),
      by = "Species"
    )
}


# =============================================================================
# LOADED
# =============================================================================

cat("FUNC_worms_lifestyle.R loaded\n")
cat("Functions:\n")
cat("  query_worms_functional_groups(species_list, cache_file, delay)\n")
cat("  classify_worms_lifestyle(worms_traits, custom_mapping)\n")
cat("  add_lifestyle_from_worms(site_level, worms_traits)\n")
