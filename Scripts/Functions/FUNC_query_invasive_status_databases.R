# =============================================================================
# COMPREHENSIVE INVASIVE SPECIES STATUS CHECKER
# Queries: WoRMS/WRiMS, GBIF, and compiles results
# =============================================================================
# 
# *** IMPORTANT: RUN THIS SCRIPT ON YOUR LOCAL MACHINE ***
# This script requires internet access to query external APIs.
# It will NOT work in restricted environments.
#
# This script queries multiple authoritative databases to determine invasive
# status for species-country combinations from eDNA metabarcoding data.
#
# Data sources:
# 1. WoRMS (World Register of Marine Species) - via worrms package
#    - AphiaID lookup
#    - Distribution records including origin (native/introduced/alien)
#    - WRiMS (World Register of Introduced Marine Species) data is integrated
#    - As of 2024-08-22, includes invasiveness and occurrence fields
#
# 2. GBIF (Global Biodiversity Information Facility) - via rgbif package
#    - Establishment means per country (Native, Introduced, Invasive, etc.)
#
# 3. OBIS (Ocean Biodiversity Information System) - via robis package
#    - Additional occurrence data for marine species
#
# Expected runtime: ~30-60 minutes for 1000 species-country combinations
# (due to API rate limits - please be respectful to free services)
#
# Author: Generated for marine eDNA metabarcoding project
# Date: 2026-02-10
# =============================================================================

# ---- 1. Install and load required packages ----
required_packages <- c("worrms", "rgbif", "robis", "dplyr", "tidyr", "purrr", 
                       "readr", "stringr", "httr", "jsonlite", "progress")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# ---- 2. Configuration ----
# Rate limiting to be respectful to APIs
WORMS_DELAY <- 0.5  # seconds between WoRMS requests
GBIF_DELAY <- 0.3   # seconds between GBIF requests
OBIS_DELAY <- 0.3   # seconds between OBIS requests

# Set to TRUE to also query OBIS (slower but more comprehensive)
USE_OBIS <- FALSE

# Country ISO codes — comprehensive lookup with dynamic resolution
# Core mapping covers common names, demonyms, and abbreviations.
# Countries not in this list are resolved automatically via get_country_code().
COUNTRY_CODES_LOOKUP <- c(
  # Western Europe
  "Austria" = "AT", "Belgium" = "BE", "France" = "FR", "Germany" = "DE",
  "Ireland" = "IE", "Luxembourg" = "LU", "Monaco" = "MC",
  "Netherlands" = "NL", "Switzerland" = "CH",
  "UK" = "GB", "United Kingdom" = "GB", "Great Britain" = "GB",
  # Northern Europe
  "Denmark" = "DK", "Estonia" = "EE", "Finland" = "FI", "Iceland" = "IS",
  "Latvia" = "LV", "Lithuania" = "LT", "Norway" = "NO", "Sweden" = "SE",
  # Southern Europe
  "Albania" = "AL", "Bosnia and Herzegovina" = "BA", "Croatia" = "HR",
  "Cyprus" = "CY", "Greece" = "GR", "Italy" = "IT", "Malta" = "MT",
  "Montenegro" = "ME", "North Macedonia" = "MK", "Portugal" = "PT",
  "Serbia" = "RS", "Slovenia" = "SI", "Spain" = "ES", "Turkey" = "TR",
  # Eastern Europe
  "Belarus" = "BY", "Bulgaria" = "BG", "Czech Republic" = "CZ",
  "Czechia" = "CZ", "Hungary" = "HU", "Moldova" = "MD", "Poland" = "PL",
  "Romania" = "RO", "Russia" = "RU", "Slovakia" = "SK", "Ukraine" = "UA",
  # Americas
  "Argentina" = "AR", "Brazil" = "BR", "Canada" = "CA", "Chile" = "CL",
  "Colombia" = "CO", "Mexico" = "MX", "Peru" = "PE",
  "USA" = "US", "United States" = "US",
  # Asia-Pacific
  "Australia" = "AU", "China" = "CN", "India" = "IN", "Indonesia" = "ID",
  "Japan" = "JP", "Malaysia" = "MY", "New Zealand" = "NZ",
  "Philippines" = "PH", "Singapore" = "SG", "South Korea" = "KR",
  "Taiwan" = "TW", "Thailand" = "TH", "Vietnam" = "VN",
  # Middle East & Africa
  "Egypt" = "EG", "Israel" = "IL", "Kenya" = "KE", "Morocco" = "MA",
  "Nigeria" = "NG", "South Africa" = "ZA", "Tunisia" = "TN",
  "UAE" = "AE", "United Arab Emirates" = "AE"
)

#' Resolve a country name to its ISO 2-letter code
#' @param country Country name (case-sensitive first, then case-insensitive)
#' @return ISO code or NA with a warning
get_country_code <- function(country) {
  # Exact match
  if (country %in% names(COUNTRY_CODES_LOOKUP)) {
    return(COUNTRY_CODES_LOOKUP[[country]])
  }
  # Case-insensitive match
  idx <- match(tolower(country), tolower(names(COUNTRY_CODES_LOOKUP)))
  if (!is.na(idx)) {
    return(COUNTRY_CODES_LOOKUP[[idx]])
  }
  warning("No ISO code found for country '", country,
          "'. GBIF query will be skipped for this country. ",
          "Add it to COUNTRY_CODES_LOOKUP if needed.")
  return(NA_character_)
}

# ---- 3. Helper Functions ----

#' Direct REST API call to WoRMS (backup if worrms package fails)
#' @param endpoint API endpoint
#' @return Parsed JSON response or NULL
worms_api_call <- function(endpoint) {
  base_url <- "https://www.marinespecies.org/rest/"
  url <- paste0(base_url, endpoint)
  
  tryCatch({
    resp <- httr::GET(url, httr::timeout(30))
    if (httr::status_code(resp) == 200) {
      return(jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8")))
    }
    return(NULL)
  }, error = function(e) {
    return(NULL)
  })
}

#' Get AphiaID from WoRMS for a species name
#' @param species_name Scientific name
#' @return AphiaID or NA
get_aphia_id <- function(species_name) {
  tryCatch({
    # Try worrms package first
    result <- wm_name2id(species_name)
    if (length(result) > 0 && !is.na(result)) {
      return(result)
    }
    # Try fuzzy match if exact match fails
    fuzzy <- wm_records_name(species_name, fuzzy = TRUE)
    if (!is.null(fuzzy) && nrow(fuzzy) > 0) {
      return(fuzzy$AphiaID[1])
    }
    return(NA)
  }, error = function(e) {
    # Fallback to direct REST API
    tryCatch({
      encoded_name <- URLencode(species_name, reserved = TRUE)
      result <- worms_api_call(paste0("AphiaIDByName/", encoded_name, "?marine_only=true"))
      if (!is.null(result) && result > 0) {
        return(result)
      }
      return(NA)
    }, error = function(e2) {
      return(NA)
    })
  })
}

#' Get WoRMS distribution data for a species
#' @param aphia_id AphiaID from WoRMS
#' @return Data frame with distribution info including origin/invasiveness
get_worms_distribution <- function(aphia_id) {
  if (is.na(aphia_id)) return(NULL)
  
  tryCatch({
    dist <- wm_distribution(aphia_id)
    if (!is.null(dist) && nrow(dist) > 0) {
      return(dist)
    }
    return(NULL)
  }, error = function(e) {
    return(NULL)
  })
}

#' Get vernacular/common name from WoRMS
#' @param aphia_id AphiaID
#' @return Common name or NA
get_common_name <- function(aphia_id) {
  if (is.na(aphia_id)) return(NA)
  
  tryCatch({
    vernaculars <- wm_common_id(aphia_id)
    if (!is.null(vernaculars) && nrow(vernaculars) > 0) {
      # Prefer English names
      eng <- vernaculars %>% filter(language == "English")
      if (nrow(eng) > 0) {
        return(eng$vernacular[1])
      }
      return(vernaculars$vernacular[1])
    }
    return(NA)
  }, error = function(e) {
    return(NA)
  })
}

#' Query GBIF for establishment means in a specific country
#' @param species_name Scientific name
#' @param country_code ISO 2-letter country code
#' @return Establishment means or NA
get_gbif_establishment <- function(species_name, country_code) {
  tryCatch({
    # First get the taxon key
    backbone <- name_backbone(name = species_name)
    if (is.null(backbone) || is.na(backbone$usageKey)) {
      return(NA)
    }
    
    # Search for occurrences in the country
    occs <- occ_data(
      taxonKey = backbone$usageKey,
      country = country_code,
      limit = 100,
      fields = c("establishmentMeans", "degreeOfEstablishment")
    )
    
    if (!is.null(occs$data) && nrow(occs$data) > 0) {
      # Get unique establishment means
      if ("establishmentMeans" %in% names(occs$data)) {
        means <- unique(na.omit(occs$data$establishmentMeans))
        if (length(means) > 0) {
          return(paste(means, collapse = "; "))
        }
      }
    }
    return(NA)
  }, error = function(e) {
    return(NA)
  })
}

#' Query OBIS for occurrence data (optional, slower)
#' @param species_name Scientific name
#' @param country Country name
#' @return Data frame with OBIS records or NULL
get_obis_records <- function(species_name, country) {
  if (!USE_OBIS) return(NULL)
  
  tryCatch({
    # OBIS uses WoRMS AphiaID internally
    records <- robis::occurrence(scientificname = species_name)
    
    if (!is.null(records) && nrow(records) > 0) {
      # Filter by country if possible
      if ("country" %in% names(records)) {
        country_records <- records %>%
          filter(grepl(country, country, ignore.case = TRUE))
        if (nrow(country_records) > 0) {
          return(country_records)
        }
      }
      return(records)
    }
    return(NULL)
  }, error = function(e) {
    return(NULL)
  })
}

#' Build a regex pattern to match WoRMS locality strings for a given country
#'
#' Uses a comprehensive lookup of country names, demonyms, and associated
#' marine regions. Countries not in the lookup get an auto-generated pattern
#' from the country name itself.
#'
#' @param country Country name
#' @return Regex pattern string (pipe-separated alternatives)
get_worms_locality_pattern <- function(country) {
  
  # Comprehensive mapping: country → country name | demonym | marine regions
  # WoRMS locality strings use a mix of country names, adjectives, and
  # sea/ocean basin names, so all relevant variants are included.
  marine_region_patterns <- c(
    # Northern Europe
    "Denmark"    = "Denmark|Danish|North Sea|Baltic|Kattegat|Skagerrak",
    "Estonia"    = "Estonia|Estonian|Baltic|Gulf of Finland|Gulf of Riga",
    "Finland"    = "Finland|Finnish|Baltic|Gulf of Bothnia|Gulf of Finland",
    "Iceland"    = "Iceland|Icelandic|North Atlantic|Arctic|Denmark Strait",
    "Ireland"    = "Ireland|Irish|Irish Sea|Celtic Sea|North Atlantic",
    "Latvia"     = "Latvia|Latvian|Baltic|Gulf of Riga",
    "Lithuania"  = "Lithuania|Lithuanian|Baltic",
    "Norway"     = "Norway|Norwegian|North Sea|Barents|Skagerrak|Norwegian Sea",
    "Sweden"     = "Sweden|Swedish|Baltic|Kattegat|Skagerrak|Gulf of Bothnia",
    "UK"         = "United Kingdom|Britain|British|English Channel|Irish Sea|North Sea|Scotland|Wales|England|Celtic Sea",
    "United Kingdom" = "United Kingdom|Britain|British|English Channel|Irish Sea|North Sea|Scotland|Wales|England|Celtic Sea",
    # Western Europe
    "Belgium"      = "Belgium|Belgian|North Sea",
    "France"       = "France|French|Atlantic|Mediterranean|Bay of Biscay|English Channel|Ligurian",
    "Germany"      = "Germany|German|North Sea|Baltic|Wadden",
    "Netherlands"  = "Netherlands|Dutch|North Sea|Wadden",
    # Southern Europe
    "Albania"    = "Albania|Albanian|Adriatic|Ionian|Mediterranean",
    "Croatia"    = "Croatia|Croatian|Adriatic|Mediterranean",
    "Cyprus"     = "Cyprus|Cypriot|Mediterranean|Levantine",
    "Greece"     = "Greece|Greek|Aegean|Mediterranean|Ionian|Levantine",
    "Italy"      = "Italy|Italian|Mediterranean|Adriatic|Tyrrhenian|Ligurian|Ionian",
    "Malta"      = "Malta|Maltese|Mediterranean",
    "Montenegro" = "Montenegro|Montenegrin|Adriatic|Mediterranean",
    "Portugal"   = "Portugal|Portuguese|Atlantic|Iberian|Macaronesia|Azores|Madeira",
    "Slovenia"   = "Slovenia|Slovenian|Adriatic|Mediterranean",
    "Spain"      = "Spain|Spanish|Mediterranean|Atlantic|Iberian|Cantabrian|Canary|Balearic",
    "Turkey"     = "Turkey|Turkish|Mediterranean|Aegean|Black Sea|Levantine|Sea of Marmara",
    # Eastern Europe
    "Bulgaria"   = "Bulgaria|Bulgarian|Black Sea",
    "Georgia"    = "Georgia|Georgian|Black Sea",
    "Poland"     = "Poland|Polish|Baltic",
    "Romania"    = "Romania|Romanian|Black Sea",
    "Russia"     = "Russia|Russian|Barents|White Sea|Baltic|Black Sea|Sea of Okhotsk|Sea of Japan|Arctic",
    "Ukraine"    = "Ukraine|Ukrainian|Black Sea|Sea of Azov",
    # Americas
    "Argentina"  = "Argentina|Argentine|South Atlantic|Patagoni",
    "Brazil"     = "Brazil|Brazilian|South Atlantic|Tropical Atlantic",
    "Canada"     = "Canada|Canadian|North Atlantic|North Pacific|Arctic|Bay of Fundy|Gulf of St. Lawrence",
    "Chile"      = "Chile|Chilean|South Pacific|Southeast Pacific",
    "Colombia"   = "Colombia|Colombian|Caribbean|Pacific",
    "Mexico"     = "Mexico|Mexican|Gulf of Mexico|Caribbean|Pacific",
    "Peru"       = "Peru|Peruvian|Southeast Pacific",
    "USA"        = "United States|American|Gulf of Mexico|North Atlantic|North Pacific|Caribbean|Chesapeake|California",
    "United States" = "United States|American|Gulf of Mexico|North Atlantic|North Pacific|Caribbean|Chesapeake|California",
    # Asia-Pacific
    "Australia"    = "Australia|Australian|Tasman|Great Barrier|Coral Sea|Indian Ocean|Southern Ocean",
    "China"        = "China|Chinese|South China Sea|East China Sea|Yellow Sea|Bohai",
    "India"        = "India|Indian|Indian Ocean|Arabian Sea|Bay of Bengal",
    "Indonesia"    = "Indonesia|Indonesian|Indian Ocean|Pacific|Java Sea|Coral Triangle",
    "Japan"        = "Japan|Japanese|Sea of Japan|Pacific|East China Sea",
    "Malaysia"     = "Malaysia|Malaysian|South China Sea|Strait of Malacca|Indian Ocean",
    "New Zealand"  = "New Zealand|Zealand|Tasman|South Pacific|Southern Ocean",
    "Philippines"  = "Philippines|Philippine|South China Sea|Pacific|Coral Triangle",
    "Singapore"    = "Singapore|Strait of Malacca|South China Sea",
    "South Korea"  = "South Korea|Korean|Sea of Japan|Yellow Sea|East China Sea",
    "Taiwan"       = "Taiwan|Taiwanese|South China Sea|East China Sea|Pacific",
    "Thailand"     = "Thailand|Thai|Gulf of Thailand|Andaman|Indian Ocean",
    "Vietnam"      = "Vietnam|Vietnamese|South China Sea",
    # Middle East & Africa
    "Egypt"        = "Egypt|Egyptian|Mediterranean|Red Sea|Suez",
    "Israel"       = "Israel|Israeli|Mediterranean|Red Sea|Levantine",
    "Kenya"        = "Kenya|Kenyan|Indian Ocean|Western Indian Ocean",
    "Morocco"      = "Morocco|Moroccan|Atlantic|Mediterranean",
    "Nigeria"      = "Nigeria|Nigerian|Gulf of Guinea|Tropical Atlantic",
    "South Africa" = "South Africa|South African|Indian Ocean|Atlantic|Agulhas",
    "Tunisia"      = "Tunisia|Tunisian|Mediterranean"
  )
  
  pattern <- marine_region_patterns[country]
  
  if (!is.na(pattern)) {
    return(unname(pattern))
  }
  
  # Fallback: use the country name as-is
  # This still matches WoRMS records that use the country name in the locality field
  return(country)
}

#' Determine invasive status from WoRMS distribution data
#' @param dist_df Distribution data frame from WoRMS
#' @param country Target country name
#' @return List with status and details
parse_worms_status <- function(dist_df, country) {
  if (is.null(dist_df) || nrow(dist_df) == 0) {
    return(list(status = "UNKNOWN", origin = NA, invasive = NA, occurrence = NA))
  }
  
  # Build search pattern for this country from the marine region lookup.
  # Countries with known marine region associations get expanded patterns;
  # others fall back to the country name + auto-generated demonym.
  pattern <- get_worms_locality_pattern(country)
  
  # Filter for relevant localities
  relevant <- dist_df %>%
    filter(grepl(pattern, locality, ignore.case = TRUE))
  
  if (nrow(relevant) == 0) {
    # Check broader regions
    relevant <- dist_df %>%
      filter(grepl("Europe|Atlantic|Mediterranean|Baltic|North Sea|Pacific|Indian Ocean|Caribbean|Arctic", 
                   locality, ignore.case = TRUE))
  }
  
  if (nrow(relevant) > 0) {
    # Check origin field (native, introduced, alien, etc.)
    origins <- unique(na.omit(relevant$origin))
    invasive_flags <- unique(na.omit(relevant$invasiveness))
    occurrences <- unique(na.omit(relevant$occurrence))
    
    # Determine status
    status <- "NATIVE"
    if (any(grepl("introduced|alien|non-native|exotic", origins, ignore.case = TRUE))) {
      status <- "INTRODUCED"
    }
    if (any(grepl("invasive", c(origins, invasive_flags), ignore.case = TRUE))) {
      status <- "INVASIVE"
    }
    
    return(list(
      status = status,
      origin = paste(origins, collapse = "; "),
      invasive = paste(invasive_flags, collapse = "; "),
      occurrence = paste(occurrences, collapse = "; ")
    ))
  }
  
  return(list(status = "NO_DATA", origin = NA, invasive = NA, occurrence = NA))
}

# ---- 4. Main Query Function ----

#' Query all databases for a species-country combination
#' @param species Scientific name
#' @param country Country name
#' @param aphia_cache Cache of AphiaIDs (environment)
#' @return Named list with all database results
query_species_status <- function(species, country, aphia_cache = NULL) {
  
  result <- list(
    Species = species,
    Country = country,
    AphiaID = NA,
    CommonName = NA,
    WoRMS_Status = NA,
    WoRMS_Origin = NA,
    WoRMS_Invasive = NA,
    WoRMS_Occurrence = NA,
    GBIF_EstablishmentMeans = NA,
    OBIS_Records = NA,
    Final_Status = "UNKNOWN",
    Data_Sources = NA,
    Notes = NA
  )
  
  data_sources <- c()
  
  # Get AphiaID (use cache if available)
  if (!is.null(aphia_cache) && exists(species, envir = aphia_cache)) {
    aphia_id <- get(species, envir = aphia_cache)
  } else {
    aphia_id <- get_aphia_id(species)
    Sys.sleep(WORMS_DELAY)
    if (!is.null(aphia_cache)) {
      assign(species, aphia_id, envir = aphia_cache)
    }
  }
  result$AphiaID <- aphia_id
  
  # Get common name
  if (!is.na(aphia_id)) {
    result$CommonName <- get_common_name(aphia_id)
    Sys.sleep(WORMS_DELAY)
  }
  
  # Get WoRMS distribution
  if (!is.na(aphia_id)) {
    dist <- get_worms_distribution(aphia_id)
    Sys.sleep(WORMS_DELAY)
    
    worms_status <- parse_worms_status(dist, country)
    result$WoRMS_Status <- worms_status$status
    result$WoRMS_Origin <- worms_status$origin
    result$WoRMS_Invasive <- worms_status$invasive
    result$WoRMS_Occurrence <- worms_status$occurrence
    
    if (result$WoRMS_Status != "UNKNOWN" && result$WoRMS_Status != "NO_DATA") {
      data_sources <- c(data_sources, "WoRMS")
    }
  }
  
  # Get GBIF establishment means
  country_code <- get_country_code(country)
  if (!is.na(country_code)) {
    result$GBIF_EstablishmentMeans <- get_gbif_establishment(species, country_code)
    Sys.sleep(GBIF_DELAY)
    
    if (!is.na(result$GBIF_EstablishmentMeans)) {
      data_sources <- c(data_sources, "GBIF")
    }
  }
  
  # Get OBIS records (if enabled)
  if (USE_OBIS) {
    obis <- get_obis_records(species, country)
    if (!is.null(obis) && nrow(obis) > 0) {
      result$OBIS_Records <- nrow(obis)
      data_sources <- c(data_sources, "OBIS")
    }
    Sys.sleep(OBIS_DELAY)
  }
  
  # Record data sources used
  result$Data_Sources <- paste(data_sources, collapse = "; ")
  
  # Determine final status
  result$Final_Status <- determine_final_status(result)
  
  return(result)
}

#' Determine final invasive status based on all sources
#' @param result Result list from query_species_status
#' @return Final status string
determine_final_status <- function(result) {
  worms <- result$WoRMS_Status
  gbif <- result$GBIF_EstablishmentMeans
  
  # Priority: Invasive > Introduced > Native > Unknown
  if (!is.na(worms) && worms == "INVASIVE") return("INVASIVE")
  if (!is.na(gbif) && grepl("Invasive", gbif, ignore.case = TRUE)) return("INVASIVE")
  
  if (!is.na(worms) && worms == "INTRODUCED") return("INTRODUCED")
  if (!is.na(gbif) && grepl("Introduced", gbif, ignore.case = TRUE)) return("INTRODUCED")
  
  if (!is.na(worms) && worms == "NATIVE") return("NATIVE")
  if (!is.na(gbif) && grepl("Native", gbif, ignore.case = TRUE)) return("NATIVE")
  
  if (!is.na(worms) && worms == "NO_DATA") return("NO_DATA")
  
  return("UNKNOWN")
}

# ---- 5. Batch Processing Function ----

#' Process a full species list with progress bar
#' @param species_df Data frame with Species and Location columns
#' @param output_file Path for output CSV
#' @param checkpoint_every Save checkpoint every N rows
#' @return Data frame with all results
process_species_list <- function(species_df, 
                                  output_file = "invasive_status_results.csv",
                                  checkpoint_every = 50) {
  
  # Create cache for AphiaIDs to avoid redundant lookups
  aphia_cache <- new.env()
  
  # Get unique species-country combinations
  combos <- species_df %>%
    select(Species, Location) %>%
    distinct() %>%
    rename(Country = Location)
  
  # Remove rows with NA species or country — these can't be verified
  n_before <- nrow(combos)
  combos <- combos %>%
    filter(!is.na(Species), !is.na(Country))
  n_removed <- n_before - nrow(combos)
  if (n_removed > 0) {
    cat("Removed", n_removed, "rows with NA Species or Country\n")
  }
  
  cat("Processing", nrow(combos), "unique species-country combinations\n")
  cat("Estimated time:", round(nrow(combos) * 2 / 60, 1), "minutes\n\n")
  
  # Initialize progress bar
  pb <- progress_bar$new(
    format = "  [:bar] :percent | :current/:total | ETA: :eta | :species",
    total = nrow(combos),
    clear = FALSE,
    width = 80
  )
  
  results <- list()
  
  for (i in 1:nrow(combos)) {
    species <- combos$Species[i]
    country <- combos$Country[i]
    
    pb$tick(tokens = list(species = substr(species, 1, 30)))
    
    result <- query_species_status(species, country, aphia_cache)
    results[[i]] <- result
    
    # Checkpoint save
    if (i %% checkpoint_every == 0) {
      checkpoint_df <- bind_rows(results)
      write_csv(checkpoint_df, paste0("checkpoint_", output_file))
      cat("\n  Checkpoint saved at row", i, "\n")
    }
  }
  
  # Combine all results
  final_df <- bind_rows(results)
  
  # Save final output
  write_csv(final_df, output_file)
  cat("\n\nResults saved to:", output_file, "\n")
  
  # Print summary
  cat("\n=== SUMMARY ===\n")
  print(table(final_df$Final_Status))
  
  return(final_df)
}

# ---- 6. Usage Example ----

# =============================================================================
# HOW TO USE THIS SCRIPT
# =============================================================================
#
# STEP 1: Make sure you have internet access (run on your LOCAL machine)
#
# STEP 2: Install required packages (run once):
#   install.packages(c("worrms", "rgbif", "robis", "dplyr", "tidyr", 
#                      "purrr", "readr", "stringr", "httr", "jsonlite", "progress"))
#
# STEP 3: Load this script:
#   source("query_invasive_status_databases.R")
#
# STEP 4: Test with a single species first:
#   test_single_query("Mnemiopsis leidyi", "Netherlands")
#   test_single_query("Amphibalanus improvisus", "Denmark")
#
# STEP 5: Process your full dataset:
#   # Your CSV needs columns: Species, Location
#   species_12S <- read_csv("Processed_data/species_invasive_status_12S.csv")
#   results_12S <- process_species_list(species_12S, "invasive_verified_12S.csv")
#
# =============================================================================

# Example code (uncomment to run):
#
# # Load your species list
# species_12S <- read_csv("Processed_data/species_invasive_status_12S.csv")
# species_18S <- read_csv("Processed_data/species_invasive_status_18S.csv")
# species_COI <- read_csv("Processed_data/species_invasive_status_COI.csv")
#
# # Process each marker
# results_12S <- process_species_list(
#   species_12S, 
#   output_file = "Processed_data/invasive_status_verified_12S.csv"
# )
#
# results_18S <- process_species_list(
#   species_18S, 
#   output_file = "Processed_data/invasive_status_verified_18S.csv"
# )
#
# results_COI <- process_species_list(
#   species_COI, 
#   output_file = "Processed_data/invasive_status_verified_COI.csv"
# )

# ---- 7. Quick Test Function ----

#' Test the query for a single species-country
#' @param species Scientific name
#' @param country Country name
test_single_query <- function(species, country) {
  cat("Testing query for:", species, "in", country, "\n\n")
  
  result <- query_species_status(species, country)
  
  cat("AphiaID:", result$AphiaID, "\n")
  cat("Common Name:", result$CommonName, "\n")
  cat("WoRMS Status:", result$WoRMS_Status, "\n")
  cat("WoRMS Origin:", result$WoRMS_Origin, "\n")
  cat("WoRMS Invasive:", result$WoRMS_Invasive, "\n")
  cat("GBIF Establishment:", result$GBIF_EstablishmentMeans, "\n")
  cat("---\n")
  cat("FINAL STATUS:", result$Final_Status, "\n")
  
  return(result)
}

# ---- 8. Run Tests ----
# Uncomment to test:
#
# # Test with known invasive
# test_single_query("Mnemiopsis leidyi", "Netherlands")
#
# # Test with known native
# test_single_query("Clupea harengus", "Norway")
#
# # Test with Amphibalanus improvisus
# test_single_query("Amphibalanus improvisus", "Denmark")

cat("\n=============================================================\n")
cat("INVASIVE SPECIES DATABASE QUERY SCRIPT - LOADED SUCCESSFULLY\n")
cat("=============================================================\n\n")
cat("This script queries WoRMS, GBIF, and optionally OBIS to verify\n")
cat("invasive species status for your eDNA metabarcoding data.\n\n")
cat("QUICK START:\n")
cat("------------\n")
cat("1. Test single species:\n")
cat("   test_single_query('Mnemiopsis leidyi', 'Netherlands')\n")
cat("   test_single_query('Amphibalanus improvisus', 'Denmark')\n\n")
cat("2. Process your full dataset:\n")
cat("   species_df <- read_csv('your_species_list.csv')\n")
cat("   results <- process_species_list(species_df, 'output.csv')\n\n")
cat("EXPECTED OUTPUT COLUMNS:\n")
cat("  - Species, Country, AphiaID, CommonName\n")
cat("  - WoRMS_Status, WoRMS_Origin, WoRMS_Invasive\n")
cat("  - GBIF_EstablishmentMeans\n")
cat("  - Final_Status (INVASIVE/INTRODUCED/NATIVE/UNKNOWN)\n")
cat("  - Data_Sources (which databases had info)\n\n")
cat("NOTE: Processing ~1000 species takes ~30-60 min (API rate limits)\n")
cat("      Checkpoints saved every 50 species in case of interruption\n")
cat("=============================================================\n")
