# =============================================================================
# INVASIVE SPECIES ANALYSIS (Steps 3–4)
# Species extraction, GISD/EASIN queries, and flagging
# =============================================================================
#
# Author: BGE Marine Invasive Species Project
# Date: 2026-02-26
#
# Functions:
#   extract_species_by_country() - Species-country pairs from phyloseq (user-defined column)
#   gisd_query_single()          - Query one species on GISD
#   query_gisd_batch()           - Batch GISD queries with caching
#   fetch_easin_list()           - Fetch EASIN list for one country
#   query_easin_all()            - Batch EASIN queries with caching
#   check_easin_status()         - Check one species against EASIN
#   process_invasives()          - Run Steps 3–4 for one marker
#   process_all_invasives()      - Run Steps 3–4 for all markers
#   parse_claude_assessment()    - Parse markdown assessment tables
#
# Requires: a cleaned phyloseq from FUNC_process_metadata.R
#
# Usage:
#   source("Scripts/Final/FUNC_process_metadata.R")
#   source("Scripts/Final/FUNC_process_invasives.R")
#   meta_12S <- process_metadata("12S", ps_input = bge12s_cs)
#   inv_12S  <- process_invasives("12S", meta_12S$phyloseq, country_col = "Location")
# =============================================================================

# ---- Load Required Packages ----
library(phyloseq)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(readr)
library(httr)
library(jsonlite)

# ---- Shared Configuration (set once) ----
if (!exists("CONFIG")) {
  CONFIG <- list(
    raw_data_dir       = "Raw_data",
    processed_data_dir = "Processed_data",
    scripts_dir        = "Scripts/Final",
    ntnu_metadata_file = "BGE_MNIS_CScodes.csv",
    project_name       = "BGE Marine Invasive Species Citizen Science",
    use_cache          = TRUE,
    gisd_cache_file    = "gisd_cache.rds",
    easin_cache_file   = "easin_cache.rds",
    rate_limit_delay   = 1.0,
    max_retries        = 3,
    verbose            = TRUE
  )
}

if (!exists("dbg", mode = "function")) {
  dbg <- function(...) {
    if (CONFIG$verbose) message(paste(...))
  }
}

# =============================================================================
# SECTION 3: SPECIES BY COUNTRY EXTRACTION
# =============================================================================

#' Extract species-country combinations from phyloseq
#' @param ps Phyloseq object
#' @param country_col Character: name of the sample_data column containing
#'   country/location information (e.g. "Location", "Country", "Site_country")
#' @return Data frame with Species and Location columns
extract_species_by_country <- function(ps, country_col) {

  dbg("Extracting species by country using column:", country_col)

  samp <- as(sample_data(ps), "data.frame")

  if (!country_col %in% colnames(samp)) {
    stop("Column '", country_col, "' not found in sample data. ",
         "Available columns: ", paste(colnames(samp), collapse = ", "))
  }

  otu <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) otu <- t(otu)

  tax <- as(tax_table(ps), "matrix") %>% as.data.frame()

  cn_lower <- tolower(colnames(tax))
  sp_col <- which(cn_lower == "species")

  if (length(sp_col) == 0) {
    stop("Species column not found in taxonomy table")
  }

  otu_pa <- (otu > 0) * 1

  results <- list()

  for (loc in unique(samp[[country_col]])) {
    loc_samples <- rownames(samp)[samp[[country_col]] == loc]
    loc_samples <- loc_samples[loc_samples %in% colnames(otu_pa)]

    if (length(loc_samples) == 0) next

    taxa_present <- rowSums(otu_pa[, loc_samples, drop = FALSE]) > 0
    present_taxa <- names(taxa_present)[taxa_present]

    species_names <- tax[present_taxa, sp_col]
    species_names <- species_names[!is.na(species_names) & species_names != ""]
    species_names <- unique(species_names)

    if (length(species_names) > 0) {
      results[[loc]] <- tibble(
        Species = species_names,
        Location = loc
      )
    }
  }

  df_species_country <- bind_rows(results) %>%
    distinct() %>%
    arrange(Location, Species)

  dbg("  Found", nrow(df_species_country), "unique species-country combinations")
  dbg("  Unique species:", n_distinct(df_species_country$Species))
  dbg("  Countries:", paste(unique(df_species_country$Location), collapse = ", "))

  return(df_species_country)
}

# =============================================================================
# SECTION 4: INVASIVE SPECIES DATABASE QUERIES
# =============================================================================

# Country code mappings
COUNTRY_MAP <- c(
  "Spain" = "ES", "Portugal" = "PT", "Norway" = "NO", "Denmark" = "DK",
  "Finland" = "FI", "France" = "FR", "Germany" = "DE", "Italy" = "IT",
  "Netherlands" = "NL", "Poland" = "PL", "UK" = "GB", "United Kingdom" = "GB",
  "Greece" = "EL", "Sweden" = "SE", "Belgium" = "BE", "Ireland" = "IE",
  "Croatia" = "HR", "Slovenia" = "SI", "Estonia" = "EE", "Latvia" = "LV",
  "Lithuania" = "LT"
)

#' Query single species from GISD
#' @param species Species name
#' @return GISD status string or NA
gisd_query_single <- function(species) {
  tryCatch({
    url <- paste0("http://www.iucngisd.org/gisd/species.php?search_species=",
                  URLencode(species, reserved = TRUE))
    res <- GET(url, timeout(15))

    if (status_code(res) != 200) return(NA_character_)

    txt <- content(res, as = "text", encoding = "UTF-8")

    if (grepl("Invasive Species", txt, ignore.case = TRUE)) {
      return("Invasive")
    } else if (grepl("species not found|no results", txt, ignore.case = TRUE)) {
      return(NA_character_)
    } else if (grepl(species, txt, ignore.case = TRUE)) {
      return("Listed")
    }

    return(NA_character_)
  }, error = function(e) {
    return(NA_character_)
  })
}

#' Query GISD for multiple species with caching
#' @param species_list Vector of species names
#' @return Data frame with Species and GISD_status
query_gisd_batch <- function(species_list) {

  cache_file <- file.path(CONFIG$processed_data_dir, CONFIG$gisd_cache_file)

  cached <- if (CONFIG$use_cache && file.exists(cache_file)) {
    dbg("Loading GISD cache...")
    readRDS(cache_file)
  } else {
    tibble(Species = character(), GISD_status = character(), queried_date = as.Date(character()))
  }

  to_query <- setdiff(species_list, cached$Species)
  dbg("GISD:", length(species_list), "species,", length(to_query), "need querying")

  if (length(to_query) > 0) {
    new_results <- tibble(Species = character(), GISD_status = character(),
                          queried_date = as.Date(character()))

    for (i in seq_along(to_query)) {
      sp <- to_query[i]

      if (i %% 10 == 1 || i == length(to_query)) {
        dbg(sprintf("  GISD: %d/%d - %s", i, length(to_query), sp))
      }

      status <- gisd_query_single(sp)
      new_results <- bind_rows(new_results, tibble(
        Species = sp,
        GISD_status = status,
        queried_date = Sys.Date()
      ))

      Sys.sleep(CONFIG$rate_limit_delay)
    }

    cached <- bind_rows(cached, new_results)

    if (CONFIG$use_cache) {
      saveRDS(cached, cache_file)
      dbg("  Saved GISD cache")
    }
  }

  cached %>%
    filter(Species %in% species_list) %>%
    select(Species, GISD_status)
}

#' Fetch EASIN species list for a country
#' @param country_code ISO country code
#' @return Vector of species names
fetch_easin_list <- function(country_code) {

  urls <- c(
    paste0("https://easin.jrc.ec.europa.eu/apixg/catxg/concernedms/", country_code, "/skip/0/take/500"),
    paste0("https://easin.jrc.ec.europa.eu/apixg/catxg/establishedms/", country_code, "/skip/0/take/500")
  )

  all_species <- character(0)

  for (url in urls) {
    for (attempt in seq_len(CONFIG$max_retries)) {
      result <- tryCatch({
        res <- GET(url, timeout(30))

        if (status_code(res) != 200) {
          Sys.sleep(1 * attempt)
          next
        }

        txt <- content(res, as = "text", encoding = "UTF-8")
        parsed <- fromJSON(txt, flatten = TRUE)

        if (is.data.frame(parsed) && "Name" %in% names(parsed)) {
          return(unique(trimws(as.character(parsed$Name))))
        }

        return(character(0))
      }, error = function(e) {
        Sys.sleep(1 * attempt)
        NULL
      })

      if (!is.null(result)) {
        all_species <- unique(c(all_species, result))
        break
      }
    }
    Sys.sleep(0.5)
  }

  all_species
}

#' Query EASIN for all countries in dataset
#' @param countries Vector of country names
#' @return List of species per country code
query_easin_all <- function(countries) {

  cache_file <- file.path(CONFIG$processed_data_dir, CONFIG$easin_cache_file)

  cached <- if (CONFIG$use_cache && file.exists(cache_file)) {
    dbg("Loading EASIN cache...")
    readRDS(cache_file)
  } else {
    list()
  }

  codes_needed <- unique(na.omit(COUNTRY_MAP[countries]))
  codes_to_query <- setdiff(codes_needed, names(cached))

  dbg("EASIN:", length(codes_needed), "countries,", length(codes_to_query), "need querying")

  if (length(codes_to_query) > 0) {
    for (code in codes_to_query) {
      dbg("  EASIN:", code, "...")
      species_list <- fetch_easin_list(code)
      cached[[code]] <- list(species = species_list, queried_date = Sys.Date())
      dbg("    Found", length(species_list), "species")
      Sys.sleep(1)
    }

    if (CONFIG$use_cache) {
      saveRDS(cached, cache_file)
      dbg("  Saved EASIN cache")
    }
  }

  lapply(cached, function(x) x$species)
}

#' Check EASIN status for a species in a country
#' @param species Species name
#' @param location Country name
#' @param easin_lists Result from query_easin_all
#' @return Status string
check_easin_status <- function(species, location, easin_lists) {
  code <- COUNTRY_MAP[location]
  if (is.na(code)) return("Country_not_mapped")

  invasives <- easin_lists[[code]]
  if (is.null(invasives) || length(invasives) == 0) return("No_EASIN_data")

  if (species %in% invasives) "EASIN_listed" else "Not_listed"
}

# =============================================================================
# PIPELINE WRAPPERS
# =============================================================================

#' Analyse invasive species status for a single marker (Steps 3–4)
#'
#' Takes a cleaned phyloseq object (output of process_metadata), extracts
#' species-by-country combinations, queries GISD and EASIN databases, and
#' flags putative invasive species.
#'
#' @param marker Character: "12S", "18S", or "COI"
#' @param ps A phyloseq object with cleaned metadata (e.g. from
#'   process_metadata()$phyloseq). Must have a country column and Species in
#'   the taxonomy table.
#' @param country_col Character: name of the sample_data column containing
#'   country/location information (e.g. "Location", "Country")
#' @param save_outputs Logical: whether to save output files
#' @return List with: marker, species_country, invasive_status
process_invasives <- function(marker, ps, country_col, save_outputs = TRUE) {

  dbg("\n", paste(rep("=", 60), collapse = ""))
  dbg("INVASIVE SPECIES ANALYSIS:", marker)
  dbg(paste(rep("=", 60), collapse = ""), "\n")

  # Step 3: Species by country
  dbg("STEP 3: Extracting species by country...")
  df_species_country <- extract_species_by_country(ps, country_col)

  # Step 4: Query invasive databases
  dbg("\nSTEP 4: Querying invasive species databases...")

  dbg("  Querying GISD...")
  gisd_res <- query_gisd_batch(unique(df_species_country$Species))

  dbg("  Querying EASIN...")
  easin_lists <- query_easin_all(unique(df_species_country$Location))

  df_species_country <- df_species_country %>%
    rowwise() %>%
    mutate(EASIN_status = check_easin_status(Species, Location, easin_lists)) %>%
    ungroup()

  df_flagged <- df_species_country %>%
    left_join(gisd_res, by = "Species") %>%
    mutate(
      invasive_flag = case_when(
        EASIN_status == "EASIN_listed" ~ "Invasive_in_country",
        !is.na(GISD_status) & GISD_status != "" ~ "GISD_listed_globally",
        TRUE ~ "Not_flagged"
      )
    ) %>%
    arrange(Location, desc(invasive_flag == "Invasive_in_country"), Species)

  dbg("\n=== SUMMARY ===")
  dbg("Total species-country pairs:", nrow(df_flagged))

  flag_summary <- df_flagged %>% count(invasive_flag)
  print(flag_summary)

  if (save_outputs) {
    out_file <- file.path(CONFIG$processed_data_dir,
                          paste0("species_invasive_status_", marker, ".csv"))
    write_csv(df_flagged, out_file)
    dbg("Saved:", out_file)

    manual_review <- df_flagged %>%
      filter(invasive_flag != "Not_flagged") %>%
      select(Species, Location, GISD_status, EASIN_status, invasive_flag) %>%
      arrange(Species, Location)

    review_file <- file.path(CONFIG$processed_data_dir,
                             paste0("species_manual_review_", marker, ".csv"))
    write_csv(manual_review, review_file)
    dbg("Saved:", review_file)
  }

  list(
    marker = marker,
    species_country = df_species_country,
    invasive_status = df_flagged
  )
}

#' Analyse invasive species for all markers in sequence
#' @param markers Vector of marker names (default: c("12S", "18S", "COI"))
#' @param meta_list Named list of process_metadata results (each must contain $phyloseq)
#' @param country_col Character: name of the sample_data column containing
#'   country/location information
#' @return Named list of process_invasives results for each marker
process_all_invasives <- function(markers = c("12S", "18S", "COI"),
                                  meta_list,
                                  country_col) {

  results <- list()

  for (m in markers) {
    if (is.null(meta_list[[m]]$phyloseq)) {
      stop("No phyloseq found in meta_list for marker: ", m)
    }
    results[[m]] <- process_invasives(m, meta_list[[m]]$phyloseq, country_col)
  }

  dbg("\n", paste(rep("=", 60), collapse = ""))
  dbg("ALL INVASIVE ANALYSES COMPLETE")
  dbg(paste(rep("=", 60), collapse = ""))

  dbg("\nCOMBINED SUMMARY:")
  for (m in markers) {
    dbg(sprintf("  %s: %d species-country pairs, %d flagged as invasive",
                m,
                nrow(results[[m]]$invasive_status),
                sum(results[[m]]$invasive_status$invasive_flag != "Not_flagged")))
  }

  return(results)
}

# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

#' Parse Claude's markdown assessment back into R
#' @param md_file Path to markdown file
#' @param marker Marker name (for output naming)
#' @return Data frame with invasive status
parse_claude_assessment <- function(md_file, marker = NULL) {

  if (is.null(marker)) {
    marker <- str_extract(basename(md_file), "12S|18S|COI")
    if (is.na(marker)) marker <- "unknown"
  }

  lines <- readLines(md_file, warn = FALSE)

  country_lines <- grep("^##\\s+", lines)
  section_bounds <- c(country_lines, length(lines) + 1)

  extract_table <- function(start_idx, end_idx) {
    section <- lines[start_idx:(end_idx - 1)]
    country <- str_replace(section[1], "^##\\s+", "") %>% str_trim()

    tbl_lines <- grep("^\\|", section, value = TRUE)
    if (length(tbl_lines) == 0) return(NULL)

    tbl_lines <- tbl_lines[!grepl("^\\|[-:]+\\|", tbl_lines)]

    tbl_split <- str_split(tbl_lines, "\\|")
    tbl_split <- lapply(tbl_split, function(x) str_trim(x))
    tbl_split <- lapply(tbl_split, function(x) x[x != ""])

    header <- tbl_split[[1]]
    rows <- tbl_split[-1]

    df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
    colnames(df) <- header
    df$Country <- country
    df <- df %>% select(Country, everything())

    return(df)
  }

  df_list <- purrr::map2(
    section_bounds[-length(section_bounds)],
    section_bounds[-1],
    extract_table
  )

  result <- bind_rows(df_list)

  dbg("Parsed", marker, "assessment:", nrow(result), "rows")

  return(result)
}

# =============================================================================
cat("\n")
cat("=============================================================\n")
cat("INVASIVE SPECIES FUNCTIONS LOADED (Steps 3-4)\n")
cat("=============================================================\n\n")
cat("  process_invasives('12S', ps, 'Location')  - Species extraction + DB queries\n")
cat("  process_all_invasives(meta_list, 'Location') - All markers at once\n")
cat("  parse_claude_assessment(f)         - Parse Claude's MD assessment\n")
cat("=============================================================\n")
