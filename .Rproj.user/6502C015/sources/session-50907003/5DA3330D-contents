### Check on if detections are novel or not considering GBIF records

# =============================================================================
# Check all MIS detections against GBIF for novel records
# =============================================================================

library(rgbif)
library(tidyverse)
library(geosphere)  # for distance calculations

# 1) Expand MIS summary to one row per species-site
# -----------------------------------------------------------------

# First, create the site coordinates lookup
site_coords <- as(sample_data(ps_cs_12S_original), "data.frame") %>%
  distinct(Sampling.event.ID, .keep_all = TRUE) %>%
  select(Site_Code = Sampling.event.ID,
         Latitude  = latitude_full,
         Longitude = longitude_full)

# Expand mis_summary to one row per species-site
mis_expanded <- mis_summary %>%
  separate_rows(Sites, sep = ", ") %>%
  rename(Site_Code = Sites) %>%
  left_join(site_coords, by = "Site_Code") %>%
  filter(!is.na(Latitude))  # Remove any without coords

cat("Total species-site combinations to check:", nrow(mis_expanded), "\n")

# 2) Function to query GBIF for nearest record
# -----------------------------------------------------------------

get_nearest_gbif_record <- function(species_name, det_lat, det_lon, search_radius_deg = 5) {
  
  # Try to get GBIF records
  tryCatch({
    # First, get the species key
    species_key <- name_backbone(name = species_name)$usageKey
    
    if (is.null(species_key)) {
      return(tibble(
        n_gbif_records = 0,
        nearest_distance_km = NA,
        nearest_year = NA,
        nearest_country = NA,
        same_country_records = 0,
        assessment = "Species not found in GBIF"
      ))
    }
    
    # Search for occurrences in wider area
    records <- occ_search(
      taxonKey = species_key,
      decimalLatitude = paste0(det_lat - search_radius_deg, ",", det_lat + search_radius_deg),
      decimalLongitude = paste0(det_lon - search_radius_deg, ",", det_lon + search_radius_deg),
      hasCoordinate = TRUE,
      limit = 500
    )
    
    if (is.null(records$data) || nrow(records$data) == 0) {
      # No records within search radius - expand search
      records_global <- occ_search(
        taxonKey = species_key,
        hasCoordinate = TRUE,
        limit = 1000
      )
      
      if (is.null(records_global$data) || nrow(records_global$data) == 0) {
        return(tibble(
          n_gbif_records = 0,
          nearest_distance_km = NA,
          nearest_year = NA,
          nearest_country = NA,
          same_country_records = 0,
          assessment = "No GBIF records with coordinates"
        ))
      }
      
      # Calculate distances to all global records
      records_df <- records_global$data %>%
        filter(!is.na(decimalLatitude), !is.na(decimalLongitude))
      
      if (nrow(records_df) == 0) {
        return(tibble(
          n_gbif_records = 0,
          nearest_distance_km = NA,
          nearest_year = NA,
          nearest_country = NA,
          same_country_records = 0,
          assessment = "No GBIF records with coordinates"
        ))
      }
      
      distances <- distHaversine(
        cbind(det_lon, det_lat),
        cbind(records_df$decimalLongitude, records_df$decimalLatitude)
      ) / 1000  # Convert to km
      
      nearest_idx <- which.min(distances)
      
      return(tibble(
        n_gbif_records = nrow(records_df),
        nearest_distance_km = round(min(distances), 1),
        nearest_year = records_df$year[nearest_idx],
        nearest_country = records_df$country[nearest_idx],
        same_country_records = 0,  # None within search radius
        assessment = ifelse(min(distances) > 500, "POTENTIALLY NOVEL - nearest record >500km", 
                            ifelse(min(distances) > 200, "Range edge - nearest record >200km",
                                   "Known in region"))
      ))
    }
    
    # Calculate distances to all records
    records_df <- records$data %>%
      filter(!is.na(decimalLatitude), !is.na(decimalLongitude))
    
    if (nrow(records_df) == 0) {
      return(tibble(
        n_gbif_records = 0,
        nearest_distance_km = NA,
        nearest_year = NA,
        nearest_country = NA,
        same_country_records = 0,
        assessment = "No GBIF records with coordinates"
      ))
    }
    
    distances <- distHaversine(
      cbind(det_lon, det_lat),
      cbind(records_df$decimalLongitude, records_df$decimalLatitude)
    ) / 1000  # Convert to km
    
    nearest_idx <- which.min(distances)
    min_dist <- min(distances)
    
    # Count records within 50km (same area)
    records_nearby <- sum(distances < 50)
    
    # Determine assessment
    if (min_dist < 50) {
      assessment <- "Well documented"
    } else if (min_dist < 100) {
      assessment <- "Known nearby (<100km)"
    } else if (min_dist < 200) {
      assessment <- "Range edge (100-200km)"
    } else if (min_dist < 500) {
      assessment <- "Possible range expansion (200-500km)"
    } else {
      assessment <- "POTENTIALLY NOVEL (>500km from nearest record)"
    }
    
    return(tibble(
      n_gbif_records = nrow(records_df),
      nearest_distance_km = round(min_dist, 1),
      nearest_year = records_df$year[nearest_idx],
      nearest_country = records_df$country[nearest_idx],
      same_country_records = records_nearby,
      assessment = assessment
    ))
    
  }, error = function(e) {
    return(tibble(
      n_gbif_records = NA,
      nearest_distance_km = NA,
      nearest_year = NA,
      nearest_country = NA,
      same_country_records = NA,
      assessment = paste("Error:", e$message)
    ))
  })
}

# 3) Run GBIF queries for all species-site combinations
# -----------------------------------------------------------------

cat("Querying GBIF for", nrow(mis_expanded), "species-site combinations...\n")
cat("This may take a few minutes...\n\n")

# Add progress tracking
gbif_results <- mis_expanded %>%
  mutate(row_id = row_number()) %>%
  rowwise() %>%
  mutate(
    gbif_check = {
      if (row_id %% 10 == 0) cat("Processing", row_id, "of", nrow(mis_expanded), "\n")
      Sys.sleep(0.5)  # Rate limiting
      list(get_nearest_gbif_record(Species, Latitude, Longitude))
    }
  ) %>%
  ungroup() %>%
  unnest(gbif_check)

# 4) Create summary table
# -----------------------------------------------------------------

# Full results table
mis_gbif_full <- gbif_results %>%
  select(
    Species, Phylum, Class, 
    Site_Code, Latitude, Longitude,
    Markers,
    n_gbif_records, nearest_distance_km, nearest_year, nearest_country,
    same_country_records, assessment
  ) %>%
  arrange(desc(nearest_distance_km))

cat("\n=============================================================================\n")
cat("FULL RESULTS: MIS detections vs GBIF records\n")
cat("=============================================================================\n\n")
print(mis_gbif_full, n = Inf)

# 5) Highlight potentially novel detections
# -----------------------------------------------------------------

novel_detections <- mis_gbif_full %>%
  filter(grepl("NOVEL|expansion|edge", assessment, ignore.case = TRUE)) %>%
  arrange(desc(nearest_distance_km))

cat("\n=============================================================================\n")
cat("POTENTIALLY NOVEL OR RANGE EXPANSION DETECTIONS\n")
cat("=============================================================================\n\n")

if (nrow(novel_detections) > 0) {
  print(novel_detections, n = Inf)
} else {
  cat("No potentially novel detections found - all species have nearby GBIF records.\n")
}

# 6) Summary by species
# -----------------------------------------------------------------

species_summary <- mis_gbif_full %>%
  group_by(Species, Phylum, Class) %>%
  summarise(
    n_sites_detected = n(),
    markers = paste(unique(Markers), collapse = ", "),
    min_distance_to_gbif = min(nearest_distance_km, na.rm = TRUE),
    max_distance_to_gbif = max(nearest_distance_km, na.rm = TRUE),
    any_novel = any(grepl("NOVEL", assessment)),
    .groups = "drop"
  ) %>%
  arrange(desc(max_distance_to_gbif))

cat("\n=============================================================================\n")
cat("SUMMARY BY SPECIES\n")
cat("=============================================================================\n\n")
print(species_summary, n = Inf)

# 7) Save results
# -----------------------------------------------------------------

write_csv(mis_gbif_full, "Processed_data/MIS_GBIF_validation.csv")
write_csv(novel_detections, "Processed_data/MIS_potentially_novel.csv")
write_csv(species_summary, "Processed_data/MIS_species_summary.csv")

cat("\n\nResults saved to Processed_data/\n")



# =============================================================================
# Add detection metrics (n_biosamples, n_pcrs) to MIS GBIF results
# =============================================================================
# Get detection metrics from combined_all$site_level - use correct columns
detection_metrics <- combined_all$site_level %>%
  filter(is_invasive) %>%
  select(Species, Site, Marker, 
         n_biosamples_positive, total_biosamples_at_site,
         n_pcr_positive, total_pcr_at_site,
         theta, p_empirical, reliability) %>%
  rename(Site_Code = Site)

# Check structure
head(detection_metrics)

# Join to mis_gbif_full - first remove old columns if they exist
mis_gbif_full <- mis_gbif_full %>%
  select(-any_of(c("n_biosamples_positive", "n_pcr_positive", "n_pcr_total", 
                   "theta", "reliability", "total_biosamples_at_site", 
                   "total_pcr_at_site", "p_empirical"))) %>%
  left_join(detection_metrics, by = c("Species", "Site_Code", "Markers" = "Marker"))

# Reorder columns for clarity
mis_gbif_full <- mis_gbif_full %>%
  select(
    Species, Phylum, Class, 
    Site_Code, Latitude, Longitude,
    Markers,
    n_biosamples_positive, total_biosamples_at_site,
    n_pcr_positive, total_pcr_at_site, 
    theta, reliability,
    n_gbif_records, nearest_distance_km, nearest_year, nearest_country,
    same_country_records, assessment
  ) %>%
  arrange(desc(nearest_distance_km))

cat("=============================================================================\n")
cat("FULL RESULTS with corrected detection metrics\n")
cat("=============================================================================\n\n")
print(mis_gbif_full, n = 30)

# Update novel_detections
novel_detections <- mis_gbif_full %>%
  filter(grepl("NOVEL|expansion|edge", assessment, ignore.case = TRUE)) %>%
  arrange(desc(nearest_distance_km))

cat("\n=============================================================================\n")
cat("POTENTIALLY NOVEL DETECTIONS\n")
cat("=============================================================================\n\n")
print(novel_detections, n = Inf)


# Full table of all MIS detections with GBIF validation

cat("=============================================================================\n")
cat("ALL MIS DETECTIONS with GBIF validation\n")
cat("=============================================================================\n\n")

print(mis_gbif_full %>% 
        mutate(Country = str_extract(Site_Code, "^[A-Za-z]+")) %>%
        select(Country, Site_Code, Species, Markers, 
               n_biosamples_positive, n_pcr_positive, total_pcr_at_site,
               theta, reliability, 
               n_gbif_records, nearest_distance_km, nearest_country, assessment) %>%
        arrange(Country, Site_Code, Species), 
      n = Inf, width = Inf)

# Export to CSV for easier viewing
write_csv(
  mis_gbif_full %>%
    mutate(Country = str_extract(Site_Code, "^[A-Za-z]+")) %>%
    select(
      Country, Site_Code, Species, Marker = Markers, Latitude, Longitude,
      n_biosamples_positive, total_biosamples_at_site,
      n_pcr_positive, total_pcr_at_site,
      theta, reliability,
      n_gbif_records, nearest_distance_km, nearest_year, nearest_country, assessment
    ) %>%
    arrange(Country, Site_Code, Species),
  "Processed_data/All_MIS_detections_GBIF_validation.csv"
)

cat("\nSaved all", nrow(mis_gbif_full), "MIS detections to:\n")
cat("Processed_data/All_MIS_detections_GBIF_validation.csv\n")





##### To get the read props..

# =============================================================================
# Calculate read metrics across ALL PCRs at site (including zeros)
# =============================================================================

# Get total PCRs per site per marker
total_pcrs_per_site <- pcr_all %>%
  distinct(Sample, Site_Code, Marker) %>%
  count(Site_Code, Marker, name = "total_pcrs_at_site")

# Join with invasive species list
invasive_species <- combined_all$site_level %>%
  filter(is_invasive) %>%
  distinct(Species, Marker)

# For each invasive species, get ALL PCRs at sites where it was detected
# First get the sites where each species was detected
species_sites <- pcr_all %>%
  inner_join(invasive_species, by = c("Species", "Marker")) %>%
  filter(reads > 0) %>%
  distinct(Species, Site_Code, Marker)

# Now get ALL PCRs at those sites (including zeros for this species)
pcr_invasive_all <- pcr_all %>%
  inner_join(species_sites, by = c("Species", "Site_Code", "Marker"))

# Calculate stats across ALL PCRs at each site
mean_prop_per_pcr <- pcr_invasive_all %>%
  group_by(Species, Site_Code, Marker) %>%
  summarise(
    total_pcrs_at_site = n(),
    n_pcrs_positive = sum(reads > 0),
    total_reads = sum(reads),
    mean_reads_per_pcr = round(mean(reads), 1),
    min_reads_in_pcr = min(reads),
    max_reads_in_pcr = max(reads),
    # Min/max among POSITIVE PCRs only
    min_reads_per_pos_pcr = ifelse(sum(reads > 0) > 0, min(reads[reads > 0]), NA),
    max_reads_per_pos_pcr = ifelse(sum(reads > 0) > 0, max(reads[reads > 0]), NA),
    mean_prop_per_pcr = mean(prop_reads_pcr),
    min_prop_in_pcr = min(prop_reads_pcr),
    max_prop_in_pcr = max(prop_reads_pcr),
    .groups = "drop"
  )

# Check Neogobius at UK_1 and Norway_4
print(mean_prop_per_pcr %>% 
        filter(Species == "Neogobius melanostomus") %>%
        arrange(n_pcrs_positive),
      n = Inf)

# Join to mis_gbif_full
mis_gbif_full <- mis_gbif_full %>%
  select(-any_of(c("n_pcrs_with_reads", "mean_reads_per_pcr", "mean_prop_reads_per_pcr",
                   "n_pcrs_positive", "total_pcrs_at_site", "total_reads",
                   "min_reads_in_pcr", "max_reads_in_pcr", 
                   "min_reads_per_pos_pcr", "max_reads_per_pos_pcr",
                   "mean_prop_per_pcr", "min_prop_in_pcr", "max_prop_in_pcr",
                   "sum_prop_positive_pcrs"))) %>%
  left_join(mean_prop_per_pcr, by = c("Species", "Site_Code", "Markers" = "Marker"))

# Update novel_detections - sort by reliability then distance
novel_detections <- mis_gbif_full %>%
  filter(grepl("NOVEL|expansion|edge", assessment, ignore.case = TRUE)) %>%
  mutate(
    reliability = factor(reliability, levels = c("Reliable", "Marginal", "Unreliable"))
  ) %>%
  arrange(reliability, desc(nearest_distance_km))

cat("=============================================================================\n")
cat("POTENTIALLY NOVEL DETECTIONS - sorted by reliability, then distance\n")
cat("=============================================================================\n\n")

print(novel_detections %>% 
        select(Species, Markers, Site_Code, n_biosamples_positive, n_pcrs_positive, total_pcrs_at_site,
               mean_reads_per_pcr, min_reads_in_pcr, max_reads_in_pcr,
               min_reads_per_pos_pcr, max_reads_per_pos_pcr,
               mean_prop_per_pcr, min_prop_in_pcr, max_prop_in_pcr,
               theta, reliability, nearest_distance_km) %>%
        mutate(
          mean_prop_pct = round(mean_prop_per_pcr * 100, 3),
          min_prop_pct = round(min_prop_in_pcr * 100, 3),
          max_prop_pct = round(max_prop_in_pcr * 100, 3)
        ) %>%
        select(-mean_prop_per_pcr, -min_prop_in_pcr, -max_prop_in_pcr), 
      n = Inf, width = Inf)

# Select and rename columns for export
novel_detections_export <- novel_detections %>%
  select(
    Species, Marker = Markers, Site_Code, Latitude, Longitude,
    n_biosamples_positive, total_biosamples_at_site,
    n_pcrs_positive, total_pcrs_at_site,
    total_reads, mean_reads_per_pcr, 
    min_reads_in_pcr, max_reads_in_pcr,
    `min_reads_per_+_pcr` = min_reads_per_pos_pcr,
    `max_reads_per_+_pcr` = max_reads_per_pos_pcr,
    mean_prop_per_pcr, min_prop_in_pcr, max_prop_in_pcr,
    theta, reliability, 
    n_gbif_records, nearest_distance_km, nearest_year, nearest_country, assessment
  )

# Save novel_detections to CSV
write_csv(novel_detections_export, "Processed_data/Potentially_novel_MISdetectionsover_100kmGBIF.csv")

cat("\nSaved", nrow(novel_detections_export), "potentially novel detections to:\n")
cat("Processed_data/Potentially_novel_MISdetectionsover_100kmGBIF.csv\n")

# =============================================================================
# PAPER TABLE: Formatted summary for notable MIS detections (>=100km from GBIF)
# =============================================================================

# Simplify assessment labels to match paper format
simplify_assessment <- function(x) {
  case_when(
    grepl("NOVEL|>500", x, ignore.case = TRUE)        ~ "NOVEL (>500km)",
    grepl("200.500|200-500|expansion", x, ignore.case = TRUE) ~ "Possible range expansion (200-500km)",
    grepl("100.200|100-200|Range edge", x, ignore.case = TRUE) ~ "Range edge (100-200km)",
    grepl("Known nearby|<100", x, ignore.case = TRUE)  ~ "Known nearby (<100km)",
    grepl("Well documented", x, ignore.case = TRUE)    ~ "Well documented",
    TRUE ~ x
  )
}

paper_table <- mis_gbif_full %>%
  # Filter to detections ≥100km from nearest GBIF record (notable range/novel detections)
  filter(!is.na(nearest_distance_km) & nearest_distance_km >= 100) %>%
  mutate(
    Country = str_extract(Site_Code, "^[A-Za-z]+"),
    `Site Code (Approx. Latitude, Longitude)` = paste0(Site_Code, " (", Latitude, ", ", Longitude, ")"),
    `Species (marker)` = paste0(Species, " (", Markers, ")"),
    `Samples +'ve (total)` = paste0(n_biosamples_positive, " (", total_biosamples_at_site, ")"),
    `PCRs +'ve (total)` = paste0(n_pcr_positive, " (", total_pcr_at_site, ")"),
    `Total reads` = total_reads,
    `Reliability Category` = reliability,
    `Min. GBIF distance (km)` = round(nearest_distance_km),
    `Min. GBIF yr` = nearest_year,
    `GBIF Assessment (closest record)` = simplify_assessment(assessment)
  ) %>%
  arrange(Country, Site_Code, Species) %>%
  select(
    Country,
    `Site Code (Approx. Latitude, Longitude)`,
    `Species (marker)`,
    `Samples +'ve (total)`,
    `PCRs +'ve (total)`,
    `Total reads`,
    `Reliability Category`,
    `Min. GBIF distance (km)`,
    `Min. GBIF yr`,
    `GBIF Assessment (closest record)`
  )

cat("\n=============================================================================\n")
cat("PAPER TABLE: Notable MIS detections (≥100km from nearest GBIF record)\n")
cat("=============================================================================\n\n")
print(paper_table, n = Inf, width = Inf)

write_csv(paper_table, "Processed_data/MIS_paper_table_GBIF_notable.csv")
cat("\nSaved", nrow(paper_table), "rows to: Processed_data/MIS_paper_table_GBIF_notable.csv\n")
