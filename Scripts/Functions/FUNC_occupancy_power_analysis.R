# =============================================================================
# SITE-LEVEL DETECTION RELIABILITY FOR eDNA METABARCODING
# =============================================================================
#
# Design: 3 biological samples × 3 PCR replicates per site (9 total)
#
# KEY DESIGN PRINCIPLES:
#   - All detection parameters (θ, p) are computed PER OTU × SITE
#   - No cross-site pooling: sites span different biogeographic regions
#   - Occupancy (ψ) requires multiple sites — not applicable per-site
#   - Detection probability (p) is estimated per site from replicates
#   - Power analysis uses Beta-posterior bootstrapping to propagate
#     uncertainty from small sample sizes into design recommendations
#
# This script provides:
#   1. Site-level species detection summaries (θ, p per OTU × site)
#   2. Detection reliability scoring (confidence, p_detect_site)
#   3. Bootstrap power analysis: P(detect) with 95% CIs across designs
#   4. Invasive status merging (EASIN/GISD + WoRMS/GBIF verification)
#   5. Visualization functions for site-level and country-level summaries
#
# Author: BGE Marine Invasive Species Project
# =============================================================================

library(phyloseq)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(purrr)

# scales is optional for nice axis labels
if (!requireNamespace("scales", quietly = TRUE)) {
  message("Note: Install 'scales' package for nicer axis formatting")
}

# =============================================================================
# SECTION 1: HIERARCHICAL DETECTION ANALYSIS
# =============================================================================
#
# Design hierarchy:
#   Site (Sampling.area.Name)
#   └── Biological sample (Name, e.g., BGEMIS0001)
#       └── PCR replicate (Replicate: R1, R2, R3)
#
# Detection parameters:
#   p     = P(PCR detection | DNA in biological sample)
#   θ     = P(DNA in sample | species present at site)  
#   p*    = P(detect in ≥1 PCR | DNA in sample) = 1-(1-p)^n_pcr
#   θ*    = P(detect in ≥1 sample | species at site) = 1-(1-θ*p*)^n_samples
# =============================================================================

#' Get hierarchical structure of samples
#' @param ps Phyloseq object
#' @return Data frame with site, biological sample, and PCR replicate info
get_sample_hierarchy <- function(ps,
                                  site_var = "Sampling.area.Name",
                                  biosample_var = "Name",
                                  pcr_var = "Replicate",
                                  country_var = "Location") {
  
  samp <- as(sample_data(ps), "data.frame")
  samp$sample_id <- rownames(samp)
  
  # Handle country - try specified var first, then fallback to Sampling.area.Country
  if (country_var %in% colnames(samp) && !all(is.na(samp[[country_var]]))) {
    samp$Country <- samp[[country_var]]
  } else if ("Sampling.area.Country" %in% colnames(samp)) {
    samp$Country <- samp[["Sampling.area.Country"]]
    message("Using Sampling.area.Country for country information")
  } else {
    samp$Country <- NA
    warning("No country information found")
  }
  
  # Fill NA countries with Sampling.area.Country if available
  if ("Sampling.area.Country" %in% colnames(samp)) {
    na_countries <- is.na(samp$Country)
    if (any(na_countries)) {
      samp$Country[na_countries] <- samp[["Sampling.area.Country"]][na_countries]
    }
  }
  
  samp %>%
    transmute(
      sample_id = sample_id,
      Country = Country,
      Site = .data[[site_var]],
      BioSample = .data[[biosample_var]],
      PCR_Rep = .data[[pcr_var]]
    )
}

#' Summarize the experimental design
#' @param hierarchy Output from get_sample_hierarchy
#' @return Summary statistics
summarize_design <- function(hierarchy) {
  
  site_summary <- hierarchy %>%
    group_by(Country, Site) %>%
    summarise(
      n_biosamples = n_distinct(BioSample),
      n_pcr_total = n(),
      pcr_per_biosample = n() / n_distinct(BioSample),
      .groups = "drop"
    )
  
  list(
    n_countries = n_distinct(hierarchy$Country, na.rm = TRUE),
    n_sites = n_distinct(hierarchy$Site),
    n_biosamples = n_distinct(hierarchy$BioSample),
    n_pcr_total = nrow(hierarchy),
    site_details = site_summary,
    
    # Typical design
    median_biosamples_per_site = median(site_summary$n_biosamples),
    median_pcr_per_biosample = median(site_summary$pcr_per_biosample)
  )
}

#' Calculate hierarchical detection rates for each OTU at each site
#' @param ps Phyloseq object
#' @param hierarchy Output from get_sample_hierarchy
#' @return Data frame with detection at PCR and biological sample levels (OTU level)
calc_hierarchical_detection <- function(ps, hierarchy) {
  
  otu <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) otu <- t(otu)
  
  # Get unique sites
  sites <- unique(hierarchy$Site)
  
  results <- list()
  
  for (site in sites) {
    site_hierarchy <- hierarchy %>% filter(Site == !!site)
    country <- site_hierarchy$Country[1]
    biosamples <- unique(site_hierarchy$BioSample)
    
    # For each biological sample, get PCR replicates
    for (biosample in biosamples) {
      biosample_data <- site_hierarchy %>% filter(BioSample == !!biosample)
      pcr_samples <- biosample_data$sample_id
      pcr_samples <- pcr_samples[pcr_samples %in% colnames(otu)]
      
      if (length(pcr_samples) == 0) next
      
      # Get detection across PCR replicates
      otu_pcr <- otu[, pcr_samples, drop = FALSE]
      pa_pcr <- (otu_pcr > 0) * 1
      
      # OTUs detected in at least one PCR
      detected <- rowSums(pa_pcr) > 0
      if (sum(detected) == 0) next
      
      for (otu_id in names(detected)[detected]) {
        n_pcr_pos <- sum(pa_pcr[otu_id, ])
        n_pcr_total <- length(pcr_samples)
        
        results[[length(results) + 1]] <- data.frame(
          Country = country,
          Site = site,
          BioSample = biosample,
          OTU = otu_id,
          n_pcr_positive = n_pcr_pos,
          n_pcr_total = n_pcr_total,
          p_pcr = n_pcr_pos / n_pcr_total,
          total_reads = sum(otu_pcr[otu_id, ]),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  bind_rows(results)
}

#' Aggregate to site level with proper hierarchical estimates (OTU level)
#' @param pcr_detection Output from calc_hierarchical_detection
#' @return Site-level detection summary per OTU
calc_site_level_detection <- function(pcr_detection) {
  
  pcr_detection %>%
    group_by(Country, Site, OTU) %>%
    summarise(
      # Biological sample level
      n_biosamples_positive = n(),
      
      # PCR level (across all biosamples at site)
      n_pcr_positive = sum(n_pcr_positive),
      n_pcr_total = sum(n_pcr_total),
      
      # Within-biosample PCR detection rates
      mean_p_pcr = mean(p_pcr),
      min_p_pcr = min(p_pcr),
      max_p_pcr = max(p_pcr),
      
      # PCR distribution across biosamples (to distinguish 3+1 vs 2+2 patterns)
      min_pcr_per_pos_biosample = min(n_pcr_positive),
      max_pcr_per_pos_biosample = max(n_pcr_positive),
      pcr_distribution = paste(sort(n_pcr_positive, decreasing = TRUE), collapse = "+"),
      
      # Total reads
      total_reads = sum(total_reads),
      
      .groups = "drop"
    ) %>%
    mutate(
      p_pcr_pooled = n_pcr_positive / n_pcr_total,
      # Is detection evenly spread across biosamples?
      pcr_spread_ratio = min_pcr_per_pos_biosample / max_pcr_per_pos_biosample
    )
}

#' Add total biosamples per site (including non-detections)
#' @param site_detection Output from calc_site_level_detection
#' @param hierarchy Output from get_sample_hierarchy
#' @return Enhanced site detection with theta estimates
add_biosample_totals <- function(site_detection, hierarchy) {
  
  # Get total biosamples per site
  site_totals <- hierarchy %>%
    group_by(Site) %>%
    summarise(
      total_biosamples_at_site = n_distinct(BioSample),
      total_pcr_at_site = n(),
      .groups = "drop"
    )
  
  site_detection %>%
    left_join(site_totals, by = "Site") %>%
    mutate(
      # === RAW EMPIRICAL METRICS (no model assumptions) ===
      
      # Simple detection rate: PCR detections / total PCR opportunities
      p_empirical = n_pcr_positive / total_pcr_at_site,
      
      # Theta = proportion of biological samples with detection
      theta = n_biosamples_positive / total_biosamples_at_site,
      
      # === MODEL-BASED METRICS (assumes OTU is truly present) ===
      
      # Cumulative detection probabilities
      p_star = 1 - (1 - mean_p_pcr)^3,  # P(>=1 PCR | DNA in sample)
      
      # P(miss all 3 PCRs | DNA in sample)
      p_miss_pcr = (1 - mean_p_pcr)^3,
      # P(miss biosample) = P(no DNA) + P(DNA but PCR fails)
      p_miss_biosample = (1 - theta) + theta * p_miss_pcr,
      # P(miss site) = P(miss all biosamples)
      p_miss_site = p_miss_biosample^total_biosamples_at_site,
      # P(detect at site | OTU present) - MODEL BASED
      p_detect_site = 1 - p_miss_site,
      
      # === CONFIDENCE SCORE (weights biological > PCR replication) ===
      # Now incorporates PCR distribution across biosamples
      
      confidence_score = case_when(
        # Detected in all biosamples with high PCR rate
        theta == 1 & mean_p_pcr >= 0.67 ~ 1.0,
        # Detected in all biosamples
        theta == 1 ~ 0.9,
        # Detected in majority of biosamples with even spread
        theta >= 0.67 & pcr_spread_ratio >= 0.5 ~ 0.75 + 0.1 * mean_p_pcr,
        # Detected in majority of biosamples with uneven spread (e.g., 3+1)
        theta >= 0.67 ~ 0.65 + 0.1 * mean_p_pcr,
        # Detected in 2 biosamples with even spread (e.g., 2+2)
        n_biosamples_positive >= 2 & pcr_spread_ratio >= 0.5 ~ 0.55 + 0.15 * theta,
        # Detected in 2 biosamples with uneven spread (e.g., 3+1)
        n_biosamples_positive >= 2 ~ 0.45 + 0.15 * theta,
        # Single biosample but ALL PCRs positive
        n_biosamples_positive == 1 & mean_p_pcr == 1 ~ 0.35,
        # Single biosample, multiple PCRs
        n_biosamples_positive == 1 & n_pcr_positive >= 2 ~ 0.25,
        # Single PCR detection only
        TRUE ~ 0.1
      ),
      
      # === RELIABILITY ASSESSMENT ===
      # Three-tier system: Reliable, Marginal, Unreliable
      reliability = case_when(
        # RELIABLE: Multiple biosamples AND strong overall detection
        n_biosamples_positive >= 2 & p_empirical >= 0.5 ~ "Reliable",
        
        # RELIABLE: All biosamples positive (regardless of PCR rate)
        n_biosamples_positive == total_biosamples_at_site & n_biosamples_positive >= 2 ~ "Reliable",
        
        # MARGINAL: Multiple biosamples with decent PCR spread (e.g., 2+2 not 3+1)
        n_biosamples_positive >= 2 & p_empirical >= 0.33 & pcr_spread_ratio >= 0.5 ~ "Marginal",
        
        # MARGINAL: Multiple biosamples but uneven spread - needs more scrutiny
        n_biosamples_positive >= 2 & p_empirical >= 0.33 ~ "Marginal",
        
        # MARGINAL: Single biosample but ALL PCRs positive (strong technical signal)
        n_biosamples_positive == 1 & mean_p_pcr == 1 & n_pcr_positive >= 3 ~ "Marginal",
        
        # UNRELIABLE: Everything else
        TRUE ~ "Unreliable"
      ),
      
      # Detection confidence category (simplified)
      confidence = case_when(
        n_biosamples_positive >= 3 & mean_p_pcr >= 0.5 ~ "Very high",
        n_biosamples_positive >= 2 & mean_p_pcr >= 0.3 ~ "High",
        n_biosamples_positive >= 2 ~ "Moderate",
        n_biosamples_positive == 1 & n_pcr_positive >= 2 ~ "Low",
        TRUE ~ "Single detection - verify"
      )
    ) %>%
    select(-p_miss_pcr, -p_miss_biosample, -p_miss_site)
}

# =============================================================================
# SECTION 2: ESTIMATE P AND THETA PER-OTU (DEPRECATED — pools across sites)
# See deprecation notes on individual functions below.
# Site-level values from Section 1 (add_biosample_totals) should be used instead.
# =============================================================================

#' Estimate PCR-level detection probability (p) per OTU
#' Uses only biological samples where OTU was detected
#' @param pcr_detection Output from calc_hierarchical_detection
#' @return OTU-level p estimates
#' 
#' @note DEPRECATED: This function pools across all sites where an OTU was
#'   detected, which is ecologically inappropriate when sites span different
#'   biogeographic regions. Site-level p values in add_biosample_totals()
#'   output (mean_p_pcr column) should be used instead. Retained for
#'   backwards compatibility only.
estimate_p_pcr <- function(pcr_detection) {
  
  pcr_detection %>%
    group_by(OTU) %>%
    summarise(
      # Number of biological samples where detected
      n_biosamples = n(),
      n_sites = n_distinct(Site),
      n_countries = n_distinct(Country),
      
      # Pooled PCR detection rate
      total_pcr_pos = sum(n_pcr_positive),
      total_pcr = sum(n_pcr_total),
      p_pcr_pooled = total_pcr_pos / total_pcr,
      
      # Variability in p across biosamples
      mean_p_pcr = mean(p_pcr),
      sd_p_pcr = sd(p_pcr),
      min_p_pcr = min(p_pcr),
      max_p_pcr = max(p_pcr),
      
      # How often detected in all PCRs of a biosample?
      pct_perfect_pcr = mean(p_pcr == 1) * 100,
      
      .groups = "drop"
    ) %>%
    mutate(
      # Cumulative P(detect | DNA in sample) with 3 PCRs
      p_star_3pcr = 1 - (1 - p_pcr_pooled)^3,
      
      # Confidence in p estimate
      p_confidence = case_when(
        n_biosamples >= 10 ~ "High (n>=10)",
        n_biosamples >= 5 ~ "Moderate (n=5-9)",
        n_biosamples >= 3 ~ "Low (n=3-4)",
        TRUE ~ "Very low (n<3)"
      )
    ) %>%
    arrange(desc(n_biosamples))
}

#' Estimate sample-level detection probability (theta) per OTU
#' What fraction of biological samples at a site detected the OTU?
#' @param site_detection Output from add_biosample_totals
#' @param p_estimates Output from estimate_p_pcr (to get p values)
#' @return OTU-level theta estimates with combined site detection
#' 
#' @note DEPRECATED: This function averages theta across all sites where an
#'   OTU was detected, which is ecologically inappropriate when sites span
#'   different biogeographic regions. Site-level theta values in
#'   add_biosample_totals() output should be used instead. Retained for
#'   backwards compatibility only.
estimate_theta <- function(site_detection, p_estimates = NULL) {
  
  theta_df <- site_detection %>%
    group_by(OTU) %>%
    summarise(
      n_sites = n(),
      n_countries = n_distinct(Country),
      
      # Theta statistics across sites
      mean_theta = mean(theta),
      sd_theta = sd(theta),
      min_theta = min(theta),
      max_theta = max(theta),
      
      # How often detected in all biosamples at a site?
      pct_all_biosamples = mean(theta == 1) * 100,
      
      # Average biosamples per site
      mean_biosamples_per_site = mean(total_biosamples_at_site),
      
      # Also grab mean p from site-level data
      mean_p_pcr_across_sites = mean(mean_p_pcr),
      
      .groups = "drop"
    )
  
  # If p_estimates provided, join to get better p estimate
  if (!is.null(p_estimates)) {
    theta_df <- theta_df %>%
      left_join(
        p_estimates %>% select(OTU, p_pcr_pooled),
        by = "OTU"
      ) %>%
      mutate(
        # Use pooled p if available, otherwise site-level mean
        p_for_calc = coalesce(p_pcr_pooled, mean_p_pcr_across_sites)
      )
  } else {
    theta_df <- theta_df %>%
      mutate(p_for_calc = mean_p_pcr_across_sites)
  }
  
  theta_df %>%
    mutate(
      # Cumulative P(detect | OTU at site) with 3 biosamples (naive)
      theta_star_3samples = 1 - (1 - mean_theta)^3,
      
      # === THE PROPER COMBINED SITE DETECTION ===
      # P(miss all 3 PCRs | DNA in sample)
      p_miss_pcr = (1 - p_for_calc)^3,
      
      # P(miss one biosample) = P(no DNA) + P(DNA but PCR fails)
      p_miss_biosample = (1 - mean_theta) + mean_theta * p_miss_pcr,
      
      # P(miss site entirely) = P(miss all 3 biosamples)
      p_miss_site = p_miss_biosample^3,
      
      # P(detect at site | OTU truly present)
      p_detect_site = 1 - p_miss_site,
      
      # Confidence in estimate
      theta_confidence = case_when(
        n_sites >= 10 ~ "High (n>=10)",
        n_sites >= 5 ~ "Moderate (n=5-9)",
        n_sites >= 3 ~ "Low (n=3-4)",
        TRUE ~ "Very low (n<3)"
      )
    ) %>%
    select(-p_for_calc, -mean_p_pcr_across_sites, -p_miss_pcr, -p_miss_biosample, -p_miss_site) %>%
    arrange(desc(n_sites))
}

# =============================================================================
# SECTION 3: OVERALL SUMMARIES & TAXONOMY JOIN
# NOTE: summarize_detection_probs() is deprecated (relies on cross-site pooling).
# join_taxonomy() remains in active use.
# =============================================================================

#' Join taxonomy to OTU-level results
#' @param otu_results Data frame with OTU column
#' @param ps Phyloseq object (for taxonomy)
#' @return Data frame with taxonomy columns added
join_taxonomy <- function(otu_results, ps) {
  
  tax <- as(tax_table(ps), "matrix") %>% 
    as.data.frame() %>%
    rownames_to_column("OTU")
  
  # Select key taxonomy columns if they exist
  tax_cols <- intersect(c("OTU", "Kingdom", "Phylum", "Class", "Order", 
                          "Family", "Genus", "Species"), colnames(tax))
  tax <- tax[, tax_cols, drop = FALSE]
  
  left_join(otu_results, tax, by = "OTU")
}

#' Get overall detection probability summary
#' @param p_estimates Output from estimate_p_pcr
#' @param theta_estimates Output from estimate_theta
#' @return Combined summary
#' 
#' @note DEPRECATED: Relies on cross-site pooled estimates from estimate_p_pcr()
#'   and estimate_theta(). Use site-level distributions and bootstrap power
#'   analysis (calc_design_power_bootstrap) instead. Retained for backwards
#'   compatibility only.
summarize_detection_probs <- function(p_estimates, theta_estimates) {
  
  # Join p and theta estimates at OTU level
  combined <- p_estimates %>%
    select(OTU, n_biosamples, p_pcr_pooled, p_star_3pcr) %>%
    left_join(
      theta_estimates %>% 
        select(OTU, n_sites, mean_theta, theta_star_3samples, p_detect_site),
      by = "OTU"
    )
  
  list(
    combined = combined,
    
    # Overall summaries (weighted by sample size)
    mean_p_pcr = weighted.mean(p_estimates$p_pcr_pooled, p_estimates$n_biosamples),
    median_p_pcr = median(p_estimates$p_pcr_pooled),
    
    mean_theta = weighted.mean(theta_estimates$mean_theta, theta_estimates$n_sites),
    median_theta = median(theta_estimates$mean_theta),
    
    # Distribution of cumulative detection
    median_p_star = median(p_estimates$p_star_3pcr),
    median_theta_star = median(theta_estimates$theta_star_3samples),
    
    # THE KEY METRIC: combined site-level detection
    median_p_detect_site = median(theta_estimates$p_detect_site),
    mean_p_detect_site = mean(theta_estimates$p_detect_site)
  )
}

#' Join invasive status to OTU-level results
#' @param otu_results Data frame with OTU column (and ideally Species after taxonomy join)
#' @param invasive_df Data frame with Species and invasive status columns
#' @return Data frame with invasive status joined
join_invasive_status <- function(otu_results, invasive_df) {
  
 if (!"Species" %in% colnames(otu_results)) {
    warning("No Species column found - join taxonomy first")
    return(otu_results)
  }
  
  # Expect invasive_df to have Species and status columns
  # Adjust column names as needed based on your invasive status data
  invasive_cols <- intersect(colnames(invasive_df), 
                              c("Species", "EASIN_status", "GISD_status", 
                                "WoRMS_Status", "Final_Status", "is_invasive"))
  
  if (length(invasive_cols) < 2) {
    warning("Invasive status dataframe doesn't have expected columns")
    return(otu_results)
  }
  
  invasive_subset <- invasive_df[, invasive_cols, drop = FALSE] %>%
    distinct()
  
  left_join(otu_results, invasive_subset, by = "Species")
}

# =============================================================================
# SECTION 4: COMBINE DETECTION + INVASIVE STATUS
# =============================================================================

#' Combine detection results with invasive species verification
#' 
#' This is the main function to create a unified analysis-ready dataset
#' 
#' @param detection_results Output from run_detection_analysis()
#' @param invasive_status Output from results_12S$invasive_status (species-country)
#' @param verified_status Output from process_species_list() (WoRMS/GBIF verified)
#' @return List with combined data at OTU and site levels
combine_detection_invasive <- function(detection_results, 
                                        invasive_status = NULL, 
                                        verified_status = NULL) {
  
  cat("Combining detection results with invasive species data...\n")
  
  # Start with site-level detection (has taxonomy)
  site_df <- detection_results$site_detection_tax
  
  # --- Join invasive_status (species × country from EASIN/GISD) ---
  if (!is.null(invasive_status)) {
    cat("  Joining EASIN/GISD status by Species × Country...\n")
    
    # invasive_status has: Species, Location (country), EASIN_status, GISD_status, etc.
    inv_cols <- intersect(colnames(invasive_status),
                          c("Species", "Location", "EASIN_status", "GISD_status", 
                            "is_EASIN_listed", "is_GISD_listed", "is_potentially_invasive"))
    
    if (length(inv_cols) >= 2 && "Species" %in% inv_cols) {
      inv_subset <- invasive_status[, inv_cols, drop = FALSE] %>%
        distinct()
      
      # Rename Location to Country if it exists (use base R to avoid NSE issues)
      if ("Location" %in% colnames(inv_subset)) {
        names(inv_subset)[names(inv_subset) == "Location"] <- "Country"
      }
      
      # Join - handle case where Country column may or may not exist
      if ("Country" %in% colnames(inv_subset) && "Country" %in% colnames(site_df)) {
        site_df <- site_df %>%
          left_join(inv_subset, by = c("Species", "Country"))
      } else {
        # Join by Species only
        site_df <- site_df %>%
          left_join(inv_subset, by = "Species")
      }
      
      cat("    Matched:", sum(!is.na(site_df$EASIN_status) | !is.na(site_df$GISD_status)), 
          "OTU-site combinations\n")
    }
  }
  
  # --- Join verified_status (WoRMS/GBIF verification) ---
  if (!is.null(verified_status)) {
    cat("  Joining WoRMS/GBIF verified status...\n")
    
    ver_cols <- intersect(colnames(verified_status),
                          c("Species", "Location", "AphiaID", "WoRMS_Status", 
                            "WoRMS_Origin", "WoRMS_Invasive", "GBIF_EstablishmentMeans",
                            "Final_Status", "Data_Sources"))
    
    if (length(ver_cols) >= 2 && "Species" %in% ver_cols) {
      ver_subset <- verified_status[, ver_cols, drop = FALSE] %>%
        distinct()
      
      # Rename Location to Country if present (use base R to avoid NSE issues)
      if ("Location" %in% colnames(ver_subset)) {
        names(ver_subset)[names(ver_subset) == "Location"] <- "Country"
        
        # Join by Species + Country (preferred - avoids many-to-many)
        site_df <- site_df %>%
          left_join(ver_subset, by = c("Species", "Country"), 
                    suffix = c("", ".verified"),
                    relationship = "many-to-many")  # Allow if needed
      } else {
        # Join by Species only - aggregate to one row per species first
        ver_agg <- ver_subset %>%
          group_by(Species) %>%
          summarise(
            AphiaID = dplyr::first(na.omit(AphiaID)),
            WoRMS_Status = paste(unique(na.omit(WoRMS_Status)), collapse = "; "),
            WoRMS_Origin = paste(unique(na.omit(WoRMS_Origin)), collapse = "; "),
            WoRMS_Invasive = any(WoRMS_Invasive == TRUE, na.rm = TRUE),
            GBIF_EstablishmentMeans = paste(unique(na.omit(GBIF_EstablishmentMeans)), collapse = "; "),
            Final_Status = case_when(
              any(Final_Status == "INVASIVE", na.rm = TRUE) ~ "INVASIVE",
              any(Final_Status == "INTRODUCED", na.rm = TRUE) ~ "INTRODUCED",
              any(Final_Status == "NATIVE", na.rm = TRUE) ~ "NATIVE",
              TRUE ~ "UNKNOWN"
            ),
            .groups = "drop"
          )
        
        site_df <- site_df %>%
          left_join(ver_agg, by = "Species", suffix = c("", ".verified"))
      }
      
      cat("    Matched:", sum(!is.na(site_df$Final_Status)), "OTU-site combinations\n")
    }
  }
  
  # --- Create unified invasive flag ---
  # Initialize is_invasive as FALSE
  site_df$is_invasive <- FALSE
  
  # Check Final_Status from verified data (most reliable)
  if ("Final_Status" %in% colnames(site_df)) {
    site_df$is_invasive <- site_df$is_invasive | 
      (!is.na(site_df$Final_Status) & site_df$Final_Status == "INVASIVE")
  }
  
 # Check EASIN_status - must have actual invasive status, not "No_EASIN_data" or "Not_listed"
  if ("EASIN_status" %in% colnames(site_df)) {
    # List of values that mean "no data" - not invasive
    easin_no_data <- c("No_EASIN_data", "Not_listed", "Not listed", 
                       "No data", "NA", "", "no_data", "not_listed")
    
    site_df$is_invasive <- site_df$is_invasive | 
      (!is.na(site_df$EASIN_status) & 
       !site_df$EASIN_status %in% easin_no_data &
       !grepl("No_EASIN|Not_listed|No data", site_df$EASIN_status, ignore.case = TRUE))
  }
  
  # Check GISD_status - look for "; Invasive" pattern (not "; Not in GISD")
  if ("GISD_status" %in% colnames(site_df)) {
    site_df$is_invasive <- site_df$is_invasive | 
      (!is.na(site_df$GISD_status) & 
       grepl("; Invasive", site_df$GISD_status, fixed = TRUE))
  }
  
  # Check WoRMS_Invasive flag if present
  if ("WoRMS_Invasive" %in% colnames(site_df)) {
    site_df$is_invasive <- site_df$is_invasive | 
      (!is.na(site_df$WoRMS_Invasive) & site_df$WoRMS_Invasive == TRUE)
  }
  
  # Check is_potentially_invasive if present (from earlier processing)
  if ("is_potentially_invasive" %in% colnames(site_df)) {
    site_df$is_invasive <- site_df$is_invasive | 
      (!is.na(site_df$is_potentially_invasive) & site_df$is_potentially_invasive == TRUE)
  }
  
  # Report what we found
  cat("  Invasive classification summary:\n")
  cat("    Total OTU-site combinations:", nrow(site_df), "\n")
  cat("    Flagged as invasive:", sum(site_df$is_invasive), "\n")
  cat("    Unique invasive species:", n_distinct(site_df$Species[site_df$is_invasive]), "\n")
  
  # Add detection categories
  site_df <- site_df %>%
    mutate(
      # Detection category
      detection_category = case_when(
        confidence %in% c("Very high", "High") ~ "Reliable",
        confidence %in% c("Moderate", "Low") ~ "Moderate",
        TRUE ~ "Unreliable"
      ),
      
      # Combined category for plotting
      detection_invasive_category = paste(
        detection_category,
        ifelse(is_invasive, "Invasive", "Non-invasive"),
        sep = " - "
      )
    )
  
  # --- Summarize at OTU level ---
  cat("  Creating OTU-level summary...\n")
  
  otu_summary <- site_df %>%
    group_by(OTU, Species, Genus, Family, Order, Class, Phylum) %>%
    summarise(
      n_sites = n(),
      n_countries = n_distinct(Country),
      
      # Detection metrics
      mean_theta = mean(theta, na.rm = TRUE),
      mean_p_pcr = mean(mean_p_pcr, na.rm = TRUE),
      mean_p_detect_site = mean(p_detect_site, na.rm = TRUE),
      
      # Confidence distribution
      n_reliable = sum(detection_category == "Reliable"),
      n_moderate = sum(detection_category == "Moderate"),
      n_unreliable = sum(detection_category == "Unreliable"),
      pct_reliable = n_reliable / n() * 100,
      
      # Invasive status
      is_invasive = any(is_invasive, na.rm = TRUE),
      n_sites_invasive = sum(is_invasive, na.rm = TRUE),
      
      # Read depth
      total_reads = sum(total_reads, na.rm = TRUE),
      
      .groups = "drop"
    ) %>%
    mutate(
      overall_reliability = case_when(
        pct_reliable >= 50 ~ "Mostly reliable",
        n_reliable >= 1 ~ "Sometimes reliable",
        TRUE ~ "Unreliable"
      )
    ) %>%
    arrange(desc(is_invasive), desc(n_sites))
  
  # --- Summary stats ---
  cat("\n  === SUMMARY ===\n")
  cat("  Total OTU × site detections:", nrow(site_df), "\n")
  cat("  Unique OTUs:", n_distinct(site_df$OTU), "\n")
  cat("  Invasive OTU-site detections:", sum(site_df$is_invasive), "\n")
  cat("  Unique invasive OTUs:", sum(otu_summary$is_invasive), "\n")
  
  cat("\n  Detection reliability:\n")
  print(table(site_df$detection_category))
  
  cat("\n  Invasive × Reliability:\n
")
  print(table(site_df$detection_category, site_df$is_invasive))
  
  return(list(
    site_level = site_df,
    otu_summary = otu_summary,
    # Quick access to key subsets
    invasive_detections = site_df %>% filter(is_invasive),
    reliable_invasive = site_df %>% filter(is_invasive & detection_category == "Reliable"),
    unreliable_invasive = site_df %>% filter(is_invasive & detection_category == "Unreliable")
  ))
}

# =============================================================================
# SECTION 5: VISUALIZATION FUNCTIONS
# =============================================================================

# -----------------------------------------------------------------------------
# COMBINING MULTIPLE MARKERS
# -----------------------------------------------------------------------------

#' Combine detection results from multiple markers
#' 
#' @param ... Named list of combined objects (e.g., `12S = combined_12S, 18S = combined_18S`)
#' @return Combined list with marker column added
#' @examples
#' combined_all <- combine_markers(
#'   `12S` = combined_12S_filtered,
#'   `18S` = combined_18S_filtered,
#'   COI = combined_COI_filtered
#' )
combine_markers <- function(...) {
  
  marker_list <- list(...)
  
  if (length(marker_list) == 0) {
    stop("Provide at least one combined object")
  }
  
  # Get marker names
  marker_names <- names(marker_list)
  if (is.null(marker_names) || any(marker_names == "")) {
    stop("All arguments must be named (e.g., `12S` = combined_12S)")
  }
  
  # Combine site_level data with marker column
  site_level_combined <- bind_rows(
    lapply(marker_names, function(m) {
      df <- marker_list[[m]]$site_level
      if (!is.null(df) && nrow(df) > 0) {
        df$Marker <- m
        return(df)
      }
      return(NULL)
    })
  )
  
  # Combine otu_summary if present
  otu_summary_combined <- bind_rows(
    lapply(marker_names, function(m) {
      df <- marker_list[[m]]$otu_summary
      if (!is.null(df) && nrow(df) > 0) {
        df$Marker <- m
        return(df)
      }
      return(NULL)
    })
  )
  
  cat("Combined", length(marker_names), "markers:", paste(marker_names, collapse = ", "), "\n")
  cat("Total OTU-site combinations:", nrow(site_level_combined), "\n")
  cat("By marker:\n")
  print(table(site_level_combined$Marker))
  
  list(
    site_level = site_level_combined,
    otu_summary = otu_summary_combined,
    markers = marker_names
  )
}

#' Plot site detection with marker differentiation
#' 
#' @param combined_markers Output from combine_markers()
#' @param site Site name
#' @param x_metric Metric for x-axis
#' @param y_metric Metric for y-axis
#' @return ggplot object
plot_site_by_marker <- function(combined_markers, site, 
                                 x_metric = "total_reads", y_metric = "p_empirical") {
  
  site_data <- combined_markers$site_level %>%
    filter(Site == !!site)
  
  if (nrow(site_data) == 0) {
    message("No data for site: ", site)
    return(NULL)
  }
  
  country <- unique(site_data$Country)[1]
  
  # Metric labels
  get_label <- function(m) {
    switch(m,
           "p_empirical" = "Detection rate",
           "theta" = "θ (biosample detection)",
           "confidence_score" = "Confidence score",
           "total_reads" = "Read count",
           m)
  }
  
  # Add labels and sort so invasive plotted last
  site_data <- site_data %>%
    mutate(
      invasive_label = ifelse(is_invasive, "Invasive", "Non-invasive"),
      label = ifelse(is_invasive, Species, NA_character_)
    ) %>%
    arrange(is_invasive)
  
  n_otus <- nrow(site_data)
  n_invasive <- sum(site_data$is_invasive)
  n_markers <- n_distinct(site_data$Marker)
  
  # Build plot with marker shapes
  p <- ggplot(site_data, aes(x = .data[[x_metric]], y = .data[[y_metric]], 
                              color = invasive_label, shape = Marker)) +
    geom_point(aes(alpha = invasive_label), size = 3) +
    scale_color_manual(values = c("Invasive" = "#e74c3c", "Non-invasive" = "#3498db")) +
    scale_alpha_manual(values = c("Invasive" = 1, "Non-invasive" = 0.15), guide = "none") +
    labs(x = get_label(x_metric),
         y = get_label(y_metric),
         title = site,
         subtitle = paste0(country, " | ", n_otus, " OTUs, ", n_invasive, " invasive, ",
                           n_markers, " markers"),
         color = "", shape = "Marker") +
    theme_minimal()
  
  # Handle axis scales
  if (x_metric == "total_reads") {
    p <- p + scale_x_log10(labels = if(requireNamespace("scales", quietly = TRUE)) scales::comma else waiver())
  } else {
    p <- p + coord_cartesian(xlim = c(0, 1))
  }
  
  if (y_metric != "total_reads" && x_metric != "total_reads") {
    p <- p + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
  } else if (y_metric != "total_reads") {
    p <- p + coord_cartesian(ylim = c(0, 1))
  }
  
  # Add labels for invasive
  if (requireNamespace("ggrepel", quietly = TRUE) && n_invasive > 0) {
    p <- p + ggrepel::geom_text_repel(
      aes(label = label),
      size = 3,
      color = "#e74c3c",
      fontface = "bold",
      max.overlaps = Inf,
      na.rm = TRUE
    )
  }
  
  p
}

#' Faceted plot by marker for a single site
#' 
#' @param combined_markers Output from combine_markers()
#' @param site Site name
#' @param metric Metric for y-axis
#' @return ggplot object
plot_site_facet_markers <- function(combined_markers, site, metric = "p_empirical") {
  
  site_data <- combined_markers$site_level %>%
    filter(Site == !!site)
  
  if (nrow(site_data) == 0) {
    message("No data for site: ", site)
    return(NULL)
  }
  
  country <- unique(site_data$Country)[1]
  
  metric_lab <- switch(metric,
                       "p_empirical" = "Detection rate",
                       "theta" = "θ (biosample detection)",
                       "confidence_score" = "Confidence score",
                       metric)
  
  site_data <- site_data %>%
    mutate(
      invasive_label = ifelse(is_invasive, "Invasive", "Non-invasive"),
      plot_value = .data[[metric]]
    ) %>%
    arrange(is_invasive)
  
  n_otus <- nrow(site_data)
  n_invasive <- sum(site_data$is_invasive)
  
  ggplot(site_data, aes(x = total_reads, y = plot_value, color = invasive_label)) +
    geom_point(aes(alpha = invasive_label), size = 2) +
    facet_wrap(~Marker, ncol = 3) +
    scale_x_log10(labels = if(requireNamespace("scales", quietly = TRUE)) scales::comma else waiver()) +
    scale_color_manual(values = c("Invasive" = "#e74c3c", "Non-invasive" = "#3498db")) +
    scale_alpha_manual(values = c("Invasive" = 1, "Non-invasive" = 0.15), guide = "none") +
    labs(x = "Read count (log scale)",
         y = metric_lab,
         title = site,
         subtitle = paste0(country, " | ", n_otus, " OTUs, ", n_invasive, " invasive"),
         color = "") +
    theme_minimal() +
    theme(strip.text = element_text(face = "bold")) +
    coord_cartesian(ylim = c(0, 1))
}

#' Plot p vs theta for a single country
#' @param combined Output from combine_detection_invasive()
#' @param country Country to plot
#' @param label_invasive Label invasive species
#' @param show_errorbars Show error bars (SD across sites)
#' @return ggplot object
plot_p_vs_theta_country <- function(combined, country, label_invasive = TRUE, 
                                     show_errorbars = TRUE) {
  
  country_data <- combined$site_level %>%
    filter(Country == !!country)
  
  if (nrow(country_data) == 0) {
    message("No data for country: ", country)
    return(NULL)
  }
  
  # Aggregate to OTU level within this country
  otu_country <- country_data %>%
    group_by(OTU, Species) %>%
    summarise(
      n_sites = n(),
      mean_theta = mean(theta, na.rm = TRUE),
      sd_theta = sd(theta, na.rm = TRUE),
      mean_p_pcr = mean(mean_p_pcr, na.rm = TRUE),
      sd_p_pcr = sd(mean_p_pcr, na.rm = TRUE),
      mean_p_detect = mean(p_detect_site, na.rm = TRUE),
      is_invasive = any(is_invasive, na.rm = TRUE),
      total_reads = sum(total_reads, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      # Replace NA SD with 0 (single observations)
      sd_theta = replace_na(sd_theta, 0),
      sd_p_pcr = replace_na(sd_p_pcr, 0),
      invasive_label = ifelse(is_invasive, "Invasive", "Non-invasive"),
      # Label invasive species OR species at multiple sites
      label = ifelse(is_invasive | n_sites >= 2, Species, NA)
    )
  
  n_sites_country <- n_distinct(country_data$Site)
  n_invasive <- sum(otu_country$is_invasive)
  
  # Base plot
  p <- ggplot(otu_country, aes(x = mean_p_pcr, y = mean_theta))
  
  # Add error bars if requested and if there's variation
  if (show_errorbars) {
    # Only show error bars for OTUs with >1 site (otherwise SD=0)
    otu_multi <- otu_country %>% filter(n_sites > 1)
    
    if (nrow(otu_multi) > 0) {
      p <- p + 
        geom_errorbar(data = otu_multi,
                      aes(ymin = pmax(0, mean_theta - sd_theta), 
                          ymax = pmin(1, mean_theta + sd_theta),
                          color = invasive_label),
                      width = 0.02, alpha = 0.4) +
        geom_errorbarh(data = otu_multi,
                       aes(xmin = pmax(0, mean_p_pcr - sd_p_pcr), 
                           xmax = pmin(1, mean_p_pcr + sd_p_pcr),
                           color = invasive_label),
                       height = 0.02, alpha = 0.4)
    }
  }
  
  # Add points and labels using ggrepel for invasive species
  # Check if ggrepel is available
  use_repel <- requireNamespace("ggrepel", quietly = TRUE)
  
  p <- p +
    geom_point(aes(size = n_sites, color = invasive_label, alpha = invasive_label))
  
  if (use_repel) {
    p <- p +
      # Labels for non-invasive (with overlap checking, regular geom_text)
      geom_text(data = otu_country %>% filter(!is_invasive),
                aes(label = label), hjust = -0.1, vjust = 0.5, size = 2.5, 
                check_overlap = TRUE, na.rm = TRUE, color = "#7f8c8d") +
      # Labels for invasive using ggrepel (always show, no overlap)
      ggrepel::geom_text_repel(
        data = otu_country %>% filter(is_invasive),
        aes(label = Species), 
        size = 3, 
        color = "#e74c3c", 
        fontface = "bold",
        box.padding = 0.5,
        point.padding = 0.3,
        max.overlaps = Inf,
        na.rm = TRUE
      )
  } else {
    p <- p +
      # Fallback without ggrepel
      geom_text(data = otu_country %>% filter(!is_invasive),
                aes(label = label), hjust = -0.1, vjust = 0.5, size = 2.5, 
                check_overlap = TRUE, na.rm = TRUE, color = "#7f8c8d") +
      geom_text(data = otu_country %>% filter(is_invasive),
                aes(label = Species), hjust = -0.1, vjust = 0.5, size = 3, 
                check_overlap = FALSE, na.rm = TRUE, color = "#e74c3c", fontface = "bold")
  }
  
  p <- p +
    scale_color_manual(values = c("Invasive" = "#e74c3c", "Non-invasive" = "#7f8c8d")) +
    scale_alpha_manual(values = c("Invasive" = 1, "Non-invasive" = 0.4)) +
    scale_size_continuous(range = c(2, 8), name = "N sites", 
                          breaks = seq(1, n_sites_country, by = max(1, n_sites_country %/% 4))) +
    geom_hline(yintercept = 0.5, linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 0.5, linetype = "dashed", alpha = 0.5) +
    annotate("rect", xmin = 0.5, xmax = 1, ymin = 0.5, ymax = 1, 
             fill = "green", alpha = 0.1) +
    labs(x = "Mean p (PCR detection rate)",
         y = "Mean θ (sample detection rate)",
         title = paste("Detection Probability -", country),
         subtitle = paste(n_sites_country, "sites,", n_distinct(otu_country$OTU), "OTUs,",
                          n_invasive, "invasive species"),
         color = "") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_minimal() +
    guides(alpha = "none")
  
  p
}

#' Plot invasive species for a single country
#' @param combined Output from combine_detection_invasive()
#' @param country Country to plot
#' @return ggplot object
plot_invasive_country <- function(combined, country) {
  
  inv_data <- combined$site_level %>%
    filter(Country == !!country, is_invasive == TRUE)
  
  if (nrow(inv_data) == 0) {
    message("No invasive species detected in: ", country)
    return(NULL)
  }
  
  # Count sites per species with detection confidence
  inv_summary <- inv_data %>%
    group_by(Species) %>%
    summarise(
      n_sites = n_distinct(Site),
      n_reliable = sum(detection_category == "Reliable"),
      n_moderate = sum(detection_category == "Moderate"),
      n_unreliable = sum(detection_category == "Unreliable"),
      mean_p_detect = mean(p_detect_site, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = c(n_reliable, n_moderate, n_unreliable),
                 names_to = "confidence", values_to = "count") %>%
    mutate(
      confidence = factor(confidence, 
                          levels = c("n_reliable", "n_moderate", "n_unreliable"),
                          labels = c("Reliable", "Moderate", "Unreliable"))
    ) %>%
    filter(count > 0)
  
  n_sites_total <- n_distinct(combined$site_level$Site[combined$site_level$Country == country])
  
  ggplot(inv_summary, aes(x = reorder(Species, -n_sites), y = count, fill = confidence)) +
    geom_col() +
    scale_fill_manual(values = c("Reliable" = "#2ecc71", 
                                  "Moderate" = "#f39c12", 
                                  "Unreliable" = "#e74c3c")) +
    labs(x = "Species", 
         y = "Number of sites detected",
         title = paste("Invasive Species -", country),
         subtitle = paste("Total sites in country:", n_sites_total),
         fill = "Detection\nConfidence") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#' Summary plot: Detection reliability across all countries
#' Shows proportion reliable/moderate/unreliable per country
#' @param combined Output from combine_detection_invasive()
#' @return ggplot object
plot_country_summary <- function(combined) {
  
  country_summary <- combined$site_level %>%
    group_by(Country) %>%
    summarise(
      n_sites = n_distinct(Site),
      n_otus = n_distinct(OTU),
      n_invasive_otus = n_distinct(OTU[is_invasive]),
      n_reliable = sum(detection_category == "Reliable"),
      n_moderate = sum(detection_category == "Moderate"),
      n_unreliable = sum(detection_category == "Unreliable"),
      pct_reliable = n_reliable / n() * 100,
      mean_p_detect = mean(p_detect_site, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_sites))
  
  # Pivot for stacked bar
  country_long <- country_summary %>%
    select(Country, n_sites, n_reliable, n_moderate, n_unreliable) %>%
    pivot_longer(cols = c(n_reliable, n_moderate, n_unreliable),
                 names_to = "category", values_to = "count") %>%
    mutate(
      category = factor(category,
                        levels = c("n_reliable", "n_moderate", "n_unreliable"),
                        labels = c("Reliable", "Moderate", "Unreliable"))
    )
  
  ggplot(country_long, aes(x = reorder(Country, -n_sites), y = count, fill = category)) +
    geom_col() +
    geom_text(data = country_summary, 
              aes(x = Country, y = n_reliable + n_moderate + n_unreliable + 5, 
                  label = paste0(n_sites, " sites")),
              inherit.aes = FALSE, size = 3, vjust = 0) +
    scale_fill_manual(values = c("Reliable" = "#2ecc71", 
                                  "Moderate" = "#f39c12", 
                                  "Unreliable" = "#e74c3c")) +
    labs(x = "Country", 
         y = "Number of OTU × site detections",
         title = "Detection Reliability Summary by Country",
         fill = "Detection\nConfidence") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#' Plot all invasive species across dataset with site counts
#' @param combined Output from combine_detection_invasive()
#' @param min_sites Minimum sites to include (default 1)
#' @return ggplot object
plot_all_invasives <- function(combined, min_sites = 1) {
  
  inv_summary <- combined$invasive_detections %>%
    group_by(Species) %>%
    summarise(
      n_countries = n_distinct(Country),
      n_sites = n_distinct(Site),
      countries = paste(unique(Country), collapse = ", "),
      n_reliable = sum(detection_category == "Reliable"),
      n_moderate = sum(detection_category == "Moderate"),
      n_unreliable = sum(detection_category == "Unreliable"),
      mean_p_detect = mean(p_detect_site, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(n_sites >= min_sites) %>%
    arrange(desc(n_sites))
  
  # Take top 30 for readability
  inv_top <- head(inv_summary, 30)
  
  inv_long <- inv_top %>%
    select(Species, n_reliable, n_moderate, n_unreliable) %>%
    pivot_longer(cols = c(n_reliable, n_moderate, n_unreliable),
                 names_to = "category", values_to = "count") %>%
    mutate(
      category = factor(category,
                        levels = c("n_reliable", "n_moderate", "n_unreliable"),
                        labels = c("Reliable", "Moderate", "Unreliable"))
    )
  
  ggplot(inv_long, aes(x = reorder(Species, -count), y = count, fill = category)) +
    geom_col() +
    scale_fill_manual(values = c("Reliable" = "#2ecc71", 
                                  "Moderate" = "#f39c12", 
                                  "Unreliable" = "#e74c3c")) +
    labs(x = "Species", 
         y = "Number of site detections",
         title = "Top Invasive Species Detections",
         subtitle = paste("Showing species detected at >=", min_sites, "sites"),
         fill = "Detection\nConfidence") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
}

#' Create heatmap of invasive species × country
#' @param combined Output from combine_detection_invasive()
#' @param min_detections Minimum total detections to include species
#' @return ggplot object
plot_invasive_heatmap <- function(combined, min_detections = 2) {
  
  inv_matrix <- combined$invasive_detections %>%
    group_by(Species, Country) %>%
    summarise(
      n_sites = n_distinct(Site),
      mean_confidence = mean(case_when(
        detection_category == "Reliable" ~ 3,
        detection_category == "Moderate" ~ 2,
        TRUE ~ 1
      )),
      .groups = "drop"
    )
  
  # Filter to species with enough detections
  species_keep <- inv_matrix %>%
    group_by(Species) %>%
    summarise(total = sum(n_sites), .groups = "drop") %>%
    filter(total >= min_detections) %>%
    pull(Species)
  
  inv_matrix <- inv_matrix %>%
    filter(Species %in% species_keep)
  
  if (nrow(inv_matrix) == 0) {
    message("No species with >= ", min_detections, " detections")
    return(NULL)
  }
  
  ggplot(inv_matrix, aes(x = Country, y = reorder(Species, n_sites), fill = n_sites)) +
    geom_tile(color = "white") +
    geom_text(aes(label = n_sites), size = 3) +
    scale_fill_gradient(low = "#fee8c8", high = "#e34a33", name = "N sites") +
    labs(x = "Country", y = "Species",
         title = "Invasive Species Detection Matrix",
         subtitle = "Number of sites where each species was detected") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 8))
}

#' Get list of countries in the dataset
#' @param combined Output from combine_detection_invasive()
#' @return Vector of country names
list_countries <- function(combined) {
  sort(unique(combined$site_level$Country))
}

#' Get list of sites in a country
#' @param combined Output from combine_detection_invasive()
#' @param country Country name
#' @return Vector of site names
list_sites <- function(combined, country) {
  combined$site_level %>%
    filter(Country == !!country) %>%
    pull(Site) %>%
    unique() %>%
    sort()
}

# =============================================================================
# SECTION 6: SITE-LEVEL DETECTION PLOTS
# =============================================================================

#' Plot two detection metrics against each other for a single site
#' 
#' @param combined Output from combine_detection_invasive()
#' @param site Site name
#' @param x_metric Metric for x-axis: "p_empirical", "theta", "confidence_score", "reads"
#' @param y_metric Metric for y-axis: "p_empirical", "theta", "confidence_score", "reads"
#' @return ggplot object
plot_site_xy <- function(combined, site, x_metric = "total_reads", y_metric = "p_empirical") {
 
  site_data <- combined$site_level %>%
    filter(Site == !!site)
  
  if (nrow(site_data) == 0) {
    message("No data for site: ", site)
    return(NULL)
  }
  
  country <- unique(site_data$Country)[1]
  n_biosamples <- unique(site_data$total_biosamples_at_site)[1]
  n_pcr <- unique(site_data$total_pcr_at_site)[1]
  
  # Metric labels
  get_label <- function(m) {
    switch(m,
           "p_empirical" = paste0("Detection rate (PCR+ / ", n_pcr, ")"),
           "theta" = paste0("θ (biosamples+ / ", n_biosamples, ")"),
           "confidence_score" = "Confidence score",
           "total_reads" = "Read count",
           m)
  }
  
  # Add invasive label and sort so invasive plotted LAST (on top)
  site_data <- site_data %>%
    mutate(
      invasive_label = ifelse(is_invasive, "Invasive", "Non-invasive"),
      label = ifelse(is_invasive, Species, NA_character_)
    ) %>%
    arrange(is_invasive)  # FALSE first, TRUE last = invasive on top
  
  n_otus <- nrow(site_data)
  n_invasive <- sum(site_data$is_invasive)
  
  # Build plot - non-invasive first (transparent), then invasive (solid)
  p <- ggplot(site_data, aes(x = .data[[x_metric]], y = .data[[y_metric]], color = invasive_label)) +
    geom_point(aes(alpha = invasive_label), size = 3) +
    scale_color_manual(values = c("Invasive" = "#e74c3c", "Non-invasive" = "#3498db")) +
    scale_alpha_manual(values = c("Invasive" = 1, "Non-invasive" = 0.15), guide = "none") +
    labs(x = get_label(x_metric),
         y = get_label(y_metric),
         title = site,
         subtitle = paste0(country, " | ", n_otus, " OTUs, ", n_invasive, " invasive"),
         color = "") +
    theme_minimal()
  
  # Log scale for reads, fixed limits for proportions
  x_is_reads <- x_metric == "total_reads"
  y_is_reads <- y_metric == "total_reads"
  
  if (x_is_reads) {
    p <- p + scale_x_log10(labels = if(requireNamespace("scales", quietly = TRUE)) scales::comma else waiver())
  }
  if (y_is_reads) {
    p <- p + scale_y_log10(labels = if(requireNamespace("scales", quietly = TRUE)) scales::comma else waiver())
  }
  
  # Set axis limits for proportion metrics
  if (!x_is_reads && !y_is_reads) {
    p <- p + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
  } else if (!x_is_reads) {
    p <- p + coord_cartesian(xlim = c(0, 1))
  } else if (!y_is_reads) {
    p <- p + coord_cartesian(ylim = c(0, 1))
  }
  
  # Add labels for invasive species using ggrepel
  if (requireNamespace("ggrepel", quietly = TRUE) && n_invasive > 0) {
    p <- p + ggrepel::geom_text_repel(
      aes(label = label),
      size = 3,
      color = "#e74c3c",
      fontface = "bold",
      max.overlaps = Inf,
      na.rm = TRUE
    )
  }
  
  p
}

#' Plot detection metrics for a single site
#' Each point is one OTU at that site
#' 
#' @param combined Output from combine_detection_invasive()
#' @param site Site name
#' @param metric Which metric for y-axis: "p_empirical", "theta", "confidence_score"
#' @return ggplot object
plot_site_detection <- function(combined, site, metric = "p_empirical") {
 
  site_data <- combined$site_level %>%
    filter(Site == !!site)
  
  if (nrow(site_data) == 0) {
    message("No data for site: ", site)
    return(NULL)
  }
  
  country <- unique(site_data$Country)[1]
  n_biosamples <- unique(site_data$total_biosamples_at_site)[1]
  n_pcr <- unique(site_data$total_pcr_at_site)[1]
  
  # Select metric
  metric_col <- switch(metric,
                       "p_empirical" = "p_empirical",
                       "theta" = "theta",
                       "confidence_score" = "confidence_score",
                       "p_empirical")
  
  metric_lab <- switch(metric,
                       "p_empirical" = paste0("Detection rate (", "PCR+ / ", n_pcr, ")"),
                       "theta" = paste0("θ (biosamples+ / ", n_biosamples, ")"),
                       "confidence_score" = "Confidence score",
                       "Detection rate")
  
  site_data$plot_value <- site_data[[metric_col]]
  
  # Add labels and sort so invasive plotted LAST (on top)
  site_data <- site_data %>%
    mutate(
      invasive_label = ifelse(is_invasive, "Invasive", "Non-invasive"),
      label = ifelse(is_invasive, Species, NA_character_)
    ) %>%
    arrange(is_invasive)  # FALSE first, TRUE last = invasive on top
  
  n_otus <- nrow(site_data)
  n_invasive <- sum(site_data$is_invasive)
  
  # Summary stats
  if (n_invasive > 0) {
    inv_median <- median(site_data$plot_value[site_data$is_invasive], na.rm = TRUE)
    noninv_median <- median(site_data$plot_value[!site_data$is_invasive], na.rm = TRUE)
    subtitle <- paste0(n_otus, " OTUs | ", n_invasive, " invasive (median ", metric, "=", 
                       round(inv_median, 2), ") vs non-invasive (", 
                       round(noninv_median, 2), ")")
  } else {
    noninv_median <- median(site_data$plot_value, na.rm = TRUE)
    subtitle <- paste0(n_otus, " OTUs | No invasive species | median ", metric, "=", 
                       round(noninv_median, 2))
  }
  
  # Create XY plot: reads (x) vs detection metric (y)
  # Non-invasive very transparent, invasive solid, invasive plotted on top
  p <- ggplot(site_data, aes(x = total_reads, y = plot_value, color = invasive_label)) +
    geom_point(aes(alpha = invasive_label), size = 3) +
    scale_x_log10(labels = scales::comma) +
    scale_color_manual(values = c("Invasive" = "#e74c3c", "Non-invasive" = "#3498db")) +
    scale_alpha_manual(values = c("Invasive" = 1, "Non-invasive" = 0.15), guide = "none") +
    labs(x = "Read count (log scale)",
         y = metric_lab,
         title = site,
         subtitle = subtitle,
         color = "") +
    theme_minimal() +
    coord_cartesian(ylim = c(0, 1))
  
  # Add labels for invasive species using ggrepel
  if (requireNamespace("ggrepel", quietly = TRUE) && n_invasive > 0) {
    p <- p + ggrepel::geom_text_repel(
      aes(label = label),
      size = 3,
      color = "#e74c3c",
      fontface = "bold",
      max.overlaps = Inf,
      na.rm = TRUE
    )
  }
  
  p
}

#' Multi-panel plot showing all metrics for a single site
#' 
#' @param combined Output from combine_detection_invasive()
#' @param site Site name
#' @return ggplot object (faceted)
plot_site_all_metrics <- function(combined, site) {
  
  site_data <- combined$site_level %>%
    filter(Site == !!site)
  
  if (nrow(site_data) == 0) {
    message("No data for site: ", site)
    return(NULL)
  }
  
  country <- unique(site_data$Country)[1]
  n_biosamples <- unique(site_data$total_biosamples_at_site)[1]
  n_pcr <- unique(site_data$total_pcr_at_site)[1]
  n_otus <- nrow(site_data)
  n_invasive <- sum(site_data$is_invasive)
  
  # Pivot to long format and sort so invasive plotted on top
  site_long <- site_data %>%
    arrange(is_invasive) %>%  # Non-invasive first, invasive last
    select(OTU, Species, total_reads, is_invasive,
           p_empirical, theta, confidence_score) %>%
    pivot_longer(cols = c(p_empirical, theta, confidence_score),
                 names_to = "metric", values_to = "value") %>%
    mutate(
      invasive_label = ifelse(is_invasive, "Invasive", "Non-invasive"),
      metric = factor(metric,
                      levels = c("p_empirical", "theta", "confidence_score"),
                      labels = c(
                        paste0("Detection rate (PCR+ / ", n_pcr, ")"),
                        paste0("θ (biosamples+ / ", n_biosamples, ")"),
                        "Confidence score"
                      ))
    )
  
  ggplot(site_long, aes(x = total_reads, y = value, color = invasive_label)) +
    geom_point(aes(alpha = invasive_label), size = 2) +
    facet_wrap(~metric, ncol = 3) +
    scale_x_log10(labels = scales::comma) +
    scale_color_manual(values = c("Invasive" = "#e74c3c", "Non-invasive" = "#3498db")) +
    scale_alpha_manual(values = c("Invasive" = 1, "Non-invasive" = 0.15), guide = "none") +
    labs(x = "Read count (log scale)",
         y = "Value",
         title = site,
         subtitle = paste0(country, " | ", n_otus, " OTUs, ", n_invasive, " invasive"),
         color = "") +
    theme_minimal() +
    theme(strip.text = element_text(face = "bold")) +
    coord_cartesian(ylim = c(0, 1))
}

#' Summary table for a single site
#' 
#' @param combined Output from combine_detection_invasive()
#' @param site Site name
#' @return Data frame with OTU detection summary
get_site_summary <- function(combined, site) {
  
  # Get available columns
  available_cols <- colnames(combined$site_level)
  
  # Base columns
  base_cols <- c("OTU", "Species", "Class", "Order", "Family",
                 "n_biosamples_positive", "total_biosamples_at_site",
                 "n_pcr_positive", "total_pcr_at_site",
                 "theta", "p_empirical", "mean_p_pcr", "confidence_score", "confidence",
                 "total_reads", "is_invasive", "EASIN_status", "GISD_status")
  
  # New columns (may not exist in older data)
  new_cols <- c("pcr_distribution", "pcr_spread_ratio", "reliability")
  
  # Select columns that exist
  select_cols <- c(base_cols[base_cols %in% available_cols], 
                   new_cols[new_cols %in% available_cols])
  
  site_data <- combined$site_level %>%
    filter(Site == !!site) %>%
    select(any_of(select_cols)) %>%
    mutate(
      detection_summary = paste0(n_biosamples_positive, "/", total_biosamples_at_site, 
                                  " samples, ", n_pcr_positive, "/", total_pcr_at_site, " PCRs")
    ) %>%
    arrange(desc(is_invasive), desc(confidence_score), desc(total_reads))
  
  site_data
}

#' Loop through all sites in a country and create plots
#' 
#' @param combined Output from combine_detection_invasive()
#' @param country Country name
#' @param output_dir Directory to save plots (NULL to just display)
#' @param metric Which metric to plot
#' @return List of plots (invisibly)
plot_all_sites_in_country <- function(combined, country, output_dir = NULL, 
                                       metric = "p_empirical") {
  
  sites <- list_sites(combined, country)
  cat("Plotting", length(sites), "sites in", country, "\n")
  
  plots <- list()
  
  for (site in sites) {
    p <- plot_site_detection(combined, site, metric)
    
    if (!is.null(p)) {
      print(p)
      plots[[site]] <- p
      
      if (!is.null(output_dir)) {
        safe_name <- gsub("[^A-Za-z0-9_-]", "_", site)
        ggsave(file.path(output_dir, paste0(safe_name, "_", metric, ".png")),
               p, width = 10, height = 6)
      }
    }
  }
  
  invisible(plots)
}

#' Compare invasive detection across all sites (one dot per site)
#' 
#' @param combined Output from combine_detection_invasive()
#' @param metric Which metric to compare
#' @return ggplot object
plot_invasive_across_sites <- function(combined, metric = "p_empirical") {
  
  # Get invasive species detections only
  inv_data <- combined$site_level %>%
    filter(is_invasive)
  
  if (nrow(inv_data) == 0) {
    message("No invasive species detected")
    return(NULL)
  }
  
  metric_col <- switch(metric,
                       "p_empirical" = "p_empirical",
                       "theta" = "theta", 
                       "confidence_score" = "confidence_score",
                       "p_empirical")
  
  metric_lab <- switch(metric,
                       "p_empirical" = "Detection rate",
                       "theta" = "θ (biosample detection)",
                       "confidence_score" = "Confidence score",
                       "Detection rate")
  
  inv_data$plot_value <- inv_data[[metric_col]]
  
  # One row per species × site
  ggplot(inv_data, aes(x = total_reads, y = plot_value, color = Species)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_x_log10(labels = scales::comma) +
    labs(x = "Read count (log scale)",
         y = metric_lab,
         title = "Invasive Species Detection Across Sites",
         subtitle = paste("Each dot = one species at one site |", 
                          n_distinct(inv_data$Species), "species,",
                          n_distinct(inv_data$Site), "sites"),
         color = "Species") +
    theme_minimal() +
    coord_cartesian(ylim = c(0, 1)) +
    theme(legend.position = "right")
}

# =============================================================================
# SECTION 7: AGGREGATED OVERVIEW PLOTS (use with caution)
# =============================================================================
# NOTE: These plots aggregate across sites/countries which loses ecological
# meaning. Use for quick data exploration only. For analysis, use the 
# site-level plots in Section 6.
# =============================================================================

#' Strip/dot plot of detection probability by invasive status
#' Shows where invasive species fall in the overall detection distribution
#' NOTE: Aggregates across sites - use for overview only
#' 
#' @param combined Output from combine_detection_invasive()
#' @param country Country to plot (NULL for all data)
#' @param metric Which metric: "p_empirical", "theta", "confidence_score", "reads", "n_biosamples"
#' @param log_reads Log-transform reads (only for metric = "reads")
#' @return ggplot object
plot_detection_distribution <- function(combined, country = NULL, 
                                         metric = "p_empirical",
                                         log_reads = TRUE) {
  
  # Filter by country if specified
  if (!is.null(country)) {
    site_data <- combined$site_level %>% filter(Country == !!country)
    title_suffix <- paste("-", country)
  } else {
    site_data <- combined$site_level
    title_suffix <- "- All Countries"
  }
  
  if (nrow(site_data) == 0) {
    message("No data available")
    return(NULL)
  }
  
  # Aggregate to OTU level
  otu_data <- site_data %>%
    group_by(OTU, Species) %>%
    summarise(
      n_sites = n(),
      mean_theta = mean(theta, na.rm = TRUE),
      mean_p_empirical = mean(p_empirical, na.rm = TRUE),
      mean_confidence_score = mean(confidence_score, na.rm = TRUE),
      total_biosamples_positive = sum(n_biosamples_positive, na.rm = TRUE),
      total_reads = sum(total_reads, na.rm = TRUE),
      prop_reads = total_reads / sum(site_data$total_reads, na.rm = TRUE),
      is_invasive = any(is_invasive, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      invasive_label = ifelse(is_invasive, "Invasive", "Non-invasive"),
      invasive_label = factor(invasive_label, levels = c("Invasive", "Non-invasive"))
    )
  
  # Select metric
  metric_map <- list(
    "p_empirical" = list(col = "mean_p_empirical", 
                         lab = "Empirical detection rate (PCR detections / total PCRs)",
                         title = "Empirical Detection Rate"),
    "theta" = list(col = "mean_theta", 
                   lab = "θ (proportion of biosamples with detection)",
                   title = "Biological Sample Detection"),
    "confidence_score" = list(col = "mean_confidence_score",
                              lab = "Confidence score (weights biological > PCR replication)",
                              title = "Detection Confidence Score"),
    "n_biosamples" = list(col = "total_biosamples_positive",
                          lab = "Total biological samples with detection",
                          title = "Biological Sample Detections"),
    "reads" = list(col = "total_reads",
                   lab = ifelse(log_reads, "Total reads (log scale)", "Total reads"),
                   title = "Read Abundance")
  )
  
  if (!metric %in% names(metric_map)) {
    stop("metric must be one of: p_empirical, theta, confidence_score, n_biosamples, reads")
  }
  
  m <- metric_map[[metric]]
  otu_data$plot_value <- otu_data[[m$col]]
  
  # Summary stats
  inv_median <- median(otu_data$plot_value[otu_data$is_invasive], na.rm = TRUE)
  noninv_median <- median(otu_data$plot_value[!otu_data$is_invasive], na.rm = TRUE)
  n_inv <- sum(otu_data$is_invasive)
  n_noninv <- sum(!otu_data$is_invasive)
  
  # Wilcoxon test for difference
  if (n_inv >= 3 && n_noninv >= 3) {
    wtest <- wilcox.test(plot_value ~ is_invasive, data = otu_data)
    pval_text <- paste0("Wilcoxon p = ", format.pval(wtest$p.value, digits = 3))
  } else {
    pval_text <- ""
  }
  
  # Create plot
  p <- ggplot(otu_data, aes(x = plot_value, y = invasive_label, color = invasive_label)) +
    geom_jitter(aes(size = n_sites), height = 0.2, alpha = 0.6) +
    geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.4) +
    scale_color_manual(values = c("Invasive" = "#e74c3c", "Non-invasive" = "#3498db")) +
    scale_size_continuous(range = c(1, 6), name = "N sites") +
    labs(x = m$lab,
         y = "",
         title = paste(m$title, title_suffix),
         subtitle = paste0("Invasive: n=", n_inv, ", median=", round(inv_median, 3),
                          " | Non-invasive: n=", n_noninv, ", median=", round(noninv_median, 3),
                          ifelse(pval_text != "", paste0("\n", pval_text), "")),
         color = "") +
    theme_minimal() +
    theme(legend.position = "right")
  
  if (metric == "reads" && log_reads) {
    p <- p + scale_x_log10()
  }
  
  p
}

#' Ranked dot plot showing all OTUs ordered by detection probability
#' Highlights invasive species in the ranking
#' 
#' @param combined Output from combine_detection_invasive()
#' @param country Country to plot (NULL for all data)
#' @param metric Which metric to rank by: "p_empirical", "theta", "confidence_score", "reads"
#' @param top_n Show top N OTUs (NULL for all)
#' @param label_invasive Add labels for invasive species
#' @return ggplot object
plot_detection_ranking <- function(combined, country = NULL, 
                                    metric = "p_empirical",
                                    top_n = NULL,
                                    label_invasive = TRUE) {
  
  # Filter by country if specified
  if (!is.null(country)) {
    site_data <- combined$site_level %>% filter(Country == !!country)
    title_suffix <- paste("-", country)
  } else {
    site_data <- combined$site_level
    title_suffix <- "- All Countries"
  }
  
  if (nrow(site_data) == 0) {
    message("No data available")
    return(NULL)
  }
  
  # Aggregate to OTU level
  otu_data <- site_data %>%
    group_by(OTU, Species) %>%
    summarise(
      n_sites = n(),
      mean_theta = mean(theta, na.rm = TRUE),
      mean_p_empirical = mean(p_empirical, na.rm = TRUE),
      mean_confidence_score = mean(confidence_score, na.rm = TRUE),
      total_reads = sum(total_reads, na.rm = TRUE),
      n_biosamples = sum(n_biosamples_positive, na.rm = TRUE),
      is_invasive = any(is_invasive, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      invasive_label = ifelse(is_invasive, "Invasive", "Non-invasive")
    )
  
  # Select metric column
  metric_col <- switch(metric,
                       "p_empirical" = "mean_p_empirical",
                       "theta" = "mean_theta",
                       "confidence_score" = "mean_confidence_score",
                       "reads" = "total_reads",
                       "n_biosamples" = "n_biosamples",
                       "mean_p_empirical")
  
  metric_lab <- switch(metric,
                       "p_empirical" = "Empirical detection rate",
                       "theta" = "θ (biosample detection)",
                       "confidence_score" = "Confidence score",
                       "reads" = "Total reads",
                       "n_biosamples" = "N biosamples positive",
                       "Empirical detection rate")
  
  otu_data$plot_value <- otu_data[[metric_col]]
  
  # Rank OTUs
  otu_data <- otu_data %>%
    arrange(desc(plot_value)) %>%
    mutate(rank = row_number())
  
  # Subset if requested
  if (!is.null(top_n)) {
    otu_data <- head(otu_data, top_n)
  }
  
  # Create plot
  p <- ggplot(otu_data, aes(x = rank, y = plot_value, color = invasive_label)) +
    geom_point(aes(size = n_sites), alpha = 0.7) +
    scale_color_manual(values = c("Invasive" = "#e74c3c", "Non-invasive" = "#3498db")) +
    scale_size_continuous(range = c(1, 5), name = "N sites") +
    labs(x = "OTU rank",
         y = metric_lab,
         title = paste("Detection Ranking", title_suffix),
         subtitle = paste("Red = Invasive species,", sum(otu_data$is_invasive), "of", nrow(otu_data), "OTUs"),
         color = "") +
    theme_minimal()
  
  # Add labels for invasive species
  if (label_invasive && requireNamespace("ggrepel", quietly = TRUE)) {
    inv_data <- otu_data %>% filter(is_invasive)
    if (nrow(inv_data) > 0) {
      p <- p + 
        ggrepel::geom_text_repel(
          data = inv_data,
          aes(label = Species),
          size = 2.5,
          max.overlaps = 20,
          na.rm = TRUE
        )
    }
  }
  
  p
}

#' Three-panel comparison: reads, theta, and p_empirical
#' 
#' @param combined Output from combine_detection_invasive()
#' @param country Country to plot (NULL for all data)
#' @return ggplot object (patchwork if available, otherwise list)
plot_detection_comparison <- function(combined, country = NULL) {
  
  # Filter by country if specified
  if (!is.null(country)) {
    site_data <- combined$site_level %>% filter(Country == !!country)
    title_suffix <- country
  } else {
    site_data <- combined$site_level
    title_suffix <- "All Countries"
  }
  
  if (nrow(site_data) == 0) {
    message("No data available")
    return(NULL)
  }
  
  # Aggregate to OTU level
  otu_data <- site_data %>%
    group_by(OTU, Species) %>%
    summarise(
      n_sites = n(),
      theta = mean(theta, na.rm = TRUE),
      p_empirical = mean(p_empirical, na.rm = TRUE),
      confidence_score = mean(confidence_score, na.rm = TRUE),
      total_reads = sum(total_reads, na.rm = TRUE),
      n_biosamples = sum(n_biosamples_positive, na.rm = TRUE),
      is_invasive = any(is_invasive, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      invasive_label = factor(ifelse(is_invasive, "Invasive", "Non-invasive"),
                              levels = c("Invasive", "Non-invasive")),
      log_reads = log10(total_reads + 1)
    )
  
  # Pivot to long format for faceting
  otu_long <- otu_data %>%
    select(OTU, Species, n_sites, is_invasive, invasive_label,
           theta, p_empirical, log_reads) %>%
    pivot_longer(cols = c(theta, p_empirical, log_reads),
                 names_to = "metric", values_to = "value") %>%
    mutate(
      metric = factor(metric, 
                      levels = c("log_reads", "theta", "p_empirical"),
                      labels = c("Log10(Reads)", "θ (biosample detection)", "Empirical detection rate"))
    )
  
  # Summary stats by metric
  summary_stats <- otu_long %>%
    group_by(metric, invasive_label) %>%
    summarise(
      median = median(value, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  
  # Create faceted plot
  p <- ggplot(otu_long, aes(x = value, y = invasive_label, color = invasive_label)) +
    geom_jitter(aes(size = n_sites), height = 0.2, alpha = 0.5) +
    geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.4) +
    facet_wrap(~metric, scales = "free_x", ncol = 3) +
    scale_color_manual(values = c("Invasive" = "#e74c3c", "Non-invasive" = "#3498db")) +
    scale_size_continuous(range = c(0.5, 4), name = "N sites") +
    labs(x = "Value",
         y = "",
         title = paste("Detection Metrics: Invasive vs Non-invasive -", title_suffix),
         subtitle = paste("n =", sum(otu_data$is_invasive), "invasive,", 
                          sum(!otu_data$is_invasive), "non-invasive OTUs"),
         color = "") +
    theme_minimal() +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
  
  p
}

#' Create summary table of invasive species
#' @param combined Output from combine_detection_invasive()
#' @return Data frame suitable for reporting
create_invasive_summary_table <- function(combined) {
  
  combined$invasive_detections %>%
    group_by(Species, OTU) %>%
    summarise(
      Countries = paste(unique(Country), collapse = ", "),
      n_countries = n_distinct(Country),
      n_sites = n_distinct(Site),  # Count unique sites
      sites = paste(unique(Site), collapse = "; "),
      mean_theta = round(mean(theta), 2),
      mean_p_pcr = round(mean(mean_p_pcr), 2),
      mean_p_detect = round(mean(p_detect_site), 2),
      n_reliable = sum(detection_category == "Reliable"),
      n_moderate = sum(detection_category == "Moderate"),
      n_unreliable = sum(detection_category == "Unreliable"),
      total_reads = sum(total_reads),
      EASIN_status = paste(unique(na.omit(EASIN_status)), collapse = "; "),
      GISD_status = paste(unique(na.omit(GISD_status)), collapse = "; "),
      .groups = "drop"
    ) %>%
    mutate(
      detection_summary = paste0(n_reliable, " reliable, ", 
                                  n_moderate, " moderate, ", 
                                  n_unreliable, " unreliable")
    ) %>%
    arrange(desc(n_sites), desc(n_countries))
}

# =============================================================================
# SECTION 8: THEORETICAL POWER ANALYSIS (HYPOTHETICAL PARAMETERS)
# NOTE: These functions take user-supplied p, theta, psi values and compute
# theoretical detection probabilities. They do NOT use empirical data.
# For empirical power analysis with uncertainty propagation, use
# calc_design_power_bootstrap() and the bootstrap output from
# run_detection_analysis() instead.
# =============================================================================

#' Calculate detection probability for nested design
#' @param p_pcr Detection probability per PCR
#' @param n_pcr Number of PCR replicates per sample
#' @param theta Detection probability per sample (DNA capture)
#' @param n_samples Number of biological samples per site
#' @return Cumulative detection probability
calc_cumulative_detection <- function(p_pcr, n_pcr = 3, theta = 1, n_samples = 3) {
  # P(detect in at least one PCR | DNA in sample)
  p_sample <- 1 - (1 - p_pcr)^n_pcr
  
  # P(detect in at least one sample | species present)
  # Incorporating sample-level detection (theta)
  p_site <- 1 - (1 - theta * p_sample)^n_samples
  
  return(p_site)
}

#' Power analysis: probability of detecting species at occupied sites
#' @param psi True occupancy
#' @param p_pcr PCR detection probability
#' @param theta Sample detection probability (default 1 = perfect DNA capture)
#' @param n_samples Number of biological samples
#' @param n_pcr Number of PCR replicates
#' @param n_sites Number of sites surveyed
#' @return Data frame with power metrics
power_analysis <- function(psi = 0.3, 
                           p_pcr = 0.7, 
                           theta = 0.9,
                           n_samples = 3, 
                           n_pcr = 3, 
                           n_sites = 30) {
  
  # Cumulative detection at site level
  p_site <- calc_cumulative_detection(p_pcr, n_pcr, theta, n_samples)
  
  # Expected number of sites where species is detected
  expected_sites_occupied <- n_sites * psi
  expected_detections <- expected_sites_occupied * p_site
  
  # Probability of detecting species at least once (across all sites)
  p_detect_somewhere <- 1 - (1 - psi * p_site)^n_sites
  
  # Power to estimate occupancy with reasonable precision
  # Rule of thumb: need ~10 detections for stable estimates
  p_10_detections <- 1 - pbinom(9, n_sites, psi * p_site)
  
  data.frame(
    true_occupancy = psi,
    p_pcr = p_pcr,
    theta = theta,
    n_samples = n_samples,
    n_pcr = n_pcr,
    n_sites = n_sites,
    p_site = p_site,
    expected_occupied = expected_sites_occupied,
    expected_detections = expected_detections,
    p_detect_anywhere = p_detect_somewhere,
    power_10_detections = p_10_detections
  )
}

#' Compare different sampling designs
#' @param psi_range Range of occupancy values to test
#' @param p_pcr_range Range of PCR detection probabilities
#' @param designs List of sampling designs (n_samples × n_pcr combinations)
#' @return Data frame comparing designs
compare_designs <- function(psi_range = seq(0.1, 0.5, 0.1),
                            p_pcr = 0.7,
                            theta = 0.9,
                            n_sites = 30,
                            designs = list(
                              "3×3 (current)" = c(3, 3),
                              "4×2" = c(4, 2),
                              "2×4" = c(2, 4),
                              "5×2" = c(5, 2),
                              "3×4" = c(3, 4)
                            )) {
  
  results <- list()
  
  for (psi in psi_range) {
    for (design_name in names(designs)) {
      d <- designs[[design_name]]
      res <- power_analysis(
        psi = psi,
        p_pcr = p_pcr,
        theta = theta,
        n_samples = d[1],
        n_pcr = d[2],
        n_sites = n_sites
      )
      res$design <- design_name
      res$total_effort <- d[1] * d[2]
      results[[length(results) + 1]] <- res
    }
  }
  
  bind_rows(results)
}

#' Minimum sample size calculation
#' @param psi Target occupancy to detect
#' @param p_site Site-level detection probability
#' @param power Desired power (probability of getting 10+ detections)
#' @param min_detections Minimum detections needed
#' @return Required number of sites
calc_min_sites <- function(psi = 0.3, p_site = 0.9, power = 0.8, min_detections = 10) {
  
  # Search for n_sites that achieves desired power
  for (n in 10:500) {
    p_success <- 1 - pbinom(min_detections - 1, n, psi * p_site)
    if (p_success >= power) {
      return(n)
    }
  }
  return(NA)
}

# =============================================================================
# SECTION 9: THEORETICAL POWER VISUALIZATION FUNCTIONS
# (Companion plots for Section 8)
# =============================================================================

#' Plot detection probability vs sampling effort
plot_detection_effort <- function(p_pcr_range = seq(0.3, 0.9, 0.1),
                                   designs = list(
                                     "3×3" = c(3, 3),
                                     "4×2" = c(4, 2),
                                     "2×4" = c(2, 4),
                                     "5×2" = c(5, 2)
                                   )) {
  
  results <- expand.grid(
    p_pcr = p_pcr_range,
    design = names(designs)
  )
  
  results$p_site <- mapply(function(p, d) {
    calc_cumulative_detection(p, designs[[d]][2], theta = 0.9, designs[[d]][1])
  }, results$p_pcr, results$design)
  
  ggplot(results, aes(x = p_pcr, y = p_site, color = design)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2) +
    labs(
      x = "PCR Detection Probability (p)",
      y = "Site-Level Detection Probability",
      title = "Cumulative Detection Probability by Sampling Design",
      subtitle = "Assuming θ = 0.9 (sample-level detection)",
      color = "Design\n(samples × PCR)"
    ) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    scale_x_continuous(labels = scales::percent) +
    theme_minimal() +
    theme(legend.position = "right")
}

#' Plot power curves for different occupancy levels
plot_power_curves <- function(n_sites_range = seq(10, 100, 5),
                               psi_values = c(0.1, 0.2, 0.3, 0.5),
                               p_pcr = 0.7,
                               n_samples = 3,
                               n_pcr = 3) {
  
  results <- expand.grid(
    n_sites = n_sites_range,
    psi = psi_values
  )
  
  p_site <- calc_cumulative_detection(p_pcr, n_pcr, theta = 0.9, n_samples)
  
  results$power <- mapply(function(n, psi) {
    1 - pbinom(9, n, psi * p_site)
  }, results$n_sites, results$psi)
  
  results$psi_label <- paste0("ψ = ", results$psi)
  
  ggplot(results, aes(x = n_sites, y = power, color = psi_label)) +
    geom_line(linewidth = 1.2) +
    geom_hline(yintercept = 0.8, linetype = "dashed", alpha = 0.5) +
    annotate("text", x = 90, y = 0.82, label = "80% power", size = 3) +
    labs(
      x = "Number of Sites",
      y = "Power (P ≥ 10 detections)",
      title = "Statistical Power by Number of Sites",
      subtitle = sprintf("Design: %d samples × %d PCR, p = %.0f%%", n_samples, n_pcr, p_pcr*100),
      color = "True Occupancy"
    ) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    theme_minimal()
}

#' Plot naive vs corrected occupancy
plot_occupancy_comparison <- function(results_df) {
  
  df <- results_df %>%
    filter(!is.na(psi_estimate)) %>%
    mutate(
      Species_short = substr(Species, 1, 25)
    )
  
  ggplot(df, aes(x = naive_occupancy, y = psi_estimate)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_errorbar(aes(ymin = psi_lower, ymax = psi_upper), 
                  alpha = 0.3, width = 0.01) +
    geom_point(aes(size = n_detections, color = p_estimate)) +
    scale_color_viridis_c(name = "Detection\nProbability (p)", 
                          limits = c(0, 1), labels = scales::percent) +
    scale_size_continuous(name = "N Detections", range = c(1, 5)) +
    labs(
      x = "Naive Occupancy (uncorrected)",
      y = "Estimated Occupancy (ψ)",
      title = "Naive vs. Detection-Corrected Occupancy",
      subtitle = "Points above line = occupancy underestimated by naive approach"
    ) +
    theme_minimal() +
    coord_equal()
}

# =============================================================================
# SECTION 10: BOOTSTRAP POWER ANALYSIS FUNCTIONS
# =============================================================================

#' Calculate detection probability with uncertainty in θ and p
#'
#' Instead of treating observed θ and p as point estimates, samples from
#' their Beta posterior distributions given the observed data, then computes
#' the theoretical detection probability analytically for each posterior draw.
#'
#' @param n_bio_pos Number of biosamples where OTU was detected at this site
#' @param n_bio_total Total number of biosamples at this site
#' @param n_pcr_pos Number of positive PCRs across positive biosamples at this site
#' @param n_pcr_total Total PCRs across positive biosamples at this site
#' @param n_bio_design Number of biosamples in the design to evaluate
#' @param n_pcr_design Number of PCRs per biosample in the design to evaluate
#' @param n_sims Number of posterior draws
#' @param prior_weight Weight for Beta prior (1 = uniform Beta(1,1))
#' @return Named vector: mean, median, lower (2.5%), upper (97.5%), sd
simulate_detection_with_uncertainty <- function(
    n_bio_pos, n_bio_total,
    n_pcr_pos, n_pcr_total,
    n_bio_design, n_pcr_design,
    n_sims = 5000,
    prior_weight = 1
) {
  
  # Sample θ from posterior: Beta(successes + prior, failures + prior)
  theta_samples <- rbeta(n_sims,
                         n_bio_pos + prior_weight,
                         n_bio_total - n_bio_pos + prior_weight)
  
  # Sample p from posterior
  p_samples <- rbeta(n_sims,
                     n_pcr_pos + prior_weight,
                     n_pcr_total - n_pcr_pos + prior_weight)
  
  # Vectorised analytical detection probability for each posterior draw
  # P(>=1 positive PCR | DNA in sample) = 1 - (1-p)^n_pcr
  # P(biosample contributes detection) = θ * [1 - (1-p)^n_pcr]
  # P(no detection in biosample) = 1 - θ*[1-(1-p)^n_pcr]
  # P(detect at site) = 1 - P(no detection)^n_bio
  p_pcr_success <- 1 - (1 - p_samples)^n_pcr_design
  p_biosample_positive <- theta_samples * p_pcr_success
  p_detect_samples <- 1 - (1 - p_biosample_positive)^n_bio_design
  
  c(
    mean   = mean(p_detect_samples),
    median = median(p_detect_samples),
    lower  = unname(quantile(p_detect_samples, 0.025)),
    upper  = unname(quantile(p_detect_samples, 0.975)),
    sd     = sd(p_detect_samples)
  )
}


#' Bootstrap design power analysis from site-level detection data
#'
#' For each unique detection pattern (combination of observed biosample/PCR
#' counts), simulates detection probability under Beta posterior uncertainty
#' across a grid of sampling designs. Results are joined back to the full
#' site-level detection data.
#'
#' @param site_detection Site-level detection data frame (from add_biosample_totals).
#'   Must contain: n_biosamples_positive, total_biosamples_at_site,
#'   n_pcr_positive, n_pcr_total (PCRs in positive biosamples only).
#' @param designs Data frame with columns n_biosamples, n_pcr, design, effort.
#' @param n_sims Number of posterior draws per pattern per design.
#' @return The site-level data joined with design power columns.
calc_design_power_bootstrap <- function(site_detection, designs, n_sims = 5000) {
  
  # Identify unique observed count patterns (avoids redundant simulation)
  unique_patterns <- site_detection %>%
    distinct(n_biosamples_positive, total_biosamples_at_site,
             n_pcr_positive, n_pcr_total)
  
  cat("    Unique detection patterns:", nrow(unique_patterns), "\n")
  cat("    Designs to evaluate:", nrow(designs), "\n")
  cat("    Running posterior simulations...\n")
  
  pb <- txtProgressBar(min = 0, max = nrow(unique_patterns), style = 3)
  
  pattern_results <- list()
  
  for (i in seq_len(nrow(unique_patterns))) {
    pat <- unique_patterns[i, ]
    
    design_results <- lapply(seq_len(nrow(designs)), function(d) {
      result <- simulate_detection_with_uncertainty(
        n_bio_pos   = pat$n_biosamples_positive,
        n_bio_total = pat$total_biosamples_at_site,
        n_pcr_pos   = pat$n_pcr_positive,
        n_pcr_total = pat$n_pcr_total,
        n_bio_design = designs$n_biosamples[d],
        n_pcr_design = designs$n_pcr[d],
        n_sims = n_sims
      )
      
      tibble(
        design        = designs$design[d],
        n_bio_design  = designs$n_biosamples[d],
        n_pcr_design  = designs$n_pcr[d],
        effort        = designs$effort[d],
        p_detect_mean  = result["mean"],
        p_detect_median = result["median"],
        p_detect_lower = result["lower"],
        p_detect_upper = result["upper"],
        p_detect_sd    = result["sd"]
      )
    }) %>% bind_rows()
    
    pattern_results[[i]] <- design_results %>%
      mutate(
        n_biosamples_positive     = pat$n_biosamples_positive,
        total_biosamples_at_site  = pat$total_biosamples_at_site,
        n_pcr_positive            = pat$n_pcr_positive,
        n_pcr_total               = pat$n_pcr_total
      )
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  # Combine and join back to site-level data
  power_lookup <- bind_rows(pattern_results)
  
  site_detection %>%
    left_join(
      power_lookup,
      by = c("n_biosamples_positive", "total_biosamples_at_site",
             "n_pcr_positive", "n_pcr_total"),
      relationship = "many-to-many"
    )
}


# =============================================================================
# SECTION 11: MAIN WORKFLOW
# =============================================================================

#' Run complete hierarchical detection analysis
#' 
#' Computes site-level detection probabilities for every OTU at every site,
#' then evaluates sampling design power using Beta-posterior bootstrapping.
#' All detection parameters (θ, p) are site-specific — no cross-site pooling.
#'
#' IMPORTANT: Subset your phyloseq to relevant samples BEFORE running:
#'   ps_cs <- subset_samples(ps, Sampling.area.Project == "BGE Marine Invasive Species Citizen Science")
#'   det_12S <- run_detection_analysis(ps_cs, "12S")
#'
#' @param ps Phyloseq object (should be pre-subsetted to relevant samples)
#' @param marker Marker name for output files
#' @param output_dir Output directory
#' @param site_var Column for site
#' @param biosample_var Column for biological sample
#' @param pcr_var Column for PCR replicate
#' @param country_var Column for country
#' @param n_sims Number of posterior draws for bootstrap power analysis
#' @return List with site-level detection results and design power analysis
run_detection_analysis <- function(ps, 
                                    marker = "12S", 
                                    output_dir = "Processed_data",
                                    site_var = "Sampling.area.Name",
                                    biosample_var = "Name",
                                    pcr_var = "Replicate",
                                    country_var = "Location",
                                    n_sims = 5000) {
  
  cat("\n", strrep("=", 70), "\n")
  cat("HIERARCHICAL DETECTION ANALYSIS FOR", marker, "\n")
  cat("(All metrics are site-specific — no cross-site pooling)\n")
  cat(strrep("=", 70), "\n\n")
  
  # ---- 1. Get sample hierarchy ----
  cat("STEP 1: Parsing sample hierarchy...\n")
  cat("  Site variable:", site_var, "\n")
  cat("  Biological sample variable:", biosample_var, "\n")
  cat("  PCR replicate variable:", pcr_var, "\n\n")
  
  hierarchy <- get_sample_hierarchy(ps, site_var, biosample_var, pcr_var, country_var)
  design <- summarize_design(hierarchy)
  
  cat("  Design summary:\n")
  cat("    Countries:", design$n_countries, "\n")
  cat("    Sites:", design$n_sites, "\n")
  cat("    Biological samples:", design$n_biosamples, "\n")
  cat("    Total PCR replicates:", design$n_pcr_total, "\n")
  cat("    Median biosamples/site:", design$median_biosamples_per_site, "\n")
  cat("    Median PCRs/biosample:", design$median_pcr_per_biosample, "\n\n")
  
  # Check for expected 3x3 design
  if (design$median_biosamples_per_site != 3 || design$median_pcr_per_biosample != 3) {
    warning("Design is not 3x3 - did you subset to citizen science samples?")
  }
  
  # ---- 2. Calculate PCR-level detection (OTU level) ----
  cat("STEP 2: Calculating PCR-level detection per biological sample...\n")
  pcr_detection <- calc_hierarchical_detection(ps, hierarchy)
  cat("  OTU × biosample detections:", nrow(pcr_detection), "\n")
  cat("  Unique OTUs:", n_distinct(pcr_detection$OTU), "\n\n")
  
  # ---- 3. Aggregate to site level (OTU level) ----
  cat("STEP 3: Aggregating to site level...\n")
  site_detection <- calc_site_level_detection(pcr_detection)
  site_detection <- add_biosample_totals(site_detection, hierarchy)
  cat("  OTU × site combinations:", nrow(site_detection), "\n\n")
  
  # ---- 4. Join taxonomy ----
  cat("STEP 4: Joining taxonomy to results...\n")
  site_detection_tax <- join_taxonomy(site_detection, ps)
  cat("  Taxonomy columns added\n\n")
  
  # ---- 5. Site-level detection parameter summary ----
  cat("STEP 5: Site-level detection parameter distributions...\n")
  cat("  (Each value is from a single OTU × site combination — no cross-site pooling)\n\n")
  
  cat("  θ (proportion of biosamples detecting OTU, per site):\n")
  cat("    Median:", round(median(site_detection$theta), 3), "\n")
  cat("    IQR:   ", round(quantile(site_detection$theta, 0.25), 3), "–",
                     round(quantile(site_detection$theta, 0.75), 3), "\n")
  cat("    Range: ", round(min(site_detection$theta), 3), "–",
                     round(max(site_detection$theta), 3), "\n")
  cat("    % with θ = 1 (all biosamples positive):",
      round(mean(site_detection$theta == 1) * 100, 1), "%\n\n")
  
  cat("  p (PCR detection rate within positive biosamples, per site):\n")
  cat("    Median:", round(median(site_detection$mean_p_pcr), 3), "\n")
  cat("    IQR:   ", round(quantile(site_detection$mean_p_pcr, 0.25), 3), "–",
                     round(quantile(site_detection$mean_p_pcr, 0.75), 3), "\n")
  cat("    Range: ", round(min(site_detection$mean_p_pcr), 3), "–",
                     round(max(site_detection$mean_p_pcr), 3), "\n")
  cat("    % with p = 1 (all PCRs positive):",
      round(mean(site_detection$mean_p_pcr == 1) * 100, 1), "%\n\n")
  
  cat("  p_detect_site (empirical site-level detection, 3×3 design):\n")
  cat("    Median:", round(median(site_detection$p_detect_site), 3), "\n")
  cat("    IQR:   ", round(quantile(site_detection$p_detect_site, 0.25), 3), "–",
                     round(quantile(site_detection$p_detect_site, 0.75), 3), "\n")
  cat("    % with P(detect) ≥ 0.80:",
      round(mean(site_detection$p_detect_site >= 0.80) * 100, 1), "%\n")
  cat("    % with P(detect) < 0.50:",
      round(mean(site_detection$p_detect_site < 0.50) * 100, 1), "%\n\n")
  
  # Show top OTUs by detection breadth
  cat("  Top 10 OTUs by number of sites detected:\n")
  top_otus <- site_detection_tax %>%
    group_by(OTU, Species) %>%
    summarise(
      n_sites = n(),
      median_theta = round(median(theta), 3),
      median_p_pcr = round(median(mean_p_pcr), 3),
      median_p_detect = round(median(p_detect_site), 3),
      .groups = "drop"
    ) %>%
    arrange(desc(n_sites)) %>%
    head(10)
  print(top_otus)
  cat("\n")
  
  # ---- 6. Bootstrap design power analysis ----
  cat("STEP 6: Bootstrap design power analysis (Beta-posterior uncertainty)...\n\n")
  
  # Design grid
  design_grid <- expand.grid(
    n_biosamples = 2:8,
    n_pcr = 2:6
  ) %>%
    mutate(
      design = paste0(n_biosamples, "\u00d7", n_pcr),
      effort = n_biosamples * n_pcr
    )
  
  # Run bootstrap across all unique detection patterns × designs
  empirical_with_power <- calc_design_power_bootstrap(
    site_detection, design_grid, n_sims = n_sims
  )
  
  # Summarise: % of OTU × site detections achieving ≥80% for each design
  design_performance <- empirical_with_power %>%
    group_by(design, n_bio_design, n_pcr_design, effort) %>%
    summarise(
      n = n(),
      pct_above_80_mean         = round(mean(p_detect_mean >= 0.80) * 100, 1),
      pct_above_80_conservative = round(mean(p_detect_lower >= 0.80) * 100, 1),
      pct_above_80_optimistic   = round(mean(p_detect_upper >= 0.80) * 100, 1),
      pct_above_50_mean         = round(mean(p_detect_mean >= 0.50) * 100, 1),
      pct_above_95_mean         = round(mean(p_detect_mean >= 0.95) * 100, 1),
      mean_p_detect             = round(mean(p_detect_mean), 3),
      mean_CI_width             = round(mean(p_detect_upper - p_detect_lower), 3),
      .groups = "drop"
    ) %>%
    arrange(effort, desc(pct_above_80_mean))
  
  # Print key designs
  cat("\n  Design comparison (% of OTU × site detections with P(detect) ≥ 80%):\n")
  cat("  'Mean' = using posterior mean; 'Conservative' = lower 95% CI ≥ 80%\n\n")
  key_designs <- design_performance %>%
    filter(design %in% c("3\u00d73", "4\u00d73", "3\u00d74", "5\u00d73",
                          "6\u00d72", "6\u00d73", "8\u00d73")) %>%
    select(design, effort, pct_above_80_mean, pct_above_80_conservative,
           mean_p_detect, mean_CI_width) %>%
    arrange(effort)
  print(key_designs)
  cat("\n")
  
  # Extract 3×3 results for convenience
  empirical_3x3 <- empirical_with_power %>%
    filter(design == "3\u00d73") %>%
    select(-design, -n_bio_design, -n_pcr_design, -effort)
  
  # ---- 7. Save outputs ----
  cat("STEP 7: Saving outputs...\n")
  
  f1 <- file.path(output_dir, paste0("pcr_detection_", marker, ".csv"))
  write.csv(pcr_detection, f1, row.names = FALSE)
  cat("  Saved:", f1, "\n")
  
  f2 <- file.path(output_dir, paste0("site_detection_", marker, ".csv"))
  write.csv(site_detection_tax, f2, row.names = FALSE)
  cat("  Saved:", f2, "\n")
  
  f3 <- file.path(output_dir, paste0("site_design_", marker, ".csv"))
  write.csv(design$site_details, f3, row.names = FALSE)
  cat("  Saved:", f3, "\n")
  
  f4 <- file.path(output_dir, paste0("design_power_", marker, ".csv"))
  write.csv(design_performance, f4, row.names = FALSE)
  cat("  Saved:", f4, "\n")
  
  cat("\n", strrep("=", 70), "\n")
  cat("ANALYSIS COMPLETE\n")
  cat(strrep("=", 70), "\n")
  
  return(list(
    marker = marker,
    hierarchy = hierarchy,
    design = design,
    # Site-level results (one row per OTU × site)
    pcr_detection = pcr_detection,
    site_detection = site_detection,
    site_detection_tax = site_detection_tax,
    # Bootstrap power analysis
    empirical_3x3 = empirical_3x3,
    empirical_with_power = empirical_with_power,
    design_performance = design_performance
  ))
}

#' Generate detection report for a specific site
#' @param results Output from run_detection_analysis
#' @param site Site name
#' @param country Country name (optional, for filtering)
#' @return Data frame with OTU list and hierarchical detection info
get_site_report <- function(results, site, country = NULL) {
  
  site_data <- results$site_detection_tax %>%
    filter(Site == !!site)
  
  if (!is.null(country)) {
    site_data <- site_data %>% filter(Country == !!country)
  }
  
  if (nrow(site_data) == 0) {
    message("No detections found for site: ", site)
    return(NULL)
  }
  
  site_data %>%
    select(OTU, Species, n_biosamples_positive, total_biosamples_at_site, theta,
           n_pcr_positive, n_pcr_total, mean_p_pcr, confidence, total_reads) %>%
    arrange(desc(theta), desc(mean_p_pcr))
}

# =============================================================================
# USAGE
# =============================================================================

cat("\n", strrep("=", 70), "\n")
cat("HIERARCHICAL DETECTION ANALYSIS SCRIPT LOADED\n")
cat(strrep("=", 70), "\n\n")

cat("=== COMPLETE WORKFLOW ===\n\n")

cat("# 1. Load and process data (already done)\n")
cat("source('Scripts/FUNC_process_metabarcoding_data.R')\n")
cat("results_12S <- process_marker('12S', ps_input = bge12s_cs)\n\n")

cat("# 2. Verify invasive status via databases (already done)\n")
cat("source('Scripts/FUNC_query_invasive_status_databases.R')\n")
cat("verified_12S <- process_species_list(results_12S$invasive_status, ...)\n\n")

cat("# 3. Run detection analysis on CITIZEN SCIENCE subset\n")
cat("source('Scripts/occupancy_power_analysis.R')\n")
cat("ps_cs <- subset_samples(results_12S$phyloseq,\n")
cat("           Sampling.area.Project == 'BGE Marine Invasive Species Citizen Science')\n")
cat("det_12S <- run_detection_analysis(ps_cs, '12S')\n\n")

cat("# 4. Combine detection + invasive status\n")
cat("combined_12S <- combine_detection_invasive(\n")
cat("  detection_results = det_12S,\n")
cat("  invasive_status = results_12S$invasive_status,\n")
cat("  verified_status = verified_12S\n")
cat(")\n\n")

cat("# 5. Explore the data\n")
cat("list_countries(combined_12S)              # See available countries\n")
cat("list_sites(combined_12S, 'Norway')        # See sites in Norway\n\n")

cat("# 6. Create COUNTRY-SPECIFIC visualizations\n")
cat("plot_p_vs_theta_country(combined_12S, 'Greece')    # p vs theta for Greece\n")
cat("plot_invasive_country(combined_12S, 'Spain')       # Invasives in Spain\n")
cat("plot_site_detection(combined_12S, 'Port of Oslo')  # Single site summary\n\n")

cat("# 7. Create SUMMARY visualizations\n")
cat("plot_country_summary(combined_12S)        # Detection reliability by country\n")
cat("plot_all_invasives(combined_12S)          # Top invasive species overall\n")
cat("plot_invasive_heatmap(combined_12S)       # Species × country heatmap\n\n")

cat("# 8. SITE-LEVEL PLOTS (the only level that makes sense)\n")
cat("# List sites in a country\n")
cat("list_sites(combined_12S_filtered, 'Norway')\n\n")

cat("# Single site - reads (x) vs detection metric (y)\n")
cat("plot_site_detection(combined_12S_filtered, 'Site Name Here')                 # p_empirical\n")
cat("plot_site_detection(combined_12S_filtered, 'Site Name Here', 'theta')        # theta\n")
cat("plot_site_detection(combined_12S_filtered, 'Site Name Here', 'confidence_score')\n\n")

cat("# Single site - all three metrics in facets\n")
cat("plot_site_all_metrics(combined_12S_filtered, 'Site Name Here')\n\n")

cat("# Full data table for a site\n")
cat("site_summary <- get_site_summary(combined_12S_filtered, 'Site Name Here')\n\n")

cat("# Loop through all sites in a country\n")
cat("plot_all_sites_in_country(combined_12S_filtered, 'Norway')\n")
cat("plot_all_sites_in_country(combined_12S_filtered, 'Greece', output_dir = 'Figures/')\n\n")

cat("# Compare invasive species across all sites\n")
cat("plot_invasive_across_sites(combined_12S_filtered, 'p_empirical')\n\n")

cat("# 9. Loop through all countries\n")
cat("for (country in list_countries(combined_12S)) {\n")
cat("  print(plot_p_vs_theta_country(combined_12S, country))\n")
cat("}\n\n")

cat("# 10. Create summary table for reporting\n")
cat("invasive_table <- create_invasive_summary_table(combined_12S)\n")
cat(strrep("=", 70), "\n")
