# =============================================================================
# METABARCODING METADATA PROCESSING (Steps 1–2)
# Loading, cleaning, and summarising phyloseq objects
# =============================================================================
#
# Author: BGE Marine Invasive Species Project
# Date: 2026-02-26
#
# Functions:
#   load_phyloseq()           - Load & subset RDS to project
#   add_ntnu_metadata()       - Join NTNU sample metadata
#   fix_citsci_metadata()     - Apply all citizen-science metadata fixes (inlined)
#   summarize_otus_by_event() - OTU summary table per sampling event
#   process_metadata()        - Run Steps 1–2 for one marker
#   process_all_metadata()    - Run Steps 1–2 for all markers
#
# Usage:
#   source("Scripts/Final/FUNC_process_metadata.R")
#   meta_12S <- process_metadata("12S", ps_input = bge12s_cs)
# =============================================================================

# ---- Load Required Packages ----
library(phyloseq)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(readr)

# ---- Shared Configuration (set once) ----
if (!exists("CONFIG")) {
  CONFIG <- list(
    raw_data_dir       = "Raw_data",
    processed_data_dir = "Processed_data",
    scripts_dir        = "Scripts/Final",
    ntnu_metadata_file = "BGE_MNIS_CScodes.csv",
    fixed_metadata     = "fixed_sample_metadata.csv",
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
# SECTION 1: DATA LOADING AND PREPROCESSING
# =============================================================================

#' Load and subset phyloseq object for a specific marker
#' @param marker Character: "12S", "18S", or "COI"
#' @param rds_file Optional: full path to RDS file. If NULL, uses default naming
#' @return Subsetted phyloseq object
load_phyloseq <- function(marker, rds_file = NULL) {

  if (is.null(rds_file)) {
    file_patterns <- list(
      "12S" = "12S_OTUs_NCBI.RDS",
      "18S" = "18S_OTUs_syntax.RDS",
      "COI" = "COI_OTUs_NCBI.RDS"
    )
    rds_file <- file.path(CONFIG$raw_data_dir, file_patterns[[marker]])
  }

  dbg("Loading", marker, "data from:", rds_file)

  if (!file.exists(rds_file)) {
    stop("File not found: ", rds_file)
  }

  ps <- readRDS(rds_file)
  dbg("  Loaded phyloseq with", nsamples(ps), "samples and", ntaxa(ps), "taxa")

  dbg("  Subsetting to:", CONFIG$project_name)
  ps_subset <- subset_samples(ps, Sampling.area.Project == CONFIG$project_name)
  ps_subset <- prune_taxa(taxa_sums(ps_subset) > 0, ps_subset)

  dbg("  After subsetting:", nsamples(ps_subset), "samples,", ntaxa(ps_subset), "taxa")

  return(ps_subset)
}

#' Add NTNU metadata to phyloseq object
#' @param ps Phyloseq object
#' @return Phyloseq object with updated sample_data
add_ntnu_metadata <- function(ps) {

  ntnu_file <- file.path(CONFIG$processed_data_dir, CONFIG$ntnu_metadata_file)

  if (!file.exists(ntnu_file)) {
    warning("NTNU metadata file not found: ", ntnu_file)
    return(ps)
  }

  dbg("Adding NTNU metadata from:", ntnu_file)
  ntnudat <- read.csv(ntnu_file, header = TRUE)

  sample_data(ps) <- sample_data(ps) %>%
    as("data.frame") %>%
    left_join(ntnudat, by = c("Name" = "Codes1")) %>%
    tibble::column_to_rownames("Library_ID") %>%
    sample_data()

  return(ps)
}

#' Fix BGE Citizen Science Phyloseq Metadata
#'
#' Applies all metadata corrections to a BGE Marine Invasive Species
#' Citizen Science phyloseq object, regardless of marker (COI/12S/18S) or
#' number of PCR replicates. Corrections target the "Name" field
#' (bare BGEMIS code), so they are agnostic to fragment identifiers
#' and replicate suffixes in rownames.
#'
#' Corrections applied:
#'   1. Whitespace trimming in Site.name, coordinate rounding (2 dp),
#'      preserve full-precision coordinates
#'   2. Copenhagen: disambiguate two sampling sites
#'   3. Portrush (UK): disambiguate two sampling events
#'   4. Heraklion (Greece): disambiguate two Dermatas bay events
#'   5. Calabria: standardize Sampling.area.Name across replicates
#'   6. Portugal (Vila do Conde): fill missing Sampling.area.Name and
#'      Sampling.area.Parent.project
#'   7. Generate unique Sampling.event.ID per sampling event
#'   8. Replace sample-level metadata with curated spreadsheet
#'      (fixed_sample_metadata.csv), matched by Name (BGEMIS code).
#'
#' @param ps A phyloseq object (subset to Citizen Science samples).
#' @param fixed_metadata Path to curated metadata CSV. Default uses CONFIG.
#' @return The same phyloseq object with corrected sample_data.
fix_citsci_metadata <- function(ps,
                                fixed_metadata = file.path(CONFIG$processed_data_dir,
                                                           CONFIG$fixed_metadata)) {

  # Extract sample_data, keeping phyloseq rownames
  sd_df <- as(sample_data(ps), "data.frame")
  sd_df$Library_ID <- rownames(sd_df)

  # --- 1. General cleanup ---
  sd_df <- sd_df %>%
    mutate(
      Site.name               = str_squish(Site.name),
      latitude_full           = Sampling.area.Latitude,
      longitude_full          = Sampling.area.Longitude,
      Sampling.area.Latitude  = round(Sampling.area.Latitude, 2),
      Sampling.area.Longitude = round(Sampling.area.Longitude, 2)
    )

  # --- 2. Copenhagen: disambiguate two sampling sites ---
  sd_df <- sd_df %>%
    mutate(Sampling.area.Name = case_when(
      Name %in% c("BGEMIS0661", "BGEMIS0662", "BGEMIS0663") ~
        "Denmark, Copenhagen, Copenhagen, Floating jetty for harbor bus",
      Name %in% c("BGEMIS0664", "BGEMIS0665") ~
        "Denmark, Copenhagen, Copenhagen, Nyhavn Harbor",
      TRUE ~ Sampling.area.Name
    ))

  # --- 3. Portrush: disambiguate two sampling events ---
  sd_df <- sd_df %>%
    mutate(Sampling.area.Name = case_when(
      Name %in% c("BGEMIS0643", "BGEMIS0644", "BGEMIS0645") ~
        "United Kingdom, North Coast, Portrush1",
      TRUE ~ Sampling.area.Name
    ))

  # --- 4. Heraklion: disambiguate Dermatas bay events ---
  sd_df <- sd_df %>%
    mutate(
      Sampling.area.Name = case_when(
        Name %in% c("BGEMIS0713", "BGEMIS0714", "BGEMIS0715") ~
          "Greece, Heraklion, Dermatas bay1",
        TRUE ~ Sampling.area.Name
      ),
      Sampling.event.Event.description = case_when(
        Name %in% c("BGEMIS0713", "BGEMIS0714", "BGEMIS0715") ~
          "Dermatas bay1",
        TRUE ~ Sampling.event.Event.description
      )
    )

  # --- 5. Calabria: standardize Sampling.area.Name ---
  sd_df <- sd_df %>%
    mutate(Sampling.area.Name = case_when(
      Name %in% c("BGEMIS0655", "BGEMIS0656", "BGEMIS0657") ~
        "Calabria, Province of Vibo Valentia, Vibo Valentia, Vibo Marina",
      TRUE ~ Sampling.area.Name
    ))

  # --- 6. Portugal (Vila do Conde): fill missing metadata ---
  sd_df <- sd_df %>%
    mutate(
      Sampling.area.Name = case_when(
        Name %in% c("BGEMIS0670", "BGEMIS0671", "BGEMIS0672") ~
          "Póvoa De Varzim",
        Name == "BGEMIS0673" ~
          "Vila Praia de Ancora",
        TRUE ~ Sampling.area.Name
      ),
      Sampling.area.Parent.project = case_when(
        Name == "BGEMIS0673" & is.na(Sampling.area.Parent.project) ~
          "BGE eDNA and Marine Invasive Species",
        TRUE ~ Sampling.area.Parent.project
      )
    )

  # --- 7. Build unique Sampling.event.ID ---
  sd_df <- sd_df %>%
    group_by(Location, Site.name, Sampling.area.Name,
             Sampling.area.Latitude, Sampling.area.Longitude) %>%
    mutate(Event_within_country = cur_group_id()) %>%
    ungroup() %>%
    group_by(Location) %>%
    mutate(
      Sampling.event.ID = paste0(Location, "_", dense_rank(Event_within_country))
    ) %>%
    ungroup()

  # --- 8. Replace sample-level metadata from curated spreadsheet ---
  if (!is.null(fixed_metadata) && file.exists(fixed_metadata)) {
    dbg("Applying curated metadata from:", fixed_metadata)

    fix_df <- read.csv(fixed_metadata, header = TRUE, check.names = FALSE) %>%
      select(-any_of(c("Sample_ID", "Extraction_Batch", "DNA_Plate",
                        "Index_Plate", "Sample_Well", "Fragment", "Replicate"))) %>%
      distinct(Name, .keep_all = TRUE)

    n_match <- sum(sd_df$Name %in% fix_df$Name)
    dbg("  Matched", n_match, "/", nrow(sd_df),
        "samples (", n_distinct(fix_df$Name), "unique Names in spreadsheet)")

    # Overwrite matching columns for matching Names
    matched_rows <- which(sd_df$Name %in% fix_df$Name)
    common_cols <- setdiff(intersect(names(sd_df), names(fix_df)), "Name")
    fix_lookup <- fix_df %>% tibble::column_to_rownames("Name")

    for (col in common_cols) {
      sd_df[[col]][matched_rows] <- fix_lookup[sd_df$Name[matched_rows], col]
    }

    # Add any new columns only present in the spreadsheet
    new_cols <- setdiff(names(fix_df), names(sd_df))
    if (length(new_cols) > 0) {
      sd_df <- sd_df %>%
        left_join(fix_df %>% select(Name, all_of(new_cols)), by = "Name")
      dbg("  Added new columns:", paste(new_cols, collapse = ", "))
    }

  } else if (!is.null(fixed_metadata)) {
    warning("Curated metadata file not found: ", fixed_metadata,
            " \u2014 skipping fix 8")
  }

  # --- Put back into phyloseq ---
  sd_df <- as.data.frame(sd_df)
  rownames(sd_df) <- sd_df$Library_ID
  sd_df <- dplyr::select(sd_df, -Library_ID)
  sample_data(ps) <- sample_data(sd_df)

  return(ps)
}

# =============================================================================
# SECTION 2: OTU SUMMARY
# =============================================================================

#' Generate OTU summary table by sampling event
#' @param ps Phyloseq object
#' @return Data frame with OTU counts and reads per event
summarize_otus_by_event <- function(ps) {

  dbg("Generating OTU summary by sampling event...")

  otu <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) otu <- t(otu)

  samp <- as(sample_data(ps), "data.frame")
  tax <- as(tax_table(ps), "matrix") %>% as.data.frame()

  if (!"Sampling.event.ID" %in% colnames(samp)) {
    stop("Sampling.event.ID column not found in sample data")
  }

  keep_samples <- rownames(samp)[!is.na(samp$Sampling.event.ID) & samp$Sampling.event.ID != ""]
  otu <- otu[, colnames(otu) %in% keep_samples, drop = FALSE]
  samp <- samp[keep_samples, , drop = FALSE]
  samp <- samp[colnames(otu), , drop = FALSE]

  otu_pa <- (otu > 0) * 1

  event_vec <- samp$Sampling.event.ID
  names(event_vec) <- colnames(otu_pa)
  grp <- split(seq_along(event_vec), event_vec)

  event_pa <- sapply(grp, function(idx) as.integer(rowSums(otu_pa[, idx, drop = FALSE]) > 0))
  event_pa <- as.matrix(event_pa)
  rownames(event_pa) <- rownames(otu_pa)

  event_reads <- sapply(grp, function(idx) rowSums(otu[, idx, drop = FALSE]))
  event_reads <- as.matrix(event_reads)
  rownames(event_reads) <- rownames(otu)

  cn_lower <- tolower(colnames(tax))
  get_col <- function(nm) {
    i <- which(cn_lower == tolower(nm))
    if (length(i)) tax[[i[1]]] else rep(NA_character_, nrow(tax))
  }

  Species <- get_col("Species")
  Genus <- get_col("Genus")
  Family <- get_col("Family")

  res_vec <- ifelse(!is.na(Species) & Species != "", "Species",
                    ifelse(!is.na(Genus) & Genus != "", "Genus",
                           ifelse(!is.na(Family) & Family != "", "Family", "Unclassified")))
  names(res_vec) <- rownames(tax)
  res_vec <- res_vec[rownames(event_pa)]

  count_if <- function(idx) {
    if (length(idx)) colSums(event_pa[idx, , drop = FALSE]) else rep(0L, ncol(event_pa))
  }
  reads_if <- function(idx) {
    if (length(idx)) colSums(event_reads[idx, , drop = FALSE]) else rep(0L, ncol(event_reads))
  }

  species_idx <- which(res_vec == "Species")
  genus_idx   <- which(res_vec == "Genus")
  family_idx  <- which(res_vec == "Family")
  uncl_idx    <- which(res_vec == "Unclassified")

  event_area <- samp %>%
    select(Sampling.event.ID, Sampling.area.Name) %>%
    distinct()

  summary_table <- tibble(
    Sampling.event.ID = colnames(event_pa),
    n_OTUs = as.integer(colSums(event_pa)),
    Species = as.integer(count_if(species_idx)),
    Genus = as.integer(count_if(genus_idx)),
    Family = as.integer(count_if(family_idx)),
    `>Family` = as.integer(count_if(uncl_idx)),
    total_reads = as.integer(colSums(event_reads)),
    species_reads_only = as.integer(reads_if(species_idx)),
    genus_reads_only = as.integer(reads_if(genus_idx)),
    family_reads_only = as.integer(reads_if(family_idx)),
    `>Family_reads_only` = as.integer(reads_if(uncl_idx))
  ) %>%
    left_join(event_area, by = "Sampling.event.ID") %>%
    relocate(Sampling.area.Name, .after = Sampling.event.ID) %>%
    arrange(Sampling.event.ID)

  dbg("  Generated summary for", nrow(summary_table), "events")

  return(summary_table)
}

# =============================================================================
# PIPELINE WRAPPERS
# =============================================================================

#' Process metadata for a single marker (Steps 1-2)
#'
#' Loads/accepts a phyloseq object, applies NTNU metadata join and citizen
#' science metadata fixes, then generates the OTU summary by sampling event.
#'
#' @param marker Character: "12S", "18S", or "COI"
#' @param ps_input Either: a phyloseq object (already loaded), OR path to RDS
#'   file, OR NULL to use defaults
#' @param save_outputs Logical: whether to save output files
#' @param skip_preprocessing Logical: if TRUE and ps_input is phyloseq, skip
#'   metadata fixes (already done)
#' @return List with: marker, phyloseq (cleaned), otu_summary
process_metadata <- function(marker, ps_input = NULL, save_outputs = TRUE,
                             skip_preprocessing = FALSE) {

  dbg("\n", paste(rep("=", 60), collapse = ""))
  dbg("METADATA PROCESSING:", marker)
  dbg(paste(rep("=", 60), collapse = ""), "\n")

  # Step 1: Load and preprocess
  dbg("STEP 1: Loading data...")

  if (inherits(ps_input, "phyloseq")) {
    dbg("  Using provided phyloseq object")
    ps <- ps_input

    if (!skip_preprocessing) {
      ps <- add_ntnu_metadata(ps)
      ps <- fix_citsci_metadata(ps)
    } else {
      dbg("  Skipping preprocessing (skip_preprocessing = TRUE)")
    }

  } else if (is.character(ps_input) && length(ps_input) == 1) {
    ps <- load_phyloseq(marker, rds_file = ps_input)
    ps <- add_ntnu_metadata(ps)
    ps <- fix_citsci_metadata(ps)

  } else {
    ps <- load_phyloseq(marker, rds_file = NULL)
    ps <- add_ntnu_metadata(ps)
    ps <- fix_citsci_metadata(ps)
  }

  if (save_outputs) {
    out_file <- file.path(CONFIG$processed_data_dir, paste0("bge_cs_", marker, ".rds"))
    saveRDS(ps, out_file)
    dbg("  Saved:", out_file)
  }

  # Step 2: OTU summary
  dbg("\nSTEP 2: Generating OTU summary...")
  otu_summary <- summarize_otus_by_event(ps)

  if (save_outputs) {
    out_file <- file.path(CONFIG$processed_data_dir, paste0("OTU_summary_by_event_", marker, ".csv"))
    write_csv(otu_summary, out_file)
    dbg("  Saved:", out_file)
  }

  list(
    marker = marker,
    phyloseq = ps,
    otu_summary = otu_summary
  )
}

#' Process metadata for all markers in sequence
#' @param markers Vector of marker names (default: c("12S", "18S", "COI"))
#' @param ps_list Optional named list of phyloseq objects OR file paths
#' @param skip_preprocessing Logical: if TRUE and using phyloseq objects, skip metadata fixes
#' @return Named list of process_metadata results for each marker
process_all_metadata <- function(markers = c("12S", "18S", "COI"),
                                 ps_list = NULL,
                                 skip_preprocessing = FALSE) {

  results <- list()

  for (m in markers) {
    ps_input <- if (!is.null(ps_list)) ps_list[[m]] else NULL
    results[[m]] <- process_metadata(m, ps_input, skip_preprocessing = skip_preprocessing)
  }

  dbg("\n", paste(rep("=", 60), collapse = ""))
  dbg("ALL METADATA PROCESSED SUCCESSFULLY")
  dbg(paste(rep("=", 60), collapse = ""))

  return(results)
}

# =============================================================================
cat("\n")
cat("=============================================================\n")
cat("METADATA PROCESSING FUNCTIONS LOADED (Steps 1-2)\n")
cat("=============================================================\n\n")
cat("  process_metadata('12S', ps_input)  - Clean metadata + OTU summary\n")
cat("  process_all_metadata()             - All markers at once\n")
cat("=============================================================\n")
