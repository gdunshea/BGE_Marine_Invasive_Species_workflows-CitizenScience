# =============================================================================
# Wheel of Life - R Wrapper for Python Generator
# =============================================================================

library(phyloseq)
library(dplyr)

# =============================================================================
# DATA EXTRACTION
# =============================================================================

extract_wheel_data <- function(ps, site_var = "Site", label_rank = "Species", 
                                min_reads = 1, chordate_by_class = TRUE,
                                collector_var = NULL, location_var = NULL, 
                                date_var = NULL, include_genus = TRUE) {
  
  stopifnot(site_var %in% colnames(sample_data(ps)))
  tax_cols <- colnames(tax_table(ps))
  sample_cols <- colnames(sample_data(ps))
  
  cat("Extracting data with full taxonomy...\n")
  df <- psmelt(ps)
  
  df$Site <- df[[site_var]]
  
  # Create label: use Species if available, otherwise "Genus sp." if genus available
  has_genus <- "Genus" %in% tax_cols
  has_species <- "Species" %in% tax_cols
  
  if (has_species && has_genus && include_genus) {
    df <- df %>%
      mutate(label = case_when(
        !is.na(Species) & Species != "" & !grepl("^Unknown", Species) ~ Species,
        !is.na(Genus) & Genus != "" & !grepl("^Unknown", Genus) ~ paste0(Genus, " sp."),
        TRUE ~ NA_character_
      ))
  } else if (has_species) {
    df$label <- df$Species
  } else if (has_genus) {
    df$label <- paste0(df$Genus, " sp.")
  } else {
    df$label <- df[[label_rank]]
  }
  
  has_phylum <- "Phylum" %in% tax_cols
  has_class <- "Class" %in% tax_cols
  
  if (chordate_by_class && has_phylum && has_class) {
    df <- df %>%
      mutate(display_group = case_when(
        Phylum == "Chordata" & !is.na(Class) ~ Class,
        !is.na(Phylum) ~ Phylum,
        !is.na(Class) ~ Class,
        TRUE ~ "Unknown"
      ))
  } else if (has_phylum) {
    df$display_group <- df$Phylum
  } else {
    df$display_group <- df$Class
  }
  
  # Filter to minimum reads and valid sites
  df <- df %>%
    dplyr::filter(Abundance >= min_reads, !is.na(Site) & Site != "")
  
  # Count TOTAL OTUs per site (before label filter) 
  otu_counts <- df %>%
    group_by(Site) %>%
    summarise(
      n_otus = n_distinct(OTU),
      n_otu_groups = n_distinct(display_group),
      .groups = 'drop'
    )
  
  # Now filter to only taxa with valid labels (species or genus sp.)
  df <- df %>%
    dplyr::filter(!is.na(label) & label != "")
  
  tax_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  available_ranks <- intersect(tax_ranks, tax_cols)
  
  # Base columns to select
  select_cols <- c("Site", "label", "display_group", available_ranks)
  
  # Add metadata columns if specified
  cat("  Checking metadata columns...\n")
  cat("    collector_var:", ifelse(is.null(collector_var), "NULL", collector_var), "\n")
  cat("    location_var:", ifelse(is.null(location_var), "NULL", location_var), "\n")
  cat("    date_var:", ifelse(is.null(date_var), "NULL", date_var), "\n")
  cat("    Available sample_cols:", paste(head(sample_cols, 10), collapse = ", "), "...\n")
  cat("    Available df cols:", paste(head(colnames(df), 15), collapse = ", "), "...\n")
  
  if (!is.null(collector_var)) {
    if (collector_var %in% colnames(df)) {
      df$Collector <- df[[collector_var]]
      select_cols <- c(select_cols, "Collector")
      cat("    -> Added Collector from:", collector_var, "\n")
    } else {
      cat("    -> WARNING: collector_var '", collector_var, "' not found in df\n")
    }
  }
  if (!is.null(location_var)) {
    if (location_var %in% colnames(df)) {
      df$Location <- df[[location_var]]
      select_cols <- c(select_cols, "Location")
      cat("    -> Added Location from:", location_var, "\n")
    } else {
      cat("    -> WARNING: location_var '", location_var, "' not found in df\n")
    }
  }
  if (!is.null(date_var)) {
    if (date_var %in% colnames(df)) {
      df$Date <- as.character(df[[date_var]])
      select_cols <- c(select_cols, "Date")
      cat("    -> Added Date from:", date_var, "\n")
    } else {
      cat("    -> WARNING: date_var '", date_var, "' not found in df\n")
    }
  }
  
  result <- df %>% 
    dplyr::select(all_of(select_cols)) %>% 
    distinct()
  
  # Join OTU counts back
  result <- result %>%
    left_join(otu_counts, by = "Site")
  
  for (rank in available_ranks) {
    result[[rank]] <- ifelse(is.na(result[[rank]]) | result[[rank]] == "", 
                              paste0("Unknown_", rank), 
                              result[[rank]])
  }
  
  cat("  Sites:", length(unique(result$Site)), "\n")
  cat("  Taxa records:", nrow(result), "\n")
  cat("  Columns exported:", paste(colnames(result), collapse = ", "), "\n")
  
  # Count genus-level vs species-level
  genus_only <- sum(grepl(" sp\\.$", result$label))
  species_level <- nrow(result) - genus_only
  cat("  Species-level:", species_level, ", Genus-level:", genus_only, "\n")
  
  # Show sample of metadata if present
  if ("Collector" %in% colnames(result)) {
    sample_collectors <- head(unique(result$Collector[!is.na(result$Collector)]), 3)
    cat("  Sample Collectors:", paste(sample_collectors, collapse = ", "), "\n")
  }
  
  return(result)
}


# =============================================================================
# MAIN FUNCTION: Generate Wheel of Life plots
# =============================================================================

#' Generate Wheel of Life plots from phyloseq objects
#'
#' @param ... Named phyloseq objects or data frames (e.g., `12S` = ps_12s, `18S` = ps_18s)
#' @param output_dir Output directory for PDF files
#' @param prefix Prefix for output filenames
#' @param site_var Sample variable for grouping/titles (determines how many plots are generated)
#' @param collector_var Sample variable containing collector name
#' @param location_var Sample variable containing location name
#' @param date_var Sample variable containing date
#' @param metadata_from Which marker to extract metadata from (default: first one). 
#'        Should match one of the names you give to the phyloseq objects.
#' @param center_image Path to image file to display in center of wheel (default: NULL, no image)
#' @param center_image_zoom Zoom level for center image (default: 0.04, smaller = smaller image)
#' @param include_genus Include genus-level IDs as "Genus sp." (default: TRUE)
#' @param min_group_size Minimum taxa per group to show separately (default: 3).
#'        Groups with fewer taxa are combined into "Other Phyla" (colors preserved).
#' @param min_species Minimum species per site to generate a wheel
#' @param script_dir Directory containing mm_wheel_of_life.py
#' @param python_path Path to Python executable (default: "python3")
#' @param temp_dir Directory for temporary CSV files (default: tempdir())
#'
#' @examples
#' generate_wheels(
#'   `12S` = ps_12s,
#'   `18S` = ps_18s,
#'   output_dir = 'Figures/wheels',
#'   site_var = 'Sampling.event.ID',
#'   collector_var = 'Collector.name',
#'   location_var = 'Site.name',
#'   date_var = 'Collection.date',
#'   metadata_from = '12S',
#'   center_image = 'path/to/logo.png',
#'   center_image_zoom = 0.04,
#'   include_genus = TRUE,
#'   min_group_size = 3,
#'   script_dir = 'Scripts/'
#' )

generate_wheels <- function(..., 
                            output_dir = "Figures/wheels",
                            prefix = "wheel",
                            site_var = "Sampling.event.ID",
                            collector_var = NULL,
                            location_var = NULL,
                            date_var = NULL,
                            metadata_from = NULL,
                            center_image = NULL,
                            center_image_zoom = 0.04,
                            include_genus = TRUE,
                            min_group_size = 3,
                            min_species = 5,
                            script_dir = ".",
                            python_path = NULL,
                            temp_dir = tempdir()) {
  
  # Find Python executable
  if (is.null(python_path)) {
    # Try to find Python with pandas
    python_candidates <- c(
      Sys.which("python3"),
      "/usr/local/bin/python3",
      "/opt/homebrew/bin/python3",
      path.expand("~/miniconda3/bin/python"),
      path.expand("~/anaconda3/bin/python"),
      "/usr/bin/python3"
    )
    
    python_path <- "python3"  # default fallback
    for (candidate in python_candidates) {
      if (nchar(candidate) > 0 && file.exists(candidate)) {
        # Test if pandas is available
        test <- suppressWarnings(
          system2(candidate, args = c("-c", "import pandas"), 
                  stdout = FALSE, stderr = FALSE)
        )
        if (test == 0) {
          python_path <- candidate
          break
        }
      }
    }
  }
  
  cat("Using Python:", python_path, "\n")
  
  # Build path to Python script
  python_script <- normalizePath(file.path(script_dir, "mm_wheel_of_life.py"), mustWork = FALSE)
  
  if (!file.exists(python_script)) {
    stop(paste("Python script not found:", python_script,
               "\nPlease set script_dir to the directory containing mm_wheel_of_life.py"))
  }
  
  # Normalize output directory path
  output_dir <- normalizePath(output_dir, mustWork = FALSE)
  
  # Get named arguments
  inputs <- list(...)
  marker_names <- names(inputs)
  
  if (is.null(marker_names) || any(marker_names == "")) {
    stop("All inputs must be named (e.g., `12S` = ps_12s, `18S` = ps_18s)")
  }
  
  cat("\n========================================\n")
  cat("  WHEEL OF LIFE GENERATOR\n")
  cat("========================================\n\n")
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
  }
  
  # Determine which marker to get metadata from (default: first one)
  if (is.null(metadata_from)) {
    metadata_from <- marker_names[1]
  }
  if (!metadata_from %in% marker_names) {
    stop(paste("metadata_from must be one of:", paste(marker_names, collapse = ", ")))
  }
  cat("Metadata (Collector, Location, Date) will be extracted from:", metadata_from, "\n")
  
  # Extract and save data for each marker
  csv_files <- list()
  for (name in marker_names) {
    obj <- inputs[[name]]
    # Sanitize filename - replace spaces and special chars
    safe_name <- gsub("[^a-zA-Z0-9]", "_", name)
    csv_path <- file.path(temp_dir, paste0("wheel_data_", safe_name, ".csv"))
    
    cat("\nProcessing", name, "marker...\n")
    
    # Only extract metadata columns from the specified marker
    if (name == metadata_from) {
      use_collector <- collector_var
      use_location <- location_var
      use_date <- date_var
    } else {
      use_collector <- NULL
      use_location <- NULL
      use_date <- NULL
    }
    
    if (inherits(obj, "phyloseq")) {
      df <- extract_wheel_data(obj, site_var = site_var,
                               collector_var = use_collector,
                               location_var = use_location,
                               date_var = use_date,
                               include_genus = include_genus)
    } else if (is.data.frame(obj)) {
      df <- obj
    } else {
      stop(paste("Input", name, "must be a phyloseq object or data frame"))
    }
    
    write.csv(df, csv_path, row.names = FALSE)
    csv_files[[name]] <- csv_path
    cat("  Saved to:", csv_path, "\n")
  }
  
  # Build Python command
  csv_dict <- paste0("'", marker_names, "': '", unlist(csv_files), "'", collapse = ", ")
  
  # Handle center_image parameter
  if (!is.null(center_image) && file.exists(center_image)) {
    center_image_path <- normalizePath(center_image)
    center_image_arg <- sprintf('center_image="%s", center_image_zoom=%f', center_image_path, center_image_zoom)
    cat("Center image:", center_image_path, "(zoom:", center_image_zoom, ")\n")
  } else {
    center_image_arg <- "center_image=None"
    if (!is.null(center_image)) {
      cat("Warning: center_image file not found:", center_image, "\n")
    }
  }
  
  # Handle include_genus parameter
  include_genus_arg <- ifelse(include_genus, "include_genus=True", "include_genus=False")
  
  # Handle min_group_size parameter
  min_group_size_arg <- sprintf("min_group_size=%d", min_group_size)
  
  python_code <- sprintf('
import pandas as pd
import sys
import traceback

try:
    exec(open("%s").read())
    
    combined_df, marker_names = load_and_combine_markers({%s})
    
    generate_site_wheels(
        combined_df, 
        output_dir="%s", 
        prefix="%s",
        min_species=%d,
        marker_names=marker_names,
        %s,
        %s,
        %s
    )
except Exception as e:
    print("ERROR:", str(e))
    traceback.print_exc()
    sys.exit(1)
', python_script, csv_dict, output_dir, prefix, min_species, center_image_arg, include_genus_arg, min_group_size_arg)
  
  # Write Python code to temp file
  py_script <- file.path(temp_dir, "run_wheel.py")
  writeLines(python_code, py_script)
  
  cat("\n--- Running Python generator ---\n")
  cat("Script:", py_script, "\n\n")
  
  # Build command string
  cmd <- paste(shQuote(python_path), shQuote(py_script), "2>&1")
  
  # Run Python
  result <- system(cmd, intern = TRUE, ignore.stderr = FALSE)
  
  # Check for errors
  exit_status <- attr(result, "status")
  if (!is.null(exit_status) && exit_status != 0) {
    cat("Python exited with status:", exit_status, "\n")
  }
  
  if (length(result) > 0) {
    cat(paste(result, collapse = "\n"), "\n")
  }
  
  cat("\n========================================\n")
  cat("  COMPLETE\n")
  cat("  Output:", output_dir, "\n")
  cat("========================================\n")
  
  invisible(output_dir)
}


# =============================================================================
# SINGLE MARKER SHORTCUT
# =============================================================================

#' Generate Wheel of Life plots from a single phyloseq object
#'
#' @param ps Phyloseq object
#' @param marker_name Name for this marker (e.g., "12S")
#' @param script_dir Directory containing mm_wheel_of_life.py
#' @param python_path Path to Python executable
#' @param ... Additional arguments passed to generate_wheels()

generate_wheels_single <- function(ps, marker_name = "marker", script_dir = ".", python_path = NULL, ...) {
  args <- list(ps)
  names(args) <- marker_name
  do.call(generate_wheels, c(args, list(script_dir = script_dir, python_path = python_path, ...)))
}


# =============================================================================
# USAGE EXAMPLES
# =============================================================================

cat("
========================================
  WHEEL OF LIFE - R Interface
========================================

USAGE:

  # Load this script
  source('wheel_of_life_R.R')

  # Basic usage:
  generate_wheels(
    `12S` = ps_12s,
    output_dir = 'Figures/wheels',
    prefix = '12S',
    site_var = 'Sampling.event.ID',
    script_dir = 'path/to/scripts'
  )

  # Multiple markers with metadata and center image:
  generate_wheels(
    `Vertebrates 12S` = ps_12s,
    `Eukaryotes 18S` = ps_18s,
    output_dir = 'Figures/wheels',
    prefix = 'combined',
    site_var = 'Sampling.event.ID',
    collector_var = 'Collector.name',
    location_var = 'Site.name',
    date_var = 'Collection.date',
    metadata_from = 'Vertebrates 12S',
    center_image = 'path/to/logo.png',
    center_image_zoom = 0.04,
    script_dir = 'path/to/scripts'
  )

  # If Python is not found automatically, specify path:
  generate_wheels(
    `12S` = ps_12s,
    script_dir = 'path/to/scripts',
    python_path = '/opt/homebrew/bin/python3'
  )

PARAMETERS:
  site_var          - Column for plot grouping/titles (must exist in ALL phyloseq objects)
  collector_var     - Column for collector name (from metadata_from phyloseq only)
  location_var      - Column for location name (from metadata_from phyloseq only)
  date_var          - Column for date (from metadata_from phyloseq only)
  metadata_from     - Which marker's phyloseq to get Collector/Location/Date from
                      (default: first marker). Must match one of the names you give.
  center_image      - Path to image file (PNG, JPG, etc.) to display in center of wheel
                      (default: NULL, no image)
  center_image_zoom - Zoom level for center image (default: 0.04)
                      Smaller number = smaller image. Try 0.02-0.1
  include_genus     - Include genus-level IDs as 'Genus sp.' (default: TRUE)
                      Set to FALSE for species-only
  min_group_size    - Minimum taxa per group to display separately (default: 3)
                      Groups with fewer taxa are combined into 'Other Phyla'
                      (individual phylum colors are preserved)

  # Find your Python path in terminal: which python3

REQUIREMENTS:
  - Python 3 with pandas, numpy, matplotlib
  - mm_wheel_of_life.py (specify location with script_dir)

========================================
")
