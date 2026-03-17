# =============================================================================
# FUNC_combine_lifestyle.R
# =============================================================================
#
# Combines the manual taxonomy-based lifestyle mapping (FUNC_assign_lifestyle.R)
# with WoRMS functional group data (FUNC_worms_lifestyle.R) into a single
# validated lifestyle column.
#
# Strategy:
#   1. Manual mapping is PRIMARY (100% coverage via Phylum/Class/Order lookup)
#   2. WoRMS functional groups are mapped to the same 8 categories
#   3. Where both exist: flag disagreements for review
#   4. Final lifestyle = manual, unless WoRMS provides a more specific
#      correction (user decides)
#
# USAGE:
#   source("Scripts/Final/FUNC_combine_lifestyle.R")
#
#   # After running assign_lifestyle() and loading WoRMS cache:
#   worms_traits <- read.csv("Processed_data/worms_functional_groups.csv")
#   validated <- combine_lifestyle(combined_all$site_level, worms_traits)
#   combined_all$site_level <- validated$site_level
#
#   # Check disagreements:
#   validated$disagreements
#
# =============================================================================

library(dplyr)
library(stringr)


#' Map WoRMS functional group strings to the 8 standard lifestyle categories
#'
#' @param fg_string Character vector of WoRMS functional_group values
#'   (may contain "; " separated multiple groups per species)
#' @param mapping Named character vector mapping WoRMS terms to lifestyle.
#'   Default handles all observed WoRMS functional groups.
#' @return Character vector of lifestyle categories
map_worms_to_lifestyle <- function(
    fg_string,
    mapping = c(
      "plankton > phytoplankton"       = "Phytoplankton",
      "plankton > zooplankton"         = "Zooplankton",
      "plankton > microplankton"       = "Microeukaryote",
      "plankton > macroplankton"       = "Zooplankton",
      "plankton > mesoplankton"        = "Zooplankton",
      "plankton > megaplankton"        = "Zooplankton",
      "plankton > nanoplankton"        = "Phytoplankton",
      "plankton > mixoplankton > Constitutive Mixoplankton"  = "Phytoplankton",
      "plankton > mixoplankton > General Non-Constitutive Mixoplankton" = "Microeukaryote",
      "plankton > mixoplankton > endosymbiotic Specialist Non-Constitutive Mixoplankton" = "Microeukaryote",
      "plankton > mixoplankton > plastid Specialists Non-Constitutive Mixoplankton" = "Microeukaryote",
      "benthos"                        = "Sessile invertebrate",
      "benthos > macrobenthos"         = "Sessile invertebrate",
      "benthos > megabenthos"          = "Mobile invertebrate",
      "benthos > meiobenthos"          = "Mobile invertebrate",
      "benthos > microbenthos"         = "Microeukaryote"
    )
) {
  
  sapply(fg_string, function(fg) {
    if (is.na(fg) || fg == "" || fg == "NA") return(NA_character_)
    
    # Take first functional group if multiple
    groups <- str_split(fg, "; ")[[1]]
    
    # Try to match each group in order
    for (g in groups) {
      g <- str_trim(g)
      if (g %in% names(mapping)) return(mapping[g])
    }
    
    # Partial matching fallback
    g1 <- groups[1]
    if (grepl("phytoplankton", g1, ignore.case = TRUE)) return("Phytoplankton")
    if (grepl("zooplankton", g1, ignore.case = TRUE)) return("Zooplankton")
    if (grepl("benthos", g1, ignore.case = TRUE)) return("Sessile invertebrate")
    if (grepl("plankton", g1, ignore.case = TRUE)) return("Zooplankton")
    
    NA_character_
  }, USE.NAMES = FALSE)
}


#' Combine manual and WoRMS lifestyle assignments
#'
#' @param site_level Data frame with `lifestyle` column from assign_lifestyle()
#' @param worms_traits Data frame from read.csv("worms_functional_groups.csv")
#'   with columns: Species, AphiaID, functional_group
#' @param worms_mapping Optional custom mapping for WoRMS terms (see
#'   map_worms_to_lifestyle for default)
#' @param prefer Character: "manual" (default) or "worms". When they disagree,
#'   which source wins? "manual" keeps the taxonomy-based assignment,
#'   "worms" uses the WoRMS functional group where available.
#' @return List with:
#'   $site_level     - Updated data frame with lifestyle_worms, lifestyle_final
#'   $disagreements  - Species where manual != worms (for review)
#'   $coverage       - Summary of WoRMS coverage and agreement
combine_lifestyle <- function(
    site_level,
    worms_traits,
    worms_mapping = NULL,
    prefer = "manual"
) {
  
  # Map WoRMS functional groups to lifestyle categories
  worms_mapped <- worms_traits %>%
    mutate(
      lifestyle_worms = if (is.null(worms_mapping)) {
        map_worms_to_lifestyle(functional_group)
      } else {
        map_worms_to_lifestyle(functional_group, mapping = worms_mapping)
      }
    ) %>%
    select(Species, functional_group, lifestyle_worms) %>%
    distinct(Species, .keep_all = TRUE)
  
  # Join to site_level
  result <- site_level %>%
    left_join(worms_mapped, by = "Species")
  
  # Identify disagreements (where both exist and differ)
  disagreements <- result %>%
    filter(!is.na(lifestyle_worms) & !is.na(lifestyle) &
             as.character(lifestyle) != lifestyle_worms) %>%
    distinct(Species, Phylum, Class, Order, lifestyle, lifestyle_worms, functional_group) %>%
    arrange(Species)
  
  # Set final lifestyle
  if (prefer == "worms") {
    result <- result %>%
      mutate(
        lifestyle_final = ifelse(!is.na(lifestyle_worms), lifestyle_worms,
                                  as.character(lifestyle)),
        lifestyle_source = ifelse(!is.na(lifestyle_worms), "worms", "manual")
      )
  } else {
    # Manual is primary; keep it unless user overrides
    result <- result %>%
      mutate(
        lifestyle_final = as.character(lifestyle),
        lifestyle_source = ifelse(!is.na(lifestyle_worms), "both", "manual")
      )
  }
  
  # Coverage summary
  n_total <- length(unique(na.omit(result$Species)))
  n_worms <- result %>%
    filter(!is.na(lifestyle_worms) & !is.na(Species)) %>%
    pull(Species) %>% unique() %>% length()
  n_agree <- result %>%
    filter(!is.na(lifestyle_worms) & as.character(lifestyle) == lifestyle_worms &
             !is.na(Species)) %>%
    pull(Species) %>% unique() %>% length()
  n_disagree <- nrow(disagreements)
  
  coverage <- tibble(
    metric = c("Total species", "WoRMS coverage", "Agreement", "Disagreement"),
    n = c(n_total, n_worms, n_agree, n_disagree),
    pct = round(c(100, n_worms/n_total*100, n_agree/max(n_worms,1)*100,
                   n_disagree/max(n_worms,1)*100), 1)
  )
  
  # Report
  cat("=== Lifestyle validation: manual vs WoRMS ===\n")
  print(coverage)
  
  if (n_disagree > 0) {
    cat("\n=== Disagreements (", n_disagree, "species) ===\n")
    print(disagreements, n = Inf, width = Inf)
    cat("\nPreference set to '", prefer, "'. lifestyle_final uses ",
        ifelse(prefer == "manual", "manual mapping", "WoRMS"), ".\n")
    cat("To override specific species, edit lifestyle_final directly.\n")
  } else {
    cat("\nNo disagreements found — manual and WoRMS agree where both exist.\n")
  }
  
  list(
    site_level = result,
    disagreements = disagreements,
    coverage = coverage
  )
}


# =============================================================================
# LOADED
# =============================================================================

cat("FUNC_combine_lifestyle.R loaded\n")
cat("Functions:\n")
cat("  map_worms_to_lifestyle(fg_string)  - Map WoRMS terms to 8 categories\n")
cat("  combine_lifestyle(site_level, worms_traits)  - Validate & combine\n")
