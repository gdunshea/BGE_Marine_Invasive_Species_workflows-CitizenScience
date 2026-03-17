# =============================================================================
# FUNC_assign_lifestyle.R
# =============================================================================
#
# Assigns ecological lifestyle categories to OTU × site detections based on
# taxonomy (Phylum, Class, Order). Uses a hierarchical lookup: Order-level
# matches take priority, then Class, then Phylum. Unmatched taxa get "Other".
#
# Categories:
#   Fish, Sessile invertebrate, Mobile invertebrate, Phytoplankton,
#   Zooplankton, Macroalgae, Microeukaryote, Other
#
# USAGE:
#   source("Scripts/Final/FUNC_assign_lifestyle.R")
#
#   # Default mapping:
#   combined_all$site_level <- assign_lifestyle(combined_all$site_level)
#
#   # Custom mapping — override or extend:
#   my_lookup <- default_lifestyle_lookup()
#   my_lookup <- bind_rows(my_lookup,
#     tibble(level = "Class", taxon = "MyNewClass", lifestyle = "Zooplankton")
#   )
#   combined_all$site_level <- assign_lifestyle(combined_all$site_level, 
#                                                lookup = my_lookup)
#
#   # Inspect unmatched taxa:
#   check_unmatched_lifestyle(combined_all$site_level)
#
# =============================================================================

library(dplyr)
library(tibble)


#' Generate the default lifestyle lookup table
#'
#' Returns a tibble with columns: level (Order/Class/Phylum), taxon, lifestyle.
#' Order-level matches take priority over Class, Class over Phylum.
#' Edit the returned table to customise assignments.
#'
#' @return tibble with lifestyle mapping rules
#'
#' @examples
#' lookup <- default_lifestyle_lookup()
#' print(lookup, n = Inf)
#'
#' # Modify: reclassify Scyphozoa from Zooplankton to Mobile invertebrate
#' lookup$lifestyle[lookup$taxon == "Scyphozoa"] <- "Mobile invertebrate"
default_lifestyle_lookup <- function() {
  
  tribble(
    ~level,   ~taxon,                              ~lifestyle,
    
    # =====================================================================
    # ORDER-LEVEL OVERRIDES (take priority over Class/Phylum)
    # =====================================================================
    # These handle classes with mixed lifestyles
    
    # Hexanauplia — mostly zooplankton, except barnacles
    "Order",  "Calanoida",                          "Zooplankton",
    "Order",  "Cyclopoida",                         "Zooplankton",
    "Order",  "Harpacticoida",                      "Zooplankton",
    "Order",  "Poecilostomatoida",                  "Zooplankton",
    "Order",  "Siphonostomatoida",                  "Zooplankton",
    "Order",  "Sessilia",                           "Sessile invertebrate",
    "Order",  "Lithoglyptida",                      "Sessile invertebrate",
    "Order",  "Podocopida",                         "Mobile invertebrate",
    
    # Hydrozoa — colonial/sessile vs pelagic
    "Order",  "Leptothecata",                       "Sessile invertebrate",
    "Order",  "Anthoathecata",                      "Sessile invertebrate",
    "Order",  "Siphonophorae",                      "Zooplankton",
    "Order",  "Trachymedusae",                      "Zooplankton",
    "Order",  "Narcomedusae",                       "Zooplankton",
    "Order",  "Limnomedusae",                       "Zooplankton",
    
    # Gastropoda — mostly benthic mobile, except pteropods
    "Order",  "Pteropoda",                          "Zooplankton",
    
    # =====================================================================
    # CLASS-LEVEL ASSIGNMENTS
    # =====================================================================
    
    # --- Fish ---
    "Class",  "Actinopteri",                        "Fish",
    "Class",  "Chondrichthyes",                     "Fish",
    
    # --- Sessile invertebrates ---
    "Class",  "Gymnolaemata",                       "Sessile invertebrate",
    "Class",  "Stenolaemata",                       "Sessile invertebrate",
    "Class",  "Phylactolaemata",                    "Sessile invertebrate",
    "Class",  "Thecostraca",                        "Sessile invertebrate",
    "Class",  "Demospongiae",                       "Sessile invertebrate",
    "Class",  "Calcarea",                           "Sessile invertebrate",
    "Class",  "Homoscleromorpha",                   "Sessile invertebrate",
    "Class",  "Anthozoa",                           "Sessile invertebrate",
    "Class",  "Staurozoa",                          "Sessile invertebrate",
    "Class",  "Ascidiacea",                         "Sessile invertebrate",
    "Class",  "Bivalvia",                           "Sessile invertebrate",
    "Class",  "Polychaeta",                         "Sessile invertebrate",
    "Class",  "Tentaculata",                        "Sessile invertebrate",
    "Class",  "Phoroniformea",                      "Sessile invertebrate",
    
    # --- Mobile invertebrates ---
    "Class",  "Malacostraca",                       "Mobile invertebrate",
    "Class",  "Gastropoda",                         "Mobile invertebrate",
    "Class",  "Polyplacophora",                     "Mobile invertebrate",
    "Class",  "Cephalopoda",                        "Mobile invertebrate",
    "Class",  "Asteroidea",                         "Mobile invertebrate",
    "Class",  "Echinoidea",                         "Mobile invertebrate",
    "Class",  "Holothuroidea",                      "Mobile invertebrate",
    "Class",  "Ophiuroidea",                        "Mobile invertebrate",
    "Class",  "Sipuncula",                          "Mobile invertebrate",
    "Class",  "Chromadorea",                        "Mobile invertebrate",
    "Class",  "Enoplea",                            "Mobile invertebrate",
    "Class",  "Gastrotrichea",                      "Mobile invertebrate",
    "Class",  "Rhabditophora",                      "Mobile invertebrate",
    "Class",  "Clitellata",                         "Mobile invertebrate",
    "Class",  "Flabellinia",                        "Mobile invertebrate",
    "Class",  "Acoelomorpha",                       "Mobile invertebrate",
    
    # --- Zooplankton ---
    "Class",  "Hexanauplia",                        "Zooplankton",
    "Class",  "Branchiopoda",                       "Zooplankton",
    "Class",  "Scyphozoa",                          "Zooplankton",
    "Class",  "Appendicularia",                     "Zooplankton",
    "Class",  "Eurotatoria",                        "Zooplankton",
    "Class",  "Monogononta",                        "Zooplankton",
    "Class",  "Sagittoidea",                        "Zooplankton",
    "Class",  "Hydrozoa",                           "Sessile invertebrate",
    
    # --- Phytoplankton ---
    "Class",  "Dinophyceae",                        "Phytoplankton",
    "Class",  "Mediophyceae",                       "Phytoplankton",
    "Class",  "Bacillariophyceae",                  "Phytoplankton",
    "Class",  "Coscinodiscophyceae",                "Phytoplankton",
    "Class",  "Fragilariophyceae",                  "Phytoplankton",
    "Class",  "Bolidophyceae",                      "Phytoplankton",
    "Class",  "Cryptophyceae",                      "Phytoplankton",
    "Class",  "Mamiellophyceae",                    "Phytoplankton",
    "Class",  "Pyramimonadophyceae",                "Phytoplankton",
    "Class",  "Chlorodendrophyceae",                "Phytoplankton",
    "Class",  "Chlorophyceae",                      "Phytoplankton",
    "Class",  "Trebouxiophyceae",                   "Phytoplankton",
    "Class",  "Chloropicophyceae",                  "Phytoplankton",
    "Class",  "Nephroselmidophyceae",               "Phytoplankton",
    "Class",  "Pedinophyceae",                      "Phytoplankton",
    "Class",  "Pseudoscourfieldiophyceae",          "Phytoplankton",
    "Class",  "Palmophyllophyceae",                 "Phytoplankton",
    "Class",  "Picomonadea",                        "Phytoplankton",
    "Class",  "Eustigmatophyceae",                  "Phytoplankton",
    "Class",  "Pelagophyceae",                      "Phytoplankton",
    "Class",  "Dictyochophyceae",                   "Phytoplankton",
    "Class",  "Raphidophyceae",                     "Phytoplankton",
    "Class",  "Chrysophyceae",                      "Phytoplankton",
    "Class",  "Synurophyceae",                      "Phytoplankton",
    "Class",  "Xanthophyceae",                      "Phytoplankton",
    "Class",  "Ellobiophyceae",                     "Phytoplankton",
    "Class",  "Urgorriphyceae",                     "Phytoplankton",
    "Class",  "Compsopogonophyceae",                "Phytoplankton",
    "Class",  "Bangiophyceae",                      "Phytoplankton",
    "Class",  "Chlorarachniophyceae",               "Phytoplankton",
    
    # --- Macroalgae ---
    "Class",  "Florideophyceae",                    "Macroalgae",
    "Class",  "Phaeophyceae",                       "Macroalgae",
    "Class",  "Ulvophyceae",                        "Macroalgae",
    
    # --- Microeukaryote (heterotrophic protists, MAST, ciliates, fungi, etc.) ---
    "Class",  "Spirotrichea",                       "Microeukaryote",
    "Class",  "Phyllopharyngea",                    "Microeukaryote",
    "Class",  "Litostomatea",                       "Microeukaryote",
    "Class",  "Oligohymenophorea",                  "Microeukaryote",
    "Class",  "Prostomatea",                        "Microeukaryote",
    "Class",  "Karyorelictea",                      "Microeukaryote",
    "Class",  "Nassophorea",                        "Microeukaryote",
    "Class",  "Thecofilosea",                       "Microeukaryote",
    "Class",  "Imbricatea",                         "Microeukaryote",
    "Class",  "Sarcomonadea",                       "Microeukaryote",
    "Class",  "Granofilosea",                       "Microeukaryote",
    "Class",  "Phaeodarea",                         "Microeukaryote",
    "Class",  "Choanoflagellatea",                  "Microeukaryote",
    "Class",  "Choanoflagellata",                   "Microeukaryote",
    "Class",  "Labyrinthulomycetes",                "Microeukaryote",
    "Class",  "Acantharea",                         "Microeukaryote",
    "Class",  "Polycystinea",                       "Microeukaryote",
    "Class",  "Centrohelea",                        "Microeukaryote",
    "Class",  "Nanomonadea",                        "Microeukaryote",
    "Class",  "Telonemidea",                        "Microeukaryote",
    "Class",  "Katablepharidophceae",               "Microeukaryote",
    "Class",  "Conoidasida",                        "Microeukaryote",
    "Class",  "Perkinsea",                          "Microeukaryote",
    "Class",  "Apusomonadea",                       "Microeukaryote",
    "Class",  "Colpodellidea",                      "Microeukaryote",
    "Class",  "Longamoebia",                        "Microeukaryote",
    "Class",  "Taxopodea",                          "Microeukaryote",
    "Class",  "Aquavolonidea",                      "Microeukaryote",
    "Class",  "Nuclearidea",                        "Microeukaryote",
    "Class",  "Ichthyosporea",                      "Microeukaryote",
    "Class",  "Pseudophyllomitidea",                "Microeukaryote",
    "Class",  "Bicoecea",                           "Microeukaryote",
    "Class",  "Bigyra",                             "Microeukaryote",
    "Class",  "Phytomyxea",                         "Microeukaryote",
    
    # Fungi
    "Class",  "Peronosporomycetes",                 "Microeukaryote",
    "Class",  "Saprolegniomycetes",                 "Microeukaryote",
    "Class",  "Blastocladiomycetes",                "Microeukaryote",
    "Class",  "Rhizophydiomycetes",                 "Microeukaryote",
    "Class",  "Lobulomycetes",                      "Microeukaryote",
    "Class",  "Pichiomycetes",                      "Microeukaryote",
    "Class",  "Agaricomycetes",                     "Microeukaryote",
    "Class",  "Microbotryomycetes",                 "Microeukaryote",
    "Class",  "Tremellomycetes",                    "Microeukaryote",
    "Class",  "Exobasidiomycetes",                  "Microeukaryote",
    "Class",  "Dothideomycetes",                    "Microeukaryote",
    "Class",  "Leotiomycetes",                      "Microeukaryote",
    "Class",  "Sordariomycetes",                    "Microeukaryote",
    "Class",  "Lecanoromycetes",                    "Microeukaryote",
    "Class",  "Entomophthoromycetes",               "Microeukaryote",
    "Class",  "Aphelidiomycetes",                   "Microeukaryote",
    
    # MAST groups (all microeukaryote)
    "Class",  "MAST-1A",                            "Microeukaryote",
    "Class",  "MAST-1B",                            "Microeukaryote",
    "Class",  "MAST-1C",                            "Microeukaryote",
    "Class",  "MAST-1D",                            "Microeukaryote",
    "Class",  "MAST-2A",                            "Microeukaryote",
    "Class",  "MAST-2B",                            "Microeukaryote",
    "Class",  "MAST-2D",                            "Microeukaryote",
    "Class",  "MAST-4A",                            "Microeukaryote",
    "Class",  "MAST-4C",                            "Microeukaryote",
    "Class",  "MAST-4D",                            "Microeukaryote",
    "Class",  "MAST-4E",                            "Microeukaryote",
    "Class",  "MAST-7A",                            "Microeukaryote",
    "Class",  "MAST-7B",                            "Microeukaryote",
    "Class",  "MAST-9A",                            "Microeukaryote",
    "Class",  "MAST-9D",                            "Microeukaryote",
    "Class",  "MAST-10",                            "Microeukaryote",
    "Class",  "MAST-11",                            "Microeukaryote",
    "Class",  "MAST-12A",                           "Microeukaryote",
    "Class",  "MAST-12B",                           "Microeukaryote",
    "Class",  "MAST-12D",                           "Microeukaryote",
    "Class",  "MAST-12E",                           "Microeukaryote",
    "Class",  "MAST-12_cl_incertae_sedis",          "Microeukaryote",
    "Class",  "MOCH-2",                             "Microeukaryote",
    "Class",  "MOCH-5",                             "Microeukaryote",
    "Class",  "GS07",                               "Microeukaryote",
    "Class",  "GS08",                               "Microeukaryote",
    "Class",  "GS55",                               "Microeukaryote",
    "Class",  "GS136",                              "Microeukaryote",
    "Class",  "GS141",                              "Microeukaryote",
    "Class",  "RAD-A",                              "Microeukaryote",
    "Class",  "RAD-B",                              "Microeukaryote",
    
    # --- Other ---
    "Class",  "Insecta",                            "Other",
    "Class",  "Arachnida",                          "Other",
    "Class",  "Collembola",                         "Other",
    "Class",  "Amphibia",                           "Other",
    "Class",  "Myxozoa",                            "Other",
    "Class",  "Enteropneusta",                      "Other",
    "Class",  "Pinopsida",                          "Other",
    "Class",  "Angiospermae",                       "Other",
    "Class",  "Pilidiophora",                       "Mobile invertebrate",
    "Class",  "Palaeonemertea",                     "Mobile invertebrate",
    "Class",  "Enopla",                             "Mobile invertebrate",
    
    # Bacteria (from 12S off-target)
    "Class",  "Gammaproteobacteria",                "Other",
    "Class",  "Alphaproteobacteria",                "Other",
    "Class",  "Betaproteobacteria",                 "Other",
    "Class",  "Flavobacteriia",                     "Other",
    
    # =====================================================================
    # PHYLUM-LEVEL FALLBACKS
    # =====================================================================
    # Catch anything not matched at Class or Order level
    
    "Phylum", "Chordata",                           "Fish",
    "Phylum", "Bacillariophyta",                    "Phytoplankton",
    "Phylum", "Dinoflagellata",                     "Phytoplankton",
    "Phylum", "cf.Dinoflagellata",                  "Phytoplankton",
    "Phylum", "Cryptophyta",                        "Phytoplankton",
    "Phylum", "Chlorophyta",                        "Phytoplankton",
    "Phylum", "Haptophyta",                         "Phytoplankton",
    "Phylum", "Eustigmatophyta",                    "Phytoplankton",
    "Phylum", "Pelagophyta",                        "Phytoplankton",
    "Phylum", "Dictyochophyta",                     "Phytoplankton",
    "Phylum", "Raphidophyta",                       "Phytoplankton",
    "Phylum", "Chrysophyta",                        "Phytoplankton",
    "Phylum", "Synchromophyta",                     "Phytoplankton",
    "Phylum", "Prasinodermophyta",                  "Phytoplankton",
    "Phylum", "Xanthophyta",                        "Phytoplankton",
    "Phylum", "Rhodophyta",                         "Macroalgae",
    "Phylum", "Phaeophyta",                         "Macroalgae",
    "Phylum", "Bryozoa",                            "Sessile invertebrate",
    "Phylum", "Porifera",                           "Sessile invertebrate",
    "Phylum", "Entoprocta",                         "Sessile invertebrate",
    "Phylum", "Phoronida",                          "Sessile invertebrate",
    "Phylum", "Cnidaria",                           "Sessile invertebrate",
    "Phylum", "Echinodermata",                      "Mobile invertebrate",
    "Phylum", "Mollusca",                           "Mobile invertebrate",
    "Phylum", "Annelida",                           "Sessile invertebrate",
    "Phylum", "Nemertea",                           "Mobile invertebrate",
    "Phylum", "Nematoda",                           "Mobile invertebrate",
    "Phylum", "Xenacoelomorpha",                    "Mobile invertebrate",
    "Phylum", "Gastrotricha",                       "Mobile invertebrate",
    "Phylum", "Platyhelminthes",                    "Mobile invertebrate",
    "Phylum", "Chaetognatha",                       "Zooplankton",
    "Phylum", "Ctenophora",                         "Zooplankton",
    "Phylum", "Rotifera",                           "Zooplankton",
    "Phylum", "Arthropoda",                         "Mobile invertebrate",
    "Phylum", "Cercozoa",                           "Microeukaryote",
    "Phylum", "Ciliophora",                         "Microeukaryote",
    "Phylum", "Radiolaria",                         "Microeukaryote",
    "Phylum", "Oomycota",                           "Microeukaryote",
    "Phylum", "Chytridiomycota",                    "Microeukaryote",
    "Phylum", "Ascomycota",                         "Microeukaryote",
    "Phylum", "Basidiomycota",                      "Microeukaryote",
    "Phylum", "Zoopagomycota",                      "Microeukaryote",
    "Phylum", "Blastocladiomycota",                 "Microeukaryote",
    "Phylum", "Rozellomycota",                      "Microeukaryote",
    "Phylum", "Aphelidiomycota",                    "Microeukaryote",
    "Phylum", "Labyrinthulidia",                    "Microeukaryote",
    "Phylum", "Picomonada",                         "Phytoplankton",
    "Phylum", "Apicomplexa",                        "Microeukaryote",
    "Phylum", "Perkinsia",                          "Microeukaryote",
    "Phylum", "Choanoflagellida",                   "Microeukaryote",
    "Phylum", "Centrohelida",                       "Microeukaryote",
    "Phylum", "Telonemidia",                        "Microeukaryote",
    "Phylum", "Katablepharidophyta",                "Microeukaryote",
    "Phylum", "Discosea",                           "Microeukaryote",
    "Phylum", "Apusomonada",                        "Microeukaryote",
    "Phylum", "Ichthyosporidia",                    "Microeukaryote",
    "Phylum", "Bicosoecida",                        "Microeukaryote",
    "Phylum", "Endomyxa",                           "Microeukaryote",
    "Phylum", "Nuclearida",                         "Microeukaryote",
    "Phylum", "Solenicolida",                       "Microeukaryote",
    "Phylum", "Pirsonida",                          "Microeukaryote",
    "Phylum", "Chromerida",                         "Microeukaryote",
    "Phylum", "Fornicatia",                         "Microeukaryote",
    "Phylum", "Gromida",                            "Microeukaryote",
    "Phylum", "Hemichordata",                       "Mobile invertebrate",
    "Phylum", "Orthonectida",                       "Microeukaryote",
    "Phylum", "Pseudomonadota",                     "Other",
    "Phylum", "Bacteroidota",                       "Other",
    "Phylum", "Tracheophyta",                       "Other",
    
    # MAST phyla (all microeukaryote)
    "Phylum", "MAST-1",                             "Microeukaryote",
    "Phylum", "MAST-2",                             "Microeukaryote",
    "Phylum", "MAST-4",                             "Microeukaryote",
    "Phylum", "MAST-6",                             "Microeukaryote",
    "Phylum", "MAST-7",                             "Microeukaryote",
    "Phylum", "MAST-9",                             "Microeukaryote",
    "Phylum", "MAST-10",                            "Microeukaryote",
    "Phylum", "MAST-11",                            "Microeukaryote",
    "Phylum", "MAST-12",                            "Microeukaryote",
    "Phylum", "MOCH-1",                             "Microeukaryote",
    "Phylum", "MOCH-2",                             "Microeukaryote",
    "Phylum", "MOCH-5",                             "Microeukaryote",
    "Phylum", "Straminipila_phy_8",                 "Microeukaryote",
    "Phylum", "Rhizaria_phy_incertae_sedis",        "Microeukaryote"
  )
}


#' Assign lifestyle categories to site-level detection data
#'
#' Uses a hierarchical lookup (Order → Class → Phylum) to assign one of
#' eight lifestyle categories to each detection. Order-level matches take
#' priority, allowing mixed-lifestyle classes to be split correctly.
#'
#' @param site_level Data frame with columns Phylum, Class, Order.
#' @param lookup Lookup table (output of default_lifestyle_lookup or custom).
#'   Must have columns: level, taxon, lifestyle.
#' @param na_lifestyle Label for detections with no Phylum/Class/Order.
#'   Default "Unclassified".
#' @return Input data frame with `lifestyle` column added (factor).
#'
#' @examples
#' # Default:
#' combined_all$site_level <- assign_lifestyle(combined_all$site_level)
#'
#' # View assignments:
#' table(combined_all$site_level$lifestyle, combined_all$site_level$Marker)
#'
#' # Check what got "Other" or "Unclassified":
#' check_unmatched_lifestyle(combined_all$site_level)
assign_lifestyle <- function(
    site_level,
    lookup = default_lifestyle_lookup(),
    na_lifestyle = "Unclassified"
) {
  
  lifestyle_levels <- c("Fish", "Sessile invertebrate", "Mobile invertebrate",
                         "Phytoplankton", "Zooplankton", "Macroalgae",
                         "Microeukaryote", "Other", na_lifestyle)
  
  # Build lookup hashes for fast matching
  order_map  <- lookup %>% filter(level == "Order") %>% 
    {setNames(.$lifestyle, .$taxon)}
  class_map  <- lookup %>% filter(level == "Class") %>% 
    {setNames(.$lifestyle, .$taxon)}
  phylum_map <- lookup %>% filter(level == "Phylum") %>% 
    {setNames(.$lifestyle, .$taxon)}
  
  # Assign hierarchically: Order > Class > Phylum
  result <- site_level %>%
    mutate(
      lifestyle = case_when(
        !is.na(Order) & Order %in% names(order_map)   ~ order_map[Order],
        !is.na(Class) & Class %in% names(class_map)   ~ class_map[Class],
        !is.na(Phylum) & Phylum %in% names(phylum_map) ~ phylum_map[Phylum],
        is.na(Phylum) & is.na(Class)                   ~ na_lifestyle,
        TRUE                                           ~ "Other"
      ),
      lifestyle = factor(lifestyle, levels = lifestyle_levels)
    )
  
  # Report
  cat("Lifestyle assignments:\n")
  cat("  Total:", nrow(result), "\n")
  counts <- table(result$lifestyle)
  for (nm in names(counts)) {
    cat("  ", nm, ":", counts[nm], 
        paste0("(", round(counts[nm] / nrow(result) * 100, 1), "%)"), "\n")
  }
  
  result
}


#' Check for unmatched or unexpected lifestyle assignments
#'
#' Reports taxa assigned to "Other" or "Unclassified" so you can extend
#' the lookup table.
#'
#' @param site_level Data frame with lifestyle column (from assign_lifestyle)
#' @return Tibble of unmatched Phylum × Class × Order combinations
check_unmatched_lifestyle <- function(site_level) {
  
  unmatched <- site_level %>%
    filter(lifestyle %in% c("Other", "Unclassified")) %>%
    group_by(Phylum, Class, Order) %>%
    summarise(n = n(), .groups = "drop") %>%
    arrange(desc(n))
  
  if (nrow(unmatched) == 0) {
    cat("All taxa matched a lifestyle category.\n")
  } else {
    cat("Taxa assigned 'Other' or 'Unclassified':\n")
    print(unmatched, n = Inf)
  }
  
  invisible(unmatched)
}


# =============================================================================
# LOADED
# =============================================================================

cat("FUNC_assign_lifestyle.R loaded\n")
cat("Functions:\n")
cat("  default_lifestyle_lookup()  - View/edit the mapping table\n")
cat("  assign_lifestyle(site_level, lookup)  - Add lifestyle column\n")
cat("  check_unmatched_lifestyle(site_level) - Find unmapped taxa\n")
