source("Scripts/Functions/FUNC_occupancy_power_analysis_FIXED.R")

# =============================================================================
# STEP 1: Subset to Citizen Science samples
# =============================================================================

ps_cs_12S <- subset_samples(meta_12S$phyloseq, 
                            Sampling.area.Project == "BGE Marine Invasive Species Citizen Science")
ps_cs_18S <- subset_samples(meta_18S$phyloseq, 
                            Sampling.area.Project == "BGE Marine Invasive Species Citizen Science")
ps_cs_COI <- subset_samples(meta_COI$phyloseq, 
                            Sampling.area.Project == "BGE Marine Invasive Species Citizen Science")

# Store pre-filtering totals
pre_filter_seqs <- c(sum(sample_sums(ps_cs_12S)), sum(sample_sums(ps_cs_18S)), sum(sample_sums(ps_cs_COI)))
pre_filter_otus <- c(sum(taxa_sums(ps_cs_12S) > 0), sum(taxa_sums(ps_cs_18S) > 0), sum(taxa_sums(ps_cs_COI) > 0))

# =============================================================================
# STEP 2: Remove mammal/bird OTUs at phyloseq level - but leave cteaceans amd pinnipeds
# =============================================================================

CETACEAN_FAMILIES <- c("Delphinidae", "Phocoenidae", "Balaenidae", "Balaenopteridae",
                       "Ziphiidae", "Physeteridae", "Kogiidae", "Monodontidae",
                       "Eschrichtiidae", "Platanistidae", "Iniidae", "Pontoporiidae")

PINNIPED_FAMILIES <- c("Phocidae", "Otariidae", "Odobenidae")

MARINE_MAMMAL_FAMILIES <- c(CETACEAN_FAMILIES, PINNIPED_FAMILIES)

filter_mambird_ps <- function(ps, marker) {
  tax <- as.data.frame(tax_table(ps))
  
  is_mammal_or_bird <- !is.na(tax$Class) & tax$Class %in% c("Mammalia", "Aves")
  is_marine_mammal  <- !is.na(tax$Family) & tax$Family %in% MARINE_MAMMAL_FAMILIES
  
  keep <- rownames(tax)[!is_mammal_or_bird | is_marine_mammal]
  
  seqs_before <- sum(sample_sums(ps))
  ps_filt <- prune_taxa(keep, ps)
  ps_filt <- prune_taxa(taxa_sums(ps_filt) > 0, ps_filt)
  seqs_after <- sum(sample_sums(ps_filt))
  
  n_kept_mm <- sum(is_mammal_or_bird & is_marine_mammal)
  
  cat(marker, ": removed", ntaxa(ps) - ntaxa(ps_filt), "mammal/bird OTUs (",
      ntaxa(ps), "->", ntaxa(ps_filt), ") |",
      formatC(seqs_before - seqs_after, format = "d", big.mark = ","), "seqs removed (",
      formatC(seqs_before, format = "d", big.mark = ","), "->",
      formatC(seqs_after, format = "d", big.mark = ","), ") |",
      n_kept_mm, "cetacean/pinniped OTUs retained\n")
  ps_filt
}

ps_cs_12S <- filter_mambird_ps(ps_cs_12S, "12S")
ps_cs_18S <- filter_mambird_ps(ps_cs_18S, "18S")
ps_cs_COI <- filter_mambird_ps(ps_cs_COI, "COI")

# =============================================================================
# STEP 3: Remove non-marine taxa using WoRMS environment flags
# =============================================================================

library(worrms)

get_species_from_ps <- function(ps) {
  tax <- as.data.frame(tax_table(ps))
  tax$OTU <- rownames(tax)
  tax %>%
    filter(!is.na(Species) & Species != "") %>%
    filter(!grepl("uncultured", Species, ignore.case = TRUE)) %>%
    filter(!is.na(Genus) & Genus != "") %>%
    distinct(Species, .keep_all = TRUE)
}

query_worms_environment <- function(species_list, delay = 0.3) {
  results <- list()
  n <- length(species_list)
  
  for (i in seq_along(species_list)) {
    sp <- species_list[i]
    cat(i, "/", n, ":", sp, "...")
    
    rec <- tryCatch(
      wm_records_name(sp, fuzzy = FALSE),
      error = function(e) NULL
    )
    
    if (is.null(rec) || nrow(rec) == 0) {
      cat(" not found\n")
      results[[i]] <- tibble(
        Species = sp, AphiaID = NA,
        isMarine = NA, isBrackish = NA, 
        isFreshwater = NA, isTerrestrial = NA
      )
    } else {
      accepted <- rec %>% filter(status == "accepted")
      if (nrow(accepted) == 0) accepted <- rec[1, , drop = FALSE] else accepted <- accepted[1, , drop = FALSE]
      
      cat(" AphiaID:", accepted$AphiaID, 
          " M:", accepted$isMarine, 
          " B:", accepted$isBrackish,
          " F:", accepted$isFreshwater, "\n")
      
      results[[i]] <- tibble(
        Species = sp, 
        AphiaID = accepted$AphiaID,
        isMarine = as.logical(accepted$isMarine),
        isBrackish = as.logical(accepted$isBrackish),
        isFreshwater = as.logical(accepted$isFreshwater),
        isTerrestrial = as.logical(accepted$isTerrestrial)
      )
    }
    
    Sys.sleep(delay)
  }
  
  bind_rows(results)
}

classify_non_marine <- function(env_df) {
  env_df %>%
    mutate(
      is_non_marine = case_when(
        is.na(isMarine) & is.na(isBrackish) & is.na(isFreshwater) ~ NA,
        isMarine == FALSE ~ TRUE,
        is.na(isMarine) & isFreshwater == TRUE & (isBrackish == FALSE | is.na(isBrackish)) ~ TRUE,
        TRUE ~ FALSE
      )
    )
}

filter_non_marine_ps <- function(ps, marker, env_data) {
  tax <- as.data.frame(tax_table(ps))
  tax$OTU <- rownames(tax)
  
  tax_env <- tax %>%
    left_join(env_data %>% select(Species, is_non_marine), by = "Species")
  
  nm_otus <- tax_env %>%
    filter(is_non_marine == TRUE)
  
  if (nrow(nm_otus) > 0) {
    cat("\n", marker, ": Removing", nrow(nm_otus), "non-marine OTUs:\n")
    nm_otus %>%
      select(OTU, Species, Genus, Family, Class, Phylum) %>%
      as_tibble() %>%
      print(n = Inf)
    
    keep <- rownames(tax)[!rownames(tax) %in% nm_otus$OTU]
    ps_filt <- prune_taxa(keep, ps)
    ps_filt <- prune_taxa(taxa_sums(ps_filt) > 0, ps_filt)
    
    cat("  Before:", ntaxa(ps), "OTUs | After:", ntaxa(ps_filt), "OTUs\n")
  } else {
    cat(marker, ": No non-marine OTUs found\n")
    ps_filt <- ps
  }
  
  return(list(ps = ps_filt, removed = nm_otus))
}

# Collect all species across all three markers
all_species <- unique(c(
  get_species_from_ps(ps_cs_12S)$Species,
  get_species_from_ps(ps_cs_18S)$Species,
  get_species_from_ps(ps_cs_COI)$Species
))

cat("Querying WoRMS for", length(all_species), "species...\n")
worms_env <- query_worms_environment(all_species, delay = 0.2)
worms_env <- classify_non_marine(worms_env)

# Manual fix: Corbicula fluminea is freshwater but not in WoRMS
worms_env <- worms_env %>%
  mutate(
    isMarine = ifelse(Species == "Corbicula fluminea", FALSE, isMarine),
    isFreshwater = ifelse(Species == "Corbicula fluminea", TRUE, isFreshwater),
    is_non_marine = ifelse(Species == "Corbicula fluminea", TRUE, is_non_marine)
  )

# Summary
cat("\n=== WoRMS Environment Summary ===\n")
cat("Total species queried:", nrow(worms_env), "\n")
cat("Marine:", sum(worms_env$isMarine == TRUE, na.rm = TRUE), "\n")
cat("Brackish:", sum(worms_env$isBrackish == TRUE, na.rm = TRUE), "\n")
cat("Freshwater:", sum(worms_env$isFreshwater == TRUE, na.rm = TRUE), "\n")
cat("Non-marine:", sum(worms_env$is_non_marine == TRUE, na.rm = TRUE), "\n")

cat("\n=== Non-marine species to be removed ===\n")
worms_env %>% filter(is_non_marine == TRUE) %>% print(n = Inf)

nm_species <- worms_env %>% filter(is_non_marine == TRUE) %>% pull(Species)

get_nm_sites <- function(ps, marker, nm_species) {
  tax <- as.data.frame(tax_table(ps))
  tax$OTU <- rownames(tax)
  
  nm_otus <- tax %>% filter(Species %in% nm_species)
  
  if (nrow(nm_otus) == 0) {
    return(tibble(Species = character(), Marker = character(), Sites = character()))
  }
  
  otu_mat <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) otu_mat <- t(otu_mat)
  
  sdata <- as(sample_data(ps), "data.frame")
  
  result <- lapply(seq_len(nrow(nm_otus)), function(i) {
    otu_id <- nm_otus$OTU[i]
    sp <- nm_otus$Species[i]
    present_samples <- names(which(otu_mat[otu_id, ] > 0))
    if (length(present_samples) == 0) return(NULL)
    sites <- unique(sdata[present_samples, "Sampling.event.ID"])
    tibble(Species = sp, Marker = marker, Sites = paste(sites, collapse = ", "))
  })
  
  bind_rows(result)
}

# Get site/marker info for non-marine species from pre-filtered phyloseq objects
ps_orig_12S <- subset_samples(meta_12S$phyloseq, 
                              Sampling.area.Project == "BGE Marine Invasive Species Citizen Science")
ps_orig_18S <- subset_samples(meta_18S$phyloseq, 
                              Sampling.area.Project == "BGE Marine Invasive Species Citizen Science")
ps_orig_COI <- subset_samples(meta_COI$phyloseq, 
                              Sampling.area.Project == "BGE Marine Invasive Species Citizen Science")

nm_sites <- bind_rows(
  get_nm_sites(ps_orig_12S, "12S", nm_species),
  get_nm_sites(ps_orig_18S, "18S", nm_species),
  get_nm_sites(ps_orig_COI, "COI", nm_species)
) %>%
  group_by(Species) %>%
  summarise(
    Markers = paste(unique(Marker), collapse = ", "),
    Sites = paste(unique(Sites), collapse = ", "),
    .groups = "drop"
  )

# Get taxonomy from original phyloseq objects
nm_tax <- bind_rows(
  as.data.frame(tax_table(ps_orig_12S)) %>% mutate(OTU = rownames(.)),
  as.data.frame(tax_table(ps_orig_18S)) %>% mutate(OTU = rownames(.)),
  as.data.frame(tax_table(ps_orig_COI)) %>% mutate(OTU = rownames(.))
) %>%
  filter(Species %in% nm_species) %>%
  distinct(Species, .keep_all = TRUE) %>%
  select(Species, Phylum, Class, Order, Family, Genus)

worms_env %>%
  filter(is_non_marine == TRUE) %>%
  left_join(nm_tax, by = "Species") %>%
  left_join(nm_sites, by = "Species") %>%
  select(Species, Phylum, Class, Order, Family, Genus, AphiaID, isMarine, isBrackish, isFreshwater, Markers, Sites) %>%
  as_tibble() %>%
  print(n = Inf, width = Inf)

write.csv(worms_env, "Processed_data/worms_environment_flags.csv", row.names = FALSE)

# Filter non-marine from each phyloseq
fw_12S <- filter_non_marine_ps(ps_cs_12S, "12S", worms_env)
fw_18S <- filter_non_marine_ps(ps_cs_18S, "18S", worms_env)
fw_COI <- filter_non_marine_ps(ps_cs_COI, "COI", worms_env)

# Replace phyloseq objects
ps_cs_12S <- fw_12S$ps
ps_cs_18S <- fw_18S$ps
ps_cs_COI <- fw_COI$ps

# Log removed taxa
non_marine_removed <- bind_rows(
  if (nrow(fw_12S$removed) > 0) fw_12S$removed %>% mutate(Marker = "12S") else NULL,
  if (nrow(fw_18S$removed) > 0) fw_18S$removed %>% mutate(Marker = "18S") else NULL,
  if (nrow(fw_COI$removed) > 0) fw_COI$removed %>% mutate(Marker = "COI") else NULL
)

cat("\n=== Total non-marine OTUs removed ===\n")
cat("12S:", nrow(fw_12S$removed), "\n")
cat("18S:", nrow(fw_18S$removed), "\n")
cat("COI:", nrow(fw_COI$removed), "\n")

cat("\n=== Unique non-marine species removed ===\n")
if (!is.null(non_marine_removed) && nrow(non_marine_removed) > 0) {
  non_marine_removed %>%
    distinct(Species, Phylum, Class, Marker) %>%
    arrange(Phylum, Class, Species) %>%
    as_tibble() %>%
    print(n = Inf)
}

write.csv(non_marine_removed, "Processed_data/non_marine_taxa_removed.csv", row.names = FALSE)

# =============================================================================
# STEP 4: Filter PCR replicates with <1000 reads
# =============================================================================

pcr_read_threshold <- 1000

filter_low_read_pcrs <- function(ps, marker, threshold) {
  reads <- sample_sums(ps)
  keep <- names(reads[reads >= threshold])
  removed <- length(reads) - length(keep)
  cat(marker, ": removed", removed, "of", length(reads), "PCR replicates (<", threshold, "reads)\n")
  
  ps_filt <- prune_samples(keep, ps)
  ps_filt <- prune_taxa(taxa_sums(ps_filt) > 0, ps_filt)
  
  cat("  Samples:", nsamples(ps_filt), " OTUs:", ntaxa(ps_filt), "\n")
  ps_filt
}

ps_cs_12S <- filter_low_read_pcrs(ps_cs_12S, "12S", pcr_read_threshold)
ps_cs_18S <- filter_low_read_pcrs(ps_cs_18S, "18S", pcr_read_threshold)
ps_cs_COI <- filter_low_read_pcrs(ps_cs_COI, "COI", pcr_read_threshold)

# =============================================================================
# STEP 5: Create COLLAPSED versions for MIS detection analysis
# =============================================================================
# Problem: Some MIS have multiple OTUs (e.g., Mnemiopsis leidyi has 2 COI OTUs)
# If OTU1 is detected in biosample A and OTU2 in biosamples B+C, the SPECIES
# was actually detected in ALL 3 biosamples (theta should = 1.0, not 0.33 each)
#
# Solution: Collapse multi-OTU species by summing reads per PCR sample
# - Keep ORIGINAL ps objects for summary statistics (OTU counts, Table 1)
# - Use COLLAPSED ps objects for detection analysis (theta, Sites, etc.)
# =============================================================================

collapse_multi_otu_species <- function(ps, marker) {
  
  # Get taxonomy
  tax <- as.data.frame(tax_table(ps))
  tax$OTU <- rownames(tax)
  
  # Find species with multiple OTUs (exclude uncultured and species without genus)
  multi_otu <- tax %>%
    filter(!is.na(Species) & Species != "") %>%
    filter(!grepl("uncultured", Species, ignore.case = TRUE)) %>%
    filter(!is.na(Genus) & Genus != "") %>%
    group_by(Species) %>%
    filter(n() > 1) %>%
    ungroup()
  
  if (nrow(multi_otu) == 0) {
    cat(marker, ": No multi-OTU species to collapse\n")
    return(list(ps = ps, collapsed_species = NULL))
  }
  
  multi_otu_species <- unique(multi_otu$Species)
  cat(marker, ": Collapsing", length(multi_otu_species), "species with multiple OTUs:\n")
  
  collapse_info <- list()
  for (sp in multi_otu_species) {
    otus <- multi_otu %>% filter(Species == sp) %>% pull(OTU)
    cat("  ", sp, ":", length(otus), "OTUs (", paste(otus, collapse = ", "), ")\n")
    collapse_info[[sp]] <- otus
  }
  
  # Get OTU table
  otu_mat <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) otu_mat <- t(otu_mat)
  
  # Get reference sequences if available
  has_refseq <- !is.null(refseq(ps, errorIfNULL = FALSE))
  if (has_refseq) {
    seqs <- refseq(ps)
  }
  
  # Process each multi-OTU species
  otus_to_remove <- c()
  
  for (sp in multi_otu_species) {
    sp_otus <- multi_otu %>% filter(Species == sp) %>% pull(OTU)
    
    # Sum reads across OTUs for this species (per sample)
    sp_reads <- colSums(otu_mat[sp_otus, , drop = FALSE])
    
    # Keep the first OTU as the "representative" and add summed reads
    keep_otu <- sp_otus[1]
    remove_otus <- sp_otus[-1]
    
    # Replace the kept OTU's reads with the summed reads
    otu_mat[keep_otu, ] <- sp_reads
    
    # Track OTUs to remove
    otus_to_remove <- c(otus_to_remove, remove_otus)
    
    cat("    -> Collapsed into", keep_otu, "(summed reads across all OTUs)\n")
  }
  
  # Remove the extra OTUs
  keep_otus <- rownames(otu_mat)[!rownames(otu_mat) %in% otus_to_remove]
  otu_mat <- otu_mat[keep_otus, , drop = FALSE]
  
  # Rebuild phyloseq object
  new_otu <- otu_table(otu_mat, taxa_are_rows = TRUE)
  new_tax <- tax_table(as.matrix(tax[keep_otus, colnames(tax) != "OTU"]))
  
  if (has_refseq) {
    new_seqs <- seqs[keep_otus]
    ps_collapsed <- phyloseq(new_otu, sample_data(ps), new_tax, new_seqs)
  } else {
    ps_collapsed <- phyloseq(new_otu, sample_data(ps), new_tax)
  }
  
  cat(marker, ": OTUs after collapsing:", ntaxa(ps), "->", ntaxa(ps_collapsed), "\n\n")
  
  return(list(
    ps = ps_collapsed, 
    collapsed_species = tibble(
      Marker = marker,
      Species = names(collapse_info),
      n_OTUs = sapply(collapse_info, length),
      OTUs = sapply(collapse_info, paste, collapse = ", ")
    )
  ))
}

# Store ORIGINAL phyloseq objects (for summary statistics - Table 1, OTU counts)
ps_cs_12S_original <- ps_cs_12S
ps_cs_18S_original <- ps_cs_18S
ps_cs_COI_original <- ps_cs_COI

# Create COLLAPSED versions (for detection analysis - theta, Sites, etc.)
collapse_12S <- collapse_multi_otu_species(ps_cs_12S, "12S")
collapse_18S <- collapse_multi_otu_species(ps_cs_18S, "18S")
collapse_COI <- collapse_multi_otu_species(ps_cs_COI, "COI")

ps_cs_12S_collapsed <- collapse_12S$ps
ps_cs_18S_collapsed <- collapse_18S$ps
ps_cs_COI_collapsed <- collapse_COI$ps

# Record which species were collapsed (for documentation)
collapsed_species_info <- bind_rows(
  collapse_12S$collapsed_species,
  collapse_18S$collapsed_species,
  collapse_COI$collapsed_species
)

if (!is.null(collapsed_species_info) && nrow(collapsed_species_info) > 0) {
  cat("\n=== Species with multiple OTUs (collapsed for detection analysis) ===\n")
  print(collapsed_species_info, n = Inf)
  write.csv(collapsed_species_info, "Processed_data/collapsed_multi_otu_species.csv", row.names = FALSE)
}


cat("Post-collapse OTU totals:",
    ntaxa(ps_cs_12S_collapsed), "(12S),",
    ntaxa(ps_cs_18S_collapsed), "(18S),",
    ntaxa(ps_cs_COI_collapsed), "(COI)\n")

# =============================================================================
# STEP 6: Run detection analysis on COLLAPSED data
# =============================================================================

det_12S <- run_detection_analysis(ps_cs_12S_collapsed, "12S", site_var = "Sampling.event.ID")
combined_12S <- combine_detection_invasive(
  detection_results = det_12S,
  invasive_status = inv_12S$invasive_status,
  verified_status = verified_12S
)

det_18S <- run_detection_analysis(ps_cs_18S_collapsed, "18S", site_var = "Sampling.event.ID")
combined_18S <- combine_detection_invasive(
  detection_results = det_18S,
  invasive_status = inv_18S$invasive_status,
  verified_status = verified_18S
)

det_COI <- run_detection_analysis(ps_cs_COI_collapsed, "COI", site_var = "Sampling.event.ID")
combined_COI <- combine_detection_invasive(
  detection_results = det_COI,
  invasive_status = inv_COI$invasive_status,
  verified_status = verified_COI
)

# =============================================================================
# STEP 7: Library sizes per site (use ORIGINAL for read/OTU stats)
# =============================================================================

get_site_total_reads <- function(ps, marker) {
  sdata <- as(sample_data(ps), "data.frame")
  sdata$sample_id <- rownames(sdata)
  
  sample_reads <- sample_sums(ps)
  sdata$sample_reads <- sample_reads[rownames(sdata)]
  
  otu_mat <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) otu_mat <- t(otu_mat)
  
  site_stats <- sdata %>%
    group_by(Location, Sampling.event.ID) %>%
    summarise(
      total_reads_at_site = sum(sample_reads),
      n_biosamples = n_distinct(Name),
      n_pcr_replicates = n(),
      sample_ids = list(sample_id),
      .groups = "drop"
    )
  
  site_stats$total_otus_at_site <- sapply(site_stats$sample_ids, function(ids) {
    sum(rowSums(otu_mat[, ids, drop = FALSE]) > 0)
  })
  
  site_stats <- site_stats %>% select(-sample_ids)
  names(site_stats)[names(site_stats) == "Location"] <- "Country"
  names(site_stats)[names(site_stats) == "Sampling.event.ID"] <- "Site"
  site_stats$Marker <- marker
  site_stats
}

# Use ORIGINAL for read/OTU statistics
site_reads_12S <- get_site_total_reads(ps_cs_12S_original, "12S")
site_reads_18S <- get_site_total_reads(ps_cs_18S_original, "18S")
site_reads_COI <- get_site_total_reads(ps_cs_COI_original, "COI")

site_total_reads <- bind_rows(site_reads_12S, site_reads_18S, site_reads_COI)

print(site_total_reads, n = Inf)

# =============================================================================
# STEP 8: Combine all three markers (uses COLLAPSED detection data)
# =============================================================================

combined_all <- combine_markers(
  `12S` = combined_12S,
  `18S` = combined_18S,
  COI = combined_COI
)

# --- Pre-correction review: all species flagged as invasive before manual fixes ---
# Review this table to identify any species that are native to Europe / incorrectly flagged.
# Add any such species to the case_when corrections block below.
combined_all$site_level %>%
  filter(is_invasive) %>%
  group_by(Species) %>%
  summarise(
    Markers = paste(sort(unique(Marker)), collapse = ", "),
    Sites   = paste(sort(unique(Site)), collapse = ", "),
    .groups = "drop"
  ) %>%
  left_join(
    bind_rows(
      as.data.frame(tax_table(ps_cs_12S_original)) %>% mutate(OTU = rownames(.)),
      as.data.frame(tax_table(ps_cs_18S_original)) %>% mutate(OTU = rownames(.)),
      as.data.frame(tax_table(ps_cs_COI_original)) %>% mutate(OTU = rownames(.))
    ) %>% distinct(Species, .keep_all = TRUE) %>%
      select(Species, Phylum, Class, Order, Family),
    by = "Species"
  ) %>%
  select(Phylum, Class, Order, Family, Species, Markers, Sites) %>%
  arrange(Phylum, Class, Order, Species) %>%
  print(n = Inf, width = Inf)

# Fixing remaining mis-identifications of invasives (native European species)
combined_all$site_level <- combined_all$site_level %>%
  mutate(
    is_invasive = case_when(
      Species == "Alitta succinea" ~ FALSE,
      Species == "Amphibalanus amphitrite" ~ FALSE,
      Species == "Barentsia benedeni" ~ FALSE,
      Species == "Carcinus maenas" ~ FALSE,
      Species == "Clava multicornis" ~ FALSE,
      Species == "Gasterosteus aculeatus" ~ FALSE,
      Species == "Littorina littorea" ~ FALSE,
      Species == "Membranipora membranacea" ~ FALSE,
      Species == "Mya arenaria" ~ FALSE,
      Species == "Salmo salar" ~ FALSE,
      Species == "Sparus aurata" ~ FALSE,
      is.na(is_invasive) ~ FALSE,
      TRUE ~ is_invasive
    )
  )

# Verify
sum(is.na(combined_all$site_level$is_invasive))

combined_all$site_level %>%
  filter(is_invasive) %>%
  distinct(Species) %>%
  arrange(Species)

# =============================================================================
# STEP 8b: Assemble combined power analysis objects for visualization
# =============================================================================
# run_detection_analysis() now returns per-marker bootstrap power results.
# Combine across markers and join Species + is_invasive from combined_all
# so the visualization script (05_visualizations_detection_power_v2.R) can use them.

# --- Combine empirical_3x3 across markers ---
empirical_3x3 <- bind_rows(
  det_12S$empirical_3x3 %>% mutate(Marker = "12S"),
  det_18S$empirical_3x3 %>% mutate(Marker = "18S"),
  det_COI$empirical_3x3 %>% mutate(Marker = "COI")
)

# Join taxonomy and is_invasive from the finalised combined_all
empirical_3x3 <- empirical_3x3 %>%
  left_join(
    combined_all$site_level %>%
      select(OTU, Site, Marker, Species, Genus, Family, Order, Class, Phylum, is_invasive) %>%
      distinct(),
    by = c("OTU", "Site", "Marker")
  ) %>%
  rename(p = mean_p_pcr)

cat("empirical_3x3 assembled:", nrow(empirical_3x3), "OTU × site records\n")
cat("  with is_invasive:", sum(empirical_3x3$is_invasive, na.rm = TRUE), "invasive detections\n")

# --- Combine design_performance across markers ---
design_performance <- bind_rows(
  det_12S$design_performance %>% mutate(Marker = "12S"),
  det_18S$design_performance %>% mutate(Marker = "18S"),
  det_COI$design_performance %>% mutate(Marker = "COI")
) %>%
  rename(n_biosamples = n_bio_design, n_pcr = n_pcr_design)

cat("design_performance assembled:", nrow(design_performance), "design × marker rows\n")

# --- Save for visualization script ---
save(empirical_3x3, design_performance,
     file = "Processed_data/power_analysis_results.RData")
cat("Saved: Processed_data/power_analysis_results.RData\n")

# --- LIST SITES ---
combined_all$site_level %>%
  distinct(Country, Site) %>%
  arrange(Country, Site) %>%
  print(n = 60)

# =============================================================================
# STEP 9: Visualizations
# =============================================================================

plot_site_by_marker(combined_all, "France_1", 
                    x_metric = "total_reads", y_metric = "p_empirical")

plot_site_by_marker(combined_all, "Denmark_1",
                    x_metric = "theta", y_metric = "confidence_score")

## Plots faceted by theta
source("Scripts/Functions/FUNC_plot_site_facet_theta.R")

countries <- c("Denmark", "Finland", "France", "Greece", "Italy", 
               "Netherlands", "Norway", "Poland", "Portugal", "Spain", "UK")

for (country in countries) {
  plot_country_facet_theta(combined_all, country, output_dir = "Figures/")
}

# =============================================================================
# STEP 10: Invasive species summary table (uses COLLAPSED for detection metrics)
# =============================================================================

invasive_summary <- combined_all$site_level %>%
  filter(is_invasive) %>%
  left_join(
    site_total_reads %>% select(Country, Site, Marker, total_reads_at_site),
    by = c("Country", "Site", "Marker")
  ) %>%
  mutate(
    prop_reads = total_reads / total_reads_at_site
  ) %>%
  select(
    Country, Site, Species, Marker,
    n_biosamples_positive, total_biosamples_at_site,
    n_pcr_positive, total_pcr_at_site,
    theta, p_empirical, mean_p_pcr, confidence_score,
    total_reads, total_reads_at_site, prop_reads,
    pcr_distribution, pcr_spread_ratio
  ) %>%
  mutate(
    Reliable = case_when(
      n_biosamples_positive >= 2 & p_empirical >= 0.5 ~ "Reliable",
      n_biosamples_positive == 3 & total_biosamples_at_site == 3 & p_empirical < 0.5 ~ "Marginal",
      n_biosamples_positive == 2 & total_biosamples_at_site == 3 & pcr_spread_ratio > 0.5 ~ "Marginal",
      TRUE ~ "Unreliable"
    )
  ) %>%
  arrange(Country, Site, Reliable, Species)

table(invasive_summary$Reliable)
summary(invasive_summary$prop_reads)
View(invasive_summary)

write.csv(invasive_summary, "Processed_data/invasive_species_detectionsSUMMARY.csv", row.names = FALSE)

####### Detection summary stats
combined_all$site_level %>%
  filter(is_invasive) %>%
  summarise(
    mean_theta = mean(theta),
    median_theta = median(theta),
    mean_p_empirical = mean(p_empirical),
    median_p_empirical = median(p_empirical),
    n = n()
  )

# By reliability tier
combined_all$site_level %>%
  filter(is_invasive) %>%
  mutate(
    Reliable = case_when(
      n_biosamples_positive >= 2 & p_empirical >= 0.5 ~ "Reliable",
      n_biosamples_positive == 3 & total_biosamples_at_site == 3 & p_empirical < 0.5 ~ "Marginal",
      n_biosamples_positive == 2 & total_biosamples_at_site == 3 & pcr_spread_ratio > 0.5 ~ "Marginal",
      TRUE ~ "Unreliable"
    )
  ) %>%
  group_by(Reliable) %>%
  summarise(
    mean_theta = mean(theta),
    median_theta = median(theta),
    mean_p_empirical = mean(p_empirical),
    median_p_empirical = median(p_empirical),
    n = n()
  )

# Build combined reliability DF for downstream GBIF novelty checks
combined_all$site_level <- combined_all$site_level %>%
  mutate(
    reliability = case_when(
      n_biosamples_positive >= 2 & p_empirical >= 0.5 ~ "Reliable",
      n_biosamples_positive == 3 & total_biosamples_at_site == 3 & p_empirical < 0.5 ~ "Marginal",
      n_biosamples_positive == 2 & total_biosamples_at_site == 3 & pcr_spread_ratio > 0.5 ~ "Marginal",
      TRUE ~ "Unreliable"
    )
  )

# =============================================================================
# STEP 11: Summary tables - Use ORIGINAL phyloseq objects for OTU counts
# =============================================================================

####### TABLE 1: Summary table per marker (uses ORIGINAL data for OTU counts)

seq_summary <- tibble(
  Marker = c("12S", "18S", "COI"),
  total_seqs = pre_filter_seqs,
  total_otus_before = pre_filter_otus
)

# Count OTUs from ORIGINAL (uncollapsed) phyloseq objects
# Excludes "uncultured" labels and species without a genus assignment
get_tax_summary_original <- function(ps, marker) {
  tax <- as.data.frame(tax_table(ps))
  tibble(
    Marker = marker,
    total_seqs_after_filter = sum(sample_sums(ps)),
    total_otus_after_filter = ntaxa(ps),
    n_to_genus = sum(!is.na(tax$Genus) & tax$Genus != "" &
                     !grepl("uncultured", tax$Genus, ignore.case = TRUE)),
    n_to_species = sum(!is.na(tax$Species) & tax$Species != "" &
                       !grepl("uncultured", tax$Species, ignore.case = TRUE) &
                       !is.na(tax$Genus) & tax$Genus != "" &
                       !grepl("uncultured", tax$Genus, ignore.case = TRUE))
  )
}

tax_summary <- bind_rows(
  get_tax_summary_original(ps_cs_12S_original, "12S"),
  get_tax_summary_original(ps_cs_18S_original, "18S"),
  get_tax_summary_original(ps_cs_COI_original, "COI")
)

# Count invasive OTUs from ORIGINAL data AFTER manual corrections
# Uses combined_all$site_level which has the corrected is_invasive flags
get_invasive_otu_count <- function(ps, marker, combined_all) {
  # Get species flagged as invasive for THIS MARKER from combined_all (after manual fixes)
  invasive_spp <- combined_all$site_level %>%
    filter(Marker == marker & is_invasive) %>%
    distinct(Species) %>%
    pull(Species)
  
  # Count OTUs for those species in ORIGINAL (uncollapsed) data
  tax <- as.data.frame(tax_table(ps))
  n_invasive_otus <- sum(tax$Species %in% invasive_spp)
  
  tibble(Marker = marker, n_invasive = n_invasive_otus)
}

invasive_otu_counts <- bind_rows(
  get_invasive_otu_count(ps_cs_12S_original, "12S", combined_all),
  get_invasive_otu_count(ps_cs_18S_original, "18S", combined_all),
  get_invasive_otu_count(ps_cs_COI_original, "COI", combined_all)
)

marker_table <- seq_summary %>%
  left_join(tax_summary, by = "Marker") %>%
  left_join(invasive_otu_counts, by = "Marker") %>%
  mutate(
    prop_seqs_retained = round(total_seqs_after_filter / total_seqs, 3),
    prop_otus_retained = round(total_otus_after_filter / total_otus_before, 3),
    prop_to_genus = round(n_to_genus / total_otus_after_filter, 3),
    prop_to_species = round(n_to_species / total_otus_after_filter, 3),
    prop_invasive = round(n_invasive / total_otus_after_filter, 4)
  ) %>%
  select(
    Marker,
    total_seqs, prop_seqs_retained,
    total_otus_before, prop_otus_retained,
    prop_to_genus, prop_to_species, n_invasive, prop_invasive
  )

marker_table2 <- seq_summary %>%
  left_join(tax_summary, by = "Marker") %>%
  left_join(invasive_otu_counts, by = "Marker") %>%
  select(
    Marker,
    total_seqs, total_seqs_after_filter,
    total_otus_before, total_otus_after_filter,
    n_to_genus, n_to_species, n_invasive
  )

cat("\n=== TABLE 1: Marker Summary Statistics ===\n")
print(marker_table)
print(marker_table2)

# =============================================================================
# STEP 12: Read stats per PCR and biosample (uses ORIGINAL)
# =============================================================================

get_read_stats <- function(ps, marker) {
  sdata <- as(sample_data(ps), "data.frame")
  pcr_reads <- sample_sums(ps)
  
  biosample_reads <- data.frame(
    Name = sdata[names(pcr_reads), "Name"],
    reads = pcr_reads
  ) %>%
    group_by(Name) %>%
    summarise(biosample_reads = sum(reads), .groups = "drop")
  
  tibble(
    Marker = marker,
    mean_reads_per_pcr = round(mean(pcr_reads)),
    sd_reads_per_pcr = round(sd(pcr_reads)),
    median_reads_per_pcr = round(median(pcr_reads)),
    iqr_reads_per_pcr = round(IQR(pcr_reads)),
    min_reads_per_pcr = min(pcr_reads),
    max_reads_per_pcr = max(pcr_reads),
    mean_reads_per_biosample = round(mean(biosample_reads$biosample_reads)),
    sd_reads_per_biosample = round(sd(biosample_reads$biosample_reads)),
    median_reads_per_biosample = round(median(biosample_reads$biosample_reads)),
    iqr_reads_per_biosample = round(IQR(biosample_reads$biosample_reads)),
    min_reads_per_biosample = min(biosample_reads$biosample_reads),
    max_reads_per_biosample = max(biosample_reads$biosample_reads)
  )
}

# Use ORIGINAL for read statistics
read_stats <- bind_rows(
  get_read_stats(ps_cs_12S_original, "12S"),
  get_read_stats(ps_cs_18S_original, "18S"),
  get_read_stats(ps_cs_COI_original, "COI")
)

cat("\n=== Read Statistics ===\n")
print(read_stats)

# =============================================================================
# STEP 13: Further MIS summaries
# =============================================================================

# Check which removed non-marine species are invasive
fw_species <- non_marine_removed %>%
  distinct(Species) %>%
  left_join(
    bind_rows(
      inv_12S$invasive_status %>% mutate(Marker = "12S"),
      inv_18S$invasive_status %>% mutate(Marker = "18S"),
      inv_COI$invasive_status %>% mutate(Marker = "COI")
    ) %>% distinct(Species, .keep_all = TRUE),
    by = "Species"
  )

cat("=== Removed non-marine species - invasive status ===\n")
fw_invasive <- fw_species %>%
  filter(grepl("Invasive", GISD_status))

cat("Freshwater invasive species removed:", nrow(fw_invasive), "\n")
fw_invasive %>% select(Species, GISD_status) %>% as_tibble() %>% print(n = Inf)

# Classify the invasive species as freshwater vs marine
invasive_env <- combined_all$site_level %>%
  filter(is_invasive) %>%
  distinct(Species) %>%
  left_join(worms_env %>% select(Species, isMarine, isBrackish, isFreshwater), by = "Species")

cat("=== Invasive species by environment ===\n")
invasive_env %>% arrange(isMarine) %>% print(n = Inf)

cat("\nMarine invasive species:", sum(invasive_env$isMarine == TRUE, na.rm = TRUE), "\n")
cat("Freshwater/non-marine invasive species:", sum(invasive_env$isMarine == FALSE | is.na(invasive_env$isMarine)), "\n")

# Marine invasive species per sampling event
marine_inv_spp <- invasive_env %>% filter(isMarine == TRUE) %>% pull(Species)

site_marine_inv <- combined_all$site_level %>%
  filter(is_invasive, Species %in% marine_inv_spp) %>%
  group_by(Country, Site) %>%
  summarise(n_marine_inv_spp = n_distinct(Species), .groups = "drop")

# Include sites with zero
all_sites <- combined_all$site_level %>% distinct(Country, Site)
site_marine_inv_full <- all_sites %>%
  left_join(site_marine_inv, by = c("Country", "Site")) %>%
  mutate(n_marine_inv_spp = replace_na(n_marine_inv_spp, 0))

cat("\n=== Marine invasive species per sampling event ===\n")
cat("Mean:", round(mean(site_marine_inv_full$n_marine_inv_spp), 2), "\n")
cat("SD:", round(sd(site_marine_inv_full$n_marine_inv_spp), 2), "\n")
cat("Range:", min(site_marine_inv_full$n_marine_inv_spp), "-", max(site_marine_inv_full$n_marine_inv_spp), "\n")
cat("Sites with 0 marine invasives:", sum(site_marine_inv_full$n_marine_inv_spp == 0), "of", nrow(site_marine_inv_full), "\n")

# =============================================================================
# STEP 14: MIS Taxonomy Table (n_OTUs from ORIGINAL, Sites from COLLAPSED)
# =============================================================================

# Get n_OTUs per species from ORIGINAL data
get_species_otu_counts <- function(ps, marker) {
  tax <- as.data.frame(tax_table(ps))
  tax %>%
    filter(!is.na(Species) & Species != "") %>%
    group_by(Species) %>%
    summarise(n_OTUs = n(), .groups = "drop") %>%
    mutate(Marker = marker)
}

species_otu_counts <- bind_rows(
  get_species_otu_counts(ps_cs_12S_original, "12S"),
  get_species_otu_counts(ps_cs_18S_original, "18S"),
  get_species_otu_counts(ps_cs_COI_original, "COI")
)

# Build MIS summary: Sites from combined_all (collapsed), n_OTUs from original
mis_summary <- combined_all$site_level %>%
  filter(is_invasive) %>%
  group_by(Species) %>%
  summarise(
    Markers = paste(sort(unique(Marker)), collapse = ", "),
    Sites = paste(sort(unique(Site)), collapse = ", "),
    .groups = "drop"
  ) %>%
  # Get n_OTUs from ORIGINAL data (sum across markers if species detected by multiple)
  left_join(
    species_otu_counts %>%
      group_by(Species) %>%
      summarise(n_OTUs = sum(n_OTUs), .groups = "drop"),
    by = "Species"
  ) %>%
  # Add taxonomy from ORIGINAL
  left_join(
    bind_rows(
      as.data.frame(tax_table(ps_cs_12S_original)) %>% mutate(OTU = rownames(.)),
      as.data.frame(tax_table(ps_cs_18S_original)) %>% mutate(OTU = rownames(.)),
      as.data.frame(tax_table(ps_cs_COI_original)) %>% mutate(OTU = rownames(.))
    ) %>% distinct(Species, .keep_all = TRUE) %>%
      select(Species, Phylum, Class, Order, Family),
    by = "Species"
  ) %>%
  select(Phylum, Class, Order, Family, Species, n_OTUs, Markers, Sites) %>%
  arrange(Phylum, Class, Order, Species)

cat("\n=== MIS TAXONOMY TABLE ===\n")
mis_summary %>% as_tibble() %>% print(n = Inf, width = Inf)

# Count by phylum
cat("\n=== MIS by Phylum ===\n")
mis_summary %>% count(Phylum, sort = TRUE) %>% print(n = Inf)


#### Also collecting WoRMS habitat / lifestyle information:

source("Scripts/Functions/FUNC_worms_lifestyle.R")
all_spp <- unique(na.omit(combined_all$site_level$Species))
worms_traits <- query_worms_functional_groups(all_spp, 
                                              cache_file = "Processed_data/worms_functional_groups.csv")

### Supplementing with manual designations:
# Step 1: Manual mapping (primary — 100% coverage)
source("Scripts/Functions/FUNC_assign_lifestyle.R")
combined_all$site_level <- assign_lifestyle(combined_all$site_level)
check_unmatched_lifestyle(combined_all$site_level)

# Step 2: Combine with WoRMS (validation — 18.6% coverage)
source("Scripts/Functions/FUNC_combine_lifestyle.R")
worms_traits <- read.csv("Processed_data/worms_functional_groups.csv")
validated <- combine_lifestyle(combined_all$site_level, worms_traits, prefer = "manual")

# Step 3: Review disagreements and apply
combined_all$site_level <- validated$site_level
print(validated$disagreements, n = Inf) # check what WoRMS disagrees on

## Review and fix manually
combined_all$site_level <- validated$site_level %>%
  mutate(lifestyle_final = case_when(
    # Ctenophores are zooplankton, not sessile
    Phylum == "Ctenophora" ~ "Zooplankton",
    # Holoplanktonic hydromedusae
    Species == "Laodicea undulata" ~ "Zooplankton",
    # Heterotrophic dinos misclassified as phytoplankton
    Species %in% c("Gyrodinium rubrum", "Islandinium tricingulatum", 
                   "Noctiluca scintillans") ~ "Zooplankton",
    # Heterotrophic flagellates
    Species %in% c("Paraphysomonas imperforata", 
                   "Picomonas judraskeda") ~ "Microeukaryote",
    # Tintinnids and choanoflagellates are zooplankton
    Order == "Tintinnida" ~ "Zooplankton",
    Order == "Choreotrichida" ~ "Zooplankton",
    Order == "Acanthoecida" ~ "Zooplankton",
    Species == "Leucocryptos marina" ~ "Microeukaryote",
    TRUE ~ lifestyle_final
  ))


# =============================================================================
# STEP 14b: Build pcr_all & metadata for downstream GBIF novelty checks
# =============================================================================

build_pcr_all <- function(ps, marker) {
  otu_mat <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) otu_mat <- t(otu_mat)
  
  tax <- as.data.frame(tax_table(ps))
  sdata <- as(sample_data(ps), "data.frame")
  
  # Total reads per sample for proportion calculation
  sample_totals <- colSums(otu_mat)
  
  # Melt to long format
  df <- as.data.frame(otu_mat) %>%
    rownames_to_column("OTU") %>%
    pivot_longer(-OTU, names_to = "Sample", values_to = "reads") %>%
    filter(reads > 0) %>%
    left_join(tax %>% rownames_to_column("OTU") %>% select(OTU, Species, Genus, Family, Class, Phylum), 
              by = "OTU") %>%
    left_join(sdata %>% rownames_to_column("Sample") %>% 
                select(Sample, Site_Code = Sampling.event.ID, Country = Location),
              by = "Sample") %>%
    mutate(
      Marker = marker,
      total_reads_at_sample = sample_totals[Sample],
      prop_reads_pcr = reads / total_reads_at_sample
    )
  df
}

pcr_all <- bind_rows(
  build_pcr_all(ps_cs_12S_original, "12S"),
  build_pcr_all(ps_cs_18S_original, "18S"),
  build_pcr_all(ps_cs_COI_original, "COI")
)

cat("pcr_all built:", nrow(pcr_all), "PCR × OTU records\n")

# =============================================================================
# FINAL: Assign standard names for downstream scripts
# =============================================================================
# Other scripts (GBIF_MIS_novel_check.R, Chi/Wilcoxon tests, etc.) expect 
# ps_cs_12S, ps_cs_18S, ps_cs_COI to exist.
# Point them to COLLAPSED versions for detection-based analyses.
# Use ps_cs_*_original for OTU counts and diversity analyses.
# =============================================================================

ps_cs_12S <- ps_cs_12S_collapsed
ps_cs_18S <- ps_cs_18S_collapsed
ps_cs_COI <- ps_cs_COI_collapsed

saveRDS(ps_cs_12S, file = "Processed_data/ps_cs_12S.rds")
saveRDS(ps_cs_18S, file = "Processed_data/ps_cs_18S.rds")
saveRDS(ps_cs_COI, file = "Processed_data/ps_cs_COI.rds")

cat("\n=============================================================================\n")
cat("PROCESSING COMPLETE\n")
cat("=============================================================================\n")
cat("\nPhyloseq objects available:\n")
cat("  ps_cs_12S, ps_cs_18S, ps_cs_COI         -> COLLAPSED (for detection analysis)\n")
cat("  ps_cs_12S_original, etc.                -> ORIGINAL (for unmodified OTU counts, e.g. paper Table 1)\n")
cat("\nKey data objects:\n")
cat("  combined_all                            -> All markers combined (COLLAPSED)\n")
cat("  collapsed_species_info                  -> Species with multiple OTUs collapsed\n")
cat("  mis_summary                             -> MIS taxonomy table\n")
cat("  marker_table / marker_table2            -> Table 1 summary stats\n")
cat("=============================================================================\n")



# =============================================================================
# STEP 15: Wheel of Life Visualizations
# =============================================================================

# Here there are a set of codes and functions with an R wrapper that generates customizable wheel of life plots from multiple markers taken from each sampling event

ps_12s <- readRDS("Processed_data/ps_cs_12S.rds")

ps_18s <- readRDS("Processed_data/ps_cs_18S.rds")

ps_COI <- readRDS("Processed_data/ps_cs_COI.rds")

####


## Lets fix some place labels first 

#### Make the wheel of life plots - 

source('Scripts/Functions/wheel_of_life_R.R')

generate_wheels(
  `Vertebrates 12S` = ps_12s,
  `Eukaryotes 18S` = ps_18s,
  `Metazoa COI` = ps_COI,
  output_dir = 'Figures/wheels',
  prefix = 'Combined',
  site_var = 'Sampling.event.ID',
  collector_var = 'Sampling.event.Collected.by',
  location_var = 'Sampling.area.Name',
  date_var = 'Sampling.event.date',
  metadata_from = 'Vertebrates 12S',
  include_genus = TRUE,
  min_group_size = 4,    # Groups with <4 taxa → "Others"
  center_image = 'Figures/BGE-DNA.png',
  center_image_zoom = 0.04,
  script_dir = 'Scripts/Functions/',
  python_path = '/Users/glenndunshea/miniconda3/bin/python3'
)

