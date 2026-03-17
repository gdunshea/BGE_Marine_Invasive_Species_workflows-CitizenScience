# =============================================================================
# FUNC_detection_reliability.R
# =============================================================================
#
# Four standalone analysis functions for characterising eDNA detection
# reliability from hierarchical sampling designs. Each function runs a
# complete analysis section: classification, summaries, figures, and
# saved outputs. All arbitrary thresholds are exposed as function arguments.
#
# FUNCTIONS:
#
#   analyse_detection_strength()
#     Classify every OTU × site detection into strength tiers based on
#     biosample/PCR replication patterns. Produces summary tables and
#     a bar chart.
#
#   analyse_filtering_tradeoffs()
#     Apply confirmation thresholds (e.g. ≥2 PCRs, ≥2 biosamples) and
#     quantify what proportion of detections and invasive species are
#     retained or lost under each. Produces trade-off table and figure.
#
#   analyse_power_by_strength()
#     Break the bootstrap power analysis down by detection strength tier
#     across multiple sampling designs. Shows which detection categories
#     benefit most from increased effort.
#
#   analyse_invasive_reliability()
#     Per-detection reliability assessment for invasive species, with
#     optional side-by-side comparison of two sampling designs (e.g.
#     3×3 vs 5×3). Produces forest plot with credible intervals.
#
# DEPENDENCIES:
#   tidyverse, viridis, scales
#
# =============================================================================

library(tidyverse)
library(viridis)


# =============================================================================
# 1. analyse_detection_strength()
# =============================================================================
#
# WHAT IT DOES:
#   Assigns every OTU × site detection a strength tier based on how many
#   biosamples and PCR replicates confirmed it. Reports the distribution
#   of tiers across markers.
#
# CALL:
#   result <- analyse_detection_strength(combined_all$site_level)
#
# CUSTOMISE TIERS:
#   Default tiers (evaluated in order, first match wins):
#     Strong:   θ = 1 AND mean p ≥ 0.67
#     Good:     θ = 1
#     Moderate: ≥2 biosamples positive
#     Weak:     1 biosample, ≥2 PCRs
#     Minimal:  1 biosample, 1 PCR (everything else)
#
#   To change, pass a named list to `tier_defs`. Each element must have:
#     - label: display name (character)
#     - fn: function(df) returning logical vector, where df has columns
#           theta, mean_p_pcr, n_biosamples_positive, n_pcr_positive
#
#   Example — add a "Perfect" tier and relax "Strong":
#     analyse_detection_strength(
#       site_level,
#       tier_defs = list(
#         perfect  = list(label = "Perfect (3/3 bio, 9/9 PCR)",
#                          fn = function(df) df$theta == 1 & df$mean_p_pcr == 1),
#         strong   = list(label = "Strong (3/3 bio, most PCR)",
#                          fn = function(df) df$theta == 1 & df$mean_p_pcr >= 0.5),
#         moderate = list(label = "Moderate (2+ bio)",
#                          fn = function(df) df$n_biosamples_positive >= 2),
#         weak     = list(label = "Weak (1 bio, 2+ PCR)",
#                          fn = function(df) df$n_pcr_positive >= 2),
#         minimal  = list(label = "Minimal (1 bio, 1 PCR)",
#                          fn = function(df) rep(TRUE, nrow(df)))
#       )
#     )
#
# RETURNS:
#   List with:
#     $classified    - Input data with detection_strength column added
#     $summary       - Counts/percentages by Marker × tier
#     $overall       - Counts/percentages across all markers
#     $plot          - ggplot bar chart
# -----------------------------------------------------------------------------

analyse_detection_strength <- function(
    site_level,
    tier_defs = list(
      strong   = list(label = "Strong (all bio, most PCR)",
                       fn = function(df) df$theta == 1 & df$mean_p_pcr >= 0.67),
      good     = list(label = "Good (all bio)",
                       fn = function(df) df$theta == 1),
      moderate = list(label = "Moderate (2+ bio)",
                       fn = function(df) df$n_biosamples_positive >= 2),
      weak     = list(label = "Weak (1 bio, 2+ PCR)",
                       fn = function(df) df$n_pcr_positive >= 2),
      minimal  = list(label = "Minimal (1 bio, 1 PCR)",
                       fn = function(df) rep(TRUE, nrow(df)))
    ),
    save_csv = NULL,
    save_fig = NULL
) {
  
  tier_labels <- sapply(tier_defs, `[[`, "label")
  
  # Classify: first match wins
  strength <- rep(NA_character_, nrow(site_level))
  for (tier in tier_defs) {
    unassigned <- is.na(strength)
    if (!any(unassigned)) break
    strength[which(unassigned & tier$fn(site_level))] <- tier$label
  }
  
  classified <- site_level %>%
    mutate(
      detection_strength = factor(strength, levels = tier_labels),
      single_pcr_single_bio = n_pcr_positive == 1 & n_biosamples_positive == 1
    )
  
  # Summary by marker
  summary_df <- classified %>%
    group_by(Marker, detection_strength) %>%
    summarise(
      n = n(),
      n_invasive = sum(is_invasive, na.rm = TRUE),
      median_reads = median(total_reads),
      median_theta = round(median(theta), 3),
      median_p = round(median(mean_p_pcr), 3),
      .groups = "drop"
    ) %>%
    group_by(Marker) %>%
    mutate(pct = round(n / sum(n) * 100, 1)) %>%
    ungroup()
  
  # Overall
  overall_df <- classified %>%
    group_by(detection_strength) %>%
    summarise(
      n = n(),
      pct = round(n() / nrow(classified) * 100, 1),
      n_invasive = sum(is_invasive, na.rm = TRUE),
      median_reads = median(total_reads),
      .groups = "drop"
    )
  
  # Plot
  fig <- ggplot(classified, aes(x = detection_strength, fill = Marker)) +
    geom_bar(position = "dodge", alpha = 0.85) +
    scale_fill_brewer(palette = "Set1") +
    labs(x = NULL, y = "Number of OTU x site detections",
         title = "Detection strength distribution",
         subtitle = paste0("N = ", nrow(classified), " detections")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  # Print
  cat("=== Detection strength by marker ===\n")
  print(summary_df, n = Inf)
  cat("\n=== Overall ===\n")
  print(overall_df)
  n_single <- sum(classified$single_pcr_single_bio)
  cat("\nSingle-PCR single-biosample:", n_single, "of", nrow(classified),
      paste0("(", round(n_single / nrow(classified) * 100, 1), "%)\n"))
  
  # Save
  if (!is.null(save_csv)) write.csv(summary_df, save_csv, row.names = FALSE)
  if (!is.null(save_fig)) ggsave(save_fig, fig, width = 11, height = 6)
  
  list(classified = classified, summary = summary_df,
       overall = overall_df, plot = fig)
}


# =============================================================================
# 2. analyse_filtering_tradeoffs()
# =============================================================================
#
# WHAT IT DOES:
#   Applies one or more confirmation thresholds to the detection data and
#   quantifies what proportion of all detections (and invasive detections)
#   pass each threshold. Identifies which specific invasive species would
#   be lost.
#
# CALL:
#   result <- analyse_filtering_tradeoffs(classified)
#
#   where `classified` is $classified from analyse_detection_strength(),
#   or any data frame with columns: n_pcr_positive, n_biosamples_positive,
#   mean_p_pcr, is_invasive, Species, Marker.
#
# CUSTOMISE THRESHOLDS:
#   Default thresholds:
#     No filter:              keep everything
#     ≥2 PCRs total:          n_pcr_positive >= 2
#     ≥2 biosamples:          n_biosamples_positive >= 2
#     ≥2 PCRs OR ≥2 bio:     either condition
#     ≥2 bio AND p≥0.67:     both conditions
#
#   To change, pass a named list to `threshold_defs`. Each element:
#     - label: display name
#     - fn: function(df) returning logical vector
#
#   Example — add a read-count threshold:
#     analyse_filtering_tradeoffs(
#       classified,
#       threshold_defs = list(
#         none     = list(label = "No filter",
#                          fn = function(df) rep(TRUE, nrow(df))),
#         reads100 = list(label = ">=100 reads",
#                          fn = function(df) df$total_reads >= 100),
#         pcr2     = list(label = ">=2 PCRs",
#                          fn = function(df) df$n_pcr_positive >= 2)
#       )
#     )
#
# ARGUMENTS:
#   inspect_threshold  Name of threshold to report invasive losses for
#                      (default "bio2" — matches a key in threshold_defs)
#
# RETURNS:
#   List with:
#     $tradeoffs       - Per-marker retention/loss for each threshold
#     $tradeoffs_all   - Same but across all markers combined
#     $invasive_losses - Invasive detections that fail inspect_threshold
#     $plot            - Grouped bar chart of retention
# -----------------------------------------------------------------------------

analyse_filtering_tradeoffs <- function(
    classified,
    threshold_defs = list(
      none = list(
        label = "No filter",
        fn = function(df) rep(TRUE, nrow(df))),
      pcr2 = list(
        label = ">=2 PCRs total",
        fn = function(df) df$n_pcr_positive >= 2),
      bio2 = list(
        label = ">=2 biosamples",
        fn = function(df) df$n_biosamples_positive >= 2),
      pcr2_or_bio2 = list(
        label = ">=2 PCRs OR >=2 bio",
        fn = function(df) df$n_pcr_positive >= 2 | df$n_biosamples_positive >= 2),
      strict = list(
        label = ">=2 bio AND p>=0.67",
        fn = function(df) df$n_biosamples_positive >= 2 & df$mean_p_pcr >= 0.67)
    ),
    inspect_threshold = "bio2",
    save_csv = NULL,
    save_fig = NULL
) {
  
  threshold_labels <- setNames(
    sapply(threshold_defs, `[[`, "label"),
    names(threshold_defs)
  )
  
  # Apply all thresholds
  for (nm in names(threshold_defs)) {
    classified[[paste0("passes_", nm)]] <- threshold_defs[[nm]]$fn(classified)
  }
  
  # Summarise per marker
  .summarise_thresholds <- function(data, group_var = "Marker") {
    map_dfr(names(threshold_defs), function(nm) {
      col <- paste0("passes_", nm)
      grp <- if (is.null(group_var)) NULL else group_var
      data %>%
        { if (!is.null(grp)) group_by(., across(all_of(grp))) else (.) } %>%
        summarise(
          threshold = threshold_labels[nm],
          total = n(),
          retained = sum(.data[[col]]),
          lost = total - retained,
          pct_retained = round(retained / total * 100, 1),
          inv_total = sum(is_invasive, na.rm = TRUE),
          inv_retained = sum(is_invasive & .data[[col]], na.rm = TRUE),
          inv_lost = inv_total - inv_retained,
          pct_inv_retained = round(inv_retained / max(inv_total, 1) * 100, 1),
          spp_total = n_distinct(Species[!is.na(Species)]),
          spp_retained = n_distinct(Species[!is.na(Species) & .data[[col]]]),
          spp_lost = spp_total - spp_retained,
          .groups = "drop"
        )
    }) %>%
      mutate(threshold = factor(threshold, levels = unname(threshold_labels)))
  }
  
  tradeoffs <- .summarise_thresholds(classified, "Marker")
  tradeoffs_all <- .summarise_thresholds(classified, NULL)
  
  # Invasive losses for inspected threshold
  inspect_col <- paste0("passes_", inspect_threshold)
  inv_losses <- classified %>%
    filter(is_invasive, !.data[[inspect_col]]) %>%
    select(Marker, Species, Site, Country,
           n_biosamples_positive, total_biosamples_at_site,
           n_pcr_positive, n_pcr_total, total_reads, pcr_distribution) %>%
    arrange(Marker, Species)
  
  # Plot
  plot_data <- tradeoffs %>%
    select(Marker, threshold, pct_retained, pct_inv_retained) %>%
    pivot_longer(c(pct_retained, pct_inv_retained),
                 names_to = "category", values_to = "pct") %>%
    mutate(category = ifelse(category == "pct_retained",
                             "All detections", "Invasive detections"))
  
  fig <- ggplot(plot_data, aes(x = threshold, y = pct / 100, fill = category)) +
    geom_col(position = "dodge", alpha = 0.85, width = 0.7) +
    geom_hline(yintercept = 0.9, linetype = "dashed", alpha = 0.4) +
    facet_wrap(~Marker) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1.05)) +
    scale_fill_manual(values = c("All detections" = "steelblue",
                                  "Invasive detections" = "tomato"), name = NULL) +
    labs(x = "Confirmation threshold", y = "% of detections retained",
         title = "Detection retention under confirmation thresholds") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          legend.position = "bottom")
  
  # Print
  cat("=== Filtering trade-offs by marker ===\n")
  print(tradeoffs, n = Inf, width = Inf)
  cat("\n=== All markers ===\n")
  print(tradeoffs_all)
  cat("\n=== Invasive losses under '", threshold_labels[inspect_threshold], "' ===\n")
  if (nrow(inv_losses) > 0) print(inv_losses, n = Inf, width = Inf) else cat("  None\n")
  
  if (!is.null(save_csv)) write.csv(tradeoffs, save_csv, row.names = FALSE)
  if (!is.null(save_fig)) ggsave(save_fig, fig, width = 14, height = 6)
  
  list(tradeoffs = tradeoffs, tradeoffs_all = tradeoffs_all,
       invasive_losses = inv_losses, plot = fig)
}


# =============================================================================
# 3. analyse_power_by_strength()
# =============================================================================
#
# WHAT IT DOES:
#   Breaks the bootstrap power analysis down by detection strength tier,
#   showing what proportion of detections achieve a target P(detect) under
#   each sampling design. Reveals which detection categories benefit most
#   from increased effort.
#
# CALL:
#   result <- analyse_power_by_strength(
#     empirical_3x3 = empirical_3x3,
#     det_list = list(`12S` = det_12S, `18S` = det_18S, COI = det_COI)
#   )
#
# CUSTOMISE:
#   p_col         Column name for PCR detection rate in empirical_3x3.
#                 Default "p" (renamed from mean_p_pcr in Step 8b).
#   p_threshold   Target detection probability. Default 0.80.
#   designs       Character vector of design labels to compare.
#                 Default: c("3×3", "4×3", "5×3", "6×3", "8×3")
#   strength_defs Named list of tier definitions (same format as
#                 analyse_detection_strength tier_defs, but operating on
#                 empirical_3x3 columns). Default uses theta, p,
#                 n_biosamples_positive, n_pcr_positive.
#
#   Example — use different tiers and a 95% threshold:
#     analyse_power_by_strength(
#       empirical_3x3, det_list,
#       p_threshold = 0.95,
#       strength_defs = list(
#         high = list(label = "High", fn = function(df) df$theta >= 0.67),
#         low  = list(label = "Low",  fn = function(df) rep(TRUE, nrow(df)))
#       )
#     )
#
# RETURNS:
#   List with:
#     $summary_3x3       - P(detect) stats by strength tier at 3×3
#     $design_comparison  - % achieving threshold per strength × design
#     $plot               - Bar chart: designs × strength tiers
# -----------------------------------------------------------------------------

analyse_power_by_strength <- function(
    empirical_3x3,
    det_list,
    p_col = "p",
    p_threshold = 0.80,
    designs = c("3\u00d73", "4\u00d73", "5\u00d73", "6\u00d73", "8\u00d73"),
    strength_defs = list(
      strong   = list(label = "Strong",
                       fn = function(df) df[[p_col]] >= 0.67 & df$theta == 1),
      good     = list(label = "Good",
                       fn = function(df) df$theta == 1),
      moderate = list(label = "Moderate",
                       fn = function(df) df$n_biosamples_positive >= 2),
      weak     = list(label = "Weak",
                       fn = function(df) df$n_pcr_positive >= 2),
      minimal  = list(label = "Minimal",
                       fn = function(df) rep(TRUE, nrow(df)))
    ),
    save_csv = NULL,
    save_fig = NULL
) {
  
  tier_labels <- sapply(strength_defs, `[[`, "label")
  
  # Helper: assign strength tiers
  .assign_strength <- function(df, p_col_name) {
    strength <- rep(NA_character_, nrow(df))
    for (tier in strength_defs) {
      unassigned <- is.na(strength)
      if (!any(unassigned)) break
      strength[which(unassigned & tier$fn(df))] <- tier$label
    }
    df$detection_strength <- factor(strength, levels = tier_labels)
    df
  }
  
  # Tag empirical_3x3
  emp_strength <- .assign_strength(empirical_3x3, p_col)
  
  # Summary at 3×3
  summary_3x3 <- emp_strength %>%
    group_by(Marker, detection_strength) %>%
    summarise(
      n = n(),
      mean_p_detect = round(mean(p_detect_mean), 3),
      median_p_detect = round(median(p_detect_mean), 3),
      mean_lower_CI = round(mean(p_detect_lower), 3),
      mean_upper_CI = round(mean(p_detect_upper), 3),
      mean_CI_width = round(mean(p_detect_upper - p_detect_lower), 3),
      pct_above = round(mean(p_detect_mean >= p_threshold) * 100, 1),
      pct_above_conservative = round(mean(p_detect_lower >= p_threshold) * 100, 1),
      .groups = "drop"
    )
  
  # Assemble all designs across markers
  all_designs <- map2_dfr(det_list, names(det_list), function(det, marker) {
    det$empirical_with_power %>% mutate(Marker = marker)
  })
  
  # For the full design grid, strength is based on mean_p_pcr (not renamed p)
  # Create a temporary column matching the name expected by strength_defs
  if (p_col != "mean_p_pcr" && p_col %in% names(empirical_3x3)) {
    # strength_defs reference p_col; all_designs has mean_p_pcr
    all_designs[[p_col]] <- all_designs$mean_p_pcr
  }
  
  all_designs_strength <- .assign_strength(all_designs, p_col)
  
  design_comparison <- all_designs_strength %>%
    filter(design %in% designs) %>%
    group_by(detection_strength, design, effort) %>%
    summarise(
      n = n(),
      pct_above = round(mean(p_detect_mean >= p_threshold) * 100, 1),
      pct_above_conservative = round(mean(p_detect_lower >= p_threshold) * 100, 1),
      mean_p_detect = round(mean(p_detect_mean), 3),
      .groups = "drop"
    ) %>%
    arrange(detection_strength, effort)
  
  # Plot
  fig <- ggplot(design_comparison,
                aes(x = design, y = pct_above / 100, fill = detection_strength)) +
    geom_col(position = "dodge", alpha = 0.85, width = 0.8) +
    geom_hline(yintercept = p_threshold, linetype = "dashed", alpha = 0.4) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1.05)) +
    scale_fill_viridis_d(name = "Detection\nstrength", direction = -1) +
    labs(x = "Sampling design (biosamples x PCRs)",
         y = paste0("% of detections with P(detect) >= ",
                     round(p_threshold * 100), "%"),
         title = "Design performance by empirical detection strength") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  # Print
  cat("=== P(detect) by strength (3x3) ===\n")
  print(summary_3x3, n = Inf, width = Inf)
  cat("\n=== Design comparison by strength ===\n")
  print(design_comparison, n = Inf, width = Inf)
  
  if (!is.null(save_csv)) write.csv(design_comparison, save_csv, row.names = FALSE)
  if (!is.null(save_fig)) ggsave(save_fig, fig, width = 12, height = 6)
  
  list(summary_3x3 = summary_3x3, design_comparison = design_comparison,
       plot = fig)
}


# =============================================================================
# 4. analyse_invasive_reliability()
# =============================================================================
#
# WHAT IT DOES:
#   Builds a per-detection reliability table for all invasive species,
#   classifying each by a reliability tier derived from bootstrap P(detect)
#   and its 95% credible interval. Optionally compares the current design
#   against an alternative (e.g. 5×3) and produces a paired forest plot.
#
# CALL (basic — 3×3 assessment only):
#   result <- analyse_invasive_reliability(empirical_3x3)
#
# CALL (with design comparison):
#   result <- analyse_invasive_reliability(
#     empirical_3x3,
#     det_list = list(`12S` = det_12S, `18S` = det_18S, COI = det_COI),
#     combined_all = combined_all,
#     compare_design = "5×3"
#   )
#
# CUSTOMISE RELIABILITY TIERS:
#   Default tiers (evaluated in order, first match wins):
#     Reliable:        lower 95% CI ≥ 0.80
#     Likely reliable: mean ≥ 0.80
#     Uncertain:       mean ≥ 0.50
#     Unreliable:      everything else
#
#   To change, pass a named list to `reliability_defs`. Each element:
#     - label: display name
#     - fn: function(mean, lower) returning logical vector
#
#   Example — stricter tiers using 90% threshold:
#     analyse_invasive_reliability(
#       empirical_3x3,
#       reliability_defs = list(
#         reliable   = list(label = "Reliable",        fn = function(m, l) l >= 0.90),
#         likely     = list(label = "Likely reliable",  fn = function(m, l) m >= 0.90),
#         uncertain  = list(label = "Uncertain",        fn = function(m, l) m >= 0.50),
#         unreliable = list(label = "Unreliable",       fn = function(m, l) rep(TRUE, length(m)))
#       )
#     )
#
# CUSTOMISE DETECTION STRENGTH (for labelling only):
#   Same format as analyse_detection_strength tier_defs, applied to
#   empirical_3x3 columns (theta, p, n_biosamples_positive, n_pcr_positive).
#
# ARGUMENTS:
#   empirical_3x3      Bootstrap 3×3 results (from Step 8b). Required.
#   det_list            Named list of per-marker det results. Required for
#                       design comparison; NULL to skip.
#   combined_all        Combined markers object. Required for design
#                       comparison (for Species/is_invasive join); NULL to skip.
#   compare_design      Design label for comparison (e.g. "5×3"). NULL to skip.
#   p_col               Column name for PCR detection rate. Default "p".
#   reliability_defs    Reliability tier definitions (see above).
#   strength_defs       Detection strength tier definitions (see above).
#
# RETURNS:
#   List with:
#     $reliability     - Per-detection table with reliability + strength
#     $summary         - Counts by Marker × reliability tier
#     $comparison      - Side-by-side design comparison (NULL if skipped)
#     $upgrade_summary - Tier transition counts (NULL if skipped)
#     $forest_plot     - ggplot forest plot (paired if comparison, single if not)
# -----------------------------------------------------------------------------

analyse_invasive_reliability <- function(
    empirical_3x3,
    det_list = NULL,
    combined_all = NULL,
    compare_design = NULL,
    p_col = "p",
    reliability_defs = list(
      reliable   = list(label = "Reliable",
                         fn = function(m, l) l >= 0.80),
      likely     = list(label = "Likely reliable",
                         fn = function(m, l) m >= 0.80),
      uncertain  = list(label = "Uncertain",
                         fn = function(m, l) m >= 0.50),
      unreliable = list(label = "Unreliable",
                         fn = function(m, l) rep(TRUE, length(m)))
    ),
    strength_defs = list(
      strong   = list(label = "Strong",
                       fn = function(df) df[[p_col]] >= 0.67 & df$theta == 1),
      good     = list(label = "Good",
                       fn = function(df) df$theta == 1),
      moderate = list(label = "Moderate",
                       fn = function(df) df$n_biosamples_positive >= 2),
      weak     = list(label = "Weak",
                       fn = function(df) df$n_pcr_positive >= 2),
      minimal  = list(label = "Minimal",
                       fn = function(df) rep(TRUE, nrow(df)))
    ),
    save_csv = NULL,
    save_fig = NULL
) {
  
  # --- Helper: classify reliability ---
  .classify_reliability <- function(mean_vec, lower_vec) {
    rel_labels <- sapply(reliability_defs, `[[`, "label")
    result <- rep(NA_character_, length(mean_vec))
    for (tier in reliability_defs) {
      unassigned <- is.na(result)
      if (!any(unassigned)) break
      result[which(unassigned & tier$fn(mean_vec, lower_vec))] <- tier$label
    }
    factor(result, levels = rel_labels)
  }
  
  # --- Helper: classify strength ---
  .classify_strength <- function(df) {
    str_labels <- sapply(strength_defs, `[[`, "label")
    strength <- rep(NA_character_, nrow(df))
    for (tier in strength_defs) {
      unassigned <- is.na(strength)
      if (!any(unassigned)) break
      strength[which(unassigned & tier$fn(df))] <- tier$label
    }
    df$detection_strength <- factor(strength, levels = str_labels)
    df
  }
  
  # --- Build reliability table ---
  inv <- empirical_3x3 %>% filter(is_invasive)
  if (nrow(inv) == 0) {
    warning("No invasive detections found in empirical_3x3")
    return(list(reliability = tibble(), summary = tibble(),
                comparison = NULL, upgrade_summary = NULL, forest_plot = NULL))
  }
  
  inv <- .classify_strength(inv)
  inv$reliability <- .classify_reliability(inv$p_detect_mean, inv$p_detect_lower)
  
  inv <- inv %>%
    mutate(
      CI_width = round(p_detect_upper - p_detect_lower, 3),
      confirmed = case_when(
        n_biosamples_positive >= 2 ~ "Multi-biosample",
        n_pcr_positive >= 2       ~ "Multi-PCR only",
        TRUE                       ~ "Single PCR only"
      )
    ) %>%
    arrange(p_detect_mean)
  
  # Summary
  inv_summary <- inv %>%
    group_by(Marker, reliability) %>%
    summarise(
      n = n(), n_species = n_distinct(Species),
      mean_p_detect = round(mean(p_detect_mean), 2),
      mean_CI_width = round(mean(CI_width), 2),
      .groups = "drop"
    ) %>%
    group_by(Marker) %>%
    mutate(pct = round(n / sum(n) * 100, 1)) %>%
    ungroup()
  
  # Print
  cat("=== All invasive detections ===\n")
  print(inv %>%
          select(Marker, Species, Site, detection_strength, confirmed,
                 theta, !!sym(p_col), p_detect_mean, p_detect_lower,
                 p_detect_upper, reliability),
        n = Inf, width = Inf)
  cat("\n=== By reliability tier ===\n")
  print(inv_summary, n = Inf)
  
  # --- Design comparison (optional) ---
  comparison <- NULL
  upgrade_summary <- NULL
  
  if (!is.null(compare_design) && !is.null(det_list) && !is.null(combined_all)) {
    
    cat("\n=== Design comparison:", compare_design, "vs baseline ===\n")
    
    # Extract comparison design
    alt <- map2_dfr(det_list, names(det_list), function(det, marker) {
      det$empirical_with_power %>%
        filter(design == compare_design) %>%
        mutate(Marker = marker)
    }) %>%
      left_join(
        combined_all$site_level %>%
          select(OTU, Site, Marker, Species, is_invasive) %>%
          distinct(),
        by = c("OTU", "Site", "Marker")
      ) %>%
      filter(is_invasive) %>%
      select(Marker, OTU, Site, Species,
             p_detect_mean_alt = p_detect_mean,
             p_detect_lower_alt = p_detect_lower,
             p_detect_upper_alt = p_detect_upper)
    
    alt$reliability_alt <- .classify_reliability(
      alt$p_detect_mean_alt, alt$p_detect_lower_alt)
    
    comparison <- inv %>%
      select(Marker, OTU, Site, Species, detection_strength, confirmed,
             p_detect_mean_base = p_detect_mean,
             p_detect_lower_base = p_detect_lower,
             p_detect_upper_base = p_detect_upper,
             reliability_base = reliability) %>%
      left_join(alt, by = c("Marker", "OTU", "Site", "Species")) %>%
      mutate(
        gain = round(p_detect_mean_alt - p_detect_mean_base, 3),
        upgraded = as.character(reliability_alt) != as.character(reliability_base)
      ) %>%
      arrange(p_detect_mean_base)
    
    upgrade_summary <- comparison %>%
      group_by(reliability_base, reliability_alt) %>%
      summarise(n = n(), .groups = "drop")
    
    print(comparison %>%
            select(Marker, Species, Site, detection_strength,
                   p_detect_mean_base, p_detect_mean_alt, gain,
                   reliability_base, reliability_alt, upgraded),
          n = Inf, width = Inf)
    cat("\nUpgraded:", sum(comparison$upgraded, na.rm = TRUE),
        "of", nrow(comparison), "\n")
  }
  
  # --- Forest plot ---
  rel_colors <- c("Reliable" = "#1a9850", "Likely reliable" = "#91cf60",
                   "Uncertain" = "#fc8d59", "Unreliable" = "#d73027")
  # Use only labels that exist in reliability_defs
  rel_labels <- sapply(reliability_defs, `[[`, "label")
  rel_colors <- rel_colors[names(rel_colors) %in% rel_labels]
  
  if (!is.null(comparison)) {
    # Paired forest plot
    base_label <- "Baseline"
    alt_label <- compare_design
    
    base_data <- comparison %>%
      mutate(Species_Site = paste0(Species, " (", Site, ")")) %>%
      select(Marker, Species_Site,
             p_detect_mean = p_detect_mean_base,
             p_detect_lower = p_detect_lower_base,
             p_detect_upper = p_detect_upper_base,
             reliability = reliability_base)
    
    alt_data <- comparison %>%
      mutate(Species_Site = paste0(Species, " (", Site, ")")) %>%
      select(Marker, Species_Site,
             p_detect_mean = p_detect_mean_alt,
             p_detect_lower = p_detect_lower_alt,
             p_detect_upper = p_detect_upper_alt,
             reliability = reliability_base)  # colour by baseline tier
    
    y_order <- base_data %>% arrange(p_detect_mean) %>% pull(Species_Site)
    base_data$Species_Site <- factor(base_data$Species_Site, levels = y_order)
    alt_data$Species_Site <- factor(alt_data$Species_Site, levels = y_order)
    
    fig <- ggplot(mapping = aes(y = Species_Site, x = p_detect_mean,
                                 color = reliability)) +
      geom_errorbarh(data = base_data,
                     aes(xmin = p_detect_lower, xmax = p_detect_upper),
                     height = 0.3, linewidth = 0.5, alpha = 0.5,
                     position = position_nudge(y = 0.12)) +
      geom_point(data = base_data, size = 2.5, shape = 16,
                 position = position_nudge(y = 0.12)) +
      geom_errorbarh(data = alt_data,
                     aes(xmin = p_detect_lower, xmax = p_detect_upper),
                     height = 0.3, linewidth = 0.5, alpha = 0.5,
                     position = position_nudge(y = -0.12)) +
      geom_point(data = alt_data, size = 2.5, shape = 21,
                 position = position_nudge(y = -0.12)) +
      geom_vline(xintercept = 0.5, linetype = "dotted", alpha = 0.5) +
      geom_vline(xintercept = 0.8, linetype = "dashed", alpha = 0.5) +
      scale_x_continuous(limits = c(0, 1), labels = scales::percent) +
      scale_color_manual(values = rel_colors, name = "Reliability") +
      facet_wrap(~Marker, scales = "free_y", ncol = 3) +
      labs(x = "P(detect | present) with 95% CI", y = NULL,
           title = paste("Invasive detection reliability: baseline vs", compare_design),
           subtitle = paste("Filled = baseline, Open =", compare_design)) +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 7),
            legend.position = "bottom",
            strip.text = element_text(face = "bold"))
    
  } else {
    # Single-design forest plot
    plot_data <- inv %>%
      mutate(Species_Site = paste0(Species, " (", Site, ")"))
    
    fig <- ggplot(plot_data,
                  aes(y = reorder(Species_Site, p_detect_mean),
                      x = p_detect_mean, color = reliability)) +
      geom_errorbarh(aes(xmin = p_detect_lower, xmax = p_detect_upper),
                     height = 0.3, linewidth = 0.5, alpha = 0.7) +
      geom_point(size = 2) +
      geom_vline(xintercept = 0.5, linetype = "dotted", alpha = 0.5) +
      geom_vline(xintercept = 0.8, linetype = "dashed", alpha = 0.5) +
      scale_x_continuous(limits = c(0, 1), labels = scales::percent) +
      scale_color_manual(values = rel_colors, name = "Reliability") +
      facet_wrap(~Marker, scales = "free_y", ncol = 3) +
      labs(x = "P(detect | present) with 95% CI", y = NULL,
           title = "Invasive species detection reliability") +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 7),
            legend.position = "bottom",
            strip.text = element_text(face = "bold"))
  }
  
  if (!is.null(save_csv)) write.csv(inv, save_csv, row.names = FALSE)
  if (!is.null(save_fig)) ggsave(save_fig, fig, width = 16, height = 12)
  
  list(reliability = inv, summary = inv_summary,
       comparison = comparison, upgrade_summary = upgrade_summary,
       forest_plot = fig)
}


# =============================================================================
# LOADED
# =============================================================================

cat("FUNC_detection_reliability.R loaded\n")
cat("Functions:\n")
cat("  analyse_detection_strength(site_level, tier_defs, ...)\n")
cat("  analyse_filtering_tradeoffs(classified, threshold_defs, ...)\n")
cat("  analyse_power_by_strength(empirical_3x3, det_list, ...)\n")
cat("  analyse_invasive_reliability(empirical_3x3, det_list, combined_all, ...)\n")
