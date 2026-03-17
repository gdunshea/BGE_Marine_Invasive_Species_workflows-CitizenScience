# =============================================================================
# 03_CS_read_proportion_validation.R
# =============================================================================
#
# Run after 02b_CS_detection_reliability.R.
# Expects: combined_all, empirical_3x3,
#          ps_cs_12S, ps_cs_18S, ps_cs_COI (phyloseq objects),
#          gbif_validation (GBIF/literature validation of invasive detections)
#
# PURPOSE:
# Adds a complementary read proportion analysis to the occupancy-based
# detection reliability framework (02b). While θ and p measure detection
# *consistency* across replicates, read proportion measures signal
# *strength* — these are independent axes that together create a 2×2
# classification:
#
#   ┌─────────────────────┬─────────────────────┐
#   │     CONCERNING      │     CONFIDENT       │  High P(detect)
#   │ Consistent but weak │ Real & reliable      │
#   ├─────────────────────┼─────────────────────┤
#   │     UNCERTAIN       │   PATCHY BUT REAL   │  Low P(detect)
#   │ Weak & inconsistent │ Strong but patchy    │
#   └─────────────────────┴─────────────────────┘
#         Low read prop       High read prop
#
# Read proportion is calculated across ALL PCRs (including zeros), so it
# reflects signal strength relative to total sampling effort. This is
# validated against GBIF occurrence records for invasive species.
#
# =============================================================================

library(tidyverse)
library(phyloseq)
library(pROC)
library(patchwork)

gbif_validation <- read.csv("Processed_data/All_MIS_detections_GBIF_validation.csv")

cat("=============================================================================\n")
cat("03: READ PROPORTION & UNIFIED CLASSIFICATION\n")
cat("=============================================================================\n\n")


# =============================================================================
# PART 1: CALCULATE READ PROPORTIONS (ALL PCRs INCLUDING ZEROS)
# =============================================================================

cat("PART 1: Calculating read proportions (all PCRs including zeros)...\n")

calculate_read_proportions <- function(ps, marker_name) {

  otu <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) otu <- t(otu)

  lib_sizes <- colSums(otu)
  otu_prop <- sweep(otu, 2, lib_sizes, "/")

  samp_data <- data.frame(sample_data(ps))
  samp_data$PCR_id <- rownames(samp_data)

  site_col <- if ("Sampling.event.ID" %in% names(samp_data)) "Sampling.event.ID" else "Site"
  biosample_col <- if ("Name" %in% names(samp_data)) "Name" else "Biosample"

  otu_long <- otu_prop %>%
    as.data.frame() %>%
    rownames_to_column("OTU") %>%
    pivot_longer(-OTU, names_to = "PCR_id", values_to = "read_proportion") %>%
    mutate(Marker = marker_name) %>%
    left_join(
      samp_data %>% select(PCR_id, all_of(site_col), all_of(biosample_col)),
      by = "PCR_id"
    ) %>%
    rename(Site = all_of(site_col), Biosample = all_of(biosample_col))

  return(otu_long)
}

prop_12S <- calculate_read_proportions(ps_cs_12S, "12S")
prop_18S <- calculate_read_proportions(ps_cs_18S, "18S")
prop_COI <- calculate_read_proportions(ps_cs_COI, "COI")
all_proportions <- bind_rows(prop_12S, prop_18S, prop_COI)

cat("  Total OTU × PCR combinations:", nrow(all_proportions), "\n")
cat("  Non-zero:", sum(all_proportions$read_proportion > 0), "\n")
cat("  Zero:", sum(all_proportions$read_proportion == 0),
    "(", round(mean(all_proportions$read_proportion == 0) * 100, 1), "%)\n")


# =============================================================================
# PART 2: HIERARCHICAL BOOTSTRAP
# =============================================================================

cat("\nPART 2: Hierarchical bootstrap (1000 iterations)...\n")

hierarchical_bootstrap <- function(site_data, n_boot = 1000) {
  biosamples <- unique(site_data$Biosample)
  n_bio <- length(biosamples)

  boot_means <- numeric(n_boot)

  for (i in 1:n_boot) {
    boot_bio <- sample(biosamples, n_bio, replace = TRUE)
    biosample_means <- sapply(boot_bio, function(b) {
      pcrs <- site_data$read_proportion[site_data$Biosample == b]
      mean(sample(pcrs, length(pcrs), replace = TRUE))
    })
    boot_means[i] <- mean(biosample_means)
  }

  return(boot_means)
}

# Filter to detected OTU-Sites only (at least one positive PCR)
detected_otu_sites <- all_proportions %>%
  group_by(Marker, OTU, Site) %>%
  summarise(ever_detected = any(read_proportion > 0), .groups = "drop") %>%
  filter(ever_detected) %>%
  select(-ever_detected)

all_proportions_filtered <- all_proportions %>%
  semi_join(detected_otu_sites, by = c("Marker", "OTU", "Site"))

n_combos <- nrow(detected_otu_sites)
cat("  Processing", n_combos, "OTU-Site combinations...\n")

bootstrap_results <- list()
pb_interval <- max(1, floor(n_combos / 10))

for (i in 1:n_combos) {
  combo <- detected_otu_sites[i, ]

  site_data <- all_proportions_filtered %>%
    filter(Marker == combo$Marker, OTU == combo$OTU, Site == combo$Site)

  obs_biosample_means <- site_data %>%
    group_by(Biosample) %>%
    summarise(bio_mean = mean(read_proportion), .groups = "drop")

  obs_mean <- mean(obs_biosample_means$bio_mean)
  boot_means <- hierarchical_bootstrap(site_data, n_boot = 1000)

  bootstrap_results[[i]] <- tibble(
    Marker = combo$Marker,
    OTU = combo$OTU,
    Site = combo$Site,
    mean_prop = obs_mean,
    boot_lower = quantile(boot_means, 0.025),
    boot_upper = quantile(boot_means, 0.975),
    boot_samples = list(boot_means),
    n_positive_pcrs = sum(site_data$read_proportion > 0),
    n_total_pcrs = nrow(site_data)
  )

  if (i %% pb_interval == 0) cat("    Progress:", round(i / n_combos * 100), "%\n")
}

cat("  Bootstrap complete.\n")
proportion_summary <- bind_rows(bootstrap_results)

cat("  Mean read proportion (all OTU-Sites):",
    round(mean(proportion_summary$mean_prop) * 100, 3), "%\n")
cat("  Median:", round(median(proportion_summary$mean_prop) * 100, 3), "%\n")


# =============================================================================
# PART 3: JOIN POWER ANALYSIS WITH READ PROPORTION
# =============================================================================

cat("\nPART 3: Joining P(detect) with read proportions...\n")

unified_data <- empirical_3x3 %>%
  left_join(
    proportion_summary %>% select(-boot_samples),
    by = c("Marker", "OTU", "Site")
  )

cat("  Matched:", sum(!is.na(unified_data$mean_prop)), "/", nrow(unified_data), "\n")


# =============================================================================
# PART 4: ROC-DERIVED THRESHOLD
# =============================================================================

cat("\nPART 4: ROC analysis for empirical read proportion threshold...\n")

PDETECT_THRESHOLD <- 0.80

# Use GBIF-validated invasive detections as ground truth
invasive_for_roc <- unified_data %>%
  filter(is_invasive) %>%
  left_join(
    gbif_validation %>% select(Species, Site, Lit_assessment),
    by = c("Species", "Site")
  ) %>%
  filter(!is.na(Lit_assessment), !is.na(mean_prop)) %>%
  mutate(well_documented = Lit_assessment == "Well documented")

cat("  Invasive detections with GBIF validation:", nrow(invasive_for_roc), "\n")
cat("  Well documented:", sum(invasive_for_roc$well_documented), "\n")
cat("  Other:", sum(!invasive_for_roc$well_documented), "\n")

# ROC: read proportion as predictor of "Well documented" status
roc_readprop <- roc(invasive_for_roc$well_documented,
                    invasive_for_roc$mean_prop,
                    levels = c(FALSE, TRUE),
                    direction = "<")

# Compare with P(detect) as predictor
roc_pdetect <- roc(invasive_for_roc$well_documented,
                   invasive_for_roc$p_detect_mean,
                   levels = c(FALSE, TRUE),
                   direction = "<")

# Theta alone
roc_theta <- roc(invasive_for_roc$well_documented,
                 invasive_for_roc$theta,
                 levels = c(FALSE, TRUE),
                 direction = "<")

# p alone
roc_p <- roc(invasive_for_roc$well_documented,
             invasive_for_roc$p,
             levels = c(FALSE, TRUE),
             direction = "<")

best_threshold <- coords(roc_readprop, "best",
                         ret = c("threshold", "sensitivity", "specificity"))
READPROP_THRESHOLD <- best_threshold$threshold[1]

cat("\n=== ROC Results ===\n")
roc_summary <- tibble(
  Metric = c("Read proportion", "p (PCR detection)", "P(detect)", "Theta"),
  AUC = round(c(auc(roc_readprop), auc(roc_p),
                auc(roc_pdetect), auc(roc_theta)), 3),
  CI_lower = round(c(ci.auc(roc_readprop)[1], ci.auc(roc_p)[1],
                     ci.auc(roc_pdetect)[1], ci.auc(roc_theta)[1]), 3),
  CI_upper = round(c(ci.auc(roc_readprop)[3], ci.auc(roc_p)[3],
                     ci.auc(roc_pdetect)[3], ci.auc(roc_theta)[3]), 3)
) %>%
  mutate(CI = paste0("(", CI_lower, "-", CI_upper, ")")) %>%
  select(Metric, AUC, CI) %>%
  arrange(-AUC)

print(roc_summary)

cat("\nDeLong test (Read proportion vs P(detect)):\n")
cat("  p-value:", round(roc.test(roc_readprop, roc_pdetect)$p.value, 4), "\n")

cat("\n=== Empirical Threshold ===\n")
cat("  P(detect) threshold:", PDETECT_THRESHOLD * 100, "%\n")
cat("  Read proportion threshold:", round(READPROP_THRESHOLD * 100, 3), "%\n")
cat("  Sensitivity:", round(best_threshold$sensitivity[1], 3), "\n")
cat("  Specificity:", round(best_threshold$specificity[1], 3), "\n")


# =============================================================================
# PART 5: PROBABILISTIC 2×2 CLASSIFICATION
# =============================================================================

cat("\nPART 5: Probabilistic classification...\n")

# Join bootstrap samples for P(above threshold) calculation
unified_with_boot <- unified_data %>%
  left_join(
    proportion_summary %>% select(Marker, OTU, Site, boot_samples),
    by = c("Marker", "OTU", "Site")
  )

cat("  Calculating P(Confident), P(Concerning), P(Patchy), P(Uncertain)...\n")

unified_with_probs <- unified_with_boot %>%
  rowwise() %>%
  mutate(
    # P(high read prop) from bootstrap
    P_high_readprop = if (!is.null(boot_samples)) {
      mean(boot_samples >= READPROP_THRESHOLD)
    } else { NA_real_ },

    # P(high P(detect)) — binary from point estimate
    P_high_pdetect = as.numeric(p_detect_mean >= PDETECT_THRESHOLD),

    # Joint probabilities (assuming independence)
    P_confident = P_high_pdetect * P_high_readprop,
    P_concerning = P_high_pdetect * (1 - P_high_readprop),
    P_patchy = (1 - P_high_pdetect) * P_high_readprop,
    P_uncertain = (1 - P_high_pdetect) * (1 - P_high_readprop),

    # Most likely class
    max_prob = max(c(P_confident, P_concerning, P_patchy, P_uncertain), na.rm = TRUE),
    reliability_class = case_when(
      P_confident == max_prob ~ "Confident",
      P_concerning == max_prob ~ "Concerning",
      P_patchy == max_prob ~ "Patchy but real",
      P_uncertain == max_prob ~ "Uncertain",
      TRUE ~ NA_character_
    ),
    class_confidence = max_prob
  ) %>%
  ungroup() %>%
  select(-boot_samples, -max_prob) %>%
  mutate(
    reliability_class = factor(reliability_class,
                                levels = c("Confident", "Patchy but real",
                                           "Concerning", "Uncertain"))
  )

cat("\n=== Classification Summary (All Detections) ===\n")
unified_with_probs %>%
  filter(!is.na(reliability_class)) %>%
  count(reliability_class) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  print()

# By marker
cat("\n=== Classification by Marker ===\n")
unified_with_probs %>%
  filter(!is.na(reliability_class)) %>%
  group_by(Marker, reliability_class) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Marker) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  print(n = Inf)


# =============================================================================
# PART 6: GBIF VALIDATION
# =============================================================================

cat("\n=============================================================================\n")
cat("PART 6: GBIF Validation\n")
cat("=============================================================================\n")

invasive_unified <- unified_with_probs %>%
  filter(is_invasive) %>%
  left_join(
    gbif_validation %>% select(Species, Site, Lit_assessment),
    by = c("Species", "Site")
  ) %>%
  mutate(well_documented = Lit_assessment == "Well documented")

cat("\n=== GBIF Assessment Distribution ===\n")
print(table(invasive_unified$Lit_assessment, useNA = "ifany"))

# Cross-tabulation
cat("\n=== Cross-tabulation: Unified Class vs Well documented ===\n")
cross_tab <- invasive_unified %>%
  filter(!is.na(Lit_assessment), !is.na(reliability_class)) %>%
  group_by(reliability_class, well_documented) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = well_documented, values_from = n, values_fill = 0) %>%
  rename(Other = `FALSE`, Well_documented = `TRUE`) %>%
  mutate(
    Total = Other + Well_documented,
    PPV = round(Well_documented / Total * 100, 1)
  )
print(cross_tab)

# PPV by class
cat("\n=== Positive Predictive Value by Class ===\n")
ppv_summary <- invasive_unified %>%
  filter(!is.na(Lit_assessment), !is.na(reliability_class)) %>%
  group_by(reliability_class) %>%
  summarise(
    n = n(),
    n_well_documented = sum(well_documented),
    PPV = round(mean(well_documented) * 100, 1),
    mean_class_confidence = round(mean(class_confidence, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(-PPV)
print(ppv_summary)


# =============================================================================
# PART 7: DETAILED TABLES
# =============================================================================

cat("\n=============================================================================\n")
cat("PART 7: DETAILED TABLES\n")
cat("=============================================================================\n")

# Concerning detections (high P(detect), low read proportion)
cat("\n=== Concerning Detections ===\n")
cat("High P(detect) but low read proportion — possible contamination or eDNA drift\n")
concerning <- invasive_unified %>%
  filter(reliability_class == "Concerning") %>%
  arrange(mean_prop) %>%
  mutate(
    read_prop_CI = paste0(round(mean_prop * 100, 3), "% [",
                          round(boot_lower * 100, 3), "-",
                          round(boot_upper * 100, 3), "]")
  ) %>%
  select(Marker, Species, Site, theta, p, read_prop_CI,
         P_high_readprop, Lit_assessment)
print(concerning, n = Inf)

# Patchy but real (low P(detect), high read proportion)
cat("\n=== Patchy but Real ===\n")
cat("Low P(detect) but strong read signal — spatially variable but genuine\n")
patchy <- invasive_unified %>%
  filter(reliability_class == "Patchy but real") %>%
  arrange(-mean_prop) %>%
  mutate(
    read_prop_CI = paste0(round(mean_prop * 100, 2), "% [",
                          round(boot_lower * 100, 2), "-",
                          round(boot_upper * 100, 2), "]")
  ) %>%
  select(Marker, Species, Site, theta, p, read_prop_CI,
         P_high_readprop, Lit_assessment)
print(patchy, n = Inf)

# Potentially novel detections
cat("\n=== Potentially Novel Detections ===\n")
novel <- invasive_unified %>%
  filter(Lit_assessment == "POTENTIALLY NOVEL - nearest record >500km") %>%
  mutate(
    read_prop_CI = paste0(round(mean_prop * 100, 3), "% [",
                          round(boot_lower * 100, 3), "-",
                          round(boot_upper * 100, 3), "]")
  ) %>%
  select(Marker, Species, Site, theta, p, read_prop_CI,
         P_high_readprop, reliability_class, P_confident, P_uncertain)
print(novel)

# Same θ/p demonstration
cat("\n=== Same theta=0.33, p=0.33 — Different read proportions ===\n")
cat("Shows read proportion adds information beyond occupancy parameters\n")
same_pattern <- invasive_unified %>%
  filter(abs(theta - 0.333) < 0.01, abs(p - 0.333) < 0.01) %>%
  arrange(-mean_prop) %>%
  mutate(
    read_prop_CI = paste0(round(mean_prop * 100, 2), "% [",
                          round(boot_lower * 100, 2), "-",
                          round(boot_upper * 100, 2), "]"),
    P_high_rp = round(P_high_readprop, 2)
  ) %>%
  select(Marker, Species, Site, read_prop_CI, P_high_rp,
         reliability_class, Lit_assessment)
print(same_pattern, n = Inf)

# Metrics by GBIF assessment category
cat("\n=== Mean Metrics by GBIF Assessment Category ===\n")
table_by_assessment <- invasive_unified %>%
  filter(!is.na(Lit_assessment)) %>%
  group_by(Lit_assessment) %>%
  summarise(
    n = n(),
    mean_readprop_pct = round(mean(mean_prop) * 100, 3),
    mean_P_high_readprop = round(mean(P_high_readprop, na.rm = TRUE), 3),
    mean_pdetect = round(mean(p_detect_mean) * 100, 1),
    pct_confident = round(mean(reliability_class == "Confident",
                                na.rm = TRUE) * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(-mean_readprop_pct)
print(table_by_assessment)


# =============================================================================
# PART 8: FIGURES
# =============================================================================

cat("\n=============================================================================\n")
cat("PART 8: FIGURES\n")
cat("=============================================================================\n")

class_colors <- c(
  "Confident" = "#1a9850",
  "Patchy but real" = "#91cf60",
  "Concerning" = "#fc8d59",
  "Uncertain" = "#d73027"
)

# ---- Figure 1: 2×2 framework scatter by marker ----
cat("\nFigure 1: 2x2 Framework scatter...\n")

fig_framework <- unified_with_probs %>%
  filter(!is.na(reliability_class)) %>%
  ggplot(aes(x = mean_prop, y = p_detect_mean, color = reliability_class)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_hline(yintercept = PDETECT_THRESHOLD, linetype = "dashed", alpha = 0.7) +
  geom_vline(xintercept = READPROP_THRESHOLD, linetype = "dashed", alpha = 0.7) +
  scale_x_log10(labels = scales::percent, limits = c(0.0001, 1)) +
  scale_y_continuous(labels = scales::percent, limits = c(0.5, 1)) +
  scale_color_manual(values = class_colors, name = "Classification") +
  facet_wrap(~Marker) +
  labs(
    x = "Read proportion (log scale, all PCRs including zeros)",
    y = "P(detect) at 3x3",
    title = "Unified Detection Reliability Framework",
    subtitle = paste0("Thresholds: P(detect) >= 80%, Read proportion >= ",
                      round(READPROP_THRESHOLD * 100, 2), "% (ROC-derived)")
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("Figures/unified_framework_scatter.png", fig_framework,
       width = 14, height = 6)
cat("Saved: Figures/unified_framework_scatter.png\n")

# ---- Figure 2: Invasive species with GBIF validation ----
cat("Figure 2: Invasive species with GBIF...\n")

fig_invasive_unified <- invasive_unified %>%
  filter(!is.na(Lit_assessment), !is.na(reliability_class)) %>%
  mutate(status = ifelse(well_documented, "Well documented", "Other")) %>%
  ggplot(aes(x = mean_prop, y = p_detect_mean,
             color = reliability_class, shape = status)) +
  geom_hline(yintercept = PDETECT_THRESHOLD, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = READPROP_THRESHOLD, linetype = "dashed", alpha = 0.5) +
  geom_errorbarh(aes(xmin = boot_lower, xmax = boot_upper),
                 alpha = 0.3, height = 0) +
  geom_point(size = 3, alpha = 0.8) +
  scale_x_log10(labels = scales::percent) +
  scale_y_continuous(labels = scales::percent, limits = c(0.5, 1)) +
  scale_color_manual(values = class_colors, name = "Classification") +
  scale_shape_manual(values = c("Well documented" = 16, "Other" = 1),
                     name = "GBIF Status") +
  labs(
    x = "Read proportion (log scale) with 95% bootstrap CI",
    y = "P(detect) at 3x3",
    title = "Invasive Species: Unified Classification with GBIF Validation",
    subtitle = "Filled = Well documented; horizontal bars = 95% bootstrap CI"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

ggsave("Figures/invasive_unified_with_gbif.png", fig_invasive_unified,
       width = 12, height = 8)
cat("Saved: Figures/invasive_unified_with_gbif.png\n")

# ---- Figure 3: PPV by classification ----
cat("Figure 3: PPV by classification...\n")

fig_ppv <- ppv_summary %>%
  ggplot(aes(x = reorder(reliability_class, PPV), y = PPV,
             fill = reliability_class)) +
  geom_col(alpha = 0.9) +
  geom_text(aes(label = paste0(PPV, "%\n(n=", n, ")")),
            vjust = -0.2, size = 3.5) +
  scale_y_continuous(limits = c(0, 110),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = class_colors) +
  labs(
    x = "Classification",
    y = "% Well documented (PPV)",
    title = "Validation: PPV by Unified Classification",
    subtitle = "Based on GBIF/literature validation of invasive species"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("Figures/ppv_by_classification.png", fig_ppv,
       width = 10, height = 6)
cat("Saved: Figures/ppv_by_classification.png\n")

# ---- Figure 4: P(high read prop) by GBIF assessment category ----
cat("Figure 4: P(high read prop) by GBIF category...\n")

assessment_order <- c(
  "Well documented",
  "Range edge (100-200km)",
  "Known nearby (<100km)",
  "POTENTIALLY NOVEL - nearest record >500km"
)

fig_assessment <- invasive_unified %>%
  filter(!is.na(Lit_assessment)) %>%
  mutate(Lit_assessment = factor(Lit_assessment, levels = assessment_order)) %>%
  ggplot(aes(x = Lit_assessment, y = P_high_readprop,
             fill = Lit_assessment)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_viridis_d(option = "plasma", direction = -1) +
  labs(
    x = "GBIF/Literature Assessment",
    y = "P(above read proportion threshold)",
    title = "Probability of exceeding signal threshold by assessment category"
  ) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1))

ggsave("Figures/P_high_readprop_by_category.png", fig_assessment,
       width = 12, height = 7)
cat("Saved: Figures/P_high_readprop_by_category.png\n")

# ---- Figure 5: N. melanostomus case study ----
cat("Figure 5: N. melanostomus case study...\n")

neogobius <- invasive_unified %>%
  filter(Species == "Neogobius melanostomus")

if (nrow(neogobius) > 0) {
  fig_neogobius <- neogobius %>%
    ggplot(aes(x = reorder(Site, -mean_prop), y = mean_prop,
               fill = reliability_class)) +
    geom_col(alpha = 0.9) +
    geom_errorbar(aes(ymin = boot_lower, ymax = boot_upper), width = 0.3) +
    geom_hline(yintercept = READPROP_THRESHOLD, linetype = "dashed",
               color = "red") +
    geom_text(aes(label = paste0("P=", round(P_high_readprop, 2))),
              vjust = -0.8, size = 2.5) +
    scale_y_continuous(labels = scales::percent,
                       expand = expansion(mult = c(0, 0.2))) +
    scale_fill_manual(values = class_colors, name = "Classification") +
    labs(
      x = "Site",
      y = "Read proportion with 95% bootstrap CI",
      title = expression(italic("Neogobius melanostomus") *
                           ": Same " * theta * "/p, different signal strength"),
      subtitle = "P = P(above threshold); dashed = ROC-derived threshold"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave("Figures/neogobius_unified.png", fig_neogobius,
         width = 12, height = 6)
  cat("Saved: Figures/neogobius_unified.png\n")
}


# =============================================================================
# PART 9: SAVE OUTPUTS
# =============================================================================

cat("\n=============================================================================\n")
cat("SAVING OUTPUTS\n")
cat("=============================================================================\n")

write.csv(unified_with_probs %>% select(-any_of("boot_samples")),
          "Processed_data/unified_detection_reliability.csv",
          row.names = FALSE)

write.csv(invasive_unified %>% select(-any_of("boot_samples")),
          "Processed_data/invasive_unified_reliability.csv",
          row.names = FALSE)

write.csv(ppv_summary,
          "Processed_data/ppv_by_classification.csv",
          row.names = FALSE)

thresholds <- tibble(
  threshold_name = c("P(detect)", "Read proportion"),
  value = c(PDETECT_THRESHOLD, READPROP_THRESHOLD),
  value_pct = c(paste0(PDETECT_THRESHOLD * 100, "%"),
                paste0(round(READPROP_THRESHOLD * 100, 3), "%")),
  source = c("Conventional (literature)",
             "Empirical (ROC with GBIF validation)")
)

write.csv(thresholds,
          "Processed_data/unified_thresholds.csv",
          row.names = FALSE)

cat("Outputs saved to Processed_data/\n")


# =============================================================================
# SUMMARY
# =============================================================================

cat("\n", strrep("=", 70), "\n")
cat("SUMMARY\n")
cat(strrep("=", 70), "\n\n")

cat("Thresholds:\n")
cat("  P(detect) >=", PDETECT_THRESHOLD * 100, "%\n")
cat("  Read proportion >=", round(READPROP_THRESHOLD * 100, 3),
    "% (ROC-derived from GBIF)\n\n")

cat("ROC performance (AUC for predicting Well documented status):\n")
print(roc_summary)

cat("\nClassification PPV:\n")
print(ppv_summary)

cat("\nKey insight: Read proportion adds information beyond theta/p.\n")
cat("Species with identical theta=0.33, p=0.33 can have very different\n")
cat("read proportions, and this difference predicts GBIF validation status.\n")
cat("The 'Concerning' category (high P(detect) but low read proportion)\n")
cat("flags detections that warrant closer scrutiny.\n\n")

cat(strrep("=", 70), "\n")
cat("03 COMPLETE\n")
cat(strrep("=", 70), "\n")
