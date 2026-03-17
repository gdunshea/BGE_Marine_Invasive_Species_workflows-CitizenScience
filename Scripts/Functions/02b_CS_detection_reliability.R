# =============================================================================
# 02b_CS_detection_reliability.R
# =============================================================================
# Run after 02_CS_detection_analysis_RESTRUCTURED.R + lifestyle assignment.
# Expects: combined_all (with lifestyle_final), empirical_3x3, 
#          det_12S, det_18S, det_COI
# =============================================================================

source("Scripts/Final/FUNC_detection_reliability.R")

det_list <- list(`12S` = det_12S, `18S` = det_18S, COI = det_COI)

# =============================================================================
# CONFIGURABLE PARAMETERS
# =============================================================================

n_perm <- 4999           # Permutations for statistical tests
n_boot <- 5000           # Bootstrap resamples for lifestyle-level CIs
set.seed(42)             # Reproducibility

# =============================================================================
# SECTIONS 1-4: Detection reliability (unchanged)
# =============================================================================

# ---- 1. Detection strength ----
str_result <- analyse_detection_strength(
  combined_all$site_level,
  save_csv = "Processed_data/detection_pattern_summary.csv",
  save_fig = "Figures/detection_pattern_distribution.png"
)

# ---- 2. Filtering trade-offs ----
filt_result <- analyse_filtering_tradeoffs(
  str_result$classified,
  inspect_threshold = "bio2",
  save_csv = "Processed_data/filtering_tradeoffs.csv",
  save_fig = "Figures/filtering_tradeoff.png"
)

# ---- 3. Power by detection strength ----
pow_result <- analyse_power_by_strength(
  empirical_3x3, det_list,
  p_threshold = 0.80,
  save_csv = "Processed_data/power_by_detection_strength.csv",
  save_fig = "Figures/power_by_detection_strength.png"
)

# ---- 4. Invasive reliability + 5x3 comparison ----
inv_result <- analyse_invasive_reliability(
  empirical_3x3, det_list, combined_all,
  compare_design = "5\u00d73",
  save_csv = "Processed_data/invasive_reliability_detail.csv",
  save_fig = "Figures/invasive_forest_3x3_vs_5x3.png"
)


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Bootstrap CI for a group-level statistic
#' Resamples detections within each group and computes the statistic.
#' Returns mean, 2.5%, 97.5% quantiles.
boot_ci <- function(x, stat_fn = mean, n_boot = 5000) {
  boot_vals <- replicate(n_boot, stat_fn(sample(x, replace = TRUE)))
  c(estimate = stat_fn(x),
    lower = unname(quantile(boot_vals, 0.025)),
    upper = unname(quantile(boot_vals, 0.975)))
}

#' Permutation test for difference in means between groups
#' Returns observed F-statistic analogue and permutation p-value.
#' More conservative than Kruskal-Wallis for bounded, tied data.
perm_test_groups <- function(values, groups, n_perm = 9999) {
  groups <- as.factor(groups)
  
  # Observed: ratio of between-group to within-group variance (F-like)
  obs_stat <- summary(aov(values ~ groups))[[1]][1, "F value"]
  
  # Permutation distribution
  perm_stats <- replicate(n_perm, {
    perm_groups <- sample(groups)
    summary(aov(values ~ perm_groups))[[1]][1, "F value"]
  })
  
  p_value <- (sum(perm_stats >= obs_stat) + 1) / (n_perm + 1)
  
  list(statistic = obs_stat, p.value = p_value, n_perm = n_perm)
}

#' Pairwise permutation tests with correction
perm_test_pairwise <- function(values, groups, n_perm = 9999, 
                                p_adjust = "bonferroni") {
  groups <- as.factor(droplevels(groups))
  levs <- levels(groups)
  n_pairs <- choose(length(levs), 2)
  
  results <- list()
  for (i in seq_along(levs)) {
    for (j in seq_len(i - 1)) {
      idx <- groups %in% c(levs[i], levs[j])
      v <- values[idx]
      g <- droplevels(groups[idx])
      
      # Observed difference in means
      obs_diff <- abs(diff(tapply(v, g, mean)))
      
      # Permutation
      perm_diffs <- replicate(n_perm, {
        pg <- sample(g)
        abs(diff(tapply(v, pg, mean)))
      })
      
      p_raw <- (sum(perm_diffs >= obs_diff) + 1) / (n_perm + 1)
      results[[paste(levs[i], levs[j], sep = " vs ")]] <- p_raw
    }
  }
  
  p_raw_vec <- unlist(results)
  p_adj_vec <- p.adjust(p_raw_vec, method = p_adjust)
  
  tibble(
    comparison = names(results),
    p_raw = round(p_raw_vec, 4),
    p_adjusted = round(p_adj_vec, 4)
  )
}


# =============================================================================
# SECTION 5: DETECTION PARAMETERS BY LIFESTYLE (WITHIN-MARKER)
# =============================================================================

cat("\n\n", strrep("=", 70), "\n")
cat("5. DETECTION PARAMETERS BY LIFESTYLE (WITHIN-MARKER)\n")
cat(strrep("=", 70), "\n\n")

# Exclude Unclassified (no taxonomy) and Other (terrestrial contaminants)
marine_data <- combined_all$site_level %>%
  filter(!lifestyle_final %in% c("Other", "Unclassified"))

cat("Marine detections (excl. Unclassified/Other):", nrow(marine_data), 
    "of", nrow(combined_all$site_level), "\n")
cat("\nNOTE: 12S is a vertebrate-targeted marker. Non-fish detections\n")
cat("are off-target and excluded from lifestyle comparisons.\n")
cat("Lifestyle analysis uses 18S and COI only.\n\n")

# Markers for lifestyle analysis (exclude 12S)
lifestyle_markers <- c("18S", "COI")

# ---- 5a. Detection parameters by lifestyle × marker ----
lifestyle_detection <- marine_data %>%
  filter(Marker %in% lifestyle_markers) %>%
  group_by(Marker, lifestyle_final) %>%
  summarise(
    n = n(),
    n_species = n_distinct(Species, na.rm = TRUE),
    median_theta = round(median(theta), 3),
    median_p = round(median(mean_p_pcr), 3),
    mean_theta = round(mean(theta), 3),
    mean_p = round(mean(mean_p_pcr), 3),
    pct_theta_1 = round(mean(theta == 1) * 100, 1),
    pct_p_1 = round(mean(mean_p_pcr == 1) * 100, 1),
    pct_minimal = round(mean(n_pcr_positive == 1 & 
                               n_biosamples_positive == 1) * 100, 1),
    median_reads = median(total_reads),
    .groups = "drop"
  ) %>%
  arrange(Marker, desc(n))

cat("=== Detection parameters by lifestyle (within marker) ===\n")
print(lifestyle_detection, n = Inf, width = Inf)

write.csv(lifestyle_detection, 
          "Processed_data/detection_by_lifestyle_marker.csv", row.names = FALSE)

# ---- 5b. Bootstrap CIs on theta by lifestyle × marker ----
cat("\n=== Bootstrap 95% CIs on mean theta by lifestyle (within marker) ===\n")

theta_boot <- marine_data %>%
  filter(Marker %in% lifestyle_markers) %>%
  group_by(Marker, lifestyle_final) %>%
  summarise(
    n = n(),
    boot = list(boot_ci(theta, mean, n_boot)),
    .groups = "drop"
  ) %>%
  mutate(
    mean_theta = sapply(boot, `[`, "estimate"),
    theta_lower = sapply(boot, `[`, "lower"),
    theta_upper = sapply(boot, `[`, "upper")
  ) %>%
  select(-boot) %>%
  mutate(across(c(mean_theta, theta_lower, theta_upper), ~round(., 3)))

print(theta_boot, n = Inf, width = Inf)

# ---- 5c. Bootstrap CIs on p by lifestyle × marker ----
cat("\n=== Bootstrap 95% CIs on mean p by lifestyle (within marker) ===\n")

p_boot <- marine_data %>%
  filter(Marker %in% lifestyle_markers) %>%
  group_by(Marker, lifestyle_final) %>%
  summarise(
    n = n(),
    boot = list(boot_ci(mean_p_pcr, mean, n_boot)),
    .groups = "drop"
  ) %>%
  mutate(
    mean_p = sapply(boot, `[`, "estimate"),
    p_lower = sapply(boot, `[`, "lower"),
    p_upper = sapply(boot, `[`, "upper")
  ) %>%
  select(-boot) %>%
  mutate(across(c(mean_p, p_lower, p_upper), ~round(., 3)))

print(p_boot, n = Inf, width = Inf)

# ---- 5d. Figure: theta with bootstrap CIs by lifestyle (within marker) ----
fig_theta_lifestyle <- theta_boot %>%
  ggplot(aes(x = reorder(lifestyle_final, -mean_theta), 
             y = mean_theta, ymin = theta_lower, ymax = theta_upper,
             color = Marker)) +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.7) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_brewer(palette = "Set1") +
  labs(
    x = NULL, 
    y = expression("Mean " * theta * " (95% bootstrap CI)"),
    title = expression("Biosample detection rate (" * theta * ") by organism lifestyle"),
    subtitle = "Within-marker comparison (18S and COI only)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

ggsave("Figures/lifestyle_theta_bootstrap.png", fig_theta_lifestyle,
       width = 11, height = 6)
cat("\nSaved: Figures/lifestyle_theta_bootstrap.png\n")

# ---- 5e. Figure: p with bootstrap CIs by lifestyle (within marker) ----
fig_p_lifestyle <- p_boot %>%
  ggplot(aes(x = reorder(lifestyle_final, -mean_p), 
             y = mean_p, ymin = p_lower, ymax = p_upper,
             color = Marker)) +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.7) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_brewer(palette = "Set1") +
  labs(
    x = NULL, 
    y = "Mean p (95% bootstrap CI)",
    title = "PCR detection rate (p) by organism lifestyle",
    subtitle = "Within-marker comparison (18S and COI only)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

ggsave("Figures/lifestyle_p_bootstrap.png", fig_p_lifestyle,
       width = 11, height = 6)
cat("Saved: Figures/lifestyle_p_bootstrap.png\n")

# ---- 5f. Figure: detection strength composition by lifestyle (within marker) ----
marine_classified <- str_result$classified %>%
  filter(!lifestyle_final %in% c("Other", "Unclassified"),
         Marker %in% lifestyle_markers)

lifestyle_strength_summary <- marine_classified %>%
  group_by(Marker, lifestyle_final, detection_strength) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Marker, lifestyle_final) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  ungroup()

fig_lifestyle_strength <- ggplot(
  lifestyle_strength_summary,
  aes(x = lifestyle_final, y = pct / 100, fill = detection_strength)
) +
  geom_col(alpha = 0.85, width = 0.8) +
  facet_wrap(~Marker) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(name = "Detection\nstrength", direction = -1) +
  labs(
    x = NULL, y = "Proportion of detections",
    title = "Detection strength composition by organism lifestyle",
    subtitle = "Within-marker comparison (18S and COI)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

ggsave("Figures/lifestyle_detection_strength.png", fig_lifestyle_strength,
       width = 14, height = 6)
cat("Saved: Figures/lifestyle_detection_strength.png\n")

write.csv(lifestyle_strength_summary, 
          "Processed_data/lifestyle_detection_strength.csv", row.names = FALSE)


# =============================================================================
# SECTION 6: BOOTSTRAP P(DETECT) BY LIFESTYLE (WITHIN-MARKER)
# =============================================================================

cat("\n\n", strrep("=", 70), "\n")
cat("6. BOOTSTRAP P(DETECT) BY LIFESTYLE (WITHIN-MARKER)\n")
cat(strrep("=", 70), "\n\n")

# Add lifestyle to empirical_3x3
emp_lifestyle <- empirical_3x3 %>%
  left_join(
    combined_all$site_level %>%
      select(OTU, Site, Marker, lifestyle_final) %>%
      distinct(),
    by = c("OTU", "Site", "Marker")
  ) %>%
  filter(!lifestyle_final %in% c("Other", "Unclassified", NA),
         Marker %in% lifestyle_markers)

# ---- 6a. P(detect) at 3x3 with bootstrap CIs by lifestyle × marker ----
power_by_lifestyle <- emp_lifestyle %>%
  group_by(Marker, lifestyle_final) %>%
  summarise(
    n = n(),
    boot_mean = list(boot_ci(p_detect_mean, mean, n_boot)),
    boot_pct80 = list(boot_ci(as.numeric(p_detect_mean >= 0.80), mean, n_boot)),
    .groups = "drop"
  ) %>%
  mutate(
    mean_p_detect = round(sapply(boot_mean, `[`, "estimate"), 3),
    pd_lower = round(sapply(boot_mean, `[`, "lower"), 3),
    pd_upper = round(sapply(boot_mean, `[`, "upper"), 3),
    pct_above_80 = round(sapply(boot_pct80, `[`, "estimate") * 100, 1),
    pct80_lower = round(sapply(boot_pct80, `[`, "lower") * 100, 1),
    pct80_upper = round(sapply(boot_pct80, `[`, "upper") * 100, 1)
  ) %>%
  select(-boot_mean, -boot_pct80) %>%
  arrange(Marker, desc(mean_p_detect))

cat("=== P(detect) at 3x3 by lifestyle with bootstrap 95% CIs ===\n")
print(power_by_lifestyle, n = Inf, width = Inf)

write.csv(power_by_lifestyle,
          "Processed_data/power_by_lifestyle_marker.csv", row.names = FALSE)

# ---- 6b. Design comparison by lifestyle (within marker) ----
all_designs_lifestyle <- map2_dfr(det_list, names(det_list), function(det, marker) {
  det$empirical_with_power %>% mutate(Marker = marker)
}) %>%
  left_join(
    combined_all$site_level %>%
      select(OTU, Site, Marker, lifestyle_final) %>%
      distinct(),
    by = c("OTU", "Site", "Marker")
  ) %>%
  filter(!lifestyle_final %in% c("Other", "Unclassified", NA),
         Marker %in% lifestyle_markers)

design_by_lifestyle <- all_designs_lifestyle %>%
  filter(design %in% c("3\u00d73", "5\u00d73", "8\u00d73")) %>%
  group_by(Marker, lifestyle_final, design, effort) %>%
  summarise(
    n = n(),
    boot_pct80 = list(boot_ci(as.numeric(p_detect_mean >= 0.80), mean, n_boot)),
    mean_p_detect = round(mean(p_detect_mean), 3),
    .groups = "drop"
  ) %>%
  mutate(
    pct_above_80 = round(sapply(boot_pct80, `[`, "estimate") * 100, 1),
    pct80_lower = round(sapply(boot_pct80, `[`, "lower") * 100, 1),
    pct80_upper = round(sapply(boot_pct80, `[`, "upper") * 100, 1)
  ) %>%
  select(-boot_pct80) %>%
  arrange(Marker, lifestyle_final, effort)

cat("\n=== Design comparison by lifestyle (within marker) ===\n")
print(design_by_lifestyle, n = Inf, width = Inf)

write.csv(design_by_lifestyle,
          "Processed_data/design_by_lifestyle.csv", row.names = FALSE)

# ---- 6c. Figure: mean P(detect) with bootstrap CIs by lifestyle ----
fig_pdetect_lifestyle <- power_by_lifestyle %>%
  ggplot(aes(x = reorder(lifestyle_final, -mean_p_detect),
             y = mean_p_detect, ymin = pd_lower, ymax = pd_upper,
             color = Marker)) +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.7) +
  geom_hline(yintercept = 0.80, linetype = "dashed", alpha = 0.4) +
  scale_y_continuous(limits = c(0.5, 1), labels = scales::percent) +
  scale_color_brewer(palette = "Set1") +
  labs(
    x = NULL,
    y = "Mean P(detect) at 3x3 (95% bootstrap CI)",
    title = "Detection probability by organism lifestyle",
    subtitle = "Within-marker comparison (18S and COI)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

ggsave("Figures/lifestyle_pdetect_bootstrap.png", fig_pdetect_lifestyle,
       width = 11, height = 6)
cat("\nSaved: Figures/lifestyle_pdetect_bootstrap.png\n")

# ---- 6d. Figure: design performance by lifestyle with CIs ----
fig_design_lifestyle <- design_by_lifestyle %>%
  ggplot(aes(x = design, y = pct_above_80 / 100, 
             ymin = pct80_lower / 100, ymax = pct80_upper / 100,
             fill = lifestyle_final)) +
  geom_col(position = position_dodge(width = 0.8), alpha = 0.7, width = 0.7) +
  geom_errorbar(position = position_dodge(width = 0.8), width = 0.2, 
                alpha = 0.6, linewidth = 0.4) +
  geom_hline(yintercept = 0.80, linetype = "dashed", alpha = 0.4) +
  facet_wrap(~Marker) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.05)) +
  scale_fill_brewer(palette = "Set2", name = "Lifestyle") +
  labs(
    x = "Sampling design",
    y = "% of detections with P(detect) >= 80%",
    title = "Design performance by organism lifestyle",
    subtitle = "With 95% bootstrap CIs (within-marker)"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("Figures/design_performance_by_lifestyle.png", fig_design_lifestyle,
       width = 14, height = 7)
cat("Saved: Figures/design_performance_by_lifestyle.png\n")

# ---- 6e. Figure: P(detect) violin/boxplot by lifestyle (within marker) ----
fig_pdetect_violin <- emp_lifestyle %>%
  ggplot(aes(x = reorder(lifestyle_final, -p_detect_mean), 
             y = p_detect_mean, fill = lifestyle_final)) +
  geom_violin(alpha = 0.4, linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.7) +
  geom_hline(yintercept = 0.80, linetype = "dashed", alpha = 0.4) +
  facet_wrap(~Marker) +
  scale_y_continuous(limits = c(0.3, 1), labels = scales::percent) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(
    x = NULL, y = "P(detect) at 3x3 design",
    title = "Distribution of detection probability by lifestyle",
    subtitle = "Violin + boxplot showing full bootstrap P(detect) distribution per detection"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

ggsave("Figures/lifestyle_pdetect_violin.png", fig_pdetect_violin,
       width = 14, height = 7)
cat("Saved: Figures/lifestyle_pdetect_violin.png\n")


# =============================================================================
# SECTION 7: PERMUTATION TESTS — LIFESTYLE EFFECT ON DETECTION
# =============================================================================

cat("\n\n", strrep("=", 70), "\n")
cat("7. PERMUTATION TESTS: LIFESTYLE EFFECT ON DETECTION\n")
cat(strrep("=", 70), "\n")
cat("(", n_perm, "permutations per test)\n\n")

# ---- 7a. Global permutation test: theta by lifestyle (per marker) ----
cat("=== Permutation test: theta ~ lifestyle (global) ===\n")

for (m in lifestyle_markers) {
  dat <- marine_data %>% filter(Marker == m)
  pt <- perm_test_groups(dat$theta, dat$lifestyle_final, n_perm)
  cat(m, ": F =", round(pt$statistic, 2), 
      ", p =", formatC(pt$p.value, format = "e", digits = 2), 
      " (", n_perm, " permutations)\n")
}

# ---- 7b. Global permutation test: p by lifestyle (per marker) ----
cat("\n=== Permutation test: p_pcr ~ lifestyle (global) ===\n")

for (m in lifestyle_markers) {
  dat <- marine_data %>% filter(Marker == m)
  pt <- perm_test_groups(dat$mean_p_pcr, dat$lifestyle_final, n_perm)
  cat(m, ": F =", round(pt$statistic, 2), 
      ", p =", formatC(pt$p.value, format = "e", digits = 2),
      " (", n_perm, " permutations)\n")
}

# ---- 7c. Global permutation test: P(detect) by lifestyle ----
cat("\n=== Permutation test: P(detect) ~ lifestyle (global) ===\n")

for (m in lifestyle_markers) {
  dat <- emp_lifestyle %>% filter(Marker == m)
  pt <- perm_test_groups(dat$p_detect_mean, dat$lifestyle_final, n_perm)
  cat(m, ": F =", round(pt$statistic, 2), 
      ", p =", formatC(pt$p.value, format = "e", digits = 2),
      " (", n_perm, " permutations)\n")
}

# ---- 7d. Pairwise permutation tests for theta ----
cat("\n=== Pairwise permutation tests: theta (Bonferroni-corrected) ===\n")

key_lifestyles <- c("Fish", "Sessile invertebrate", "Phytoplankton", 
                     "Zooplankton", "Macroalgae", "Microeukaryote",
                     "Mobile invertebrate")

for (m in lifestyle_markers) {
  dat <- marine_data %>% 
    filter(Marker == m, lifestyle_final %in% key_lifestyles) %>%
    mutate(lifestyle_final = droplevels(factor(lifestyle_final)))
  
  # Only test groups with n >= 10
  grp_n <- table(dat$lifestyle_final)
  keep_grps <- names(grp_n[grp_n >= 10])
  dat <- dat %>% filter(lifestyle_final %in% keep_grps) %>%
    mutate(lifestyle_final = droplevels(lifestyle_final))
  
  if (length(unique(dat$lifestyle_final)) < 2) next
  
  cat("\n", m, "- theta (groups with n >= 10):\n")
  pw <- perm_test_pairwise(dat$theta, dat$lifestyle_final, 
                            n_perm = n_perm, p_adjust = "bonferroni")
  print(pw, n = Inf)
}

# ---- 7e. Pairwise permutation tests for P(detect) ----
cat("\n=== Pairwise permutation tests: P(detect) (Bonferroni-corrected) ===\n")

for (m in lifestyle_markers) {
  dat <- emp_lifestyle %>% 
    filter(Marker == m, lifestyle_final %in% key_lifestyles) %>%
    mutate(lifestyle_final = droplevels(factor(lifestyle_final)))
  
  grp_n <- table(dat$lifestyle_final)
  keep_grps <- names(grp_n[grp_n >= 10])
  dat <- dat %>% filter(lifestyle_final %in% keep_grps) %>%
    mutate(lifestyle_final = droplevels(lifestyle_final))
  
  if (length(unique(dat$lifestyle_final)) < 2) next
  
  cat("\n", m, "- P(detect) (groups with n >= 10):\n")
  pw <- perm_test_pairwise(dat$p_detect_mean, dat$lifestyle_final, 
                            n_perm = n_perm, p_adjust = "bonferroni")
  print(pw, n = Inf)
}


# =============================================================================
# SECTION 8: FILTERING TRADE-OFFS BY LIFESTYLE (WITHIN-MARKER)
# =============================================================================

cat("\n\n", strrep("=", 70), "\n")
cat("8. FILTERING IMPACT BY LIFESTYLE (WITHIN-MARKER)\n")
cat(strrep("=", 70), "\n\n")

filter_by_lifestyle <- marine_classified %>%
  group_by(Marker, lifestyle_final) %>%
  summarise(
    n_total = n(),
    n_pass_pcr2 = sum(n_pcr_positive >= 2),
    pct_pass_pcr2 = round(n_pass_pcr2 / n_total * 100, 1),
    n_pass_bio2 = sum(n_biosamples_positive >= 2),
    pct_pass_bio2 = round(n_pass_bio2 / n_total * 100, 1),
    n_minimal = sum(n_pcr_positive == 1 & n_biosamples_positive == 1),
    pct_minimal = round(n_minimal / n_total * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(Marker, desc(pct_minimal))

cat("=== Filtering survival by lifestyle (within marker) ===\n")
print(filter_by_lifestyle, n = Inf, width = Inf)

write.csv(filter_by_lifestyle,
          "Processed_data/filtering_by_lifestyle.csv", row.names = FALSE)

# Figure
fig_filter_lifestyle <- filter_by_lifestyle %>%
  filter(Marker %in% lifestyle_markers) %>%
  select(Marker, lifestyle_final, pct_pass_pcr2, pct_pass_bio2) %>%
  pivot_longer(c(pct_pass_pcr2, pct_pass_bio2),
               names_to = "threshold", values_to = "pct") %>%
  mutate(threshold = ifelse(threshold == "pct_pass_pcr2",
                            ">=2 PCRs", ">=2 biosamples")) %>%
  ggplot(aes(x = reorder(lifestyle_final, -pct), y = pct / 100, 
             fill = threshold)) +
  geom_col(position = "dodge", alpha = 0.85, width = 0.7) +
  geom_hline(yintercept = 0.5, linetype = "dotted", alpha = 0.4) +
  facet_wrap(~Marker) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.05)) +
  scale_fill_manual(values = c(">=2 PCRs" = "steelblue",
                                ">=2 biosamples" = "tomato"), name = NULL) +
  labs(
    x = NULL, y = "% of detections retained",
    title = "Confirmation threshold impact by organism lifestyle",
    subtitle = "Within-marker comparison (18S and COI)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "bottom")

ggsave("Figures/filtering_by_lifestyle.png", fig_filter_lifestyle,
       width = 14, height = 6)
cat("Saved: Figures/filtering_by_lifestyle.png\n")


# =============================================================================
# SECTION 9: SUMMARY
# =============================================================================

cat("\n\n", strrep("=", 70), "\n")
cat("SUMMARY\n")
cat(strrep("=", 70), "\n\n")

cat("Sections 1-4: Standard detection reliability characterisation\n")
cat("  - Detection strength tiers, filtering trade-offs, power analysis\n")
cat("  - Invasive species reliability with 3x3 vs 5x3 comparison\n\n")

cat("Sections 5-8: Lifestyle-stratified analysis (18S and COI only)\n")
cat("  - 12S excluded from lifestyle comparisons (vertebrate-targeted marker)\n")
cat("  - All estimates include bootstrap 95% CIs (", n_boot, " resamples)\n")
cat("  - Statistical tests are permutation-based (", n_perm, " permutations)\n\n")

# Key findings per marker
for (m in lifestyle_markers) {
  cat("---", m, "---\n")
  
  pw_m <- power_by_lifestyle %>% filter(Marker == m)
  best <- pw_m %>% slice_max(mean_p_detect, n = 2)
  worst <- pw_m %>% slice_min(mean_p_detect, n = 1)
  
  cat("  Highest P(detect):", 
      paste(best$lifestyle_final, "(", best$mean_p_detect, 
            "[", best$pd_lower, "-", best$pd_upper, "])", 
            collapse = ", "), "\n")
  cat("  Lowest P(detect):", 
      paste(worst$lifestyle_final, "(", worst$mean_p_detect,
            "[", worst$pd_lower, "-", worst$pd_upper, "])", 
            collapse = ", "), "\n\n")
}

cat("CAVEATS:\n")
cat("  1. Analysis characterises P(re-detect) for detected species only.\n")
cat("     False negative rate for undetected species cannot be estimated.\n")
cat("  2. Lifestyle differences reflect both biological (eDNA shedding,\n")
cat("     local abundance) and methodological (primer affinity) factors.\n")
cat("  3. Within-marker comparisons control for primer/PCR effects but\n")
cat("     taxonomic composition differs between markers by design.\n\n")

cat(strrep("=", 70), "\n")
cat("02b COMPLETE\n")
cat(strrep("=", 70), "\n")

