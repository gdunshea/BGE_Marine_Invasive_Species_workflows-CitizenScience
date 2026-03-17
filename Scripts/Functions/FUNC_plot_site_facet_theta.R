# =============================================================================
# Plot detection by theta bins (1/3, 2/3, 3/3 biosamples)
# =============================================================================

library(ggplot2)
library(dplyr)

#' Plot reads vs p_empirical, faceted by theta (biosample detection)
#' 
#' @param combined_markers Output from combine_markers()
#' @param site Site name
#' @return ggplot object
plot_site_facet_theta <- function(combined_markers, site) {
  
  site_data <- combined_markers$site_level %>%
    filter(Site == !!site)
  
  if (nrow(site_data) == 0) {
    message("No data for site: ", site)
    return(NULL)
  }
  
  country <- unique(site_data$Country)[1]
  
  # Create theta bin labels
  site_data <- site_data %>%
    mutate(
      invasive_label = ifelse(is_invasive, "Invasive", "Non-invasive"),
      label = ifelse(is_invasive, Species, NA_character_),
      # Convert theta to fraction label
      theta_bin = case_when(
        theta <= 0.34 ~ "1/3 biosamples",
        theta <= 0.67 ~ "2/3 biosamples",
        TRUE ~ "3/3 biosamples"
      ),
      theta_bin = factor(theta_bin, levels = c("1/3 biosamples", "2/3 biosamples", "3/3 biosamples"))
    ) %>%
    arrange(is_invasive)  # Plot invasive on top
  
  n_otus <- nrow(site_data)
  n_invasive <- sum(site_data$is_invasive)
  n_markers <- n_distinct(site_data$Marker)
  
  p <- ggplot(site_data, aes(x = total_reads, y = p_empirical, 
                              color = invasive_label, shape = Marker)) +
    geom_point(aes(alpha = invasive_label), size = 3) +
    facet_wrap(~theta_bin, ncol = 3) +
    scale_x_log10(labels = scales::comma) +
    scale_color_manual(values = c("Invasive" = "#e74c3c", "Non-invasive" = "#3498db")) +
    scale_alpha_manual(values = c("Invasive" = 1, "Non-invasive" = 0.15), guide = "none") +
    labs(x = "Read count (log scale)",
         y = "Detection rate (p_empirical)",
         title = site,
         subtitle = paste0(country, " | ", n_otus, " OTUs, ", n_invasive, " invasive, ",
                           n_markers, " markers"),
         color = "", shape = "Marker") +
    theme_minimal() +
    theme(strip.text = element_text(face = "bold", size = 11)) +
    coord_cartesian(ylim = c(0, 1))
  
  # Add labels for invasive species
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

#' Plot all sites in a country, faceted by theta
#' Each site gets one row (3 panels for 1/3, 2/3, 3/3 biosamples)
#' 
#' @param combined_markers Output from combine_markers()
#' @param country Country name
#' @param sites_per_page Number of sites (rows) per page (default 4)
#' @param output_dir Directory to save PDFs (NULL to return plots without saving)
#' @return List of ggplot objects (one per page)
plot_country_facet_theta <- function(combined_markers, country, 
                                      sites_per_page = 4, 
                                      output_dir = NULL) {
  
  # Get all sites in this country
  country_data <- combined_markers$site_level %>%
    filter(Country == !!country)
  
  if (nrow(country_data) == 0) {
    message("No data for country: ", country)
    return(NULL)
  }
  
  sites <- sort(unique(country_data$Site))
  n_sites <- length(sites)
  n_pages <- ceiling(n_sites / sites_per_page)
  
  cat("Country:", country, "\n")
  cat("Sites:", n_sites, "\n")
  cat("Pages:", n_pages, "\n\n")
  
  # Prepare data with theta bins and short site names
  country_data <- country_data %>%
    mutate(
      invasive_label = ifelse(is_invasive, "Invasive", "Non-invasive"),
      label = ifelse(is_invasive, Species, NA_character_),
      theta_bin = case_when(
        theta <= 0.34 ~ "1/3 biosamples",
        theta <= 0.67 ~ "2/3 biosamples",
        TRUE ~ "3/3 biosamples"
      ),
      theta_bin = factor(theta_bin, levels = c("1/3 biosamples", "2/3 biosamples", "3/3 biosamples")),
      # Create shorter site label (remove country prefix)
      Site_short = gsub(paste0("^", country, ", ?"), "", Site)
    ) %>%
    arrange(is_invasive)
  
  plots <- list()
  
  for (page in 1:n_pages) {
    # Get sites for this page
    start_idx <- (page - 1) * sites_per_page + 1
    end_idx <- min(page * sites_per_page, n_sites)
    page_sites <- sites[start_idx:end_idx]
    n_sites_this_page <- length(page_sites)
    
    # Filter data for these sites
    page_data <- country_data %>%
      filter(Site %in% page_sites) %>%
      mutate(Site_short = factor(Site_short, levels = gsub(paste0("^", country, ", ?"), "", page_sites)))
    
    # Create facet grid: Site (rows) x theta_bin (columns)
    p <- ggplot(page_data, aes(x = total_reads, y = p_empirical, 
                                color = invasive_label, shape = Marker)) +
      geom_point(aes(alpha = invasive_label), size = 2.5) +
      facet_grid(Site_short ~ theta_bin, scales = "free_x") +
      scale_x_log10(labels = scales::comma) +
      scale_color_manual(values = c("Invasive" = "#e74c3c", "Non-invasive" = "#3498db")) +
      scale_alpha_manual(values = c("Invasive" = 1, "Non-invasive" = 0.15), guide = "none") +
      labs(x = "Read count (log scale)",
           y = "Detection rate (p_empirical)",
           title = country,
           subtitle = paste0("Page ", page, "/", n_pages, " | Sites ", start_idx, "-", end_idx, " of ", n_sites),
           color = "", shape = "Marker") +
      theme_minimal() +
      theme(
        strip.text = element_text(face = "bold", size = 9),
        strip.text.y = element_text(angle = 0, hjust = 0),
        panel.spacing = unit(0.3, "lines")
      ) +
      coord_cartesian(ylim = c(0, 1))
    
    # Add labels for invasive species
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      inv_data <- page_data %>% filter(is_invasive)
      if (nrow(inv_data) > 0) {
        p <- p + ggrepel::geom_text_repel(
          data = inv_data,
          aes(label = Species),
          size = 2.5,
          color = "#e74c3c",
          fontface = "bold",
          max.overlaps = 10,
          na.rm = TRUE
        )
      }
    }
    
    plots[[page]] <- p
    
    # Save if output_dir specified
    if (!is.null(output_dir)) {
      safe_country <- gsub("[^A-Za-z0-9_-]", "_", country)
      filename <- file.path(output_dir, paste0(safe_country, "_page", page, ".png"))
      ggsave(filename, p, width = 14, height = 3.5 * n_sites_this_page, limitsize = FALSE)
      cat("Saved:", filename, "\n")
    }
  }
  
  cat("\nReturning", length(plots), "plot(s). Use ggsave() to save or index with [[1]], [[2]], etc.\n")
  cat("Example: ggsave('plot.png', plots[[1]], width = 14, height = 12)\n")
  
  invisible(plots)
}

# =============================================================================
# USAGE
# =============================================================================

# Single site:
# plot_site_facet_theta(combined_all, "France, Narbona, Port de la Nouvelle")

# All sites in a country (4 sites per page):
# plot_country_facet_theta(combined_all, "Greece")
# plot_country_facet_theta(combined_all, "Norway")

# Save to files:
# plot_country_facet_theta(combined_all, "Greece", output_dir = "Figures/")

# Change sites per page:
# plot_country_facet_theta(combined_all, "Spain", sites_per_page = 3)
