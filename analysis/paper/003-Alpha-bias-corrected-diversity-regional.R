#Version that seems to add shannon
generate_alphaCoverageRichness_plots <- function(Listdf, rects, output_dir = here::here("analysis", "figures")) {

  plot_list <- list()
  results_list <- list()

  for (region_name in names(Listdf)) {
    region_data <- Listdf[[region_name]]
    region_results <- data.frame()

    for (time_name in names(region_data)) {
      sublist <- region_data[[time_name]]
      if (is.null(sublist) || nrow(sublist) == 0) next

      # Convert to matrix, drop Time column if exists
      comm_data <- as.matrix(sublist)
      if (ncol(comm_data) == 0) next

      # ---- Observed richness (raw species count) ----
      obs_rich <- tryCatch({
        sum(rowSums(comm_data) > 0)
      }, error = function(e) NA)

      # ---- Sample coverage ----
      cov <- tryCatch({
        Coverage(comm_data)
      }, error = function(e) NA)

      # ---- Shannon diversity ----
      shannon_div <- tryCatch({
        MC <- MetaCommunity(comm_data)
        diversity_obj <- AlphaDiversity(MC, q = 1, Correction = "Best")
        diversity_obj$Total
      }, error = function(e) NA)

      # ---- Append results ----
      if (!is.na(obs_rich)) {
        region_results <- dplyr::bind_rows(
          region_results,
          data.frame(Time = time_name, Metric = "Observed Richness", Value = obs_rich)
        )
      }

      if (!is.na(cov)) {
        region_results <- dplyr::bind_rows(
          region_results,
          data.frame(Time = time_name, Metric = "Coverage", Value = cov)
        )
      }

      if (!is.na(shannon_div)) {
        region_results <- dplyr::bind_rows(
          region_results,
          data.frame(Time = time_name, Metric = "Shannon Diversity", Value = shannon_div)
        )
      }
    }

    if (nrow(region_results) == 0) {
      message("No valid results for region ", region_name, " — skipping plotting.")
      next
    }

    # Convert Time to numeric for plotting
    time_numeric <- suppressWarnings(as.numeric(region_results$Time))
    region_results$TimeNumeric <- ifelse(is.na(time_numeric),
                                         seq_along(region_results$Time), time_numeric)

    # Plot range and crop rectangles
    plot_min <- min(region_results$TimeNumeric, na.rm = TRUE) - 500
    plot_max <- max(region_results$TimeNumeric, na.rm = TRUE) + 500
    rects_filtered <- rects[rects$ystart <= plot_max & rects$yend >= plot_min, ]
    rects_filtered$ystart <- pmin(rects_filtered$ystart, plot_max)
    rects_filtered$yend   <- pmax(rects_filtered$yend, plot_min)

    # Time axis breaks
    time_breaks <- pretty(region_results$TimeNumeric, n = 10)

    # ---- Plot ----
    p <- ggplot() +
      geom_rect(
        data = rects_filtered,
        inherit.aes = FALSE,
        aes(ymin = ystart, ymax = yend, xmin = -Inf, xmax = Inf, fill = col),
        alpha = 0.5
      ) +
      geom_point(data = region_results, aes(x = Value, y = TimeNumeric), color = "black", size = 3) +
      geom_path(data = region_results, aes(x = Value, y = TimeNumeric, group = Metric),
                color = "black", linetype = "dashed") +
      facet_wrap(~Metric, scales = "free_x") +
      ggsci::scale_fill_jco(guide = guide_legend(reverse = TRUE)) +
      scale_y_reverse(breaks = time_breaks, labels = time_breaks, expand = expansion(mult = c(0.01, 0.01))) +
      labs(title = paste("Alpha diversity metrics for", region_name),
           x = "Value", y = "Time", fill = "Period") +
      coord_cartesian(ylim = c(plot_max, plot_min)) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black", angle = 45, vjust = .5),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
        strip.background = element_rect(fill = "#f0f0f0", color = NA),
        strip.text.x = element_text(size = 14, face = "bold", colour = "black"),
        legend.title = element_text(size = 16, face = "bold", colour = "black"),
        legend.position = "right"
      )

    # Save plot
    safe_name <- gsub("[^A-Za-z0-9_]", "_", region_name)
    file_name <- paste0("003-", safe_name, "-alpha-diversity.jpg")
    ggsave(filename = file_name, plot = p, device = "jpg",
           path = output_dir, width = 3300, height = 4200, units = "px", dpi = 300)

    plot_list[[region_name]] <- p
    results_list[[region_name]] <- region_results
  }

  return(list(plots = plot_list, results = results_list))
}



generate_alphaCoverageRichness_plots(Listdf2, rects, output_dir = here::here("analysis", "figures"))

# Version that works for richness and Coverage
generate_alphaCoverageRichness_plots <- function(Listdf, rects, output_dir = here::here("analysis", "figures")) {

  plot_list <- list()
  results_list <- list()

  for (region_name in names(Listdf)) {
    region_data <- Listdf[[region_name]]
    region_results <- data.frame()

    for (time_name in names(region_data)) {
      sublist <- region_data[[time_name]]
      if (is.null(sublist) || nrow(sublist) == 0) next

      # Convert to matrix, drop Time column if exists
      comm_data <- as.matrix(sublist)
      if (ncol(comm_data) == 0) next

      # ---- Observed richness (raw species count) ----
      obs_rich <- tryCatch({
        sum(rowSums(comm_data) > 0)
      }, error = function(e) NA)

      # ---- Sample coverage ----
      cov <- tryCatch({
        Coverage(comm_data)
      }, error = function(e) NA)

      # ---- Append results ----
      if (!is.na(obs_rich)) {
        region_results <- dplyr::bind_rows(
          region_results,
          data.frame(Time = time_name, Metric = "Observed Richness", Value = obs_rich)
        )
      }

      if (!is.na(cov)) {
        region_results <- dplyr::bind_rows(
          region_results,
          data.frame(Time = time_name, Metric = "Coverage", Value = cov)
        )
      }
    }

    if (nrow(region_results) == 0) {
      message("No valid results for region ", region_name, " — skipping plotting.")
      next
    }

    # Convert Time to numeric for plotting
    time_numeric <- suppressWarnings(as.numeric(region_results$Time))
    region_results$TimeNumeric <- ifelse(is.na(time_numeric),
                                         seq_along(region_results$Time), time_numeric)

    # Plot range and crop rectangles
    plot_min <- min(region_results$TimeNumeric, na.rm = TRUE) - 500
    plot_max <- max(region_results$TimeNumeric, na.rm = TRUE) + 500
    rects_filtered <- rects[rects$ystart <= plot_max & rects$yend >= plot_min, ]
    rects_filtered$ystart <- pmin(rects_filtered$ystart, plot_max)
    rects_filtered$yend   <- pmax(rects_filtered$yend, plot_min)

    # Time axis breaks
    time_breaks <- pretty(region_results$TimeNumeric, n = 10)

    # ---- Plot ----
    p <- ggplot() +
      geom_rect(
        data = rects_filtered,
        inherit.aes = FALSE,
        aes(ymin = ystart, ymax = yend, xmin = -Inf, xmax = Inf, fill = col),
        alpha = 0.5
      ) +
      geom_point(data = region_results, aes(x = Value, y = TimeNumeric), color = "black", size = 3) +
      geom_path(data = region_results, aes(x = Value, y = TimeNumeric, group = Metric),
                color = "black", linetype = "dashed") +
      facet_wrap(~Metric, scales = "free_x") +
      ggsci::scale_fill_jco(guide = guide_legend(reverse = TRUE)) +
      scale_y_reverse(breaks = time_breaks, labels = time_breaks, expand = expansion(mult = c(0.01, 0.01))) +
      labs(title = paste("Alpha diversity metrics for", region_name),
           x = "Value", y = "Time", fill = "Period") +
      coord_cartesian(ylim = c(plot_max, plot_min)) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black", angle = 45, vjust = .5),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
        strip.background = element_rect(fill = "#f0f0f0", color = NA),
        strip.text.x = element_text(size = 14, face = "bold", colour = "black"),
        legend.title = element_text(size = 16, face = "bold", colour = "black"),
        legend.position = "right"
      )

    # Save plot
    safe_name <- gsub("[^A-Za-z0-9_]", "_", region_name)
    file_name <- paste0("003-", safe_name, "-alpha-diversity.jpg")
    ggsave(filename = file_name, plot = p, device = "jpg",
           path = output_dir, width = 3300, height = 4200, units = "px", dpi = 300)

    plot_list[[region_name]] <- p
    results_list[[region_name]] <- region_results
  }

  return(list(plots = plot_list, results = results_list))
}

generate_alphaCoverageRichness_plots(Listdf2, rects, output_dir = here::here("analysis", "figures"))
