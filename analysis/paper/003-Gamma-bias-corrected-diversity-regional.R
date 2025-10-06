generate_gammaCoverage_plots <- function(Listdf, rects,
                                         output_dir = here::here("analysis", "figures")) {
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  requireNamespace("ggsci")
  requireNamespace("entropart")

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  plot_list <- list()
  results_list <- list()

  for (region_name in names(Listdf)) {
    region_data <- Listdf[[region_name]]
    region_results <- data.frame(Time = character(), Metric = character(), Value = numeric(),
                                 stringsAsFactors = FALSE)

    for (time_name in names(region_data)) {
      sublist <- region_data[[time_name]]
      if (is.null(sublist) || nrow(sublist) == 0) next

      comm_mat <- tryCatch(as.matrix(sublist), error = function(e) NULL)
      if (is.null(comm_mat) || ncol(comm_mat) == 0) next

      # MetaCommunity creation
      meta <- tryCatch({
        entropart::MetaCommunity(Abundances = comm_mat)
      }, error = function(e) {
        warning(sprintf("MetaCommunity failed for %s %s: %s", region_name, time_name, e$message))
        NULL
      })
      if (is.null(meta)) next

      # Coverage (for pooled abundances)
      cov_time <- tryCatch({
        entropart::Coverage(rowSums(comm_mat))
      }, error = function(e) NA_real_)

      # Gamma metrics (whole MetaCommunity)
      richness <- tryCatch({
        entropart::GammaDiversity(meta, q = 0)
      }, error = function(e) NA_real_)

      shannon <- tryCatch({
        entropart::GammaDiversity(meta, q = 1)
      }, error = function(e) NA_real_)

      simpson <- tryCatch({
        entropart::GammaDiversity(meta, q = 2)
      }, error = function(e) NA_real_)

      # Collect results
      region_results <- dplyr::bind_rows(
        region_results,
        data.frame(Time = time_name, Metric = "Coverage", Value = cov_time),
        data.frame(Time = time_name, Metric = "Richness", Value = richness),
        data.frame(Time = time_name, Metric = "Shannon", Value = shannon),
        data.frame(Time = time_name, Metric = "Inverse Simpson", Value = simpson)
      )
    }

    if (nrow(region_results) == 0) next

    # Prepare numeric time for plotting
    region_results$TimeNumeric <- suppressWarnings(as.numeric(region_results$Time))
    region_results <- region_results[!is.na(region_results$TimeNumeric), ]

    # Ensure facet order
    region_results$Metric <- factor(
      region_results$Metric,
      levels = c("Coverage", "Richness", "Shannon", "Inverse Simpson")
    )

    # Plot limits and rect prep
    plot_min <- min(region_results$TimeNumeric, na.rm = TRUE)
    plot_max <- max(region_results$TimeNumeric, na.rm = TRUE)

    rects_filtered <- rects[rects$ystart <= plot_max & rects$yend >= plot_min, , drop = FALSE]
    if (nrow(rects_filtered) > 0) {
      rects_filtered$ystart <- pmin(rects_filtered$ystart, plot_max)
      rects_filtered$yend <- pmax(rects_filtered$yend, plot_min)
    }

    # Time axis breaks
    min_time_ticks <- floor(plot_min / 500) * 500
    max_time_ticks <- ceiling(plot_max / 500) * 500
    time_breaks <- seq(from = max_time_ticks, to = min_time_ticks, by = -500)

    # Plot
    p <- ggplot2::ggplot() +
      ggplot2::geom_rect(
        data = rects_filtered,
        ggplot2::aes(ymin = ystart, ymax = yend, xmin = -Inf, xmax = Inf, fill = col),
        alpha = 0.5
      ) +
      ggplot2::geom_point(data = region_results,
                          ggplot2::aes(x = Value, y = TimeNumeric), color = "black", size = 3) +
      ggplot2::geom_path(data = region_results,
                         ggplot2::aes(x = Value, y = TimeNumeric, group = Metric),
                         color = "black", linetype = "dashed") +
      ggplot2::facet_wrap(~Metric, scales = "free_x") +
      ggsci::scale_fill_jco(guide = ggplot2::guide_legend(reverse = TRUE)) +
      ggplot2::scale_y_reverse(breaks = time_breaks, labels = time_breaks,
                               expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
      ggplot2::labs(title = paste("Gamma diversity metrics for", region_name),
                    x = "Diversity Value", y = "Time (Years BP)", fill = "Period") +
      ggplot2::coord_cartesian(ylim = c(plot_max, plot_min)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 12, colour = "black"),
        axis.text.x = ggplot2::element_text(size = 12, colour = "black", angle = 45, vjust = .5),
        axis.title.y = ggplot2::element_text(size = 12, face = "bold", colour = "black"),
        strip.background = ggplot2::element_rect(fill = "#f0f0f0", color = NA),
        strip.text.x = ggplot2::element_text(size = 14, face = "bold", colour = "black"),
        legend.title = ggplot2::element_text(size = 16, face = "bold", colour = "black"),
        legend.position = "right"
      )

    # Save plot
    safe_name <- gsub("[^A-Za-z0-9_]", "_", region_name)
    file_name <- paste0("004-", safe_name, "gamma-diversity.jpg")
    ggplot2::ggsave(filename = file_name, plot = p, device = "jpg",
                    path = output_dir, width = 3300, height = 4200, units = "px", dpi = 300)

    plot_list[[region_name]] <- p
    results_list[[region_name]] <- region_results
  }

  return(list(plots = plot_list, results = results_list))
}


generate_gammaCoverage_plots(Listdf2, rects,
                                         output_dir = here::here("analysis", "figures"))
