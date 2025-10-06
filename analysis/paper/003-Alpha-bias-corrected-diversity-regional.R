# Period definitions from Pilotto et al. (2022)
rects <- data.frame(
  ystart = c(12000, 9500, 8000, 6500, 4500, 4000, 3500, -500),
  yend   = c(16000,12000, 9500, 8000, 6500, 4500, 4000, 3500),
  col    = factor(
    c("TM 1", "TM 2", "TM 3", "TM 4", "TM 5", "TM 6", "TM 7", "TM 8"),
    levels = c("TM 1", "TM 2", "TM 3", "TM 4", "TM 5", "TM 6", "TM 7", "TM 8")
  )
)

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
        entropart::Coverage(comm_data)
      }, error = function(e) NA)

      # ---- Build MetaCommunity once (if possible) ----
      MC <- tryCatch({
        entropart::MetaCommunity(Abundances = comm_data)
      }, error = function(e) {
        message(sprintf("MetaCommunity failed for %s - %s: %s", region_name, time_name, e$message))
        NULL
      })

      # ---- Shannon diversity (q = 1) ----
      shannon_div <- NA
      if (!is.null(MC)) {
        shannon_div <- tryCatch({
          ad <- entropart::AlphaDiversity(MC, q = 1, Correction = "Best")
          # AlphaDiversity returns an MCdiversity object; the total alpha is in $Total
          as.numeric(ad$Total)
        }, error = function(e) {
          message(sprintf("Alpha Shannon failed for %s - %s: %s", region_name, time_name, e$message))
          NA
        })
      }

      # ---- Simpson diversity (q = 2) as Inverse Simpson ----
      simpson_div <- NA
      if (!is.null(MC)) {
        simpson_div <- tryCatch({
          ad2 <- entropart::AlphaDiversity(MC, q = 2, Correction = "Best")
          as.numeric(ad2$Total)
        }, error = function(e) {
          message(sprintf("Alpha Simpson failed for %s - %s: %s", region_name, time_name, e$message))
          NA
        })
      }

      # ---- Append results ----
      if (!is.na(obs_rich)) {
        region_results <- dplyr::bind_rows(
          region_results,
          data.frame(Time = time_name, Metric = "Observed Richness", Value = obs_rich, stringsAsFactors = FALSE)
        )
      }

      if (!is.na(cov)) {
        region_results <- dplyr::bind_rows(
          region_results,
          data.frame(Time = time_name, Metric = "Coverage", Value = cov, stringsAsFactors = FALSE)
        )
      }

      if (!is.na(shannon_div)) {
        region_results <- dplyr::bind_rows(
          region_results,
          data.frame(Time = time_name, Metric = "Shannon", Value = shannon_div, stringsAsFactors = FALSE)
        )
      }

      if (!is.na(simpson_div)) {
        region_results <- dplyr::bind_rows(
          region_results,
          data.frame(Time = time_name, Metric = "Inverse Simpson", Value = simpson_div, stringsAsFactors = FALSE)
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

    # Ensure facet order: Coverage, Observed Richness, Shannon, Inverse Simpson
    region_results$Metric <- factor(region_results$Metric,
                                    levels = c("Coverage", "Observed Richness", "Shannon", "Inverse Simpson"))

    # Plot range and crop rectangles
    plot_min <- min(region_results$TimeNumeric, na.rm = TRUE) - 500
    plot_max <- max(region_results$TimeNumeric, na.rm = TRUE) + 500
    rects_filtered <- rects[rects$ystart <= plot_max & rects$yend >= plot_min, ]
    rects_filtered$ystart <- pmin(rects_filtered$ystart, plot_max)
    rects_filtered$yend   <- pmax(rects_filtered$yend, plot_min)

    # Time axis breaks
    time_breaks <- pretty(region_results$TimeNumeric, n = 10)

    # ---- Plot ----
    p <- ggplot2::ggplot() +
      ggplot2::geom_rect(
        data = rects_filtered,
        inherit.aes = FALSE,
        ggplot2::aes(ymin = ystart, ymax = yend, xmin = -Inf, xmax = Inf, fill = col),
        alpha = 0.5
      ) +
      ggplot2::geom_point(data = region_results, ggplot2::aes(x = Value, y = TimeNumeric), color = "black", size = 3) +
      ggplot2::geom_path(data = region_results, ggplot2::aes(x = Value, y = TimeNumeric, group = Metric),
                         color = "black", linetype = "dashed") +
      ggplot2::facet_wrap(~Metric, scales = "free_x") +
      ggsci::scale_fill_jco(guide = ggplot2::guide_legend(reverse = TRUE)) +
      ggplot2::scale_y_reverse(breaks = time_breaks, labels = time_breaks, expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
      ggplot2::labs(title = paste("Alpha diversity metrics for", region_name),
                    x = "Value", y = "Time", fill = "Period") +
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
    file_name <- paste0("003-", safe_name, "-alpha-diversity.jpg")
    ggplot2::ggsave(filename = file_name, plot = p, device = "jpg",
                    path = output_dir, width = 3300, height = 4200, units = "px", dpi = 300)

    plot_list[[region_name]] <- p
    results_list[[region_name]] <- region_results
  }

  return(list(plots = plot_list, results = results_list))
}

# Generate and save plots
generate_alphaCoverageRichness_plots(Listdf2, rects, output_dir = here::here("analysis", "figures"))
