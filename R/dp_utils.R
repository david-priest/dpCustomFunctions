# ==============================================================================
# dpCustomFunctions — Additional / Updated Utilities
# ==============================================================================
# Functions added here either:
#   (a) did not exist in the original CustomFunctionsSept25 build, or
#   (b) are updated replacements for older versions in functions.R
#
# plot_spill replaces the original plot_spill in functions.R (improved version
# with valid_interactions mode and rasterisation support).
# ==============================================================================

# ==============================================================================
# SCE Marker Class Utilities
# ==============================================================================

#' Set custom type and state marker classes on a SingleCellExperiment
#'
#' Sets all markers to \code{"state"} class, then promotes the named
#' \code{type_markers} to \code{"type"}. Use this before CATALYST clustering
#' to ensure only the intended markers drive FlowSOM/metaclustering.
#'
#' @param sce A \code{SingleCellExperiment} object. Must have
#'   \code{rowData} with a \code{marker_name} column.
#' @param type_markers Character vector of marker names to assign class
#'   \code{"type"}. Names not present in \code{rowData(sce)$marker_name}
#'   are silently ignored.
#' @return The input SCE with updated \code{rowData(sce)$marker_class}.
#' @export
set_custom_marker_classes <- function(sce, type_markers) {
  rd <- SummarizedExperiment::rowData(sce)
  rd$marker_class <- "state"
  type_markers_in_sce <- type_markers[type_markers %in% rd$marker_name]
  rd$marker_class[rd$marker_name %in% type_markers_in_sce] <- "type"
  SummarizedExperiment::rowData(sce) <- rd
  return(sce)
}

# ==============================================================================
# Output / Save Utilities
# ==============================================================================

#' Save a plot to PDF with optional Excel data sidecar
#'
#' Saves any ggplot, ComplexHeatmap, or cowplot object to a PDF file.
#' When the plot is a ggplot with a data frame attached, an Excel sidecar is
#' also written for downstream traceability.
#'
#' @param p A plot object. Accepts \code{gg}/\code{ggplot},
#'   \code{Heatmap} (ComplexHeatmap), or \code{gtable} (cowplot).
#' @param name Filename stem (no extension, no directory path).
#' @param w Width in inches. Default \code{10}.
#' @param h Height in inches. Default \code{10}.
#' @param out_dir Output directory. Created recursively if it does not exist.
#'   Default \code{"outputs/figures"}.
#' @param save_data If \code{TRUE}, saves an Excel sidecar when the plot
#'   carries a data frame (ggplot objects only). Default \code{TRUE}.
#' @param width Alias for \code{w}; overrides \code{w} when supplied.
#'   Default \code{NULL}.
#' @param height Alias for \code{h}; overrides \code{h} when supplied.
#'   Default \code{NULL}.
#' @return Invisibly returns the path to the saved PDF.
#' @export
save_plot_extended <- function(p, name, w = 10, h = 10,
                               out_dir = "outputs/figures",
                               save_data = TRUE,
                               width = NULL, height = NULL) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  file_latest <- file.path(out_dir, paste0(name, ".pdf"))

  pdf_width  <- if (!is.null(width))  width  else w
  pdf_height <- if (!is.null(height)) height else h
  grDevices::pdf(file_latest, width = pdf_width, height = pdf_height)
  tryCatch({
    if (inherits(p, "Heatmap")) {
      ComplexHeatmap::draw(p)
    } else if (inherits(p, "gtable")) {
      grid::grid.newpage()
      grid::grid.draw(p)
    } else {
      print(p)
    }
  }, error = function(e) message("Error saving plot: ", e$message))
  grDevices::dev.off()

  if (save_data) {
    data_latest <- file.path(out_dir, paste0(name, ".xlsx"))
    plot_data <- NULL
    if (inherits(p, "gg") && !is.null(p$data) && is.data.frame(p$data)) {
      plot_data <- p$data
    } else if (is.data.frame(p)) {
      plot_data <- p
    }
    if (!is.null(plot_data) && nrow(plot_data) > 0) {
      tryCatch({
        wb <- openxlsx::createWorkbook()
        openxlsx::modifyBaseFont(wb, fontSize = 12, fontColour = "black",
                                 fontName = "Arial")
        openxlsx::addWorksheet(wb, "Sheet 1")
        openxlsx::writeData(wb, "Sheet 1", plot_data)
        openxlsx::setColWidths(wb, "Sheet 1",
                               cols = seq_len(ncol(plot_data)), widths = "auto")
        openxlsx::freezePane(wb, "Sheet 1", firstRow = TRUE)
        openxlsx::addFilter(wb, "Sheet 1", rows = 1,
                            cols = seq_len(ncol(plot_data)))
        openxlsx::saveWorkbook(wb, data_latest, overwrite = TRUE)
      }, error = function(e) {
        message("Could not save Excel sidecar: ", e$message)
      })
    }
  }
  invisible(file_latest)
}

#' Save a cluster-merging template Excel from a plotExprHeatmap1 result
#'
#' Creates an Excel workbook with two columns: \code{old_cluster} (numeric,
#' in heatmap visual order top-to-bottom) and \code{new_cluster} (empty, for
#' the user to fill in before the next pipeline phase). The file is saved to
#' \code{out_dir/<name>_merging_template.xlsx}.
#'
#' @param row_labels_df Data frame with a \code{row_labels} column, as
#'   returned in the \code{$row_labels_df} element of
#'   \code{\link{plotExprHeatmap1}}.
#' @param name Filename stem (no extension, no directory path).
#' @param out_dir Output directory. Created recursively if it does not exist.
#'   Default \code{"outputs"}.
#' @return Path to the saved Excel file.
#' @export
save_merging_template <- function(row_labels_df, name, out_dir = "outputs") {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  df <- data.frame(
    old_cluster = as.numeric(row_labels_df$row_labels),
    new_cluster = "",
    stringsAsFactors = FALSE
  )
  file_path <- file.path(out_dir, paste0(name, "_merging_template.xlsx"))
  wb <- openxlsx::createWorkbook()
  openxlsx::modifyBaseFont(wb, fontSize = 12, fontColour = "black",
                           fontName = "Arial")
  openxlsx::addWorksheet(wb, "Sheet 1")
  openxlsx::writeData(wb, "Sheet 1", df)
  openxlsx::setColWidths(wb, "Sheet 1", cols = 1:2, widths = "auto")
  openxlsx::freezePane(wb, "Sheet 1", firstRow = TRUE)
  openxlsx::addFilter(wb, "Sheet 1", rows = 1, cols = 1:2)
  openxlsx::saveWorkbook(wb, file_path, overwrite = TRUE)
  return(file_path)
}

#' Archive all output figures to a timestamped subfolder
#'
#' Copies all PDF, PNG, and XLSX files from \code{output_dir} into a new
#' timestamped subdirectory under \code{archive_base_dir}. The originals
#' remain in \code{output_dir} so that a \code{{targets}} pipeline does not
#' consider them stale on the next run.
#'
#' @param output_dir Directory containing the current output files.
#'   Default \code{"outputs/figures"}.
#' @param archive_base_dir Parent directory for the archive subfolder.
#'   Defaults to \code{dirname(output_dir)}.
#' @return Path to the created archive directory.
#' @export
archive_figures <- function(output_dir = "outputs/figures",
                            archive_base_dir = NULL) {
  if (is.null(archive_base_dir)) archive_base_dir <- dirname(output_dir)
  timestamp   <- format(Sys.time(), "%Y-%m-%d_%H.%M.%S")
  archive_dir <- file.path(archive_base_dir, paste0(timestamp, "_figures"))
  dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)
  figure_files <- list.files(output_dir, pattern = "\\.(png|pdf|xlsx)$",
                             full.names = TRUE)
  if (length(figure_files) > 0) {
    file.copy(figure_files, archive_dir, overwrite = TRUE, recursive = TRUE)
  }
  return(archive_dir)
}

#' Write a JSON status summary for a CATALYST SingleCellExperiment
#'
#' Saves key metadata from a CATALYST SCE to a JSON file: colData column
#' names, reduced dimension names, cluster codes, cell and feature counts,
#' type and state marker lists, and a timestamp. Intended as a "flight
#' recorder" target at the end of a \code{{targets}} pipeline run.
#'
#' @param sce A \code{SingleCellExperiment} processed by CATALYST.
#'   Requires the \code{cluster_codes} entry in \code{metadata(sce)}.
#' @param path Output file path. Default \code{"status.json"}.
#' @return The \code{path} argument, invisibly.
#' @export
write_status_json <- function(sce, path = "status.json") {
  status <- list(
    colnames_colData = colnames(SingleCellExperiment::colData(sce)),
    reducedDimNames  = SingleCellExperiment::reducedDimNames(sce),
    cluster_codes    = colnames(CATALYST::cluster_codes(sce)),
    n_cells          = ncol(sce),
    n_features       = nrow(sce),
    type_markers     = CATALYST::type_markers(sce),
    state_markers    = CATALYST::state_markers(sce),
    timestamp        = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
  jsonlite::write_json(status, path, pretty = TRUE)
  cat("Status written to:", path, "\n")
  invisible(path)
}

# ==============================================================================
# Spillover QC (improved version replacing plot_spill in functions.R)
# ==============================================================================

#' Plot spillover scatter grids for compensation QC
#'
#' Generates a grid of scatter plots (one per channel pair) visualising
#' spillover relationships in a CyTOF experiment. Supports two modes:
#' \code{"spill_matrix"} (all non-zero entries) and
#' \code{"valid_interactions"} (physically meaningful isotope neighbours only).
#'
#' @param spill_matrix Numeric matrix. Spillover/compensation matrix with
#'   channel names as row and column names.
#' @param cytof_data A CATALYST \code{SingleCellExperiment} to plot cell-level
#'   expression values from.
#' @param panel_table Data frame with columns \code{fcs_colname} (channel
#'   names matching \code{spill_matrix} dimnames) and \code{antigen}.
#' @param n_cells Maximum number of cells per sample to sub-sample before
#'   plotting. Default \code{1000}.
#' @param my_seed Random seed for sub-sampling. Default \code{1234}.
#' @param exclude_channels Character vector of source channel names to skip.
#'   Default \code{NULL}.
#' @param n_col Number of columns in the output grid. Default \code{6}.
#' @param mode One of \code{"spill_matrix"} (plot non-zero entries only) or
#'   \code{"valid_interactions"} (plot physically plausible isotope pairs).
#'   Default \code{"spill_matrix"}.
#' @param valid_interactions Character vector controlling which interaction
#'   types are shown in \code{"valid_interactions"} mode. Any combination of
#'   \code{"same"} (same element), \code{"ox"} (+16 Da oxide), \code{"one_plus"}
#'   (+1 Da), \code{"one_minus"} (-1 Da). Default: all four.
#' @param rasterize Logical. If \code{TRUE}, rasterise hex layers via
#'   \code{ggrastr} to keep PDF file sizes small. Default \code{TRUE}.
#' @param dpi Resolution for rasterised layers. Default \code{150}.
#' @param title_size Font size for plot titles. Default \code{8}.
#' @param axis_size Font size for axis labels. Default \code{7}.
#' @param organize_by One of \code{"source"} or \code{"destination"}.
#'   Controls whether plots are grouped by the emitting or the receiving
#'   channel. Default \code{"source"}.
#' @return A \code{cowplot} composite plot, or \code{NULL} if no pairs pass
#'   the selection criteria.
#' @export
plot_spill <- function(spill_matrix, cytof_data, panel_table,
                       n_cells   = 1000,
                       my_seed   = 1234,
                       exclude_channels = NULL,
                       n_col     = 6,
                       mode      = c("spill_matrix", "valid_interactions"),
                       valid_interactions = c("same", "ox", "one_plus", "one_minus"),
                       rasterize = TRUE,
                       dpi       = 150,
                       title_size = 8,
                       axis_size  = 7,
                       organize_by = c("source", "destination")) {

  mode        <- match.arg(mode)
  organize_by <- match.arg(organize_by)

  # Sub-sample cytof_data
  set.seed(my_seed)
  idx <- split(seq(ncol(cytof_data)), cytof_data$sample_id)
  idx <- lapply(idx, function(.) sample(., min(n_cells, length(.))))
  cytof_data <- cytof_data[, unlist(idx)]

  # Helper: classify isotope interactions
  is_valid_interaction <- function(source, dest,
                                   allowed = valid_interactions) {
    same_element <- (gsub("[0-9]", "", dest) == gsub("[0-9]", "", source))
    ox           <- (readr::parse_number(dest) == readr::parse_number(source) + 16)
    one_plus     <- (readr::parse_number(dest) == readr::parse_number(source) + 1)
    one_minus    <- (readr::parse_number(dest) == readr::parse_number(source) - 1)
    (("same"      %in% allowed && same_element) ||
     ("ox"        %in% allowed && ox)           ||
     ("one_plus"  %in% allowed && one_plus)     ||
     ("one_minus" %in% allowed && one_minus))
  }

  # Select entries to plot
  if (mode == "spill_matrix") {
    relevant_entries <- which(spill_matrix != 0, arr.ind = TRUE)
  } else {
    relevant_entries <- which(spill_matrix >= 0 & spill_matrix <= 1,
                              arr.ind = TRUE)
  }
  if (nrow(relevant_entries) == 0) {
    warning("No entries pass the selection criteria. Returning NULL.")
    return(NULL)
  }

  # Sort by source or destination channel number
  order_channels <- if (organize_by == "source") {
    rownames(spill_matrix)[relevant_entries[, "row"]]
  } else {
    colnames(spill_matrix)[relevant_entries[, "col"]]
  }
  relevant_entries <- relevant_entries[
    order(as.numeric(gsub("[^0-9]", "", order_channels))), , drop = FALSE
  ]

  grouped_plots <- list()

  for (i in seq_len(nrow(relevant_entries))) {
    source <- rownames(spill_matrix)[relevant_entries[i, "row"]]
    dest   <- colnames(spill_matrix)[relevant_entries[i, "col"]]

    if (source == dest) next
    if (!is.null(exclude_channels) && source %in% exclude_channels) next
    if (mode == "valid_interactions" && !is_valid_interaction(source, dest)) next

    source_antigen <- panel_table$antigen[panel_table$fcs_colname == source]
    dest_antigen   <- panel_table$antigen[panel_table$fcs_colname == dest]

    if (length(source_antigen) < 1 || length(dest_antigen) < 1) {
      message("Could not map channels ", source, " and/or ", dest,
              " to antigens. Skipping.")
      next
    }
    if (!(source_antigen %in% rownames(cytof_data)) ||
        !(dest_antigen   %in% rownames(cytof_data))) {
      message("Antigens ", source_antigen, " and/or ", dest_antigen,
              " not found in cytof_data. Skipping.")
      next
    }

    spill_value <- sprintf("%.4f",
      spill_matrix[relevant_entries[i, "row"], relevant_entries[i, "col"]])

    title <- if (organize_by == "source") {
      paste0("Dest: ", dest, " (", dest_antigen, ")\nComp: ", spill_value)
    } else {
      paste0("Src: ", source, " (", source_antigen, ")\nComp: ", spill_value)
    }

    p <- tryCatch({
      bp <- suppressWarnings(
        CATALYST::plotScatter(cytof_data, c(source_antigen, dest_antigen),
                              zeros = TRUE)
      )
      bp <- bp +
        ggplot2::ggtitle(title) +
        ggplot2::theme(
          plot.title   = ggplot2::element_text(size = title_size, hjust = 0.5),
          axis.title   = ggplot2::element_text(size = axis_size),
          axis.text    = ggplot2::element_text(size = axis_size - 1),
          plot.margin  = ggplot2::margin(2, 2, 2, 2, "pt"),
          legend.position = "none"
        )
      if (rasterize && requireNamespace("ggrastr", quietly = TRUE)) {
        for (j in seq_along(bp$layers)) {
          if (inherits(bp$layers[[j]]$geom, "GeomHex")) {
            bp$layers[[j]] <- ggrastr::rasterise(bp$layers[[j]], dpi = dpi)
          }
        }
      }
      bp
    }, error = function(e) {
      message("plotScatter failed for ", source, " -> ", dest, ": ", e$message)
      NULL
    })

    if (!is.null(p) && inherits(p, "ggplot")) {
      group_key <- if (organize_by == "source") source else dest
      if (!group_key %in% names(grouped_plots)) grouped_plots[[group_key]] <- list()
      grouped_plots[[group_key]][[length(grouped_plots[[group_key]]) + 1]] <- p
    }
  }

  if (length(grouped_plots) == 0) {
    warning("No plots were generated. Returning NULL.")
    return(NULL)
  }

  all_plots <- list()
  for (channel_name in names(grouped_plots)) {
    plist          <- grouped_plots[[channel_name]]
    channel_antigen <- panel_table$antigen[panel_table$fcs_colname == channel_name]
    if (length(channel_antigen) < 1) channel_antigen <- "Unknown"
    label_text <- if (organize_by == "source") "Source:" else "Destination:"

    label_plot <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.58,
                        label = label_text, size = 5, fontface = "bold",
                        color = "darkblue") +
      ggplot2::annotate("text", x = 0.5, y = 0.50,
                        label = channel_name, size = 5, fontface = "bold",
                        color = "darkred") +
      ggplot2::annotate("text", x = 0.5, y = 0.42,
                        label = channel_antigen, size = 5, fontface = "bold",
                        color = "darkgreen") +
      ggplot2::theme_void() +
      ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0.35, 0.65)) +
      ggplot2::theme(
        plot.margin = ggplot2::margin(1, 1, 1, 1, "pt"),
        panel.background = ggplot2::element_rect(fill = "grey95",
                                                 color = "grey80"),
        aspect.ratio = 0.8
      )
    all_plots <- c(all_plots, list(label_plot), plist)
  }

  cowplot::plot_grid(plotlist = all_plots, ncol = n_col, align = "hv")
}
