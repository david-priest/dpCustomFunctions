## ----------------------------------------------------------------------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param title PARAM_DESCRIPTION, Default: 'your_title'
#' @param save_dir PARAM_DESCRIPTION, Default: '.'
#' @return OUTPUT_DESCRIPTION
#' @export 
backup <- function(title = "your_title", save_dir = ".") {
  # Get the name of the current RMD file
  library(rmdhelp)
  rmd_file_path <- get_this_rmd_file()
  rmd_file_name <- basename(rmd_file_path)

  # Check if the directory exists, and create it if it doesn't
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }
  
  # Create the timestamp for the file names
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H.%M.%S")
  
  # Save a copy of the RMD file
  rmd_filename <- file.path(save_dir, paste0(timestamp, "_", title, "_", rmd_file_name))
  file.copy(from = rmd_file_path, to = rmd_filename)

  cat("RMD file saved successfully in folder", save_dir, "\n")
}


## ----------------------------------------------------------------------------------------------------------
# Change the directory to your desired save location

# Uncomment this and run in between making the package
#backup("", save_dir = '/Users/davidpriest/My Drive/Wing Lab/Custom R functions/Custom R Functions Backups')



## ----------------------------------------------------------------------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param sce PARAM_DESCRIPTION
#' @param n_cells PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @export 
filterSCEvents <- function(sce, n_cells = NULL)
  
{
  
test = table(sce$sample_id)
test = test[test>n_cells]
test = as.data.frame(test)
sce = filterSCE(sce, sample_id %in% test$Var1)
  
return(sce)
}


## ----------------------------------------------------------------------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param sce PARAM_DESCRIPTION
#' @param n_cells PARAM_DESCRIPTION, Default: 1000
#' @param my_seed PARAM_DESCRIPTION, Default: 1234
#' @return OUTPUT_DESCRIPTION
#' @export 
subSCE <- function (sce, n_cells = 1000, my_seed = 1234)
  
{
  
  
  set.seed(my_seed) # Subsetting depends on the seed
  
  # split cell indices by sample
  idx <- split(seq(ncol(sce)), sce$sample_id)
  # downsample to at most 'n' cells per sample
  # (I would not! sample with replacement here, 
  # as you are essentially duplicating cells; 
  # instead, use at most as many cells as there are)
  idx <- lapply(idx, \(.) sample(., min(n_cells, length(.))))
  # subset the data
  sub <- sce[, unlist(idx)]
  
  return(sub)
}


## ----------------------------------------------------------------------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param sce PARAM_DESCRIPTION
#' @param n_cells PARAM_DESCRIPTION, Default: 500
#' @param my_seed PARAM_DESCRIPTION, Default: 1234
#' @param metalevel PARAM_DESCRIPTION, Default: 'merging2'
#' @return OUTPUT_DESCRIPTION
#' @export 
subSCE_cluster <- function(sce, n_cells = 500, my_seed = 1234, metalevel = "merging2") {
  
  set.seed(my_seed) # Subsetting depends on the seed
  
  # split cell indices by cluster label
  idx <- split(seq(ncol(sce)), sce[[metalevel]])
  
  # downsample to at most 'n' cells per cluster
  idx <- lapply(idx, \(.) sample(., min(n_cells, length(.))))
  
  # subset the data
  sub <- sce[, unlist(idx)]
  
  return(sub)
}



## ----------------------------------------------------------------------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param k PARAM_DESCRIPTION, Default: 'meta20'
#' @param by PARAM_DESCRIPTION, Default: c("sample_id", "cluster_id")
#' @param group_by PARAM_DESCRIPTION, Default: 'condition'
#' @param shape_by PARAM_DESCRIPTION, Default: NULL
#' @param col_clust PARAM_DESCRIPTION, Default: TRUE
#' @param n_cols PARAM_DESCRIPTION, Default: 4
#' @param log PARAM_DESCRIPTION, Default: FALSE
#' @param miny PARAM_DESCRIPTION, Default: 0.01
#' @param maxy PARAM_DESCRIPTION, Default: NA
#' @param point_size PARAM_DESCRIPTION, Default: 2
#' @param distance PARAM_DESCRIPTION, Default: c("euclidean", "maximum", "manhattan", "canberra", "binary", 
#'    "minkowski")
#' @param linkage PARAM_DESCRIPTION, Default: c("average", "ward.D", "single", "complete", "mcquitty", "median", 
#'    "centroid", "ward.D2")
#' @param k_pal PARAM_DESCRIPTION, Default: CATALYST:::.cluster_cols
#' @param clusters_order PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @export 
#' @importFrom CATALYST .cluster_cols
plotAbundances22 <- function (x, k = "meta20", by = c("sample_id", "cluster_id"), 
                              group_by = "condition", shape_by = NULL, col_clust = TRUE, n_cols = 4, log = FALSE, miny = 0.01, maxy = NA, point_size = 2,
                              distance = c("euclidean", "maximum", "manhattan", "canberra", 
                                           "binary", "minkowski"), linkage = c("average", "ward.D", 
                                                                               "single", "complete", "mcquitty", "median", "centroid", 
                                                                               "ward.D2"), k_pal = CATALYST:::.cluster_cols, clusters_order = NULL) {
  by <- match.arg(by)
  .check_sce(x, TRUE)
  k <- .check_k(x, k)
  .check_cd_factor(x, group_by)
  .check_cd_factor(x, shape_by)
  .check_pal(k_pal)
  linkage <- match.arg(linkage)
  distance <- match.arg(distance)
  stopifnot(is.logical(col_clust), length(col_clust) == 1)
  shapes <- .get_shapes(x, shape_by)
  if (is.null(shapes)) 
    shape_by <- NULL
  if (by == "sample_id") {
    nk <- nlevels(cluster_ids(x, k))
    if (length(k_pal) < nk) 
      k_pal <- colorRampPalette(k_pal)(nk)
  }
  ns <- table(cluster_id = cluster_ids(x, k), sample_id = sample_ids(x))
  fq <- prop.table(ns, 2) * 100
  df <- as.data.frame(fq)
  m <- match(df$sample_id, x$sample_id)
  for (i in c(shape_by, group_by)) df[[i]] <- x[[i]][m]
  if (by == "sample_id" && col_clust && length(unique(df$sample_id)) > 1) {
    d <- dist(t(fq), distance)
    h <- hclust(d, linkage)
    o <- colnames(fq)[h$order]
    df$sample_id <- factor(df$sample_id, o)
  }
  
  # Apply clusters_order if provided
  if (!is.null(clusters_order)) {
    df$cluster_id <- factor(df$cluster_id, levels = clusters_order)
  }
  
  # Decide whether to plot on a log scale. Add 0.02 to all values
  if (log == TRUE) {
    df$Freq <- df$Freq + 0.02
  }
  
  dfout <<- df
  
  p <- ggplot(df, aes_string(y = "Freq")) + 
    labs(x = NULL, y = "Proportion [%]") + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          strip.text = element_text(face = "bold"), 
          strip.background = element_rect(fill = NA, color = NA), 
          axis.text = element_text(color = "black"), 
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
          legend.key.height = unit(0.8, "lines"))
  
  p <- p + facet_wrap(~cluster_id, scales = "free_y", ncol = n_cols) + 
    geom_boxplot(aes_string(x = group_by, fill = group_by), color = "grey16", position = position_dodge(), alpha = 0.8, outlier.color = NA, show.legend = TRUE)
  
  if (!is.null(shape_by)) {
    p <- p + geom_quasirandom(aes_string(x = group_by, shape = shape_by), size = point_size, width = 0.2)
  } else {
    p <- p + geom_quasirandom(aes_string(x = group_by), fill = "grey84", size = point_size, width = 0.2, shape = 21)
  }
  
  if (log == TRUE) {
    p + scale_y_continuous(trans = 'log10', limits = c(miny, maxy), breaks = c(0.01, 0.1, 1, 10, 100), labels = c(0.01, 0.1, 1, 10, 100)) + 
      annotation_logticks(base = 10, sides = "l", outside = TRUE) + 
      coord_cartesian(clip = "off") + 
      scale_size_area(max_size = 15) + 
      theme(axis.text.y = element_text(margin = margin(r = 8)))
  } else {
    p + scale_y_continuous(limits = c(0, maxy)) + 
      coord_cartesian(clip = "off") + 
      scale_size_area(max_size = 15)
  }
}

environment(plotAbundances22) <- asNamespace('CATALYST')




## ----------------------------------------------------------------------------------------------------------
 #Firstly here's a custome anno_counts function to fix the size of the cluster percentages label.
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param perc PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @export 
anno_counts1 <- function (x, perc) 
 {
     ns <- table(x)
     fq <- round(ns/sum(ns) * 100, 2)
     if (perc) {
         #txt <- sprintf("%s%%(%s)", fq, names(fq))  # I removed the part where it adds the cell type name
         txt <- paste(fq, "%")
         foo <- row_anno_text(txt, just = "center", gp = gpar(fontsize = 12),  # I increased the font size
             location = unit(0.5, "npc"))
     }
     else foo <- NULL
     rowAnnotation(n_cells = row_anno_barplot(x = as.matrix(ns), 
         width = unit(2, "cm"), gp = gpar(fill = "grey", col = "white"), 
         border = FALSE, axis = TRUE, bar_width = 0.8), foo = foo)
 }


#Editing anno_clusters to simply return the merging table directly


## ----------------------------------------------------------------------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param k PARAM_DESCRIPTION
#' @param m PARAM_DESCRIPTION
#' @param k_pal PARAM_DESCRIPTION
#' @param m_pal PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @export 
merging_ids <- function (x, k, m, k_pal, m_pal) 
{
    kids <- levels(x$cluster_id)
    nk <- length(kids)
    if (nk > length(k_pal)) 
        k_pal <- colorRampPalette(k_pal)(nk)
    k_pal <- k_pal[seq_len(nk)]
    names(k_pal) <- kids
    df <- data.frame(cluster_id = kids)
    col <- list(cluster_id = k_pal)
    if (!is.null(m)) {
        i <- match(kids, cluster_codes(x)[, k])
        mids <- droplevels(cluster_codes(x)[, m][i])
        nm <- nlevels(mids)
        if (nm > length(m_pal)) 
            m_pal <- colorRampPalette(m_pal)(nm)
        m_pal <- m_pal[seq_len(nm)]
        names(m_pal) <- levels(mids)
        df$merging_id <- mids
        col$merging_id <- m_pal
    }
    df <- mutate_all(df, function(u) factor(u, unique(u)))
    #rowAnnotation(df = df, col = col, gp = gpar(col = "white"))
    return(df)
}


## ----------------------------------------------------------------------------------------------------------

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param k PARAM_DESCRIPTION
#' @param m PARAM_DESCRIPTION
#' @param k_pal PARAM_DESCRIPTION
#' @param m_pal PARAM_DESCRIPTION
#' @param named_k_pal PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @export 
anno_clusters1 <- function (x, k, m, k_pal, m_pal, named_k_pal = NULL) 
{
    kids <- levels(x$cluster_id)
    nk <- length(kids)
    if (!is.null(named_k_pal)) {
        k_pal <- named_k_pal[kids]
    } else {
        if (nk > length(k_pal)) 
            k_pal <- colorRampPalette(k_pal)(nk)
        k_pal <- k_pal[seq_len(nk)]
        names(k_pal) <- kids
    }
    df <- data.frame(cluster_id = kids)
    col <- list(cluster_id = k_pal)
    colout <<- col
    if (!is.null(m)) {
        i <- match(kids, cluster_codes(x)[, k])
        mids <- droplevels(cluster_codes(x)[, m][i])
        nm <- nlevels(mids)
        if (nm > length(m_pal)) 
            m_pal <- colorRampPalette(m_pal)(nm)
        m_pal <- m_pal[seq_len(nm)]
        names(m_pal) <- levels(mids)
        df$merging_id <- mids
        col$merging_id <- m_pal
    }
    df <- mutate_all(df, function(u) factor(u, unique(u)))
    rowAnnotation(df = df, col = col, gp = gpar(col = "white"))
}


## ----------------------------------------------------------------------------------------------------------

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param features PARAM_DESCRIPTION, Default: NULL
#' @param by PARAM_DESCRIPTION, Default: c("sample_id", "cluster_id", "both")
#' @param k PARAM_DESCRIPTION, Default: 'meta20'
#' @param m PARAM_DESCRIPTION, Default: NULL
#' @param assay PARAM_DESCRIPTION, Default: 'exprs'
#' @param fun PARAM_DESCRIPTION, Default: c("median", "mean", "sum")
#' @param scale PARAM_DESCRIPTION, Default: c("first", "last", "never")
#' @param q PARAM_DESCRIPTION, Default: 0.01
#' @param sub PARAM_DESCRIPTION, Default: FALSE
#' @param row_anno PARAM_DESCRIPTION, Default: TRUE
#' @param row_names PARAM_DESCRIPTION, Default: FALSE
#' @param col_anno PARAM_DESCRIPTION, Default: TRUE
#' @param row_clust PARAM_DESCRIPTION, Default: TRUE
#' @param col_clust PARAM_DESCRIPTION, Default: TRUE
#' @param row_dend PARAM_DESCRIPTION, Default: TRUE
#' @param col_dend PARAM_DESCRIPTION, Default: TRUE
#' @param bars PARAM_DESCRIPTION, Default: FALSE
#' @param perc PARAM_DESCRIPTION, Default: FALSE
#' @param title PARAM_DESCRIPTION, Default: ' '
#' @param bin_anno PARAM_DESCRIPTION, Default: FALSE
#' @param hm_pal PARAM_DESCRIPTION, Default: rev(brewer.pal(11, "RdYlBu"))
#' @param k_pal PARAM_DESCRIPTION, Default: CATALYST:::.cluster_cols
#' @param m_pal PARAM_DESCRIPTION, Default: k_pal
#' @param named_k_pal PARAM_DESCRIPTION, Default: NULL
#' @param distance PARAM_DESCRIPTION, Default: c("euclidean", "maximum", "manhattan", "canberra", "binary", 
#'    "minkowski")
#' @param linkage PARAM_DESCRIPTION, Default: c("average", "ward.D", "single", "complete", "mcquitty", "median", 
#'    "centroid", "ward.D2")
#' @param plot_m_clusters PARAM_DESCRIPTION, Default: T
#' @return OUTPUT_DESCRIPTION
#' @export 
#' @importFrom CATALYST .cluster_cols .get_features .check_k .agg .scale_exprs .anno_factors
plotExprHeatmap1 <- function (x, features = NULL, by = c("sample_id", "cluster_id", 
    "both"), k = "meta20", m = NULL, assay = "exprs", fun = c("median", 
    "mean", "sum"), scale = c("first", "last", "never"), q = 0.01, sub = FALSE,
    row_anno = TRUE, row_names = FALSE, col_anno = TRUE, row_clust = TRUE, col_clust = TRUE, 
    row_dend = TRUE, col_dend = TRUE, bars = FALSE, perc = FALSE, title = " ",
    bin_anno = FALSE, hm_pal = rev(brewer.pal(11, "RdYlBu")), 
    k_pal = CATALYST:::.cluster_cols, m_pal = k_pal, named_k_pal = NULL, distance = c("euclidean", 
        "maximum", "manhattan", "canberra", "binary", "minkowski"), 
    linkage = c("average", "ward.D", "single", "complete", "mcquitty", 
        "median", "centroid", "ward.D2"), plot_m_clusters = T) 
{
  
  # 
  if (sub == TRUE)
  {
    x <- subSCE(x, n_cells = 1000)
    title <- "WARNING: downsampled SCE used."
  }
  
    args <- as.list(environment())
    # CATALYST:::.check_args_plotExprHeatmap(args)
    distance <- match.arg(distance)
    linkage <- match.arg(linkage)
    scale <- match.arg(scale)
    fun <- match.arg(fun)
    by <- match.arg(by)
    x <- x[unique(CATALYST:::.get_features(x, features)), ]
    if (by != "sample_id") {
       CATALYST:::.check_k(x, k)
        x$cluster_id <- cluster_ids(x, k)
        #x$cluster_id <- x[[k]]
    }
    if (by == "both") 
        by <- c("cluster_id", "sample_id")
    .do_agg <- function() {
        z <- CATALYST:::.agg(x, by, fun, assay)
        if (length(by) == 1) 
            return(z)
        set_rownames(do.call("rbind", z), levels(x$cluster_id))
    }
    .do_scale <- function() {
        if (scale == "first") {
            z <- assay(x, assay)
            z <- CATALYST:::.scale_exprs(z, 1, q)
            assay(x, assay, FALSE) <- z
            return(x)
        }
        else CATALYST:::.scale_exprs(z, 1, q)
    }
    z <- switch(scale, first = {
        x <- .do_scale()
        .do_agg()
    }, last = {
        z <- .do_agg()
        .do_scale()
    }, never = {
        .do_agg()
    })
    if (length(by) == 1) 
        z <- t(z)
    if (scale != "never" && !(assay == "counts" && fun == "sum")) {
        qs <- round(quantile(z, c(0.01, 0.99)) * 5)/5
        lgd_aes <- list(at = seq(qs[1], qs[2], 0.2))
    }
    else lgd_aes <- list()
    lgd_aes$title_gp <- gpar(fontsize = 10, fontface = "bold", 
        lineheight = 0.8)
    if (!isFALSE(row_anno)) {
  if (plot_m_clusters) {
      left_anno <- switch(by[1], sample_id = CATALYST:::.anno_factors(x, levels(x$sample_id), row_anno, "row"), anno_clusters1(x, k, m, k_pal, m_pal, named_k_pal))
    } else {
      left_anno <- switch(by[1], sample_id = CATALYST:::.anno_factors(x, levels(x$sample_id), row_anno, "row"), anno_clusters1(x, k, NULL, k_pal, m_pal, named_k_pal))
    }
        catalyst_merge <<- merging_ids(x, k, m, k_pal, m_pal) # This will edit a global variable in the workspace with catalyst's choice for merging.
    }
    else left_anno <- NULL
    if (!isFALSE(col_anno) && length(by) == 2) {
        top_anno <- CATALYST:::.anno_factors(x, levels(x$sample_id), col_anno, 
            "colum")
    }
    else top_anno <- NULL
    if (bars) {
        right_anno <- anno_counts1(x[[by[1]]], perc)
    }
    else right_anno <- NULL
    if (bin_anno) {
        cell_fun <- function(j, i, x, y, ...) grid.text(gp = gpar(fontsize = 8), 
            sprintf("%.2f", z[i, j]), x, y)
    }
    else cell_fun <- NULL
    a <- ifelse(assay == "exprs", "expression", assay)
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param p PARAM_DESCRIPTION
#' @param w PARAM_DESCRIPTION, Default: 20
#' @param h PARAM_DESCRIPTION, Default: 20
#' @param title PARAM_DESCRIPTION, Default: 'your_title'
#' @return OUTPUT_DESCRIPTION
#' @export 
#' @importFrom tools file_path_sans_ext
    f <- switch(fun, median = "med", fun)
    hm_title <- switch(scale, first = sprintf("%s %s\n%s", fun, 
        "scaled", a), last = sprintf("%s %s\n%s", "scaled", fun, 
        a), never = paste(fun, a, sep = "\n"))
    if (length(by) == 2) {
        col_title <- features
    }
    else if (length(features) == 1 && features %in% c("type", 
        "state")) {
        col_title <- paste0(features, "_markers")
    }
    else col_title <- ""
    
        if (!is.null(left_anno)) {
      heatmap_palette <<- left_anno@anno_list$cluster_id@color_mapping@colors
    }
  
    cn <- colnames(z) # Get the colunm names (names of each channel)
    
    p <- Heatmap(matrix = z, name = hm_title, 
        col = colorRamp2(seq(min(z), 
        max(z), l = n <- 100), 
        colorRampPalette(hm_pal)(n)), 
        column_title = col_title, column_title_side = ifelse(length(by) == 2, "top", "bottom"), 
        cell_fun = cell_fun, 
        cluster_rows = row_clust, 
        cluster_columns = col_clust, 
        show_row_dend = row_dend, 
        show_column_dend = col_dend, 
        clustering_distance_rows = distance, 
        clustering_method_rows = linkage, 
        clustering_distance_columns = distance, 
        clustering_method_columns = linkage, 
        show_row_names = row_names, # whether to show row names
        row_names_side = ifelse(by[1] == "cluster_id" || isFALSE(row_anno) && !row_dend || isFALSE(row_clust), "left", "right"), #which side for row names
        top_annotation = top_anno,  # Not sure what top annotation is
        left_annotation = left_anno, # Left annotation is the colours for clusters and merging (not row names)
        right_annotation = right_anno, # Right annotation is the bar plots and %
        rect_gp = gpar(col = "white"), # How to draw heatmap body colours
        heatmap_legend_param = lgd_aes,
        row_names_max_width = unit(13, "cm"),
        row_title = title,
        column_names_rot = 45,
        row_names_gp = gpar(fontsize = 12)) # Font size for row names (clusters)

  # Rotate the merging etc annotations too, but only if they exist
  if (!is.null(left_anno)) {
    p@left_annotation@anno_list[[1]]@name_param$rot = 45
    if (plot_m_clusters && length(p@left_annotation@anno_list) > 1) {
      p@left_annotation@anno_list[[2]]@name_param$rot = 45
    }
  }
    
    # Draw the heatmap
    ht <- draw(p)

    # Get the row labels in the order they appear on the clustered heatmap
    row_labels_after_clustering <- rownames(z)[row_order(ht)]
    row_labels_df <- data.frame(row_labels = row_labels_after_clustering)

    # Return both the heatmap object and the dataframe containing row labels
    list(heatmap = p, row_labels_df = row_labels_df)
}

environment(plotExprHeatmap1) <- asNamespace('CATALYST')


## ----------------------------------------------------------------------------------------------------------

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param k PARAM_DESCRIPTION
#' @param k_pal PARAM_DESCRIPTION
#' @param named_k_pal PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @export 
anno_clusters22 <- function (x, k, k_pal, named_k_pal = NULL) 
{
    kids <- levels(x$cluster_id)
    nk <- length(kids)
    if (!is.null(named_k_pal)) {
        k_pal <- named_k_pal[kids]
    } else {
        if (nk > length(k_pal)) 
            k_pal <- colorRampPalette(k_pal)(nk)
        k_pal <- k_pal[seq_len(nk)]
        names(k_pal) <- kids
    }
    df <- data.frame(cluster_id = kids)
    col <- list(cluster_id = k_pal)
    df <- mutate_all(df, function(u) factor(u, unique(u)))
    rowAnnotation(df = df, col = col, gp = gpar(col = "white"))
}


## ----------------------------------------------------------------------------------------------------------
library(ComplexHeatmap)
library(grid)
library(dplyr)

# Custom annotation function to draw circles with the correct colors
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param index PARAM_DESCRIPTION
#' @param k_pal PARAM_DESCRIPTION
#' @param cluster_order PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @export 
circle_anno <- function(index, k_pal, cluster_order) {
  n <- length(index)
  # Use the cluster order to get the correct color for each row
  colors <- k_pal[cluster_order[index]]
  grid.points(x = rep(0.5, n), 
              y = unit((n:1 - 0.5) / n, "npc"), 
              pch = 16, 
              size = unit(8, "mm"), 
              gp = gpar(col = colors))
}



## ----------------------------------------------------------------------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param cluster_order PARAM_DESCRIPTION
#' @param k_pal PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @export 
anno_clusters2 <- function(cluster_order, k_pal) {
  # Create rowAnnotation with the correct color mapping
  rowAnnotation(
    annotation_function = function(index) {
      circle_anno(index, k_pal, cluster_order)
    },
    show_annotation_name = FALSE, 
    gp = gpar(col = "white")
  )
}


## ----------------------------------------------------------------------------------------------------------

# Plot heatmap function
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param features PARAM_DESCRIPTION, Default: NULL
#' @param by PARAM_DESCRIPTION, Default: c("sample_id", "cluster_id", "both")
#' @param k PARAM_DESCRIPTION, Default: 'meta20'
#' @param assay PARAM_DESCRIPTION, Default: 'exprs'
#' @param fun PARAM_DESCRIPTION, Default: c("median", "mean", "sum")
#' @param scale PARAM_DESCRIPTION, Default: c("first", "last", "never")
#' @param q PARAM_DESCRIPTION, Default: 0.01
#' @param sub PARAM_DESCRIPTION, Default: FALSE
#' @param row_anno PARAM_DESCRIPTION, Default: TRUE
#' @param row_names PARAM_DESCRIPTION, Default: FALSE
#' @param col_anno PARAM_DESCRIPTION, Default: TRUE
#' @param row_clust PARAM_DESCRIPTION, Default: TRUE
#' @param col_clust PARAM_DESCRIPTION, Default: TRUE
#' @param row_dend PARAM_DESCRIPTION, Default: TRUE
#' @param col_dend PARAM_DESCRIPTION, Default: TRUE
#' @param bars PARAM_DESCRIPTION, Default: FALSE
#' @param perc PARAM_DESCRIPTION, Default: FALSE
#' @param title PARAM_DESCRIPTION, Default: ' '
#' @param bin_anno PARAM_DESCRIPTION, Default: FALSE
#' @param hm_pal PARAM_DESCRIPTION, Default: rev(brewer.pal(11, "RdYlBu"))
#' @param k_pal PARAM_DESCRIPTION, Default: NULL
#' @param named_k_pal PARAM_DESCRIPTION, Default: NULL
#' @param distance PARAM_DESCRIPTION, Default: c("euclidean", "maximum", "manhattan", "canberra", "binary", 
#'    "minkowski")
#' @param linkage PARAM_DESCRIPTION, Default: c("average", "ward.D", "single", "complete", "mcquitty", "median", 
#'    "centroid", "ward.D2")
#' @param row_order PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @export 
#' @importFrom CATALYST .get_features .agg .scale_exprs .anno_factors
plotExprHeatmapCol <- function (x, features = NULL, by = c("sample_id", "cluster_id", 
    "both"), k = "meta20", assay = "exprs", fun = c("median", 
    "mean", "sum"), scale = c("first", "last", "never"), q = 0.01, sub = FALSE,
    row_anno = TRUE, row_names = FALSE, col_anno = TRUE, row_clust = TRUE, col_clust = TRUE, 
    row_dend = TRUE, col_dend = TRUE, bars = FALSE, perc = FALSE, title = " ",
    bin_anno = FALSE, hm_pal = rev(brewer.pal(11, "RdYlBu")), 
    k_pal = NULL, named_k_pal = NULL, distance = c("euclidean", 
        "maximum", "manhattan", "canberra", "binary", "minkowski"), 
    linkage = c("average", "ward.D", "single", "complete", "mcquitty", 
        "median", "centroid", "ward.D2"), row_order = NULL) 
{
    if (sub == TRUE)
    {
      x <- subSCE(x, n_cells = 1000)
      title <- "WARNING: downsampled SCE used."
    }

    args <- as.list(environment())
    distance <- match.arg(distance)
    linkage <- match.arg(linkage)
    scale <- match.arg(scale)
    fun <- match.arg(fun)
    by <- match.arg(by)
    x <- x[unique(CATALYST:::.get_features(x, features)), ]

    if (by != "sample_id") {
        x$cluster_id <- x[[k]]
    }
    if (by == "both") 
        by <- c("cluster_id", "sample_id")
    
    .do_agg <- function() {
        z <- CATALYST:::.agg(x, by, fun, assay)
        if (length(by) == 1) 
            return(z)
        set_rownames(do.call("rbind", z), levels(x$cluster_id))
    }
    .do_scale <- function() {
        if (scale == "first") {
            z <- assay(x, assay)
            z <- CATALYST:::.scale_exprs(z, 1, q)
            assay(x, assay, FALSE) <- z
            return(x)
        }
        else CATALYST:::.scale_exprs(z, 1, q)
    }
    z <- switch(scale, first = {
        x <- .do_scale()
        .do_agg()
    }, last = {
        z <- .do_agg()
        .do_scale()
    }, never = {
        .do_agg()
    })

    if (length(by) == 1) 
        z <- t(z)
    if (scale != "never" && !(assay == "counts" && fun == "sum")) {
        qs <- round(quantile(z, c(0.01, 0.99)) * 5)/5
        lgd_aes <- list(at = seq(qs[1], qs[2], 0.2))
    }
    else lgd_aes <- list()
    lgd_aes$title_gp <- gpar(fontsize = 10, fontface = "bold", 
        lineheight = 0.8)
    
    # Apply custom row order
    if (!is.null(row_order)) {
      z <- z[row_order, , drop = FALSE]
    }

    # Extract the cluster IDs in the order of the heatmap rows
    cluster_order <- rownames(z)

    # Ensure named_k_pal is provided and matches the cluster_order
    if (!is.null(named_k_pal)) {
        k_pal <- named_k_pal[cluster_order]
    } else {
        stop("Please provide a named color palette 'named_k_pal' matching your clusters.")
    }

    # Create the left annotation
    if (!isFALSE(row_anno)) {
        left_anno <- anno_clusters2(cluster_order, k_pal)
    } else {
        left_anno <- NULL
    }

    if (!isFALSE(col_anno) && length(by) == 2) {
        top_anno <- CATALYST:::.anno_factors(x, levels(x$sample_id), col_anno, 
            "column")
    } else {
        top_anno <- NULL
    }

    if (bars) {
        right_anno <- anno_counts1(x[[by[1]]], perc)
    } else {
        right_anno <- NULL
    }

    if (bin_anno) {
        cell_fun <- function(j, i, x, y, ...) grid.text(gp = gpar(fontsize = 6), 
            sprintf("%.2f", z[i, j]), x, y)
    } else {
        cell_fun <- NULL
    }

    a <- ifelse(assay == "exprs", "expression", assay)
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param p PARAM_DESCRIPTION
#' @param w PARAM_DESCRIPTION, Default: 20
#' @param h PARAM_DESCRIPTION, Default: 20
#' @param title PARAM_DESCRIPTION, Default: 'your_title'
#' @return OUTPUT_DESCRIPTION
#' @export 
#' @importFrom tools file_path_sans_ext
    f <- switch(fun, median = "med", fun)
    hm_title <- switch(scale, first = sprintf("%s %s\n%s", fun, 
        "scaled", a), last = sprintf("%s %s\n%s", "scaled", fun, 
        a), never = paste(fun, a, sep = "\n"))
    if (length(by) == 2) {
        col_title <- features
    } else if (length(features) == 1 && features %in% c("type", 
        "state")) {
        col_title <- paste0(features, "_markers")
    } else {
        col_title <- ""
    }

    p <- Heatmap(matrix = z, name = hm_title, 
        col = colorRamp2(seq(min(z), 
        max(z), length.out = 100), 
        colorRampPalette(hm_pal)(100)), 
        column_title = col_title, column_title_side = ifelse(length(by) == 2, "top", "bottom"), 
        cell_fun = cell_fun, 
        cluster_rows = row_clust, 
        cluster_columns = col_clust, 
        show_row_dend = row_dend, 
        show_column_dend = col_dend, 
        clustering_distance_rows = distance, 
        clustering_method_rows = linkage, 
        clustering_distance_columns = distance, 
        clustering_method_columns = linkage, 
        show_row_names = row_names, 
        row_names_side = ifelse(by[1] == "cluster_id" || isFALSE(row_anno) && !row_dend || isFALSE(row_clust), "left", "right"), 
        top_annotation = top_anno, 
        left_annotation = left_anno, 
        right_annotation = right_anno, 
        rect_gp = gpar(col = "white"), 
        heatmap_legend_param = lgd_aes,
        row_names_max_width = unit(13, "cm"),
        row_title = title,
        column_names_rot = 45,
        column_names_gp = gpar(fontsize = 12),
        row_names_gp = gpar(fontsize = 14))

    ht <- draw(p)

    row_labels_after_clustering <- rownames(z)[row_order(ht)]
    row_labels_df <- data.frame(row_labels = row_labels_after_clustering)

    list(heatmap = p, row_labels_df = row_labels_df)
}

environment(plotExprHeatmapCol) <- asNamespace('CATALYST')


## ----------------------------------------------------------------------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param sce PARAM_DESCRIPTION
#' @param cell_id_column PARAM_DESCRIPTION
#' @param search_terms PARAM_DESCRIPTION, Default: NULL
#' @param replacement_terms PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @export 
create_superscript_column <- function(sce, cell_id_column, search_terms = NULL, replacement_terms = NULL) {

  # Convert factor to character if necessary
  cell_ids <- as.character(sce[[cell_id_column]])
  
  # Check for any NA values and warn if found
  if (any(is.na(cell_ids))) {
    warning("There are NA values in the cell IDs. They will not be modified.")
  }

  print("Original Cell IDs:")
  print(unique(cell_ids))  # Debug: Check the original unique values

  # Convert minus and plus to superscript versions by default
  cell_ids <- gsub("-", "^'−'*", cell_ids)   # Superscript minus sign from hyphen
  cell_ids <- gsub("\\+", "^'+'*", cell_ids) # Superscript plus sign
  cell_ids <- gsub(" ", "~", cell_ids)       # Convert spaces to tildes

  print("Cell IDs after default superscripting:")
  print(unique(cell_ids))  # Debug: Check after minus and plus transformation

  # Apply additional search and replace terms if provided
  if (!is.null(search_terms)) {
    for (i in seq_along(search_terms)) {
      print(paste0("Superscripting term: ", search_terms[[i]], " to: ^'", replacement_terms[[i]], "'*"))
      cell_ids <- gsub(search_terms[[i]], paste0("^'", replacement_terms[[i]], "'*"), cell_ids)
    }
  }

  # Remove any asterisks (*) at the end of the string
  cell_ids <- gsub("\\*+$", "", cell_ids)

  print("Final Cell IDs:")
  print(unique(cell_ids))  # Debug: Check final transformed IDs

  # Add the new superscripted column to the sce object
  sce[[paste0(cell_id_column, "_superscript")]] <- cell_ids

  return(sce)
}




## ----------------------------------------------------------------------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param sce PARAM_DESCRIPTION
#' @param k PARAM_DESCRIPTION, Default: 'merging1'
#' @param meta PARAM_DESCRIPTION, Default: c("patient_id")
#' @param meta2 PARAM_DESCRIPTION, Default: 'sample_id ~ cluster_id'
#' @param ncells PARAM_DESCRIPTION, Default: NULL
#' @param merging_col PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @export 
#' @importFrom CATALYST .check_k
sceget <- function(sce, k = "merging1", meta = c("patient_id"), meta2 = "sample_id ~ cluster_id",ncells = NULL, merging_col = NULL)
  
{
  
  # Filter the sce by number of cells if desired
if(!is.null(ncells)){
  sce <- filterSCEvents(sce, n_cells = ncells)
}
  
if (!is.null(merging_col)) { 
cluster_ids <- sce[[k]] 
} else { 
k <- CATALYST:::.check_k(sce, k) 
cluster_ids <- cluster_ids(sce, k) 
}
  
# Make a data frame containing cluster abundances
ns <- table(cluster_id = cluster_ids, sample_id = sample_ids(sce))
fq <- prop.table(ns, 2) * 100
df <- as.data.frame(fq)
m <- match(df$sample_id, sce$sample_id)

# Extract any other metadata into the df data frame (replace the variables in c() with your own choices).
for (i in meta) df[[i]] <- sce[[i]][m]

#dfout <<- df

# dcast to get rows as patient ids
df = dcast(df, meta2, value.var = "Freq")

#dfout2 <<- df
# put it in the same order as ei(sce)  !! Important for diffcyt!  I'm not doing this for now!
#sampids <- ei(sce)

#df <- df[match(sampids$sample_id, df$sample_id),]

return(df)
}


## ----------------------------------------------------------------------------------------------------------

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param sce PARAM_DESCRIPTION
#' @param k PARAM_DESCRIPTION, Default: 'merging1'
#' @param meta PARAM_DESCRIPTION, Default: c("patient_id")
#' @return OUTPUT_DESCRIPTION
#' @export 
sceget2 <- function(sce, k = "merging1", meta = c("patient_id"))
  
{
  
# Make a data frame containing cluster abundances
ns <- table(cluster_id = cluster_ids(sce, k), sample_id = sample_ids(sce))
fq <- prop.table(ns, 2) * 100
df <- as.data.frame(fq)
m <- match(df$sample_id, sce$sample_id)

# Extract any other metadata into the df data frame (replace the variables in c() with your own choices).
for (i in meta) df[[i]] <- sce[[i]][m]

return(df)
}


## ----------------------------------------------------------------------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param sce PARAM_DESCRIPTION
#' @param filterColumn PARAM_DESCRIPTION
#' @param filterValues PARAM_DESCRIPTION
#' @param exclude PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @export 
filterSCE_simple <- function(sce, filterColumn, filterValues, exclude = FALSE) {
  # Get the condition
  if(exclude) {
    condition <- !sce[[filterColumn]] %in% filterValues
  } else {
    condition <- sce[[filterColumn]] %in% filterValues
  }

  # Filter the sce object
  sce_filtered <- sce[, condition]

  return(sce_filtered)
}




## ----------------------------------------------------------------------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param sce PARAM_DESCRIPTION
#' @param filterColumn PARAM_DESCRIPTION
#' @param filterValues PARAM_DESCRIPTION
#' @param exclude PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @export 
filterSCE_simple2 <- function(sce, filterColumn, filterValues, exclude = FALSE) {
  # Get the condition
  if(exclude) {
    condition <- !sce[[filterColumn]] %in% filterValues
  } else {
    condition <- sce[[filterColumn]] %in% filterValues
  }

  # Filter the sce object
  sce_filtered <- sce[, condition]
  
  # Drop unused levels in all factor columns in the colData
  is_factor <- sapply(colData(sce_filtered), is.factor)
  colData(sce_filtered)[is_factor] <- lapply(colData(sce_filtered)[is_factor], droplevels)

  return(sce_filtered)
}



## ----------------------------------------------------------------------------------------------------------

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param sce PARAM_DESCRIPTION
#' @param days_threshold PARAM_DESCRIPTION, Default: 24
#' @return OUTPUT_DESCRIPTION
#' @export 
filterSceByDays <- function(sce, days_threshold = 24) {
  # Convert the day column to numeric
  sce$day <- as.numeric(as.character(sce$day))
  
  # Identify patients who have at least one sample recorded beyond the specified day threshold
  eligible_patients <- unique(sce$patient_id[sce$day > days_threshold])
  
  # Filter the sce object to include only samples from eligible patients
  sce_filtered <- sce[, sce$patient_id %in% eligible_patients]
  
  return(sce_filtered)
}




## ----------------------------------------------------------------------------------------------------------
# Experiment info 2 function.  Kind of hacked together, but it works at generating the metadata from the base data
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param sce PARAM_DESCRIPTION
#' @param meta PARAM_DESCRIPTION, Default: c("patient_id")
#' @return OUTPUT_DESCRIPTION
#' @export 
ei2 <- function(sce, meta = c("patient_id"))
  
{
# Make a data frame of sample_ids
ns <- table(sample_id = sample_ids(sce))
df <- as.data.frame(ns)
m <- match(df$sample_id, sce$sample_id)

# Extract metadata
for (i in meta) df[[i]] <- sce[[i]][m]

return(df)
}



## ----------------------------------------------------------------------------------------------------------

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param dr PARAM_DESCRIPTION, Default: NULL
#' @param textsize PARAM_DESCRIPTION, Default: 18
#' @param legendpointsize PARAM_DESCRIPTION, Default: 7
#' @param color_by PARAM_DESCRIPTION, Default: 'condition'
#' @param facet_by PARAM_DESCRIPTION, Default: NULL
#' @param hide_axis PARAM_DESCRIPTION, Default: FALSE
#' @param border_width PARAM_DESCRIPTION, Default: 1
#' @param ncol PARAM_DESCRIPTION, Default: NULL
#' @param pointsize PARAM_DESCRIPTION, Default: 0.4
#' @param assay PARAM_DESCRIPTION, Default: 'exprs'
#' @param scale PARAM_DESCRIPTION, Default: TRUE
#' @param random_order PARAM_DESCRIPTION, Default: FALSE
#' @param q PARAM_DESCRIPTION, Default: 0.01
#' @param dims PARAM_DESCRIPTION, Default: c(1, 2)
#' @param alpha PARAM_DESCRIPTION, Default: 0.8
#' @param k_pal PARAM_DESCRIPTION, Default: CATALYST:::.cluster_cols
#' @param a_pal PARAM_DESCRIPTION, Default: hcl.colors(10, "Viridis")
#' @param rast PARAM_DESCRIPTION, Default: FALSE
#' @param panel_spacing PARAM_DESCRIPTION, Default: 1
#' @param plot_order PARAM_DESCRIPTION, Default: NULL
#' @param highlight_cluster PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @export 
#' @importFrom CATALYST .cluster_cols .check_assay .check_pal .check_cd_factor .check_sce .check_k .scale_exprs
plotDR1 <- function (x, dr = NULL, textsize = 18, legendpointsize = 7, color_by = "condition", 
                     facet_by = NULL, hide_axis = FALSE, border_width = 1, ncol = NULL, pointsize = 0.4,
                     assay = "exprs", scale = TRUE, random_order = FALSE, q = 0.01, dims = c(1,2), alpha = 0.8,
                     k_pal = CATALYST:::.cluster_cols, a_pal = hcl.colors(10, "Viridis"), rast = FALSE, 
                     panel_spacing = 1, plot_order = NULL, highlight_cluster = NULL) 
{
    stopifnot(is(x, "SingleCellExperiment"), CATALYST:::.check_assay(x, assay), 
              length(reducedDims(x)) != 0, is.logical(scale), length(scale) == 1, 
              is.numeric(q), length(q) == 1, q >= 0, q < 0.5)

    CATALYST:::.check_pal(a_pal)
    CATALYST:::.check_cd_factor(x, facet_by)

    if (!is.null(ncol)) stopifnot(is.numeric(ncol), length(ncol) == 1, ncol%%1 == 0)

    if (is.null(dr)) {
        dr <- reducedDimNames(x)[1]
    } else {
        stopifnot(is.character(dr), length(dr) == 1, dr %in% reducedDimNames(x))
    }

    stopifnot(is.numeric(dims), length(dims) == 2, dims %in% seq_len(ncol(reducedDim(x, dr))))

    if (!all(color_by %in% rownames(x))) {
        stopifnot(length(color_by) == 1)
        if (!color_by %in% names(colData(x))) {
            CATALYST:::.check_sce(x, TRUE)
            CATALYST:::.check_pal(k_pal)
            CATALYST:::.check_k(x, color_by)
            kids <- cluster_ids(x, color_by)
            nk <- nlevels(kids)
            if (length(k_pal) < nk) k_pal <- colorRampPalette(k_pal)(nk)
            
            # Save the color palette to a global variable
            plotDR1_colorpal <<- k_pal
            names(plotDR1_colorpal) <<- levels(kids)
        } else kids <- NULL
    }

    xy <- reducedDim(x, dr)[, dims]
    colnames(xy) <- c("x", "y")
    df <- data.frame(colData(x), xy, check.names = FALSE)

    if (all(color_by %in% rownames(x))) {
        es <- as.matrix(assay(x, assay))
        es <- es[color_by, , drop = FALSE]
        if (scale) es <- CATALYST:::.scale_exprs(es, 1, q)
        df <- melt(cbind(df, t(es)), id.vars = colnames(df))
        l <- switch(assay, exprs = "expression", assay)
        l <- paste0("scaled\n"[scale], l)
        scale <- scale_colour_gradientn(l, colors = a_pal)
        thm <- guide <- NULL
        color_by <- "value"
        facet <- facet_wrap("variable", ncol = ncol)
    } else if (is.numeric(df[[color_by]])) {
        if (scale) {
            vs <- as.matrix(df[[color_by]])
            df[[color_by]] <- CATALYST:::.scale_exprs(vs, 2, q)
        }
        l <- paste0("scaled\n"[scale], color_by)
        scale <- scale_colour_gradientn(l, colors = a_pal)
        color_by <- sprintf("`%s`", color_by)
        facet <- thm <- guide <- NULL
    } else {
        facet <- NULL
        if (!is.null(kids)) {
            df[[color_by]] <- kids
            scale <- scale_color_manual(values = k_pal)
        } else scale <- NULL

        n <- nlevels(droplevels(factor(df[[color_by]])))
        guide <- guides(col = guide_legend(ncol = ifelse(n > 12, 2, 1), 
                                           override.aes = list(alpha = 1, size = legendpointsize)))
        thm <- theme(legend.key.height = unit(0.8, "lines"), text = element_text(size = textsize)) 
    }

    if (dr %in% c("PCA", "MDS")) {
        asp <- coord_equal()
    } else asp <- NULL

    if (dr == "PCA") {
        labs <- paste0("PC", dims)
    } else labs <- paste(dr, "dim.", dims)
    
    df <- df[!(is.na(df$x) | is.na(df$y)), ]
    
    # Apply plotting order based on color_by grouping
    if (!is.null(plot_order)) {
        # Ensure df[[color_by]] is a factor with levels specified by plot_order
        df[[color_by]] <- factor(df[[color_by]], levels = plot_order)
        # Order df according to the factor levels of color_by
        df <- df[order(df[[color_by]]), ]
    } else if (random_order) {
        # Shuffle rows of df if random_order is TRUE
        set.seed(1)
        df <- df[sample(nrow(df)), ]
    }

    library(ggrastr)
    p <- ggplot(df, aes_string("x", "y", col = color_by))  
    
    if(rast){
      p <- p + geom_point_rast(size = pointsize, alpha = alpha, shape = 16, raster.dpi = 600)
    } else {
      p <- p + geom_point(size = pointsize, alpha = alpha, shape = 16)
    }
    
    # Add highlighted cluster with density contours
    if (!is.null(highlight_cluster)) {
        # Plot points for the highlighted cluster without mapping 'fill' to 'color_by'
        p <- p + geom_point(
            data = df[df[[color_by]] == highlight_cluster, ],
            aes(x = x, y = y),
            size = 2, alpha = 0.8, shape = 21, color = "black", fill = "yellow"
        ) +
        # Add density contours for the highlighted cluster
        stat_density_2d(
            data = df[df[[color_by]] == highlight_cluster, ],
            aes(x = x, y = y, fill = after_stat(level)),
            geom = "polygon", color = "black", alpha = 0.3
        )
    }

    p <- p + labs(x = labs[1], y = labs[2]) + facet + scale + guide + asp + 
      theme_minimal() + 
      thm + 
      theme(panel.grid.minor = element_blank(), 
            strip.text = element_text(face = "bold"), 
            panel.grid.major = element_blank(),
            axis.text = element_text(color = "black"), 
            panel.spacing = unit(panel_spacing, "lines"),
            aspect.ratio = if (is.null(asp)) 1 else NULL)
    
    p <- p + theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=border_width))

    if (hide_axis == TRUE) {
      p <- p + theme(axis.text.x = element_blank(), 
                     axis.text.y = element_blank(),
                     axis.ticks.x = element_blank(), 
                     axis.ticks.y = element_blank())
    }
    
    if (is.null(facet_by)) return(p)
    
    if (is.null(facet)) {
        p + facet_wrap(facet_by, ncol = ncol)
    } else {
        if (nlevels(df$variable) == 1) {
            p + facet_wrap(facet_by, ncol = ncol) + ggtitle(levels(df$variable))
        } else {
            fs <- c("variable", facet_by)
            ns <- vapply(df[fs], nlevels, numeric(1))
            if (ns[2] > ns[1]) fs <- rev(fs)
            p + facet_grid(reformulate(fs[1], fs[2]))
        }
    }
}




## ----------------------------------------------------------------------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param dr PARAM_DESCRIPTION, Default: NULL
#' @param textsize PARAM_DESCRIPTION, Default: 18
#' @param legendpointsize PARAM_DESCRIPTION, Default: 7
#' @param color_by PARAM_DESCRIPTION, Default: 'condition'
#' @param facet_by PARAM_DESCRIPTION, Default: NULL
#' @param hide_axis PARAM_DESCRIPTION, Default: F
#' @param alpha PARAM_DESCRIPTION, Default: 0.8
#' @param pointsize PARAM_DESCRIPTION, Default: 0.4
#' @param ncol PARAM_DESCRIPTION, Default: NULL
#' @param assay PARAM_DESCRIPTION, Default: 'exprs'
#' @param scale PARAM_DESCRIPTION, Default: TRUE
#' @param random_order PARAM_DESCRIPTION, Default: FALSE
#' @param q PARAM_DESCRIPTION, Default: 0.01
#' @param dims PARAM_DESCRIPTION, Default: c(1, 2)
#' @param k_pal PARAM_DESCRIPTION, Default: CATALYST:::.cluster_cols
#' @param a_pal PARAM_DESCRIPTION, Default: hcl.colors(10, "Viridis")
#' @param rast PARAM_DESCRIPTION, Default: FALSE
#' @param panel_spacing PARAM_DESCRIPTION, Default: 1
#' @param parse_color_by PARAM_DESCRIPTION, Default: FALSE
#' @param legend_title PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @export 
#' @importFrom CATALYST .cluster_cols .check_assay .check_pal .check_cd_factor .check_sce .check_k .scale_exprs
plotDR2 <- function (x, dr = NULL, textsize = 18, legendpointsize = 7, color_by = "condition", facet_by = NULL, hide_axis = F, alpha = 0.8, pointsize = 0.4,
                     ncol = NULL, assay = "exprs", scale = TRUE, random_order = FALSE, q = 0.01, dims = c(1,2), k_pal = CATALYST:::.cluster_cols, a_pal = hcl.colors(10, "Viridis"), rast = FALSE, panel_spacing = 1, parse_color_by = FALSE, legend_title = NULL) 
{
  stopifnot(is(x, "SingleCellExperiment"), CATALYST:::.check_assay(x, assay), length(reducedDims(x)) != 0, 
            is.logical(scale), length(scale) == 1, is.numeric(q), length(q) == 1, q >= 0, q < 0.5)
  
  CATALYST:::.check_pal(a_pal)
  CATALYST:::.check_cd_factor(x, facet_by)
  
  if (!is.null(ncol)) stopifnot(is.numeric(ncol), length(ncol) == 1, ncol%%1 == 0)
  
  if (is.null(dr)) {
    dr <- reducedDimNames(x)[1]
  } else {
    stopifnot(is.character(dr), length(dr) == 1, dr %in% reducedDimNames(x))
  }
  
  stopifnot(is.numeric(dims), length(dims) == 2, dims %in% seq_len(ncol(reducedDim(x, dr))))
  
  if (!all(color_by %in% rownames(x))) {
    stopifnot(length(color_by) == 1)
    if (!color_by %in% names(colData(x))) {
      CATALYST:::.check_sce(x, TRUE)
      CATALYST:::.check_pal(k_pal)
      CATALYST:::.check_k(x, color_by)
      kids <- cluster_ids(x, color_by)
      nk <- nlevels(kids)
      if (length(k_pal) < nk) k_pal <- colorRampPalette(k_pal)(nk)
      
      # Save the color palette to a global variable
      plotDR1_colorpal <<- k_pal
      names(plotDR1_colorpal) <<- levels(kids)
    } else {
      # Clustering is in colData
      kids <- factor(colData(x)[[color_by]])
    }
  }
  
  xy <- reducedDim(x, dr)[, dims]
  colnames(xy) <- c("x", "y")
  df <- data.frame(colData(x), xy, check.names = FALSE)
  
  if (all(color_by %in% rownames(x))) {
    es <- as.matrix(assay(x, assay))
    es <- es[color_by, , drop = FALSE]
    if (scale) es <- CATALYST:::.scale_exprs(es, 1, q)
    df <- melt(cbind(df, t(es)), id.vars = colnames(df))
    l <- switch(assay, exprs = "expression", assay)
    l <- paste0("scaled\n"[scale], l)
    scale <- scale_colour_gradientn(l, colors = a_pal)
    thm <- guide <- NULL
    color_by <- "value"
    facet <- facet_wrap("variable", ncol = ncol)
  } else if (is.numeric(df[[color_by]])) {
    if (scale) {
      vs <- as.matrix(df[[color_by]])
      df[[color_by]] <- CATALYST:::.scale_exprs(vs, 2, q)
    }
    l <- paste0("scaled\n"[scale], color_by)
    scale <- scale_colour_gradientn(l, colors = a_pal)
    color_by <- sprintf("`%s`", color_by)
    facet <- thm <- guide <- NULL
  } else {
    facet <- NULL
    if (!is.null(kids)) {
      df[[color_by]] <- kids
      scale <- scale_color_manual(values = k_pal)
    } else scale <- NULL
    
    n <- nlevels(droplevels(factor(df[[color_by]])))
    guide <- guides(col = guide_legend(ncol = ifelse(n > 12, 2, 1), 
                                       override.aes = list(alpha = 1, size = legendpointsize)))
    thm <- theme(legend.key.height = unit(0.8, "lines"), text = element_text(size = textsize)) 
  }
  
  if (dr %in% c("PCA", "MDS")) {
    asp <- coord_equal()
  } else asp <- NULL
  
  if (dr == "PCA") {
    labs <- paste0("PC", dims)
  } else labs <- paste(dr, "dim.", dims)
  
  df <- df[!(is.na(df$x) | is.na(df$y)), ]
  
  # Shuffle rows of df if random_order is TRUE
  if (random_order) {
    df <- df[sample(nrow(df)), ]
  }
  
  library(ggrastr)
  p <- ggplot(df, aes_string("x", "y", col = color_by))  
  
  if(rast){
    p <- p + geom_point_rast(size = pointsize, alpha = alpha, shape = 16, raster.dpi = 600)
  } else {
    p <- p + geom_point(size = pointsize, alpha = alpha, shape = 16)
  }
  
  p <- p + labs(x = labs[1], y = labs[2]) + facet + scale + guide + asp + 
    theme_minimal() + thm + theme(panel.grid.minor = element_blank(), 
                                  strip.text = element_text(face = "bold"), 
                                  panel.grid.major = element_blank(),
                                  axis.text = element_text(color = "black"), 
                                  panel.spacing = unit(panel_spacing, "lines"),
                                  aspect.ratio = if (is.null(asp)) 1 else NULL)
  
  #  if(!is.null(facet_by)) {
  p <- p + theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
  # }
  
  if (hide_axis == TRUE) {
    p <- p + theme(axis.text.x = element_blank(), 
                   axis.text.y = element_blank(),
                   axis.ticks.x = element_blank(), 
                   axis.ticks.y = element_blank())
  }
  
  if (!is.null(legend_title)) {
    p <- p + guides(col = guide_legend(title = legend_title, 
                                       override.aes = list(alpha = 1, size = legendpointsize)))
  }
  
  if (parse_color_by) {
    p <- p + scale_color_manual(values = k_pal, labels = sapply(levels(df[[color_by]]), function(name) as.expression(parse(text = name))))
  }
  
  if (is.null(facet_by)) return(p)
  
  if (is.null(facet)) {
    p + facet_wrap(facet_by, ncol = ncol)
  } else {
    if (nlevels(df$variable) == 1) {
      p + facet_wrap(facet_by, ncol = ncol) + ggtitle(levels(df$variable))
    } else {
      fs <- c("variable", facet_by)
      ns <- vapply(df[fs], nlevels, numeric(1))
      if (ns[2] > ns[1]) fs <- rev(fs)
      p + facet_grid(reformulate(fs[1], fs[2]))
    }
  }
}


## ----------------------------------------------------------------------------------------------------------

# metafile is the name of the file containing updated metadata fields.  It needs to have sample_id
# metafields is array containing the names of the fields in metafile to be added to the sce

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param sce PARAM_DESCRIPTION
#' @param metafields PARAM_DESCRIPTION
#' @param md_u PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @export 
sce_meta <- function(sce, metafields,md_u)
  
{

  #https://stackoverflow.com/questions/35636315/replace-values-in-a-dataframe-based-on-lookup-table
  df = as.data.frame(sce@colData$sample_id)
  lookup = md_u
  
  for (i in 1:length(metafields)){
  
  field <- metafields[i]  
  #print(paste("Now updating: ",field))
  new <- df
  new[] <- lapply(df, function(x) lookup[[field]][match(x, lookup$sample_id)])
  sce[[field]] = as.factor(new$`sce@colData$sample_id`)
    
  }

  return(sce)
}


## ----------------------------------------------------------------------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param sce PARAM_DESCRIPTION
#' @param metafields PARAM_DESCRIPTION
#' @param md_u PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @export 
sce_meta2 <- function(sce, metafields, md_u) {

  # Create a data frame from the sample_id in sce
  df = as.data.frame(sce@colData$sample_id)
  
  # This is your lookup data frame
  lookup = md_u
  
  # Loop over each metafield that needs to be updated
  for (i in 1:length(metafields)) {
    
    # Get the field name
    field <- metafields[i]
    
    # Create a copy of the original dataframe
    new <- df
    
    # Match and update only if sample_id is found in lookup
    new[] <- lapply(df, function(x) {
      matched_values <- match(x, lookup$sample_id)
      ifelse(is.na(matched_values), NA, lookup[[field]][matched_values])
    })
    
    # Check if field already exists, if not create it, otherwise update it
    if (field %in% colnames(sce@colData)) {
      sce@colData[[field]] <- as.factor(new$`sce@colData$sample_id`)
    } else {
      sce[[field]] <- as.factor(new$`sce@colData$sample_id`)
    }
  }
  
  # Return the updated sce
  return(sce)
}



## ----------------------------------------------------------------------------------------------------------

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param k PARAM_DESCRIPTION, Default: 'meta20'
#' @param features PARAM_DESCRIPTION, Default: 'type'
#' @param clust PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @export 
plotClusterExprs1 <- function (x, k = "meta20", features = "type", clust = FALSE) 
{
    .check_sce(x, TRUE)
    k <- .check_k(x, k)
    x$cluster_id <- cluster_ids(x, k)
    features <- .get_features(x, features)
    ms <- t(.agg(x[features, ], "cluster_id", "median"))
    d <- dist(ms, method = "euclidean")
    o <- hclust(d, method = "average")$order
    cd <- colData(x)
    es <- assay(x[features, ], "exprs")
    df <- data.frame(t(es), cd, check.names = FALSE)
    df <- melt(df, id.vars = names(cd), variable.name = "antigen", 
        value.name = "expression")
    df$avg <- "no"
    avg <- df
    avg$cluster_id <- "avg"
    avg$avg <- "yes"
    df <- rbind(df, avg)
    fq <- tabulate(x$cluster_id)/ncol(x)
    fq <- round(fq * 100, 2)
    names(fq) <- levels(x$cluster_id)
    if (clust) {
        df$cluster_id <- factor(df$cluster_id, levels = rev(c("avg", 
            levels(x$cluster_id)[o])), labels = rev(c("average", 
            paste0(names(fq), " (", fq, "%)")[o])))
    } else {
        df$cluster_id <- factor(df$cluster_id, levels = rev(c("avg", 
            levels(x$cluster_id))), labels = rev(c("average", 
            paste0(names(fq), " (", fq, "%)"))))
    }
    ggplot(df, aes_string(x = "expression", y = "cluster_id", 
        col = "avg", fill = "avg")) + facet_wrap(~antigen, scales = "free_x", 
        nrow = 2) + geom_density_ridges(alpha = 0.2) + theme_ridges() + 
        theme(legend.position = "none", strip.background = element_blank(), 
            strip.text = element_text(face = "bold")) + 
        scale_x_continuous(breaks = seq(min(df$expression, na.rm = TRUE), max(df$expression, na.rm = TRUE), by = 1))
}

environment(plotClusterExprs1) <- asNamespace('CATALYST')


## ----------------------------------------------------------------------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param k PARAM_DESCRIPTION, Default: 'meta20'
#' @param features PARAM_DESCRIPTION, Default: 'type'
#' @param cluster_rows PARAM_DESCRIPTION, Default: TRUE
#' @param add_border PARAM_DESCRIPTION, Default: TRUE
#' @param text_size PARAM_DESCRIPTION, Default: 12
#' @param heatmap_palette PARAM_DESCRIPTION, Default: NULL
#' @param y_shift PARAM_DESCRIPTION, Default: 0.01
#' @param y_scale PARAM_DESCRIPTION, Default: 1.1
#' @param panel_spacing PARAM_DESCRIPTION, Default: 0.4
#' @param alpha_amount PARAM_DESCRIPTION, Default: 0.6
#' @param merging_col PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @export 
plotClusterExprs2 <- function (x, k = "meta20", features = "type", cluster_rows = TRUE, add_border = TRUE, text_size = 12, heatmap_palette = NULL, y_shift = 0.01, y_scale = 1.1, panel_spacing = 0.4, alpha_amount = 0.6, merging_col = FALSE) 
{
    .check_sce(x, TRUE)
    
    # Use the merging column from colData if merging_col is TRUE
    if (merging_col) {
        if (!k %in% names(colData(x))) {
            stop("Clustering column not found in colData.")
        }
        cluster_ids <- x[[k]]
    } else {
        k <- .check_k(x, k)
        cluster_ids <- cluster_ids(x, k)
    }
    
    x$cluster_id <- cluster_ids
    features <- .get_features(x, features)
    ms <- t(.agg(x[features, ], "cluster_id", "median"))
    
    if (cluster_rows) {
        d <- dist(ms, method = "euclidean")
        o <- hclust(d, method = "average")$order
    } else {
        o <- seq_along(levels(x$cluster_id))
    }
    
    cd <- colData(x)
    es <- assay(x[features, ], "exprs")
    df <- data.frame(t(es), cd, check.names = FALSE)
    df <- melt(df, id.vars = names(cd), variable.name = "antigen", 
               value.name = "expression")
    
    fq <- tabulate(x$cluster_id)/ncol(x)
    fq <- round(fq * 100, 2)
    names(fq) <- levels(x$cluster_id)
    cluster_levels <- rev(levels(x$cluster_id)[o])
    cluster_labels <- rev(names(fq)[o])
    df$cluster_id <- factor(df$cluster_id, levels = cluster_levels, labels = cluster_labels)
    
    # If a heatmap_palette is provided, use it to fill the colors, otherwise use red
    if (!is.null(heatmap_palette)) {
        cluster_colors <- heatmap_palette[levels(df$cluster_id)]
        scale_fill <- scale_fill_manual(values = cluster_colors)
    } else {
        scale_fill <- scale_fill_manual(values = rep("red", length(levels(df$cluster_id))))
    }
    
    # Define the base plot
    plot <- ggplot(df, aes(x = expression, y = cluster_id, fill = cluster_id)) +
        facet_wrap(~antigen, scales = "free_x", nrow = 1) + 
        geom_density_ridges(alpha = alpha_amount, rel_min_height = 0.01, scale = y_scale) + 
        scale_fill +
        scale_x_continuous(breaks = seq(0, 10, by = 2), limits = c(-0.5, NA), expand = c(0, 0)) +
        scale_y_discrete(expand = expand_scale(mult = c(y_shift, 0.2))) +
        theme(
              legend.position = "none", 
              strip.background = element_blank(), 
              strip.text = element_text(face = "bold", color = "black", size = text_size),
              panel.spacing = unit(panel_spacing, "lines"),
              plot.margin = unit(c(0, 1, 0, 0), "lines"),
              panel.background = element_rect(fill = "transparent", color = NA), # Make panel background transparent
              panel.border = if (add_border) element_rect(color = "black", fill = NA, size = 0.8) else element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_line(color = "black"), # Set y major gridlines as black
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              axis.line.x = element_blank(), # Hide main x-axis line
              axis.line.y = element_blank(), # Hide main y-axis line
              axis.text.x = element_text(color = "black"),
              axis.text.y = element_text(color = "black"),
              axis.ticks.x = element_line(color = "black"),
              text = element_text(size = text_size))
    
    plot
}

environment(plotClusterExprs2) <- asNamespace('CATALYST')


## ----------------------------------------------------------------------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param spill_matrix PARAM_DESCRIPTION
#' @param panel_table PARAM_DESCRIPTION
#' @param source_marker PARAM_DESCRIPTION
#' @param dest_marker PARAM_DESCRIPTION
#' @param value PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @export 
update_spill_matrix <- function(spill_matrix, panel_table, source_marker, dest_marker, value) {
  
  # Get the corresponding channel names
  source_channel <- panel_table$fcs_colname[panel_table$antigen == source_marker]
  dest_channel <- panel_table$fcs_colname[panel_table$antigen == dest_marker]
  
  # Check if the source and destination markers are valid
  if (length(source_channel) < 1 | length(dest_channel) < 1) {
    stop(paste0("Could not map markers ", source_marker, " and/or ", dest_marker, " to channels."))
  }
  
  # Update the spill matrix
  spill_matrix[source_channel, dest_channel] <- value
  
  # Return the updated spill matrix
  return(spill_matrix)
}



## ----------------------------------------------------------------------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param spill_matrix PARAM_DESCRIPTION
#' @param panel_table PARAM_DESCRIPTION
#' @param source_channel PARAM_DESCRIPTION
#' @param dest_channel PARAM_DESCRIPTION
#' @param value PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @export 
update_spill_matrix2 <- function(spill_matrix, panel_table, source_channel, dest_channel, value) {
  
  # Check if the source and destination channels are valid
  if (!source_channel %in% panel_table$fcs_colname | !dest_channel %in% panel_table$fcs_colname) {
    stop(paste0("Could not find channels ", source_channel, " and/or ", dest_channel, " in panel_table."))
  }
  
  # Update the spill matrix
  spill_matrix[source_channel, dest_channel] <- value
  
  # Return the updated spill matrix
  return(spill_matrix)
}


## ----------------------------------------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(patchwork)
library(SingleCellExperiment)

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param pal PARAM_DESCRIPTION
#' @param title PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @export 
plot_colors <- function(pal, title) {
  if (is.null(names(pal))) {
    names(pal) <- paste0(seq_along(pal), ": ", pal)
  }
  
  df <- data.frame(
    name = names(pal),
    color = pal,
    stringsAsFactors = FALSE
  )
  
  ggplot(df, aes(x = 1, y = factor(name, levels = rev(name)), fill = color)) +
    geom_tile(width = 2) +  # Increase the width of the tiles
    scale_fill_identity() +
    geom_text(aes(label = name), color = "black", size = 4, hjust = 0.5) +
    theme_void() +
    theme(legend.position = "none") +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    coord_fixed(ratio = 1/10)  # Adjust the aspect ratio
}



## ----echo=FALSE--------------------------------------------------------------------------------------------

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param sce PARAM_DESCRIPTION
#' @param x_channel PARAM_DESCRIPTION
#' @param y_channel PARAM_DESCRIPTION
#' @param gate_name PARAM_DESCRIPTION
#' @param n_samples PARAM_DESCRIPTION, Default: NULL
#' @param plot_above_one PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @export 
getPolygonGate <- function(sce, x_channel, y_channel, gate_name, n_samples = NULL, plot_above_one = FALSE) {
 # Extract data from the SingleCellExperiment object
  library(sp)
  library(shiny)
  test <- t(sce@assays@data$exprs)
  
  # Create a data frame for gating and downsampling if requested
  test_plot <- data.frame('x' = test[, x_channel], 'y' = test[, y_channel])
  
  if (!is.null(n_samples) && n_samples < nrow(test)) {
    set.seed(42)  # For reproducibility
    indices <- sample(nrow(test_plot), n_samples)
    test_plot <- test_plot[indices, ]
    
    # If plot_above_one is TRUE, select points outside the (0,0)-(1,1) square and add them to the test_plot
    if (plot_above_one) {
      above_one <- test[!(test[, x_channel] <= 0.5 & test[, y_channel] <= 0.5), c(x_channel, y_channel)]
      above_one <- data.frame('x' = above_one[, x_channel], 'y' = above_one[, y_channel])
      test_plot <- rbind(test_plot, above_one)
    }
  }
  
  
  # Create a shiny app
  ui <- fluidPage(
    plotOutput("plot", click = "plot_click"),
    actionButton("done", "Done")
  )
  
  server <- function(input, output, session) {
    clicks <- reactiveValues(data = data.frame(x = numeric(), y = numeric()))
    
    output$plot <- renderPlot({
      p <- ggplot(test_plot, aes(x, y)) +
        geom_point(alpha = 0.4) +
        expand_limits(x = c(min(test[,x_channel]), max(test[,x_channel]) * 1.05),
                      y = c(min(test[,y_channel]), max(test[,y_channel]) * 1.05)) +
        labs(x = x_channel, y = y_channel) +
        theme_minimal()
      if(nrow(clicks$data) > 0) {
        p <- p + geom_polygon(data = clicks$data, aes(x, y), fill = NA, color = "red", linetype = "dashed", linewidth = 1)
      }
      p
    })
    
    observeEvent(input$plot_click, {
      new_point <- data.frame(x = input$plot_click$x, y = input$plot_click$y)
      isolate({
        clicks$data <- rbind(clicks$data, new_point)
      })
    })
    
    observeEvent(input$done, {
      stopApp(clicks$data)
    })
  }
  
  # Run the shiny app and save the result
  vertices <- runApp(shinyApp(ui, server), port = 5823)
  
  # Add the first vertex to the end to close the polygon
  vertices <- rbind(vertices, vertices[1, ])
  
  # Add a new column 'gate' to the test_plot
  full_test_plot <- data.frame('x' = test[, x_channel], 'y' = test[, y_channel])
  full_test_plot$gate <- point.in.polygon(full_test_plot$x, full_test_plot$y, vertices$x, vertices$y) %in% c(1, 2)
  
  # Create a ggplot of the gate, but downsampled to n_samples cells (the same as for making the gate)
  p <- ggplot(full_test_plot %>% sample_n(min(nrow(full_test_plot), n_samples)), aes(x = x, y = y)) +
    geom_point(aes(color = gate), alpha = 0.4, size = 0.5) +
    geom_polygon(data = vertices, aes(x = x, y = y), fill = NA, color = "red", linetype = "dashed") +
    labs(x = x_channel, y = y_channel) +
    scale_color_manual(values = c("FALSE" = "grey80", "TRUE" = "green4")) +
    theme_minimal() +
    ggtitle(gate_name) +
    coord_equal()
  
  # Count the number of cells that passed the gate
  gate_passed_count <- sum(full_test_plot$gate)
  
  # Calculate the proportion of total cells that passed the gate
  gate_passed_prop <- gate_passed_count / nrow(full_test_plot)
  
  # Add the count and proportion to the plot as a text label
  p <- p + geom_text(aes(label = sprintf("%s: %d (%.2f%%) of cells passed", gate_name, gate_passed_count, gate_passed_prop * 100), x = Inf, y = Inf), hjust = "right", vjust = "top")
  
  # Update the SingleCellExperiment object with the gating information
  sce[[gate_name]] <- factor(full_test_plot$gate)
  
  # Return the updated SingleCellExperiment object and the plot
  list(sce = sce, plot = p)
}





## ----------------------------------------------------------------------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param object PARAM_DESCRIPTION
#' @param features PARAM_DESCRIPTION
#' @param dims PARAM_DESCRIPTION, Default: c(1, 2)
#' @param cells PARAM_DESCRIPTION, Default: NULL
#' @param cols PARAM_DESCRIPTION, Default: if (blend) {
#'    c("lightgrey", "#ff0000", "#00ff00")
#'} else {
#'    c("lightgrey", "blue")
#'}
#' @param pt.size PARAM_DESCRIPTION, Default: 1.5
#' @param alpha PARAM_DESCRIPTION, Default: 1
#' @param order PARAM_DESCRIPTION, Default: FALSE
#' @param min.cutoff PARAM_DESCRIPTION, Default: NA
#' @param textsize PARAM_DESCRIPTION, Default: 12
#' @param max.cutoff PARAM_DESCRIPTION, Default: NA
#' @param reduction PARAM_DESCRIPTION, Default: NULL
#' @param split.by PARAM_DESCRIPTION, Default: NULL
#' @param keep.scale PARAM_DESCRIPTION, Default: 'feature'
#' @param shape.by PARAM_DESCRIPTION, Default: NULL
#' @param slot PARAM_DESCRIPTION, Default: 'data'
#' @param blend PARAM_DESCRIPTION, Default: FALSE
#' @param blend.threshold PARAM_DESCRIPTION, Default: 0.5
#' @param myTheme PARAM_DESCRIPTION, Default: TRUE
#' @param label PARAM_DESCRIPTION, Default: FALSE
#' @param label.size PARAM_DESCRIPTION, Default: 4
#' @param label.color PARAM_DESCRIPTION, Default: 'black'
#' @param repel PARAM_DESCRIPTION, Default: FALSE
#' @param ncol PARAM_DESCRIPTION, Default: NULL
#' @param coord.fixed PARAM_DESCRIPTION, Default: FALSE
#' @param by.col PARAM_DESCRIPTION, Default: TRUE
#' @param sort.cell PARAM_DESCRIPTION, Default: deprecated()
#' @param interactive PARAM_DESCRIPTION, Default: FALSE
#' @param combine PARAM_DESCRIPTION, Default: TRUE
#' @param raster PARAM_DESCRIPTION, Default: NULL
#' @param raster.dpi PARAM_DESCRIPTION, Default: c(512, 512)
#' @param panel_spacing PARAM_DESCRIPTION, Default: 1
#' @param hide_axis PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @export 
#' @importFrom tidyr pivot_longer
#' @importFrom scattermore geom_scattermore
FeaturePlot3 <- function(object, features, dims = c(1, 2), cells = NULL, 
                           cols = if (blend) { c("lightgrey", "#ff0000", "#00ff00") } else { c("lightgrey", "blue") },
                           pt.size = 1.5, alpha = 1, order = FALSE, min.cutoff = NA, textsize = 12,
                           max.cutoff = NA, reduction = NULL, split.by = NULL, keep.scale = "feature", 
                           shape.by = NULL, slot = "data", blend = FALSE, blend.threshold = 0.5, myTheme = TRUE,
                           label = FALSE, label.size = 4, label.color = "black", repel = FALSE, 
                           ncol = NULL, coord.fixed = FALSE, by.col = TRUE, sort.cell = deprecated(), 
                           interactive = FALSE, combine = TRUE, raster = NULL, raster.dpi = c(512,512),
                           panel_spacing = 1, hide_axis = FALSE) {
    if (is_present(arg = sort.cell)) {
        deprecate_stop(when = "4.9.0", what = "FeaturePlot(sort.cell = )", 
            with = "FeaturePlot(order = )")
    }
    if (isTRUE(x = interactive)) {
        return(IFeaturePlot(object = object, feature = features[1], 
            dims = dims, reduction = reduction, slot = slot))
    }
    if (!is.null(x = keep.scale)) {
        keep.scale <- arg_match0(arg = keep.scale, values = c("feature", "all"))
    }
    reduction <- reduction %||% DefaultDimReduc(object = object)
    if (!is_integerish(x = dims, n = 2L, finite = TRUE) || !all(dims > 0L)) {
        abort(message = "'dims' must be a two-length integer vector")
    }
    if (isTRUE(x = blend)) {
        abort("Blending (blend = TRUE) is not supported when using facetting")
    }
    dims <- paste0(Key(object = object[[reduction]]), dims)
    cells <- cells %||% Cells(x = object[[reduction]])
    data <- FetchData(object = object, vars = c(dims, "ident", features), cells = cells, slot = slot)
    if (ncol(x = data) < 4) {
        abort(message = paste("None of the requested features were found:", 
            paste(features, collapse = ", "), "in slot ", slot))
    } else if (!all(dims %in% colnames(x = data))) {
        abort(message = "The dimensions requested were not found")
    }
    # Replicate min.cutoff and max.cutoff if a single value is provided
    if(length(min.cutoff) == 1) {
        min.cutoff <- rep(min.cutoff, length(features))
        names(min.cutoff) <- features
    }
    if(length(max.cutoff) == 1) {
        max.cutoff <- rep(max.cutoff, length(features))
        names(max.cutoff) <- features
    }
    # Adjust min and max cutoffs for each feature
    for (f in features) {
        current_min <- min.cutoff[[f]]
        current_max <- max.cutoff[[f]]
        min_val <- if (is.na(current_min)) min(data[[f]]) else current_min
        max_val <- if (is.na(current_max)) max(data[[f]]) else current_max
        min.use <- SetQuantile(cutoff = min_val, data = data[[f]])
        max.use <- SetQuantile(cutoff = max_val, data = data[[f]])
        data[[f]] <- pmin(pmax(data[[f]], min.use), max.use)
    }
    # Create a split factor for facetting. If split.by is NULL, create a dummy factor.
    data$split <- if (is.null(x = split.by)) {
        "all"
    } else {
        switch(EXPR = split.by, 
            ident = Idents(object = object)[cells, drop = TRUE], 
            object[[split.by, drop = TRUE]][cells, drop = TRUE])
    }
    if (!is.factor(x = data$split)) {
        data$split <- factor(x = data$split)
    }
    if (!is.null(x = shape.by)) {
        data[, shape.by] <- object[[shape.by, drop = TRUE]]
    }
    if (!requireNamespace("tidyr", quietly = TRUE)) {
        stop("Package 'tidyr' is required for pivoting the data")
    }
    data_long <- tidyr::pivot_longer(data, cols = all_of(features), names_to = "feature", 
                                     values_to = "expression")
    # Set the order of facets according to the input feature list
    data_long$feature <- factor(data_long$feature, levels = features)
    
    library(ggplot2)
    # Use scattermore for rasterized plotting if raster is TRUE
    if (!is.null(raster) && isTRUE(raster)) {
        p <- ggplot(data_long, aes_string(x = dims[1], y = dims[2], color = "expression")) +
             scattermore::geom_scattermore(pointsize = pt.size, alpha = alpha) 
    } else {
        p <- ggplot(data_long, aes_string(x = dims[1], y = dims[2], color = "expression")) +
             geom_point(size = pt.size, alpha = alpha)
    }
    
    p <- p + scale_color_gradientn(colors = cols) +
        theme_cowplot() + 
        CenterTitle() +
        labs(x = dims[1], y = dims[2])
    
    if (coord.fixed) {
        p <- p + coord_fixed()
    }
    
    if (is.null(split.by)) {
        p <- p + facet_wrap(~ feature, ncol = ncol %||% 2)
    } else {
        p <- p + facet_wrap(feature ~ split, ncol = ncol %||% length(levels(data$split)))
    }
    
    if (myTheme) {
        p <- p + theme(
            panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.8),
            strip.text = element_text(face = "bold", color = "black", size = textsize),
            axis.line = element_blank(),
            aspect.ratio = 1,
            text = element_text(size = textsize),
            strip.background = element_blank(), 
            panel.spacing = unit(panel_spacing, "lines")
        )
        if (isTRUE(hide_axis)) {
            p <- p + theme(
                axis.title.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_blank(),
                axis.title.y = element_blank(),
                axis.text.x = element_blank(),
                axis.text.y = element_blank()
            )
        }
    }
    return(p)
}

environment(FeaturePlot3) <- asNamespace('Seurat')


## ----------------------------------------------------------------------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param object PARAM_DESCRIPTION
#' @param features PARAM_DESCRIPTION
#' @param dims PARAM_DESCRIPTION, Default: c(1, 2)
#' @param cells PARAM_DESCRIPTION, Default: NULL
#' @param cols PARAM_DESCRIPTION, Default: if (blend) {
#'    c("lightgrey", "#ff0000", "#00ff00")
#'} else {
#'    c("lightgrey", "blue")
#'}
#' @param pt.size PARAM_DESCRIPTION, Default: NULL
#' @param alpha PARAM_DESCRIPTION, Default: 1
#' @param order PARAM_DESCRIPTION, Default: FALSE
#' @param min.cutoff PARAM_DESCRIPTION, Default: NA
#' @param textsize PARAM_DESCRIPTION, Default: 12
#' @param max.cutoff PARAM_DESCRIPTION, Default: NA
#' @param reduction PARAM_DESCRIPTION, Default: NULL
#' @param split.by PARAM_DESCRIPTION, Default: NULL
#' @param keep.scale PARAM_DESCRIPTION, Default: 'feature'
#' @param shape.by PARAM_DESCRIPTION, Default: NULL
#' @param slot PARAM_DESCRIPTION, Default: 'data'
#' @param blend PARAM_DESCRIPTION, Default: FALSE
#' @param blend.threshold PARAM_DESCRIPTION, Default: 0.5
#' @param myTheme PARAM_DESCRIPTION, Default: T
#' @param label PARAM_DESCRIPTION, Default: FALSE
#' @param label.size PARAM_DESCRIPTION, Default: 4
#' @param label.color PARAM_DESCRIPTION, Default: 'black'
#' @param repel PARAM_DESCRIPTION, Default: FALSE
#' @param ncol PARAM_DESCRIPTION, Default: NULL
#' @param coord.fixed PARAM_DESCRIPTION, Default: FALSE
#' @param by.col PARAM_DESCRIPTION, Default: TRUE
#' @param sort.cell PARAM_DESCRIPTION, Default: deprecated()
#' @param interactive PARAM_DESCRIPTION, Default: FALSE
#' @param combine PARAM_DESCRIPTION, Default: TRUE
#' @param raster PARAM_DESCRIPTION, Default: NULL
#' @param raster.dpi PARAM_DESCRIPTION, Default: c(512, 512)
#' @return OUTPUT_DESCRIPTION
#' @export 
FeaturePlot2 <- function (object, features, dims = c(1, 2), cells = NULL, cols = if (blend) {
    c("lightgrey", "#ff0000", "#00ff00")
} else {
    c("lightgrey", "blue")
}, pt.size = NULL, alpha = 1, order = FALSE, min.cutoff = NA, textsize = 12,
    max.cutoff = NA, reduction = NULL, split.by = NULL, keep.scale = "feature", 
    shape.by = NULL, slot = "data", blend = FALSE, blend.threshold = 0.5, myTheme = T,
    label = FALSE, label.size = 4, label.color = "black", repel = FALSE, 
    ncol = NULL, coord.fixed = FALSE, by.col = TRUE, sort.cell = deprecated(), 
    interactive = FALSE, combine = TRUE, raster = NULL, raster.dpi = c(512, 
        512)) 
{
    if (is_present(arg = sort.cell)) {
        deprecate_stop(when = "4.9.0", what = "FeaturePlot(sort.cell = )", 
            with = "FeaturePlot(order = )")
    }
    if (isTRUE(x = interactive)) {
        return(IFeaturePlot(object = object, feature = features[1], 
            dims = dims, reduction = reduction, slot = slot))
    }
    if (!is.null(x = keep.scale)) {
        keep.scale <- arg_match0(arg = keep.scale, values = c("feature", 
            "all"))
    }
    no.right <- theme(axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), 
        axis.text.y.right = element_blank(), axis.title.y.right = element_text(face = "bold", 
            size = 14, margin = margin(r = 7)))
    reduction <- reduction %||% DefaultDimReduc(object = object)
    if (!is_integerish(x = dims, n = 2L, finite = TRUE) && !all(dims > 
        0L)) {
        abort(message = "'dims' must be a two-length integer vector")
    }
    if (isTRUE(x = blend) && length(x = features) != 2) {
        abort(message = "Blending feature plots only works with two features")
    }
    if (isTRUE(x = blend)) {
        default.colors <- eval(expr = formals(fun = FeaturePlot)$cols)
        cols <- switch(EXPR = as.character(x = length(x = cols)), 
            `0` = {
                warn(message = "No colors provided, using default colors")
                default.colors
            }, `1` = {
                warn(message = paste("Only one color provided, assuming", 
                  sQuote(x = cols), "is double-negative and augmenting with default colors"))
                c(cols, default.colors[2:3])
            }, `2` = {
                warn(message = paste("Only two colors provided, assuming specified are for features and agumenting with", 
                  sQuote(default.colors[1]), "for double-negatives", 
                  ))
                c(default.colors[1], cols)
            }, `3` = cols, {
                warn(message = "More than three colors provided, using only first three")
                cols[1:3]
            })
    }
    if (isTRUE(x = blend) && length(x = cols) != 3) {
        abort("Blending feature plots only works with three colors; first one for negative cells")
    }
    dims <- paste0(Key(object = object[[reduction]]), dims)
    cells <- cells %||% Cells(x = object[[reduction]])
    data <- FetchData(object = object, vars = c(dims, "ident", 
        features), cells = cells, slot = slot)
    if (ncol(x = data) < 4) {
        abort(message = paste("None of the requested features were found:", 
            paste(features, collapse = ", "), "in slot ", slot))
    }
    else if (!all(dims %in% colnames(x = data))) {
        abort(message = "The dimensions requested were not found")
    }
    features <- setdiff(x = names(x = data), y = c(dims, "ident"))
    min.cutoff <- mapply(FUN = function(cutoff, feature) {
        return(ifelse(test = is.na(x = cutoff), yes = min(data[, 
            feature]), no = cutoff))
    }, cutoff = min.cutoff, feature = features)
    max.cutoff <- mapply(FUN = function(cutoff, feature) {
        return(ifelse(test = is.na(x = cutoff), yes = max(data[, 
            feature]), no = cutoff))
    }, cutoff = max.cutoff, feature = features)
    check.lengths <- unique(x = vapply(X = list(features, min.cutoff, 
        max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
    if (length(x = check.lengths) != 1) {
        abort(message = "There must be the same number of minimum and maximum cuttoffs as there are features")
    }
    names(x = min.cutoff) <- names(x = max.cutoff) <- features
    brewer.gran <- ifelse(test = length(x = cols) == 1, yes = brewer.pal.info[cols, 
        ]$maxcolors, no = length(x = cols))
    for (i in seq_along(along.with = features)) {
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param p PARAM_DESCRIPTION
#' @param w PARAM_DESCRIPTION, Default: 20
#' @param h PARAM_DESCRIPTION, Default: 20
#' @param title PARAM_DESCRIPTION, Default: 'your_title'
#' @return OUTPUT_DESCRIPTION
#' @export 
#' @importFrom tools file_path_sans_ext
        f <- features[i]
        data.feature <- data[[f]]
        min.use <- SetQuantile(cutoff = min.cutoff[f], data = data.feature)
        max.use <- SetQuantile(cutoff = max.cutoff[f], data = data.feature)
        data.feature[data.feature < min.use] <- min.use
        data.feature[data.feature > max.use] <- max.use
        if (brewer.gran != 2) {
            data.feature <- if (all(data.feature == 0)) {
                rep_len(x = 0, length.out = length(x = data.feature))
            }
            else {
                as.numeric(x = as.factor(x = cut(x = as.numeric(x = data.feature), 
                  breaks = 2)))
            }
        }
        data[[f]] <- data.feature
    }
    data$split <- if (is.null(x = split.by)) {
        RandomName()
    }
    else {
        switch(EXPR = split.by, ident = Idents(object = object)[cells, 
            drop = TRUE], object[[split.by, drop = TRUE]][cells, 
            drop = TRUE])
    }
    if (!is.factor(x = data$split)) {
        data$split <- factor(x = data$split)
    }
    if (!is.null(x = shape.by)) {
        data[, shape.by] <- object[[shape.by, drop = TRUE]]
    }
    plots <- vector(mode = "list", length = ifelse(test = blend, 
        yes = 4, no = length(x = features) * length(x = levels(x = data$split))))
    xlims <- c(floor(x = min(data[, dims[1]])), ceiling(x = max(data[, 
        dims[1]])))
    ylims <- c(floor(min(data[, dims[2]])), ceiling(x = max(data[, 
        dims[2]])))
    if (blend) {
        ncol <- 4
        color.matrix <- BlendMatrix(two.colors = cols[2:3], col.threshold = blend.threshold, 
            negative.color = cols[1])
        cols <- cols[2:3]
        colors <- list(color.matrix[, 1], color.matrix[1, ], 
            as.vector(x = color.matrix))
    }
    for (i in 1:length(x = levels(x = data$split))) {
        ident <- levels(x = data$split)[i]
        data.plot <- data[as.character(x = data$split) == ident, 
            , drop = FALSE]
        if (isTRUE(x = blend)) {
            features <- features[1:2]
            no.expression <- features[colMeans(x = data.plot[, 
                features]) == 0]
            if (length(x = no.expression) != 0) {
                abort(message = paste("The following features have no value:", 
                  paste(no.expression, collapse = ", ")))
            }
            data.plot <- cbind(data.plot[, c(dims, "ident")], 
                BlendExpression(data = data.plot[, features[1:2]]))
            features <- colnames(x = data.plot)[4:ncol(x = data.plot)]
        }
        for (j in 1:length(x = features)) {
            feature <- features[j]
            if (isTRUE(x = blend)) {
                cols.use <- as.numeric(x = as.character(x = data.plot[, 
                  feature])) + 1
                cols.use <- colors[[j]][sort(x = unique(x = cols.use))]
            }
            else {
                cols.use <- NULL
            }
            data.single <- data.plot[, c(dims, "ident", feature, 
                shape.by)]
            plot <- SingleDimPlot(data = data.single, dims = dims, 
                col.by = feature, order = order, pt.size = pt.size, 
                alpha = alpha, cols = cols.use, shape.by = shape.by, 
                label = FALSE, raster = raster, raster.dpi = raster.dpi) + 
                scale_x_continuous(limits = xlims) + scale_y_continuous(limits = ylims) + 
                theme_cowplot() + CenterTitle()
            if (isTRUE(x = label)) {
                plot <- LabelClusters(plot = plot, id = "ident", 
                  repel = repel, size = label.size, color = label.color)
            }
            if (length(x = levels(x = data$split)) > 1) {
                plot <- plot + theme(panel.border = element_rect(fill = NA, 
                  colour = "black"), linewidth=1)
                plot <- plot + if (i == 1) {
                  labs(title = feature)
                }
                else {
                  labs(title = NULL)
                }
                if (j == length(x = features) && !blend) {
                  suppressMessages(expr = plot <- plot + scale_y_continuous(sec.axis = dup_axis(name = ident), 
                    limits = ylims) + no.right)
                }
                if (j != 1) {
                  plot <- plot + theme(axis.line.y = element_blank(), 
                    axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
                    axis.title.y.left = element_blank())
                }
                if (i != length(x = levels(x = data$split))) {
                  plot <- plot + theme(axis.line.x = element_blank(), 
                    axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
                    axis.title.x = element_blank())
                }
            }
            else {
              
              if (myTheme == TRUE){
              # David update to add border around UMAP plots as well as remove axes labels
                        # Put a border around the plot!
        plot <- plot + theme(panel.border = element_rect(fill = NA, 
                  colour = "black", linewidth=0.8), axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), aspect.ratio = 1,
        axis.title.y = element_blank(), axis.line = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), text = element_text(size = textsize))}
                
                
            }
            if (!blend) {
                plot <- plot + guides(color = NULL)
                cols.grad <- cols
                if (length(x = cols) == 1) {
                  plot <- plot + scale_color_brewer(palette = cols)
                }
                else if (length(x = cols) > 1) {
                  unique.feature.exp <- unique(data.plot[, feature])
                  if (length(unique.feature.exp) == 1) {
                    warn(message = paste0("All cells have the same value (", 
                      unique.feature.exp, ") of ", dQuote(x = feature)))
                    if (unique.feature.exp == 0) {
                      cols.grad <- cols[1]
                    }
                    else {
                      cols.grad <- cols
                    }
                  }
                  plot <- suppressMessages(expr = plot + scale_color_gradientn(colors = cols.grad, 
                    guide = "colorbar"))
                }
            }
            if (!(is.null(x = keep.scale)) && keep.scale == "feature" && 
                !blend) {
                max.feature.value <- max(data[, feature])
                min.feature.value <- min(data[, feature])
                plot <- suppressMessages(plot & scale_color_gradientn(colors = cols, 
                  limits = c(min.feature.value, max.feature.value)))
            }
            if (coord.fixed) {
                plot <- plot + coord_fixed()
            }
            plot <- plot
            plots[[(length(x = features) * (i - 1)) + j]] <- plot
        }
    }
    if (isTRUE(x = blend)) {
        blend.legend <- BlendMap(color.matrix = color.matrix)
        for (ii in 1:length(x = levels(x = data$split))) {
            suppressMessages(expr = plots <- append(x = plots, 
                values = list(blend.legend + scale_y_continuous(sec.axis = dup_axis(name = ifelse(test = length(x = levels(x = data$split)) > 
                  1, yes = levels(x = data$split)[ii], no = "")), 
                  expand = c(0, 0)) + labs(x = features[1], y = features[2], 
                  title = if (ii == 1) {
                    paste("Color threshold:", blend.threshold)
                  } else {
                    NULL
                  }) + no.right), after = 4 * ii - 1))
        }
    }
    plots <- Filter(f = Negate(f = is.null), x = plots)
    if (is.null(x = ncol)) {
        ncol <- 2
        if (length(x = features) == 1) {
            ncol <- 1
        }
        if (length(x = features) > 6) {
            ncol <- 3
        }
        if (length(x = features) > 9) {
            ncol <- 4
        }
    }
    ncol <- ifelse(test = is.null(x = split.by) || isTRUE(x = blend), 
        yes = ncol, no = length(x = features))
    legend <- if (isTRUE(x = blend)) {
        "none"
    }
    else {
        split.by %iff% "none"
    }
    if (isTRUE(x = combine)) {
        if (by.col && !is.null(x = split.by) && !blend) {
            plots <- lapply(X = plots, FUN = function(x) {
                return(suppressMessages(expr = x + theme_cowplot() + 
                  ggtitle("") + scale_y_continuous(sec.axis = dup_axis(name = ""), 
                  limits = ylims) + no.right))
            })
            nsplits <- length(x = levels(x = data$split))
            idx <- 1
            for (i in (length(x = features) * (nsplits - 1) + 
                1):(length(x = features) * nsplits)) {
                plots[[i]] <- suppressMessages(expr = plots[[i]] + 
                  scale_y_continuous(sec.axis = dup_axis(name = features[[idx]]), 
                    limits = ylims) + no.right)
                idx <- idx + 1
            }
            idx <- 1
            for (i in which(x = 1:length(x = plots)%%length(x = features) == 
                1)) {
                plots[[i]] <- plots[[i]] + ggtitle(levels(x = data$split)[[idx]]) + 
                  theme(plot.title = element_text(hjust = 0.5))
                idx <- idx + 1
            }
            idx <- 1
            if (length(x = features) == 1) {
                for (i in 1:length(x = plots)) {
                  plots[[i]] <- plots[[i]] + ggtitle(levels(x = data$split)[[idx]]) + 
                    theme(plot.title = element_text(hjust = 0.5))
                  idx <- idx + 1
                }
                ncol <- 1
                nrow <- nsplits
            }
            else {
                nrow <- split.by %iff% length(x = levels(x = data$split))
            }
            plots <- plots[c(do.call(what = rbind, args = split(x = 1:length(x = plots), 
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param p PARAM_DESCRIPTION
#' @param w PARAM_DESCRIPTION, Default: 20
#' @param h PARAM_DESCRIPTION, Default: 20
#' @param title PARAM_DESCRIPTION, Default: 'your_title'
#' @return OUTPUT_DESCRIPTION
#' @export 
#' @importFrom tools file_path_sans_ext
                f = ceiling(x = seq_along(along.with = 1:length(x = plots))/length(x = features)))))]
            plots <- wrap_plots(plots, ncol = nrow, nrow = ncol)
            if (!is.null(x = legend) && legend == "none") {
                plots <- plots & NoLegend()
            }
        }
        else {
            plots <- wrap_plots(plots, ncol = ncol, nrow = split.by %iff% 
                length(x = levels(x = data$split)))
        }
        if (!is.null(x = legend) && legend == "none") {
            plots <- plots & NoLegend()
        }
        if (!(is.null(x = keep.scale)) && keep.scale == "all" && 
            !blend) {
            max.feature.value <- max(data[, features])
            min.feature.value <- min(data[, features])
            plots <- suppressMessages(plots & scale_color_gradientn(colors = cols, 
                limits = c(min.feature.value, max.feature.value)))
        }
    }
    return(plots)
}

environment(FeaturePlot2) <- asNamespace('Seurat')



## ----------------------------------------------------------------------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param object PARAM_DESCRIPTION
#' @param dims PARAM_DESCRIPTION, Default: c(1, 2)
#' @param cells PARAM_DESCRIPTION, Default: NULL
#' @param cols PARAM_DESCRIPTION, Default: NULL
#' @param pt.size PARAM_DESCRIPTION, Default: NULL
#' @param reduction PARAM_DESCRIPTION, Default: NULL
#' @param group.by PARAM_DESCRIPTION, Default: NULL
#' @param split.by PARAM_DESCRIPTION, Default: NULL
#' @param legend_title PARAM_DESCRIPTION, Default: group.by
#' @param legendpointsize PARAM_DESCRIPTION, Default: 6
#' @param shape.by PARAM_DESCRIPTION, Default: NULL
#' @param order PARAM_DESCRIPTION, Default: NULL
#' @param shuffle PARAM_DESCRIPTION, Default: FALSE
#' @param seed PARAM_DESCRIPTION, Default: 1
#' @param textsize PARAM_DESCRIPTION, Default: 12
#' @param label PARAM_DESCRIPTION, Default: FALSE
#' @param label.size PARAM_DESCRIPTION, Default: 4
#' @param label.color PARAM_DESCRIPTION, Default: 'black'
#' @param label.box PARAM_DESCRIPTION, Default: FALSE
#' @param repel PARAM_DESCRIPTION, Default: FALSE
#' @param alpha PARAM_DESCRIPTION, Default: 1
#' @param cells.highlight PARAM_DESCRIPTION, Default: NULL
#' @param cols.highlight PARAM_DESCRIPTION, Default: '#DE2D26'
#' @param sizes.highlight PARAM_DESCRIPTION, Default: 1
#' @param na.value PARAM_DESCRIPTION, Default: 'grey50'
#' @param ncol PARAM_DESCRIPTION, Default: NULL
#' @param combine PARAM_DESCRIPTION, Default: TRUE
#' @param raster PARAM_DESCRIPTION, Default: NULL
#' @param raster.dpi PARAM_DESCRIPTION, Default: c(512, 512)
#' @return OUTPUT_DESCRIPTION
#' @export 
DimPlot2 <- function (object, dims = c(1, 2), cells = NULL, cols = NULL, 
    pt.size = NULL, reduction = NULL, group.by = NULL, split.by = NULL, legend_title = group.by, legendpointsize = 6,
    shape.by = NULL, order = NULL, shuffle = FALSE, seed = 1, textsize = 12,
    label = FALSE, label.size = 4, label.color = "black", label.box = FALSE, 
    repel = FALSE, alpha = 1, cells.highlight = NULL, cols.highlight = "#DE2D26", 
    sizes.highlight = 1, na.value = "grey50", ncol = NULL, combine = TRUE, 
    raster = NULL, raster.dpi = c(512, 512)) 
{
    if (!is_integerish(x = dims, n = 2L, finite = TRUE) || !all(dims > 
        0L)) {
        abort(message = "'dims' must be a two-length integer vector")
    }
    reduction <- reduction %||% DefaultDimReduc(object = object)
    cells <- cells %||% Cells(x = object, assay = DefaultAssay(object = object[[reduction]]))
    dims <- paste0(Key(object = object[[reduction]]), dims)
    orig.groups <- group.by
    group.by <- group.by %||% "ident"
    data <- FetchData(object = object, vars = c(dims, group.by), 
        cells = cells, clean = "project")
    group.by <- colnames(x = data)[3:ncol(x = data)]
    for (group in group.by) {
        if (!is.factor(x = data[, group])) {
            data[, group] <- factor(x = data[, group])
        }
    }
    if (!is.null(x = shape.by)) {
        data[, shape.by] <- object[[shape.by, drop = TRUE]]
    }
    if (!is.null(x = split.by)) {
        split <- FetchData(object = object, vars = split.by, 
            clean = TRUE)[split.by]
        data <- data[rownames(split), ]
        data[, split.by] <- split
    }
    if (isTRUE(x = shuffle)) {
        set.seed(seed = seed)
        data <- data[sample(x = 1:nrow(x = data)), ]
    }
    plots <- lapply(X = group.by, FUN = function(x) {
        plot <- SingleDimPlot(data = data[, c(dims, x, split.by, 
            shape.by)], dims = dims, col.by = x, cols = cols, 
            pt.size = pt.size, shape.by = shape.by, order = order, 
            alpha = alpha, label = FALSE, cells.highlight = cells.highlight, 
            cols.highlight = cols.highlight, sizes.highlight = sizes.highlight, 
            na.value = na.value, raster = raster, raster.dpi = raster.dpi)
        if (label) {
            plot <- LabelClusters(plot = plot, id = x, repel = repel, 
                size = label.size, split.by = split.by, box = label.box, 
                color = label.color)
        }
        if (!is.null(x = split.by)) {
            plot <- plot + FacetTheme() + facet_wrap(facets = vars(!!sym(x = split.by)), 
                ncol = if (length(x = group.by) > 1 || is.null(x = ncol)) {
                  length(x = unique(x = data[, split.by]))
                }
                else {
                  ncol
                })
        }
        plot <- if (is.null(x = orig.groups)) {
            plot + labs(title = NULL)
        }
        else {
            plot + CenterTitle()
        }
        

        # Put a border around the plot!
        plot <- plot + theme(panel.border = element_rect(fill = NA, colour = "black", linewidth=0.8), 
                  axis.title.x = element_blank(), 
                  axis.ticks.x = element_blank(), 
                  axis.ticks.y = element_blank(), 
                  aspect.ratio = 1,
                  axis.title.y = element_blank(), 
                  axis.line = element_blank(), 
                  axis.text.x = element_blank(), 
                  axis.text.y = element_blank(), 
                  text = element_text(size = textsize)) + 
          guides(col = guide_legend(title = legend_title, override.aes = list(alpha = 1, size = legendpointsize)))
        
    })
    if (!is.null(x = split.by)) {
        ncol <- 1
    }
    if (combine) {
        plots <- wrap_plots(plots, ncol = orig.groups %iff% ncol)
    }
    return(plots)
}

environment(DimPlot2) <- asNamespace('Seurat')


## ----------------------------------------------------------------------------------------------------------
# Not in function
'%!in%' <- function(x,y)!('%in%'(x,y))


## ----------------------------------------------------------------------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param p PARAM_DESCRIPTION
#' @param w PARAM_DESCRIPTION, Default: 20
#' @param h PARAM_DESCRIPTION, Default: 20
#' @param title PARAM_DESCRIPTION, Default: 'your_title'
#' @return OUTPUT_DESCRIPTION
#' @export 
#' @importFrom tools file_path_sans_ext
f <- function(p, w = 20, h = 20, title = "your_title") {
  # Get the name of the current RMD file
  library(rmdhelp)
  rmd_file_path <- get_this_rmd_file()
  rmd_file_name <- basename(rmd_file_path)
  rmd_file_name_noext <- tools::file_path_sans_ext(rmd_file_name)

  # Create folder name with today's date + name of the RMD file without extension
  folder_name <- paste0(Sys.Date(), "_", rmd_file_name_noext)

  # Check if the directory exists, and create it if it doesn't
  if (!dir.exists(folder_name)) {
    dir.create(folder_name)
  }
  
  # Create the timestamp for the file names
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H.%M.%S")
  
  # Create the filename for the PDF figure using file.path for compatibility
  pdf_filename <- file.path(folder_name, paste0(timestamp, "_", title, ".pdf"))
  
  # Create the PDF file
  pdf(pdf_filename, width = w, height = h)
  invisible(print(p))
  dev.off()

  # Check if an RMD or sessionInfo file has been saved within the last 5 minutes
  recent_files <- list.files(folder_name, pattern = "\\.rmd$|_sessioninfo\\.txt$", full.names = TRUE)
  latest_file_time <- max(file.info(recent_files)$mtime, na.rm = TRUE)
  current_time <- Sys.time()

  # If no RMD or sessionInfo file saved in the last 5 minutes, save a copy
  if (is.na(latest_file_time) || difftime(current_time, latest_file_time, units = "mins") > 20) {
    rmd_filename <- file.path(folder_name, paste0(timestamp, "_", rmd_file_name))
    file.copy(from = rmd_file_path, to = rmd_filename)

    # Save the output of sessionInfo() as a text file
    sessioninfo_filename <- file.path(folder_name, paste0(timestamp, "_sessioninfo.txt"))
    writeLines(capture.output(sessionInfo()), con = sessioninfo_filename)

    cat("Figure, RMD file, and session information saved successfully in folder", folder_name, "\n")
  } else {
    cat("Figure saved successfully in folder", folder_name, ". RMD file and session information were not saved as a recent copy already exists.\n")
  }
}



## ----------------------------------------------------------------------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param p PARAM_DESCRIPTION
#' @param w PARAM_DESCRIPTION, Default: 20
#' @param h PARAM_DESCRIPTION, Default: 20
#' @param title PARAM_DESCRIPTION, Default: 'your_title'
#' @param saveRData PARAM_DESCRIPTION, Default: FALSE
#' @param saveExcel PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @export 
#' @importFrom tools file_path_sans_ext
f2 <- function(p, w = 20, h = 20, title = "your_title", saveRData = FALSE, saveExcel = FALSE) {
  # Get the name of the current RMD file
  library(rmdhelp)
  rmd_file_path <- get_this_rmd_file()
  rmd_file_name <- basename(rmd_file_path)
  rmd_file_name_noext <- tools::file_path_sans_ext(rmd_file_name)

  # Create folder name with today's date + name of the RMD file without extension
  folder_name <- paste0(Sys.Date(), "_", rmd_file_name_noext)

  # Check if the directory exists, and create it if it doesn't
  if (!dir.exists(folder_name)) {
    dir.create(folder_name)
  }
  
  # Create the timestamp for the file names
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H.%M.%S")
  
  # Create the filename for the PDF figure using file.path for compatibility
  pdf_filename <- file.path(folder_name, paste0(timestamp, "_", title, ".pdf"))
  
  # Create the PDF file
  pdf(pdf_filename, width = w, height = h)
  invisible(print(p))
  dev.off()

 # Save RData
  if (saveRData) {
    rdata_filename <- file.path(folder_name, paste0(timestamp, "_", title, ".RData"))
    cat("Saving RData...\n")
    save(p, file = rdata_filename)
  }

  # Save Excel
  if (saveExcel && !is.null(p$data)) {
  library(openxlsx)
    excel_filename <- paste0(folder_name, "/", timestamp, "_", title, "_data.xlsx")
    cat("Saving Excel...\n")
    wb <- createWorkbook()
    modifyBaseFont(wb, fontSize = 12, fontColour = "black", fontName = "Arial")
    addWorksheet(wb, "Sheet 1")
    writeData(wb, "Sheet 1", p$data)
    setColWidths(wb, "Sheet 1", cols = 1:ncol(p$data), widths = "auto")
    freezePane(wb, "Sheet 1", firstRow = TRUE)
    addFilter(wb, "Sheet 1", rows = 1, cols = 1:ncol(p$data))
    saveWorkbook(wb, excel_filename, overwrite = TRUE)
  }

  # Check if an RMD or sessionInfo file has been saved within the last 5 minutes
  recent_files <- list.files(folder_name, pattern = "\\.rmd$|_sessioninfo\\.txt$", full.names = TRUE)
  latest_file_time <- max(file.info(recent_files)$mtime, na.rm = TRUE)
  current_time <- Sys.time()

  # If no RMD or sessionInfo file saved in the last 5 minutes, save a copy
  if (is.na(latest_file_time) || difftime(current_time, latest_file_time, units = "mins") > 20) {
    rmd_filename <- file.path(folder_name, paste0(timestamp, "_", rmd_file_name))
    file.copy(from = rmd_file_path, to = rmd_filename)

    # Save the output of sessionInfo() as a text file
    sessioninfo_filename <- file.path(folder_name, paste0(timestamp, "_sessioninfo.txt"))
    writeLines(capture.output(sessionInfo()), con = sessioninfo_filename)

    cat("Figure, RMD file, and session information saved successfully in folder", folder_name, "\n")
  } else {
    cat("Figure saved successfully in folder", folder_name, ". RMD file and session information were not saved as a recent copy already exists.\n")
  }
}



## ----------------------------------------------------------------------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param p_list PARAM_DESCRIPTION
#' @param w PARAM_DESCRIPTION, Default: 20
#' @param h PARAM_DESCRIPTION, Default: 20
#' @param title PARAM_DESCRIPTION, Default: 'your_title'
#' @return OUTPUT_DESCRIPTION
#' @export 
#' @importFrom tools file_path_sans_ext
fl <- function(p_list, w = 20, h = 20, title = "your_title") {
  # Get the name of the current RMD file
  library(rmdhelp)
  rmd_file_path <- get_this_rmd_file()
  rmd_file_name <- basename(rmd_file_path)
  rmd_file_name_noext <- tools::file_path_sans_ext(rmd_file_name)

  # Create folder name with today's date + name of the RMD file without extension
  folder_name <- paste0(Sys.Date(), "_", rmd_file_name_noext)

  # Check if the directory exists, and create it if it doesn't
  if (!dir.exists(folder_name)) {
    dir.create(folder_name)
  }
  
  # Create the timestamp for the file names
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H.%M.%S")
  
  # Create the filename for the PDF figure
  pdf_filename <- paste0(folder_name, "/", timestamp, "_", title, ".pdf")
  
  # Create the PDF file
  pdf(pdf_filename, width = w, height = h)
  
  # Loop through each figure in p_list and print it to a new page in the PDF
  for(p in p_list) {
    invisible(print(p))
  }
  
  dev.off()

  # Check if an RMD or sessionInfo file has been saved within the last 5 minutes
  recent_files <- list.files(folder_name, pattern = "\\.rmd$|_sessioninfo\\.txt$", full.names = TRUE)
  latest_file_time <- max(file.info(recent_files)$mtime, na.rm = TRUE)
  current_time <- Sys.time()

  # If no RMD or sessionInfo file saved in the last 5 minutes, save a copy
  if (is.na(latest_file_time) || difftime(current_time, latest_file_time, units = "mins") > 20) {
    rmd_filename <- paste0(folder_name, "/", timestamp, "_", rmd_file_name)
    file.copy(from = rmd_file_path, to = rmd_filename)

    # Save the output of sessionInfo() as a text file
    sessioninfo_filename <- paste0(folder_name, "/", timestamp, "_sessioninfo.txt")
    writeLines(capture.output(sessionInfo()), con = sessioninfo_filename)

    cat("Figures, RMD file, and session information saved successfully in folder", folder_name, "\n")
  } else {
    cat("Figures saved successfully in folder", folder_name, ". RMD file and session information were not saved as a recent copy already exists.\n")
  }
}



## ----------------------------------------------------------------------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param p_list PARAM_DESCRIPTION
#' @param w PARAM_DESCRIPTION, Default: 20
#' @param h PARAM_DESCRIPTION, Default: 20
#' @param title PARAM_DESCRIPTION, Default: 'your_title'
#' @param saveRData PARAM_DESCRIPTION, Default: FALSE
#' @param saveExcel PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @export 
#' @importFrom tools file_path_sans_ext
fl2 <- function(p_list, w = 20, h = 20, title = "your_title", saveRData = FALSE, saveExcel = FALSE) {
  # Get the name of the current RMD file
  library(rmdhelp)
  rmd_file_path <- get_this_rmd_file()
  rmd_file_name <- basename(rmd_file_path)
  rmd_file_name_noext <- tools::file_path_sans_ext(rmd_file_name)

  # Create folder name with today's date + name of the RMD file without extension
  folder_name <- paste0(Sys.Date(), "_", rmd_file_name_noext)

  # Check if the directory exists, and create it if it doesn't
  if (!dir.exists(folder_name)) {
    dir.create(folder_name)
  }
  
  # Create the timestamp for the file names
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H.%M.%S")
  
  # Create the filename for the PDF figure
  pdf_filename <- paste0(folder_name, "/", timestamp, "_", title, ".pdf")
  
  # Create the PDF file
  pdf(pdf_filename, width = w, height = h)
  
  # Loop through each figure in p_list and print it to a new page in the PDF
  for(p in p_list) {
    invisible(print(p))
  }
  
  dev.off()

  # Save RData
  if (saveRData) {
    rdata_filename <- file.path(folder_name, paste0(timestamp, "_", title, ".RData"))
    cat("Saving RData...\n")
    save(p_list, file = rdata_filename)
  }

  # Save Excel
  if (saveExcel) {
    library(openxlsx)
    excel_filename <- paste0(folder_name, "/", timestamp, "_", title, "_data.xlsx")
    cat("Saving Excel...\n")
    wb <- createWorkbook()
    modifyBaseFont(wb, fontSize = 12, fontColour = "black", fontName = "Arial")
    
    for (i in seq_along(p_list)) {
      p <- p_list[[i]]
      sheet_name <- names(p_list)[i]
      if (!is.null(p$data)) {
        addWorksheet(wb, sheet_name)
        writeData(wb, sheet_name, p$data)
        setColWidths(wb, sheet_name, cols = 1:ncol(p$data), widths = "auto")
        freezePane(wb, sheet_name, firstRow = TRUE)
        addFilter(wb, sheet_name, rows = 1, cols = 1:ncol(p$data))
      }
    }
    
    saveWorkbook(wb, excel_filename, overwrite = TRUE)
  }

  # Check if an RMD or sessionInfo file has been saved within the last 5 minutes
  recent_files <- list.files(folder_name, pattern = "\\.rmd$|_sessioninfo\\.txt$", full.names = TRUE)
  latest_file_time <- max(file.info(recent_files)$mtime, na.rm = TRUE)
  current_time <- Sys.time()

  # If no RMD or sessionInfo file saved in the last 5 minutes, save a copy
  if (is.na(latest_file_time) || difftime(current_time, latest_file_time, units = "mins") > 20) {
    rmd_filename <- paste0(folder_name, "/", timestamp, "_", rmd_file_name)
    file.copy(from = rmd_file_path, to = rmd_filename)

    # Save the output of sessionInfo() as a text file
    sessioninfo_filename <- paste0(folder_name, "/", timestamp, "_sessioninfo.txt")
    writeLines(capture.output(sessionInfo()), con = sessioninfo_filename)

    cat("Figures, RMD file, and session information saved successfully in folder", folder_name, "\n")
  } else {
    cat("Figures saved successfully in folder", folder_name, ". RMD file and session information were not saved as a recent copy already exists.\n")
  }
}


## ----------------------------------------------------------------------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param df PARAM_DESCRIPTION
#' @param title PARAM_DESCRIPTION, Default: 'your_data'
#' @return OUTPUT_DESCRIPTION
#' @export 
#' @importFrom tools file_path_sans_ext
fx <- function(df, title = "your_data") {
  
  library(rmdhelp)
  library(openxlsx)
  rmd_file_path <- get_this_rmd_file()
  rmd_file_name <- basename(rmd_file_path)
  rmd_file_name_noext <- tools::file_path_sans_ext(rmd_file_name)

  folder_name <- paste0(Sys.Date(), "_", rmd_file_name_noext)

  if (!dir.exists(folder_name)) {
    dir.create(folder_name)
  }

  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  excel_filename <- paste0(folder_name, "/", timestamp, "_", title, ".xlsx")
  

  wb <- createWorkbook()
  modifyBaseFont(wb, fontSize = 12, fontColour = "black", fontName = "Arial")
  addWorksheet(wb, "Sheet 1")
  writeData(wb, "Sheet 1", df)
  setColWidths(wb, "Sheet 1", cols = 1:ncol(df), widths = "auto")
  freezePane(wb, "Sheet 1", firstRow = TRUE)
  addFilter(wb, "Sheet 1", rows = 1, cols = 1:ncol(df))


  saveWorkbook(wb, excel_filename, overwrite = TRUE)

}


## ----------------------------------------------------------------------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param df_list PARAM_DESCRIPTION
#' @param title PARAM_DESCRIPTION, Default: 'your_data'
#' @param separate_sheets PARAM_DESCRIPTION, Default: TRUE
#' @return OUTPUT_DESCRIPTION
#' @export 
#' @importFrom tools file_path_sans_ext
fxl <- function(df_list, title = "your_data", separate_sheets = TRUE) {
  
  library(rmdhelp)
  library(openxlsx)
  rmd_file_path <- get_this_rmd_file()
  rmd_file_name <- basename(rmd_file_path)
  rmd_file_name_noext <- tools::file_path_sans_ext(rmd_file_name)

  folder_name <- paste0(Sys.Date(), "_", rmd_file_name_noext)

  if (!dir.exists(folder_name)) {
    dir.create(folder_name)
  }

  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  excel_filename <- paste0(folder_name, "/", timestamp, "_", title, ".xlsx")
  
  wb <- createWorkbook()
  modifyBaseFont(wb, fontSize = 12, fontColour = "black", fontName = "Arial")
  
  if (separate_sheets) {
    for (i in seq_along(df_list)) {
      df <- df_list[[i]]
      sheet_name <- paste0("Sheet ", i)
      addWorksheet(wb, sheet_name)
      writeData(wb, sheet_name, df)
      setColWidths(wb, sheet_name, cols = 1:ncol(df), widths = "auto")
      freezePane(wb, sheet_name, firstRow = TRUE)
      addFilter(wb, sheet_name, rows = 1, cols = 1:ncol(df))
    }
  } else {
    addWorksheet(wb, "Concatenated Data")
    start_row <- 1
    for (i in seq_along(df_list)) {
      df <- df_list[[i]]
      writeData(wb, "Concatenated Data", df, startRow = start_row)
      setColWidths(wb, "Concatenated Data", cols = 1:ncol(df), widths = "auto")
      freezePane(wb, "Concatenated Data", firstRow = TRUE)
      addFilter(wb, "Concatenated Data", rows = start_row, cols = 1:ncol(df))
      start_row <- start_row + nrow(df) + 2  # Move to the next row after the current data frame plus one empty row
    }
  }

  saveWorkbook(wb, excel_filename, overwrite = TRUE)
  cat("Data frames saved successfully in", excel_filename, "\n")
}


## ----------------------------------------------------------------------------------------------------------
# Load necessary library
library(SingleCellExperiment)
library(ggplot2)
library(ggbeeswarm)
library(scales)
library(dplyr)

# Define the function
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param k PARAM_DESCRIPTION, Default: 'meta20'
#' @param by PARAM_DESCRIPTION, Default: c("sample_id", "cluster_id")
#' @param group_by PARAM_DESCRIPTION, Default: 'condition'
#' @param shape_by PARAM_DESCRIPTION, Default: NULL
#' @param col_clust PARAM_DESCRIPTION, Default: TRUE
#' @param n_cols PARAM_DESCRIPTION, Default: 4
#' @param log PARAM_DESCRIPTION, Default: FALSE
#' @param miny PARAM_DESCRIPTION, Default: 0.01
#' @param maxy PARAM_DESCRIPTION, Default: NA
#' @param panel_spacing PARAM_DESCRIPTION, Default: 1
#' @param point_size PARAM_DESCRIPTION, Default: 2
#' @param facet_ratio PARAM_DESCRIPTION, Default: 1.5
#' @param label_parse PARAM_DESCRIPTION, Default: F
#' @param distance PARAM_DESCRIPTION, Default: c("euclidean", "maximum", "manhattan", "canberra", "binary", 
#'    "minkowski")
#' @param linkage PARAM_DESCRIPTION, Default: c("average", "ward.D", "single", "complete", "mcquitty", "median", 
#'    "centroid", "ward.D2")
#' @param k_pal PARAM_DESCRIPTION, Default: CATALYST:::.cluster_cols
#' @param clusters_order PARAM_DESCRIPTION, Default: NULL
#' @param merging_col PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @export 
#' @importFrom CATALYST .cluster_cols
#' @importFrom dplyr summarize
plotAbundances222 <- function (x, k = "meta20", by = c("sample_id", "cluster_id"), 
                              group_by = "condition", shape_by = NULL, col_clust = TRUE, n_cols = 4, 
                              log = FALSE, miny = 0.01, maxy = NA, panel_spacing = 1, point_size = 2, facet_ratio = 1.5, label_parse = F,
                              distance = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"), 
                              linkage = c("average", "ward.D", "single", "complete", "mcquitty", "median", "centroid", "ward.D2"), 
                              k_pal = CATALYST:::.cluster_cols, clusters_order = NULL, merging_col = NULL) {
  by <- match.arg(by)
  .check_sce(x, TRUE)
  
  # Use the merging column from colData if provided
  if (!is.null(merging_col)) {
    cluster_ids <- x[[k]]
  } else {
    k <- .check_k(x, k)
    cluster_ids <- cluster_ids(x, k)
  }
  
  .check_cd_factor(x, group_by)
  .check_cd_factor(x, shape_by)
  .check_pal(k_pal)
  linkage <- match.arg(linkage)
  distance <- match.arg(distance)
  stopifnot(is.logical(col_clust), length(col_clust) == 1)
  shapes <- .get_shapes(x, shape_by)
  if (is.null(shapes)) 
    shape_by <- NULL
  if (by == "sample_id") {
    nk <- nlevels(factor(cluster_ids))
    if (length(k_pal) < nk) 
      k_pal <- colorRampPalette(k_pal)(nk)
  }
  
  ns <- table(cluster_id = cluster_ids, sample_id = sample_ids(x))
  fq <- prop.table(ns, 2) * 100
  df <- as.data.frame(fq)
  m <- match(df$sample_id, x$sample_id)
  for (i in c(shape_by, group_by)) df[[i]] <- x[[i]][m]
  if (by == "sample_id" && col_clust && length(unique(df$sample_id)) > 1) {
    d <- dist(t(fq), distance)
    h <- hclust(d, linkage)
    o <- colnames(fq)[h$order]
    df$sample_id <- factor(df$sample_id, o)
  }
  
  # Apply clusters_order if provided
  if (!is.null(clusters_order)) {
    df$cluster_id <- factor(df$cluster_id, levels = clusters_order)
  }
  
  # Decide whether to plot on a log scale. Add 0.01 to all values
  if (log == TRUE) {
    df$Freq <- df$Freq + 0.01
  }
  
  dfout <<- df
  
  # Calculate the maximum y value for each cluster_id and round up to the nearest whole number
  maxy_values <- df %>%
    group_by(cluster_id) %>%
    dplyr::summarize(maxy = ceiling(max(Freq, na.rm = TRUE) * 1.1))  # Add 10% buffer
  
  p <- ggplot(df, aes_string(y = "Freq")) + 
    labs(x = NULL, y = "Proportion [%]") + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          strip.text = element_text(face = "bold"), 
          strip.background = element_rect(fill = NA, color = NA), 
          axis.text = element_text(color = "black"), 
          aspect.ratio = facet_ratio,
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
          panel.border = element_rect(color = "black", fill = NA, size = 0.3),
          panel.spacing = unit(panel_spacing, "lines"),
          legend.key.height = unit(1.5, "lines"))
  
  if (label_parse == TRUE) {
    p <- p + facet_wrap(~cluster_id, scales = "free", ncol = n_cols, labeller = label_parsed) + 
      geom_boxplot(aes_string(x = group_by, fill = group_by), color = "grey16", position = position_dodge(), size = 0.5, alpha = 0.8, outlier.color = NA, show.legend = TRUE)
  } else {
    p <- p + facet_wrap(~cluster_id, scales = "free", ncol = n_cols) + 
      geom_boxplot(aes_string(x = group_by, fill = group_by), color = "grey16", position = position_dodge(), size = 0.5, alpha = 0.8, outlier.color = NA, show.legend = TRUE)
  }
  
  if (!is.null(shape_by)) {
    p <- p + geom_quasirandom(aes_string(x = group_by, shape = shape_by), size = point_size, width = 0.2)
  } else {
    p <- p + geom_quasirandom(aes_string(x = group_by), fill = "grey84", size = point_size, width = 0.2, shape = 21, alpha = 0.8)
    #p <- p + geom_quasirandom(aes_string(x = group_by, fill = group_by), size = point_size, width = 0.2, shape = 21, alpha = 0.6)
  }
  
  if (log == TRUE) {
    p + scale_y_continuous(trans = 'log10', limits = c(miny, maxy), breaks = c(0.01, 0.1, 1, 10, 100), labels = c(0.01, 0.1, 1, 10, 100)) + 
      annotation_logticks(base = 10, sides = "l", outside = TRUE) + 
      coord_cartesian(clip = "off") + 
      scale_size_area(max_size = 15) + 
      theme(axis.text.y = element_text(margin = margin(r = 8)))
  } else {
    p + scale_y_continuous(limits = c(0, NA)) + 
      coord_cartesian(clip = "off") + 
      scale_size_area(max_size = 15) + 
      facet_wrap(~cluster_id, scales = "free_y", ncol = n_cols, labeller = label_both)
  }
  
  # Apply the calculated maxy values to each facet
  p <- p + geom_blank(data = maxy_values, aes(y = maxy))
  
  return(p)
}

environment(plotAbundances222) <- asNamespace('CATALYST')


## ----------------------------------------------------------------------------------------------------------
library(ggiraph)
library(ggplot2)
library(htmlwidgets)

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param k PARAM_DESCRIPTION, Default: 'meta20'
#' @param by PARAM_DESCRIPTION, Default: c("sample_id", "cluster_id")
#' @param group_by PARAM_DESCRIPTION, Default: 'condition'
#' @param text_size PARAM_DESCRIPTION, Default: 16
#' @param shape_by PARAM_DESCRIPTION, Default: NULL
#' @param col_clust PARAM_DESCRIPTION, Default: TRUE
#' @param n_cols PARAM_DESCRIPTION, Default: 4
#' @param log PARAM_DESCRIPTION, Default: FALSE
#' @param miny PARAM_DESCRIPTION, Default: 0.01
#' @param maxy PARAM_DESCRIPTION, Default: NA
#' @param point_size PARAM_DESCRIPTION, Default: 2
#' @param distance PARAM_DESCRIPTION, Default: c("euclidean", "maximum", "manhattan", "canberra", "binary", 
#'    "minkowski")
#' @param linkage PARAM_DESCRIPTION, Default: c("average", "ward.D", "single", "complete", "mcquitty", "median", 
#'    "centroid", "ward.D2")
#' @param k_pal PARAM_DESCRIPTION, Default: CATALYST:::.cluster_cols
#' @param clusters_order PARAM_DESCRIPTION, Default: NULL
#' @param width_svg PARAM_DESCRIPTION, Default: 10
#' @param height_svg PARAM_DESCRIPTION, Default: 8
#' @return OUTPUT_DESCRIPTION
#' @export 
#' @importFrom CATALYST .cluster_cols
plotAbundancesInteractive <- function (x, k = "meta20", by = c("sample_id", "cluster_id"), group_by = "condition", text_size = 16,shape_by = NULL, col_clust = TRUE, n_cols = 4, log = FALSE, miny = 0.01, maxy = NA, point_size = 2, distance = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"), linkage = c("average", "ward.D", "single", "complete", "mcquitty", "median", "centroid", "ward.D2"), k_pal = CATALYST:::.cluster_cols, clusters_order = NULL, width_svg = 10, height_svg = 8) { 
    by <- match.arg(by) 
    .check_sce(x, TRUE) 
    k <- .check_k(x, k) 
    .check_cd_factor(x, group_by) 
    .check_cd_factor(x, shape_by) 
    .check_pal(k_pal) 
    linkage <- match.arg(linkage) 
    distance <- match.arg(distance) 
    stopifnot(is.logical(col_clust), length(col_clust) == 1) 
    shapes <- .get_shapes(x, shape_by) 
    if (is.null(shapes)) shape_by <- NULL 
    if (by == "sample_id") { 
        nk <- nlevels(cluster_ids(x, k)) 
        if (length(k_pal) < nk) k_pal <- colorRampPalette(k_pal)(nk) 
    } 
    ns <- table(cluster_id = cluster_ids(x, k), sample_id = sample_ids(x)) 
    fq <- prop.table(ns, 2) * 100 
    df <- as.data.frame(fq) 
    m <- match(df$sample_id, x$sample_id) 
    for (i in c(shape_by, group_by)) df[[i]] <- x[[i]][m] 
    if (by == "sample_id" && col_clust && length(unique(df$sample_id)) > 1) { 
        d <- dist(t(fq), distance) 
        h <- hclust(d, linkage) 
        o <- colnames(fq)[h$order] 
        df$sample_id <- factor(df$sample_id, o) 
    }

    # Apply clusters_order if provided
    if (!is.null(clusters_order)) { 
        df$cluster_id <- factor(df$cluster_id, levels = clusters_order) 
    }

    # Decide whether to plot on a log scale. Add 0.02 to all values
    if (log == TRUE) { 
        df$Freq <- df$Freq + 0.02 
    }

    dfout <<- df

    p <- ggplot(df, aes_string(y = "Freq", tooltip = group_by)) + 
        labs(x = NULL, y = "Proportion [%]") + 
        theme_bw() + 
        theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"), strip.background = element_rect(fill = NA, color = NA), axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.key.height = unit(0.8, "lines"))

    p <- p + facet_wrap(~cluster_id, scales = "free_y", ncol = n_cols) + 
        geom_boxplot_interactive(aes_string(x = group_by, fill = group_by), color = "grey16", position = position_dodge(), alpha = 0.8, outlier.color = NA, show.legend = TRUE)

    if (!is.null(shape_by)) { 
        p <- p + geom_point_interactive(aes_string(x = group_by, shape = shape_by), size = point_size, position = position_jitter(width = 0.2)) 
    } else { 
        p <- p + geom_point_interactive(aes_string(x = group_by), fill = "grey84", size = point_size, position = position_jitter(width = 0.2), shape = 21) 
    }

    if (log == TRUE) { 
        p <- p + scale_y_continuous(trans = 'log10', limits = c(miny, maxy), breaks = c(0.01, 0.1, 1, 10, 100), labels = c(0.01, 0.1, 1, 10, 100)) + 
            annotation_logticks(base = 10, sides = "l", outside = TRUE) + 
            coord_cartesian(clip = "off") + 
            scale_size_area(max_size = 15) + 
            theme(axis.text.y = element_text(margin = margin(r = 8))) 
    } else { 
        p <- p + scale_y_continuous(limits = c(0, maxy)) + 
            coord_cartesian(clip = "off") + 
            scale_size_area(max_size = 15) 
    }
    
    p <- p + theme(text = element_text(size = text_size))


    # Create the interactive plot object with specified size
    interactive_plot <- girafe(ggobj = p, width_svg = width_svg, height_svg = height_svg)

    # Set the size of the HTML widget container
    interactive_plot <- interactive_plot %>% girafe_options(opts_sizing(rescale = TRUE))

    # Return the interactive plot object
    return(interactive_plot) 
}

environment(plotAbundancesInteractive) <- asNamespace('CATALYST')


## ----------------------------------------------------------------------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param k PARAM_DESCRIPTION, Default: 'meta20'
#' @param by PARAM_DESCRIPTION, Default: c("sample_id", "cluster_id")
#' @param sample_order PARAM_DESCRIPTION, Default: NULL
#' @param group_by PARAM_DESCRIPTION, Default: 'condition'
#' @param shape_by PARAM_DESCRIPTION, Default: NULL
#' @param col_clust PARAM_DESCRIPTION, Default: TRUE
#' @param n_cols PARAM_DESCRIPTION, Default: 4
#' @param rotang PARAM_DESCRIPTION, Default: 45
#' @param distance PARAM_DESCRIPTION, Default: c("euclidean", "maximum", "manhattan", "canberra", "binary", 
#'    "minkowski")
#' @param linkage PARAM_DESCRIPTION, Default: c("average", "ward.D", "single", "complete", "mcquitty", "median", 
#'    "centroid", "ward.D2")
#' @param k_pal PARAM_DESCRIPTION, Default: CATALYST:::.cluster_cols
#' @return OUTPUT_DESCRIPTION
#' @export 
#' @importFrom CATALYST .cluster_cols
plotAbundances23 <- function (x, k = "meta20", by = c("sample_id", "cluster_id"), sample_order = NULL,
    group_by = "condition", shape_by = NULL, col_clust = TRUE, n_cols = 4, rotang = 45,
    distance = c("euclidean", "maximum", "manhattan", "canberra", 
        "binary", "minkowski"), linkage = c("average", "ward.D", 
        "single", "complete", "mcquitty", "median", "centroid", 
        "ward.D2"), k_pal = CATALYST:::.cluster_cols) 
{
    by <- match.arg(by)
    .check_sce(x, TRUE)
    k <- .check_k(x, k)
    .check_cd_factor(x, group_by)
    .check_cd_factor(x, shape_by)
    .check_pal(k_pal)
    linkage <- match.arg(linkage)
    distance <- match.arg(distance)
    stopifnot(is.logical(col_clust), length(col_clust) == 1)
    shapes <- .get_shapes(x, shape_by)
    if (is.null(shapes)) 
        shape_by <- NULL
    if (by == "sample_id") {
        nk <- nlevels(cluster_ids(x, k))
        if (length(k_pal) < nk) 
            k_pal <- colorRampPalette(k_pal)(nk)
    }
    ns <- table(cluster_id = cluster_ids(x, k), sample_id = sample_ids(x))
    fq <- prop.table(ns, 2) * 100
    df <- as.data.frame(fq)
    m <- match(df$sample_id, x$sample_id)
    for (i in c(shape_by, group_by)) df[[i]] <- x[[i]][m]
    if (by == "sample_id" && col_clust && length(unique(df$sample_id)) > 
        1) {
        d <- dist(t(fq), distance)
        h <- hclust(d, linkage)
        o <- colnames(fq)[h$order]
        df$sample_id <- factor(df$sample_id, o)
    }
    
    # Force desired sample order
    if (!is.null(sample_order)) {
        df$sample_id <- factor(df$sample_id, ordered = T, levels = sample_order)
    }
    
    p <- ggplot(df, aes_string(y = "Freq")) + labs(x = NULL, 
        y = "Proportion [%]") + theme_bw() + theme(panel.grid = element_blank(), 
        strip.text = element_text(face = "bold"), strip.background = element_rect(fill = NA, 
            color = NA), axis.text = element_text(color = "black"), 
        axis.text.x = element_text(angle = rotang, hjust = 1, vjust = 1), 
        legend.key.height = unit(0.8, "lines"))
    switch(by, sample_id = p + (if (!is.null(group_by)) facet_wrap(group_by, 
        scales = "free_x",ncol = n_cols)) + geom_bar(aes_string(x = "sample_id", 
        fill = "cluster_id"), position = "fill", stat = "identity") + 
        scale_fill_manual("cluster_id", values = k_pal) + scale_x_discrete(expand = c(0, 
        0)) + scale_y_continuous(expand = c(0, 0), labels = seq(0, 
        100, 25)) + theme(panel.border = element_blank(), panel.spacing.x = unit(1, 
        "lines")), cluster_id = {
        p <- p + scale_shape_manual(values = shapes) + guides(col = guide_legend(order = 1, 
            override.aes = list(size = 3)), shape = guide_legend(override.aes = list(size = 3)))
        if (is.null(group_by)) {
            p + geom_boxplot(aes_string(x = "cluster_id"), alpha = 0.2, 
                position = position_dodge(), outlier.color = NA) + 
                geom_point(aes_string("cluster_id", shape = shape_by), 
                  position = position_jitter(width = 0.2))
        } else {
            p + facet_wrap("cluster_id", scales = "free_y", ncol = n_cols) + 
                geom_boxplot(aes_string(x = group_by, color = group_by, 
                  fill = group_by), position = position_dodge(), 
                  alpha = 0.2, outlier.color = NA, show.legend = FALSE) + 
                geom_point(aes_string(x = group_by, col = group_by, 
                  shape = shape_by), position = position_jitter(width = 0.2))
        }
    })
}

environment(plotAbundances23) <- asNamespace('CATALYST')



## ----------------------------------------------------------------------------------------------------------
# Make sure you have the necessary libraries loaded
# install.packages(c("ggplot2", "dplyr", "CATALYST", "scales", "ggcyto"))
library(ggplot2)
library(dplyr)
library(CATALYST)
library(scales)
library(ggcyto)


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param k_abundances PARAM_DESCRIPTION, Default: 'merging1'
#' @param sample_order PARAM_DESCRIPTION, Default: NULL
#' @param meta PARAM_DESCRIPTION, Default: c("sample_id", "patient_id", "condition")
#' @param facet_by PARAM_DESCRIPTION, Default: NULL
#' @param group_by PARAM_DESCRIPTION, Default: 'condition'
#' @param shape_by PARAM_DESCRIPTION, Default: NULL
#' @param text_size PARAM_DESCRIPTION, Default: 16
#' @param n_cols PARAM_DESCRIPTION, Default: 4
#' @param wrap_cols PARAM_DESCRIPTION, Default: 1
#' @param rotang PARAM_DESCRIPTION, Default: 45
#' @param k_pal PARAM_DESCRIPTION, Default: CATALYST:::.cluster_cols
#' @param average_across_samples PARAM_DESCRIPTION, Default: TRUE
#' @param plot_ci PARAM_DESCRIPTION, Default: FALSE
#' @param show_cell_counts PARAM_DESCRIPTION, Default: FALSE
#' @param cell_count_position PARAM_DESCRIPTION, Default: c("above", "below")
#' @param cell_count_size PARAM_DESCRIPTION, Default: 3
#' @return OUTPUT_DESCRIPTION
#' @export 
#' @importFrom CATALYST .cluster_cols
#' @importFrom dplyr group_by summarise left_join arrange mutate ungroup across all_of
#' @importFrom scales percent_format
plotAbundanceStacked <- function(x,
                                 k_abundances = "merging1",
                                 sample_order = NULL,
                                 meta = c("sample_id", "patient_id", "condition"),
                                 facet_by = NULL,
                                 group_by = "condition",
                                 shape_by = NULL,
                                 text_size = 16,
                                 n_cols = 4,
                                 wrap_cols = 1,
                                 rotang = 45,
                                 k_pal = CATALYST:::.cluster_cols,
                                 average_across_samples = TRUE,
                                 plot_ci = FALSE,
                                 show_cell_counts = FALSE,
                                 cell_count_position = c("above", "below"),
                                 cell_count_size = 3) {
  # --- Argument Validation ---
  cell_count_position <- match.arg(cell_count_position)

  if (average_across_samples && show_cell_counts) {
    warning("`show_cell_counts = TRUE` is only implemented when `average_across_samples = FALSE`.",
            "\nIgnoring `show_cell_counts`.",
            call. = FALSE)
    show_cell_counts <- FALSE
  }

  # --- Data Preparation ---
  cluster_ids_abundances <- colData(x)[[k_abundances]]
  ns <- table(cluster_id = cluster_ids_abundances, sample_id = colData(x)$sample_id)
  fq <- prop.table(ns, 2) * 100

  df_fq <- as.data.frame(fq)
  df_counts <- as.data.frame(ns)

  df <- merge(df_fq,
              df_counts,
              by = c("cluster_id", "sample_id"),
              suffixes = c("_pct", "_count"))
  names(df)[names(df) == "Freq_count"] <- "cell_count"
  names(df)[names(df) == "Freq_pct"] <- "Freq"

  m <- match(df$sample_id, colData(x)$sample_id)
  for (i in meta)
    df[[i]] <- colData(x)[[i]][m]

  dfout <<- df

  # --- Data Averaging (if requested) ---
  if (average_across_samples) {
    df_avg <- df %>%
      dplyr::group_by(cluster_id,!!sym(group_by)) %>%
      dplyr::summarise(Freq = mean(Freq), .groups = "drop")
    
    df_ci <- df %>%
      dplyr::group_by(cluster_id,!!sym(group_by)) %>%
      dplyr::summarise(
        CI_lower = mean(Freq) - 1.96 * (sd(Freq) / sqrt(n())),
        CI_upper = mean(Freq) + 1.96 * (sd(Freq) / sqrt(n())),
        .groups = "drop"
      )
    
    df2 <- dplyr::left_join(df_avg, df_ci, by = c("cluster_id", group_by))
    
    df2 <- df2 %>%
      dplyr::arrange(desc(cluster_id)) %>%
      dplyr::group_by(!!sym(group_by)) %>%
      dplyr::mutate(cumFreq = cumsum(Freq)) %>%
      dplyr::ungroup()
    
    meta_df <- ei(x)[, c(group_by, facet_by)]
    
    df2 <- dplyr::left_join(df2, meta_df, by = group_by)
    df <- df2
  }
  
  dfoutavg <<- df

  # --- Factor Ordering ---
  if (!is.null(sample_order)) {
    df$sample_id <-
      factor(df$sample_id, ordered = TRUE, levels = sample_order)
  }

  # --- Plotting ---
  if (average_across_samples) {
    p <- ggplot(df, aes_string(x = group_by, y = "Freq", fill = "cluster_id")) +
      labs(x = NULL, y = "Proportion [%]") +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.background = element_rect(fill = "grey88", color = NA),
        axis.text = element_text(color = "black"),
        text = element_text(size = text_size),
        axis.text.x = element_text(angle = rotang, hjust = 1, vjust = 1),
        legend.key.height = unit(0.8, "lines"),
        panel.spacing.x = unit(0.4, "lines")
      ) +
      geom_bar(stat = "identity", position = "stack", color = "black", size = 0.2) +
      scale_fill_manual("cluster_id", values = k_pal) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0), labels = scales::percent_format(scale = 1)) +
      theme(panel.border = element_blank(), panel.spacing.x = unit(1, "lines"))
    
    if (!is.null(facet_by)) {
      p <- p + facet_grid2(cols = vars(!!sym(facet_by)), scales = "free_x", space = "free") +
        theme(strip.placement = "outside")
    }
    
    if (plot_ci) {
      p <- p + geom_errorbar(aes_string(ymin = "cumFreq - (Freq - CI_lower)", ymax = "cumFreq + (CI_upper - Freq)"),
                             width = 0.5, position = position_dodge(width = 0.9))
    }
    
  } else {
    p <- ggplot(df, aes_string(x = "sample_id", y = "Freq", fill = "cluster_id")) +
      labs(x = NULL, y = "Proportion [%]") +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.background = element_rect(fill = "grey88", color = NA),
        axis.text = element_text(color = "black"),
        text = element_text(size = text_size),
        axis.text.x = element_text(angle = rotang, hjust = 1, vjust = 1),
        legend.key.height = unit(0.8, "lines"),
        panel.spacing.x = unit(0.1, "lines")
      ) +
      geom_bar(stat = "identity", position = "stack", color = "black", size = 0.2) +
      scale_fill_manual("cluster_id", values = k_pal) +
      scale_x_discrete(expand = c(0, 0)) +
      theme(panel.border = element_blank(), panel.spacing.x = unit(1, "lines"))

    # Robust faceting
    if (!is.null(facet_by)) {
      p <- p + facet_wrap(vars(!!sym(facet_by)), scales = "free_x", ncol = n_cols)
    }

    # --- ADD CELL COUNTS (NEW FEATURE) ---
    if (show_cell_counts) {
      # Be explicit with dplyr:: to avoid conflicts
      grouping_vars <- c("sample_id", facet_by)
      sample_counts <- df %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(grouping_vars))) %>%
        dplyr::summarise(total_cells = sum(cell_count), .groups = 'drop')

      if (cell_count_position == "above") {
        y_pos <- 100
        vjust_val <- -0.5
        y_expand <- expansion(mult = c(0, 0.1))
      } else {
        y_pos <- 0
        vjust_val <- 1.5
        y_expand <- expansion(mult = c(0.1, 0))
      }

      p <- p +
        geom_text(
          data = sample_counts,
          aes(x = .data$sample_id, y = y_pos, label = .data$total_cells),
          inherit.aes = FALSE,
          vjust = vjust_val,
          size = cell_count_size,
          color = "black"
        ) +
        scale_y_continuous(expand = y_expand, labels = scales::percent_format(scale = 1))
      
    } else {
      p <- p + scale_y_continuous(expand = c(0, 0), labels = scales::percent_format(scale = 1))
    }
  }

  return(p)
}

# It's good practice to set the environment for non-exported functions
environment(plotAbundanceStacked) <- asNamespace('CATALYST')


## ----------------------------------------------------------------------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param k_abundances PARAM_DESCRIPTION, Default: 'merging1'
#' @param sample_order PARAM_DESCRIPTION, Default: NULL
#' @param meta PARAM_DESCRIPTION, Default: c("sample_id", "patient_id", "condition")
#' @param facet_by PARAM_DESCRIPTION, Default: NULL
#' @param group_by PARAM_DESCRIPTION, Default: 'condition'
#' @param shape_by PARAM_DESCRIPTION, Default: NULL
#' @param text_size PARAM_DESCRIPTION, Default: 16
#' @param n_cols PARAM_DESCRIPTION, Default: 4
#' @param wrap_cols PARAM_DESCRIPTION, Default: 1
#' @param rotang PARAM_DESCRIPTION, Default: 45
#' @param k_pal PARAM_DESCRIPTION, Default: CATALYST:::.cluster_cols
#' @param average_across_samples PARAM_DESCRIPTION, Default: TRUE
#' @param plot_ci PARAM_DESCRIPTION, Default: FALSE
#' @param cell_threshold PARAM_DESCRIPTION, Default: 0
#' @param drop_condition_if_any_sample_below_threshold PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @export 
#' @importFrom CATALYST .cluster_cols
#' @importFrom scales percent_format
plotAbundanceStacked2 <- function (x, 
                                  k_abundances = "merging1",
                                  sample_order = NULL, 
                                  meta = c("sample_id", "patient_id", "condition"), 
                                  facet_by = NULL,
                                  group_by = "condition", 
                                  shape_by = NULL, 
                                  text_size = 16,
                                  n_cols = 4, 
                                  wrap_cols = 1,
                                  rotang = 45, 
                                  k_pal = CATALYST:::.cluster_cols, 
                                  average_across_samples = TRUE,
                                  plot_ci = FALSE,
                                  cell_threshold = 0,
                                  drop_condition_if_any_sample_below_threshold = FALSE) 
{
  # Extract cluster IDs from colData
  cluster_ids_abundances <- colData(x)[[k_abundances]]
  
  ns <- table(cluster_id = cluster_ids_abundances, sample_id = colData(x)$sample_id)
  
  # Initialize dropped samples and conditions
  dropped_samples <- character(0)
  dropped_conditions <- character(0)
  
  # Filter samples based on the cell threshold if specified
  if (cell_threshold > 0) {
    sample_totals <- colSums(ns)
    
    # When dropping conditions: if any sample in a condition has count below threshold,
    # then drop the entire condition (only when averaging).
    if (drop_condition_if_any_sample_below_threshold && average_across_samples) {
      samp_df <- data.frame(sample_id = names(sample_totals), 
                            count = as.numeric(sample_totals),
                            stringsAsFactors = FALSE)
      meta_samp <- as.data.frame(colData(x)[, c("sample_id", group_by)], 
                                  stringsAsFactors = FALSE)
      samp_df <- merge(samp_df, meta_samp, by = "sample_id")
      conditions_to_drop <- unique(samp_df[[group_by]][samp_df$count < cell_threshold])
      valid_conditions <- setdiff(unique(samp_df[[group_by]]), conditions_to_drop)
      valid_samples <- samp_df$sample_id[samp_df[[group_by]] %in% valid_conditions & samp_df$count >= cell_threshold]
      dropped_conditions <- conditions_to_drop
      dropped_samples <- setdiff(names(sample_totals), valid_samples)
    } else {
      valid_samples <- names(sample_totals)[sample_totals >= cell_threshold]
      dropped_samples <- names(sample_totals)[sample_totals < cell_threshold]
      
      # Additionally, report conditions that have no valid samples.
      meta_samp <- as.data.frame(colData(x)[, c("sample_id", group_by)], stringsAsFactors = FALSE)
      all_conditions <- unique(meta_samp[[group_by]])
      valid_cond <- unique(meta_samp[[group_by]][meta_samp$sample_id %in% valid_samples])
      dropped_conditions <- setdiff(all_conditions, valid_cond)
    }
    
    ns <- ns[, colnames(ns) %in% valid_samples, drop = FALSE]
    
    if (length(dropped_samples) > 0) {
      message("Dropped samples (total cell count < ", cell_threshold, "):")
      for (s in dropped_samples) {
        message(" - ", s)
      }
    }
    if (length(dropped_conditions) > 0) {
      message("Dropped conditions (no samples with total cell count >= ", cell_threshold, "):")
      for (c in dropped_conditions) {
        message(" - ", c)
      }
    }
  }
  
  fq <- prop.table(ns, 2) * 100
  
  # Create data frames for frequencies and counts
  df_fq <- as.data.frame(fq)
  df_counts <- as.data.frame(ns)
  
  # Merge frequency and cell count info; 
  # df_fq has Freq as percentage, df_counts has Freq as count.
  df <- merge(df_fq, df_counts, by = c("cluster_id", "sample_id"), 
              suffixes = c("_pct", "_count"))
  names(df)[names(df) == "Freq_count"] <- "cell_count"
  names(df)[names(df) == "Freq_pct"] <- "Freq"
  
  m <- match(df$sample_id, colData(x)$sample_id)
  for (i in meta) df[[i]] <- colData(x)[[i]][m]
  
  dfout <<- df
  
  if (average_across_samples) {
    # Calculate the average frequencies
    df_avg <- df %>%
      group_by(cluster_id, !!sym(group_by)) %>%
      summarise(Freq = mean(Freq), .groups = "drop")
    
    # Calculate the 95% CI
    df_ci <- df %>%
      group_by(cluster_id, !!sym(group_by)) %>%
      summarise(CI_lower = mean(Freq) - 1.96 * (sd(Freq) / sqrt(n())),
                CI_upper = mean(Freq) + 1.96 * (sd(Freq) / sqrt(n())),
                .groups = "drop")
    
    # Join the CI back to the averages
    df2 <- left_join(df_avg, df_ci, by = c("cluster_id", group_by))
    
    # Calculate cumulative sum for error bar positioning
    df2 <- df2 %>%
      arrange(desc(cluster_id)) %>%
      group_by(!!sym(group_by)) %>%
      mutate(cumFreq = cumsum(Freq)) %>%
      ungroup()
    
    meta_df <- ei2(x, meta = c(group_by, facet_by))
    
    # Use left_join to avoid duplications
    df2 <- left_join(df2, meta_df, by = group_by)
    
    df <- df2
  }
  
  dfoutavg <<- df
  
  if (!is.null(sample_order)) {
    df$sample_id <- factor(df$sample_id, ordered = TRUE, levels = sample_order)
  }
  
  if (average_across_samples) {
    p <- ggplot(df, aes_string(x = group_by, y = "Freq", fill = "cluster_id")) + 
      labs(x = NULL, y = "Proportion [%]") + 
      theme_bw() + 
      theme(panel.grid = element_blank(), 
            strip.text = element_text(face = "bold"), 
            strip.background = element_rect(fill = "grey88", color = NA), 
            axis.text = element_text(color = "black"), 
            text = element_text(size = text_size),
            axis.text.x = element_text(angle = rotang, hjust = 1, vjust = 1), 
            legend.key.height = unit(0.8, "lines"),
            panel.spacing.x = unit(0.4, "lines")) + 
      geom_bar(stat = "identity", position = "stack", color = "black", size = 0.2) + 
      scale_fill_manual("cluster_id", values = k_pal) + 
      scale_x_discrete(expand = c(0, 0)) + 
      scale_y_continuous(expand = c(0, 0), labels = scales::percent_format(scale = 1)) + 
      theme(panel.border = element_blank(), panel.spacing.x = unit(1, "lines"))
    
    if (!is.null(facet_by)) {
      p <- p + facet_grid2(cols = vars(!!sym(facet_by)), scales = "free_x", space = "free") +
              theme(strip.placement = "outside")
    }
    
    if (plot_ci) {
      p <- p + geom_errorbar(aes_string(ymin = "cumFreq - (Freq - CI_lower)", 
                                          ymax = "cumFreq + (CI_upper - Freq)"), 
                             width = 0.5, position = position_dodge(width = 0.9))
    }
  } else {
    p <- ggplot(df, aes_string(x = "sample_id", y = "Freq", fill = "cluster_id")) + 
      labs(x = NULL, y = "Proportion [%]") + 
      theme_bw() + 
      theme(panel.grid = element_blank(), 
            strip.text = element_text(face = "bold"), 
            strip.background = element_rect(fill = "grey88", color = NA), 
            axis.text = element_text(color = "black"), 
            text = element_text(size = text_size),
            axis.text.x = element_text(angle = rotang, hjust = 1, vjust = 1), 
            legend.key.height = unit(0.8, "lines"),
            panel.spacing.x = unit(0.1, "lines")) +  
      facet_wrap(facet_by, scales = "free_x", ncol = n_cols) + 
      geom_bar(stat = "identity", position = "stack", color = "black", size = 0.2) + 
      scale_fill_manual("cluster_id", values = k_pal) + 
      scale_x_discrete(expand = c(0, 0)) + 
      scale_y_continuous(expand = c(0, 0), labels = scales::percent_format(scale = 1)) + 
      theme(panel.border = element_blank(), panel.spacing.x = unit(1, "lines"))
  }
  
  p <- p + scale_fill_manual(values = k_pal)
  
  return(p)
}

environment(plotAbundanceStacked) <- asNamespace('CATALYST')


## ----------------------------------------------------------------------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param k_abundances PARAM_DESCRIPTION, Default: 'merging1'
#' @param group_by_order PARAM_DESCRIPTION, Default: NULL
#' @param sample_order PARAM_DESCRIPTION, Default: NULL
#' @param meta PARAM_DESCRIPTION, Default: c("sample_id", "patient_id", "condition")
#' @param group_by PARAM_DESCRIPTION, Default: 'condition'
#' @param shape_by PARAM_DESCRIPTION, Default: NULL
#' @param text_size PARAM_DESCRIPTION, Default: 16
#' @param n_cols PARAM_DESCRIPTION, Default: 4
#' @param rotang PARAM_DESCRIPTION, Default: 45
#' @param k_pal PARAM_DESCRIPTION, Default: CATALYST:::.cluster_cols
#' @param average_across_samples PARAM_DESCRIPTION, Default: TRUE
#' @param plot_ci PARAM_DESCRIPTION, Default: FALSE
#' @param cell_threshold PARAM_DESCRIPTION, Default: 0
#' @param drop_condition_if_any_sample_below_threshold PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @export 
#' @importFrom CATALYST .cluster_cols
plotAbundanceDonut <- function (x, 
                                k_abundances = "merging1",
                                group_by_order = NULL, 
                                sample_order = NULL,
                                meta = c("sample_id", "patient_id", "condition"), 
                                group_by = "condition", 
                                shape_by = NULL, 
                                text_size = 16,
                                n_cols = 4, 
                                rotang = 45, 
                                k_pal = CATALYST:::.cluster_cols, 
                                average_across_samples = TRUE,
                                plot_ci = FALSE,
                                cell_threshold = 0,
                                drop_condition_if_any_sample_below_threshold = FALSE) 
{
  cluster_ids_abundances <- colData(x)[[k_abundances]]
  
  ns <- table(cluster_id = cluster_ids_abundances, sample_id = colData(x)$sample_id)
  
  # Initialize dropped samples and conditions
  dropped_samples <- character(0)
  dropped_conditions <- character(0)
  
  # Filter samples based on the cell threshold if specified
  if (cell_threshold > 0) {
    sample_totals <- colSums(ns)
    if (drop_condition_if_any_sample_below_threshold && average_across_samples) {
      samp_df <- data.frame(sample_id = names(sample_totals), count = as.numeric(sample_totals),
                            stringsAsFactors = FALSE)
      meta_samp <- as.data.frame(colData(x)[, c("sample_id", group_by)], stringsAsFactors = FALSE)
      samp_df <- merge(samp_df, meta_samp, by = "sample_id")
      conditions_to_drop <- unique(samp_df[[group_by]][samp_df$count < cell_threshold])
      valid_conditions <- setdiff(unique(samp_df[[group_by]]), conditions_to_drop)
      valid_samples <- samp_df$sample_id[samp_df[[group_by]] %in% valid_conditions & samp_df$count >= cell_threshold]
      dropped_conditions <- conditions_to_drop
      dropped_samples <- setdiff(names(sample_totals), valid_samples)
    } else {
      valid_samples <- names(sample_totals)[sample_totals >= cell_threshold]
      dropped_samples <- names(sample_totals)[sample_totals < cell_threshold]
      
      # Additionally, report conditions that have no valid samples.
      meta_samp <- as.data.frame(colData(x)[, c("sample_id", group_by)], stringsAsFactors = FALSE)
      all_conditions <- unique(meta_samp[[group_by]])
      valid_cond <- unique(meta_samp[[group_by]][meta_samp$sample_id %in% valid_samples])
      dropped_conditions <- setdiff(all_conditions, valid_cond)
    }
    ns <- ns[, colnames(ns) %in% valid_samples, drop = FALSE]
    
    if (length(dropped_samples) > 0) {
      message("Dropped samples (total cell count < ", cell_threshold, "):")
      for (s in dropped_samples) {
        message(" - ", s)
      }
    }
    if (length(dropped_conditions) > 0) {
      message("Dropped conditions (no samples with total cell count >= ", cell_threshold, "):")
      for (c in dropped_conditions) {
        message(" - ", c)
      }
    }
  }
  
  fq <- prop.table(ns, 2) * 100
  df <- as.data.frame(fq)
  m <- match(as.character(df$sample_id), as.character(colData(x)$sample_id))
  for (i in meta) df[[i]] <- colData(x)[[i]][m]
  
  dfout <<- df
  
  if (average_across_samples) {
    df_avg <- df %>%
      group_by(cluster_id, !!sym(group_by)) %>%
      summarise(Freq = mean(Freq), .groups = 'drop')
    
    df_ci <- df %>%
      group_by(cluster_id, !!sym(group_by)) %>%
      summarise(CI_lower = mean(Freq) - 1.96 * (sd(Freq) / sqrt(n())),
                CI_upper = mean(Freq) + 1.96 * (sd(Freq) / sqrt(n())), .groups = 'drop')
    
    df2 <- left_join(df_avg, df_ci, by = c("cluster_id", group_by))
    
    df2 <- df2 %>%
      arrange(desc(cluster_id)) %>%
      group_by(!!sym(group_by)) %>%
      mutate(cumFreq = cumsum(Freq)) %>%
      ungroup()
    
    meta_data <- ei2(x, meta = c(group_by))
    df2 <- left_join(df2, meta_data, by = group_by)
    
    df <- df2
    if (!is.null(group_by_order)) {
      df[[group_by]] <- factor(df[[group_by]], ordered = TRUE, levels = group_by_order)
    }
  } else {
    # When not averaging, allow reordering of donuts by sample_order.
    if (!is.null(sample_order)) {
      df$sample_id <- factor(df$sample_id, ordered = TRUE, levels = sample_order)
    }
  }
  
  dfoutavg <<- df

  if (average_across_samples) {
    p <- ggplot(df, aes(x = 2, y = Freq, fill = cluster_id)) + 
      geom_bar(width = 1, stat = "identity", color = "black", size = 0.2) +
      coord_polar(theta = "y") + 
      xlim(1, 2.5) + 
      labs(x = NULL, y = "Proportion [%]") + 
      theme_void() + 
      theme(
        strip.text = element_text(face = "bold"), 
        strip.background = element_rect(fill = NA, color = NA),
        text = element_text(size = text_size),
        legend.key.height = unit(0.8, "lines")
      ) +
      scale_fill_manual("cluster_id", values = k_pal) +
      facet_wrap(vars(!!sym(group_by)), ncol = n_cols)
    
    if (plot_ci) {
      p <- p + geom_errorbar(aes(ymin = cumFreq - (Freq - CI_lower), ymax = cumFreq + (CI_upper - Freq)), 
                             width = 0.5, position = position_dodge(width = 0.9))
    }
  } else {
    # For non-averaged data, facet by sample_id for one donut per sample.
    p <- ggplot(df, aes(x = 2, y = Freq, fill = cluster_id)) + 
      geom_bar(width = 1, stat = "identity", color = "black", size = 0.2) +
      coord_polar(theta = "y") + 
      xlim(1, 2.5) + 
      labs(x = NULL, y = "Proportion [%]") + 
      theme_void() + 
      theme(
        strip.text = element_text(face = "bold"), 
        strip.background = element_rect(fill = NA, color = NA),
        text = element_text(size = text_size),
        legend.key.height = unit(0.8, "lines")
      ) +
      scale_fill_manual("cluster_id", values = k_pal) +
      facet_wrap(vars(sample_id), ncol = n_cols)
  }
  
  p <- p + scale_fill_manual(values = k_pal)
  
  return(p)
}

environment(plotAbundanceDonut) <- asNamespace('CATALYST')


## ----------------------------------------------------------------------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param k_abundances PARAM_DESCRIPTION, Default: 'merging1'
#' @param k_daughter PARAM_DESCRIPTION, Default: NULL
#' @param k_pal PARAM_DESCRIPTION, Default: CATALYST:::.cluster_cols
#' @param k_pal_daughter PARAM_DESCRIPTION, Default: NULL
#' @param group_by_order PARAM_DESCRIPTION, Default: NULL
#' @param sample_order PARAM_DESCRIPTION, Default: NULL
#' @param meta PARAM_DESCRIPTION, Default: c("sample_id", "patient_id", "condition")
#' @param group_by PARAM_DESCRIPTION, Default: 'condition'
#' @param shape_by PARAM_DESCRIPTION, Default: NULL
#' @param text_size PARAM_DESCRIPTION, Default: 16
#' @param n_cols PARAM_DESCRIPTION, Default: 4
#' @param rotang PARAM_DESCRIPTION, Default: 45
#' @param average_across_samples PARAM_DESCRIPTION, Default: TRUE
#' @param plot_ci PARAM_DESCRIPTION, Default: FALSE
#' @param cell_threshold PARAM_DESCRIPTION, Default: 0
#' @param drop_condition_if_any_sample_below_threshold PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @export 
#' @importFrom CATALYST .cluster_cols
#' @importFrom dplyr select
#' @importFrom ggnewscale new_scale_fill
plotAbundanceDonut2 <- function (x, 
                                 k_abundances = "merging1",
                                 k_daughter = NULL,
                                 k_pal = CATALYST:::.cluster_cols, 
                                 k_pal_daughter = NULL,
                                 group_by_order = NULL, 
                                 sample_order = NULL,
                                 meta = c("sample_id", "patient_id", "condition"), 
                                 group_by = "condition", 
                                 shape_by = NULL, 
                                 text_size = 16,
                                 n_cols = 4, 
                                 rotang = 45, 
                                 average_across_samples = TRUE,
                                 plot_ci = FALSE,
                                 cell_threshold = 0,
                                 drop_condition_if_any_sample_below_threshold = FALSE) 
{
  library(dplyr)
  library(ggplot2)
  library(rlang)
  library(tidyr)
  library(ggnewscale)
  
  # ----- MASTER DONUT DATA -----
  cluster_ids_abundances <- colData(x)[[k_abundances]]
  ns <- table(cluster_id = cluster_ids_abundances, sample_id = colData(x)$sample_id)
  
  dropped_samples <- character(0)
  dropped_conditions <- character(0)
  
  if (cell_threshold > 0) {
    sample_totals <- colSums(ns)
    if (drop_condition_if_any_sample_below_threshold && average_across_samples) {
      samp_df <- data.frame(sample_id = names(sample_totals), count = as.numeric(sample_totals),
                            stringsAsFactors = FALSE)
      meta_samp <- as.data.frame(colData(x)[, c("sample_id", group_by)], stringsAsFactors = FALSE)
      samp_df <- merge(samp_df, meta_samp, by = "sample_id")
      conditions_to_drop <- unique(samp_df[[group_by]][samp_df$count < cell_threshold])
      valid_conditions <- setdiff(unique(samp_df[[group_by]]), conditions_to_drop)
      valid_samples <- samp_df$sample_id[samp_df[[group_by]] %in% valid_conditions & samp_df$count >= cell_threshold]
      dropped_conditions <- conditions_to_drop
      dropped_samples <- setdiff(names(sample_totals), valid_samples)
    } else {
      valid_samples <- names(sample_totals)[sample_totals >= cell_threshold]
      dropped_samples <- names(sample_totals)[sample_totals < cell_threshold]
      meta_samp <- as.data.frame(colData(x)[, c("sample_id", group_by)], stringsAsFactors = FALSE)
      all_conditions <- unique(meta_samp[[group_by]])
      valid_cond <- unique(meta_samp[[group_by]][meta_samp$sample_id %in% valid_samples])
      dropped_conditions <- setdiff(all_conditions, valid_cond)
    }
    ns <- ns[, colnames(ns) %in% valid_samples, drop = FALSE]
    if (length(dropped_samples) > 0) {
      message("Dropped samples (total cell count < ", cell_threshold, "):")
      for (s in dropped_samples) message(" - ", s)
    }
    if (length(dropped_conditions) > 0) {
      message("Dropped conditions (no samples with total cell count >= ", cell_threshold, "):")
      for (c in dropped_conditions) message(" - ", c)
    }
  }
  
  fq <- prop.table(ns, 2) * 100
  df <- as.data.frame(fq)
  m <- match(as.character(df$sample_id), as.character(colData(x)$sample_id))
  for (i in meta) {
    df[[i]] <- colData(x)[[i]][m]
  }
  dfout <<- df
  
  if (average_across_samples) {
    df_avg <- df %>%
      group_by(cluster_id, !!sym(group_by)) %>%
      summarise(Freq = mean(Freq), .groups = "drop")
    
    df_ci <- df %>%
      group_by(cluster_id, !!sym(group_by)) %>%
      summarise(CI_lower = mean(Freq) - 1.96 * (sd(Freq) / sqrt(n())),
                CI_upper = mean(Freq) + 1.96 * (sd(Freq) / sqrt(n())), .groups = "drop")
    
    df2 <- left_join(df_avg, df_ci, by = c("cluster_id", group_by))
    
    df2 <- df2 %>%
      arrange(desc(cluster_id)) %>%
      group_by(!!sym(group_by)) %>%
      mutate(cumFreq = cumsum(Freq)) %>%
      ungroup()
    
    meta_data <- ei2(x, meta = c(group_by))
    df2 <- left_join(df2, meta_data, by = group_by)
    df <- df2
    if (!is.null(group_by_order)) {
      df[[group_by]] <- factor(df[[group_by]], ordered = TRUE, levels = group_by_order)
    }
  } else {
    if (!is.null(sample_order)) {
      df$sample_id <- factor(df$sample_id, ordered = TRUE, levels = sample_order)
    }
  }
  
  # Re-normalize so each group sums to exactly 100.
  df <- df %>%
    group_by(!!sym(group_by)) %>%
    mutate(Freq = Freq/sum(Freq)*100) %>%
    ungroup()
  
  dfoutavg <<- df
  
  # ----- DAUGHTER DONUT DATA -----
  if (!is.null(k_daughter) && !is.null(k_pal_daughter)) {
    daughter_df <- as.data.frame(colData(x)) %>%
      dplyr::select(sample_id, !!sym(group_by), daughter = !!sym(k_daughter))
    
    daughter_counts <- daughter_df %>%
      group_by(sample_id, !!sym(group_by), daughter) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(sample_id, !!sym(group_by)) %>%
      mutate(prop = n/sum(n)) %>%
      ungroup()
    
    if (average_across_samples) {
      daughter_data <- daughter_counts %>%
        group_by(!!sym(group_by), daughter) %>%
        summarise(Freq = mean(prop)*100, .groups = "drop")
    } else {
      daughter_data <- daughter_counts
      if (!is.null(sample_order))
        daughter_data$sample_id <- factor(daughter_data$sample_id, ordered = TRUE, levels = sample_order)
    }
    if (!is.null(group_by_order))
      daughter_data[[group_by]] <- factor(daughter_data[[group_by]], ordered = TRUE, levels = group_by_order)
    # Re-normalize so each group sums to 100.
    daughter_data <- daughter_data %>%
      group_by(!!sym(group_by)) %>%
      mutate(Freq = Freq/sum(Freq)*100) %>%
      ungroup()
  }
  
  # ----- PLOTTING -----
  # Master donut (inner ring) at x = 2 with width = 1 (covers approx. 1.5 to 2.5)
  if (average_across_samples) {
    p <- ggplot(df, aes(x = 2, y = Freq, fill = cluster_id)) + 
      geom_bar(width = 1, stat = "identity", color = "black", size = 0.2) +
      coord_polar(theta = "y", start = 0) +
      labs(x = NULL, y = "Proportion [%]") +
      theme_void() +
      theme(strip.text = element_text(face = "bold"),
            strip.background = element_rect(fill = NA, color = NA),
            text = element_text(size = text_size),
            legend.key.height = unit(0.8, "lines")) +
      facet_wrap(vars(!!sym(group_by)), ncol = n_cols)
  } else {
    p <- ggplot(df, aes(x = 2, y = Freq, fill = cluster_id)) + 
      geom_bar(width = 1, stat = "identity", color = "black", size = 0.2) +
      coord_polar(theta = "y", start = 0) +
      labs(x = NULL, y = "Proportion [%]") +
      theme_void() +
      theme(strip.text = element_text(face = "bold"),
            strip.background = element_rect(fill = NA, color = NA),
            text = element_text(size = text_size),
            legend.key.height = unit(0.8, "lines")) +
      facet_wrap(vars(sample_id), ncol = n_cols)
  }
  
  p <- p + scale_fill_manual("Master", values = k_pal)
  
  # Daughter donut (outer ring) at x = 3.5 with width = 1
  if (!is.null(k_daughter) && !is.null(k_pal_daughter)) {
    p <- p + ggnewscale::new_scale_fill()
    p <- p +
      geom_bar(data = daughter_data,
               mapping = aes(x = 3.2, y = Freq, fill = daughter),
               stat = "identity", color = "black", size = 0.2, width = 1,
               inherit.aes = FALSE) +
      scale_fill_manual("Daughter", values = k_pal_daughter)
  }
  
  # Use scale_x_continuous to remove extra padding and ensure full circle.
  p <- p + scale_x_continuous(limits = c(1, 4), expand = c(0,0))
  
  # Also force y-scale to [0,100] to fill the circle.
  p <- p + scale_y_continuous(limits = c(0,100), expand = c(0,0))
  
  return(p)
}

environment(plotAbundanceDonut2) <- asNamespace('CATALYST')


## ----------------------------------------------------------------------------------------------------------
library(ConsensusClusterPlus)

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param features PARAM_DESCRIPTION, Default: 'type'
#' @param xdim PARAM_DESCRIPTION, Default: 10
#' @param ydim PARAM_DESCRIPTION, Default: 10
#' @param maxK PARAM_DESCRIPTION, Default: 20
#' @param verbose PARAM_DESCRIPTION, Default: TRUE
#' @param seed PARAM_DESCRIPTION, Default: 1
#' @param alg PARAM_DESCRIPTION, Default: 'hc'
#' @return OUTPUT_DESCRIPTION
#' @export 
cluster2 <- function (x, features = "type", xdim = 10, ydim = 10, maxK = 20, 
    verbose = TRUE, seed = 1, alg = "hc") 
{
    stopifnot(is(x, "SingleCellExperiment"))
    stopifnot(is.logical(verbose), length(verbose) == 1, vapply(list(xdim, 
        ydim, maxK, seed), function(arg) is.numeric(arg) && length(arg) == 
        1, logical(1)))
    features <- .get_features(x, features)
    if (is.null(marker_classes(x))) {
        rowData(x)$marker_class <- factor(c("state", "type")[as.numeric(rownames(x) %in% 
            features) + 1], levels = c("type", "state", "none"))
    }
    rowData(x)$used_for_clustering <- rownames(x) %in% features
    if (verbose) 
        message("o running FlowSOM clustering...")
    fsom <- ReadInput(flowFrame(t(assay(x, "exprs"))))
    som <- BuildSOM(fsom, colsToUse = features, silent = TRUE, 
        xdim = xdim, ydim = ydim)
    if (verbose) 
        message("o running ConsensusClusterPlus metaclustering...")
    pdf(NULL)
    mc <- suppressWarnings(suppressMessages(ConsensusClusterPlus(t(som$map$codes), 
        maxK = maxK, reps = 100, distance = "euclidean", seed = seed, clusterAlg = alg,
        plot = NULL)))
    dev.off()
    k <- xdim * ydim
    mcs <- seq_len(maxK)[-1]
    codes <- data.frame(seq_len(k), map(mc[-1], "consensusClass"))
    codes <- mutate_all(codes, function(u) factor(u, levels = sort(unique(u))))
    colnames(codes) <- c(sprintf("som%s", k), sprintf("meta%s", 
        mcs))
    x$cluster_id <- factor(som$map$mapping[, 1])
    metadata(x)$cluster_codes <- codes
    metadata(x)$SOM_codes <- som$map$codes
    metadata(x)$delta_area <- .plot_delta_area(mc)
    return(x)
}

environment(cluster2) <- asNamespace('CATALYST')


## ----------------------------------------------------------------------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param k PARAM_DESCRIPTION, Default: 'meta20'
#' @param features PARAM_DESCRIPTION, Default: 'type'
#' @param facet_by PARAM_DESCRIPTION, Default: 'condition'
#' @param cluster_rows PARAM_DESCRIPTION, Default: TRUE
#' @param add_border PARAM_DESCRIPTION, Default: TRUE
#' @param text_size PARAM_DESCRIPTION, Default: 12
#' @param heatmap_palette PARAM_DESCRIPTION, Default: NULL
#' @param y_shift PARAM_DESCRIPTION, Default: 0.01
#' @param panel_spacing PARAM_DESCRIPTION, Default: 0.4
#' @param alpha_amount PARAM_DESCRIPTION, Default: 0.6
#' @return OUTPUT_DESCRIPTION
#' @export 
plotClusterExprs3 <- function (x, k = "meta20", features = "type", facet_by = "condition",cluster_rows = TRUE, add_border = TRUE, text_size = 12, heatmap_palette = NULL, y_shift = 0.01, panel_spacing = 0.4, alpha_amount = 0.6) 
{
    .check_sce(x, TRUE)
    k <- .check_k(x, k)
    x$cluster_id <- cluster_ids(x, k)
    features <- .get_features(x, features)
    
    ms <- t(.agg(x[features, ], "cluster_id", "median"))
    
    if (cluster_rows) {
        d <- dist(ms, method = "euclidean")
        o <- hclust(d, method = "average")$order
    } else {
        o <- seq_along(levels(x$cluster_id))
    }
    
    cd <- colData(x)
    es <- assay(x[features, ], "exprs")
    df <- data.frame(t(es), cd, check.names = FALSE)
    df <- melt(df, id.vars = names(cd), variable.name = "antigen", 
               value.name = "expression")
    
    fq <- tabulate(x$cluster_id)/ncol(x)
    fq <- round(fq * 100, 2)
    names(fq) <- levels(x$cluster_id)
    cluster_levels <- rev(levels(x$cluster_id)[o])
    cluster_labels <- rev(names(fq)[o])
    df$cluster_id <- factor(df$cluster_id, levels = cluster_levels, labels = cluster_labels)
    
        dfhist <<- df
    
    # If a heatmap_palette is provided, use it to fill the colors, otherwise use red
    if (!is.null(heatmap_palette)) {
        cluster_colors <- heatmap_palette[levels(df$cluster_id)]
        scale_fill <- scale_fill_manual(values = cluster_colors)
    } else {
        scale_fill <- scale_fill_manual(values = rep("red", length(levels(df$cluster_id))))
    }
    
    # Define the base plot
    plot <- ggplot(df, aes(x = expression, y = cluster_id, fill = cluster_id)) +
        facet_wrap(~condition_day, scales = "free_x", nrow = 1) + 
        geom_density_ridges(alpha = alpha_amount, rel_min_height = 0.01, scale = 1.1) + 
        scale_fill +
        scale_x_continuous(breaks = seq(0, 10, by = 2), limits = c(-0.5, NA), expand = c(0, 0)) +
        scale_y_discrete(expand = expand_scale(mult = c(y_shift, 0.2))) +
        theme(
              legend.position = "none", 
              strip.background = element_blank(), 
              strip.text = element_text(face = "bold", color = "black", size = text_size),
              panel.spacing = unit(panel_spacing, "lines"),
              plot.margin = unit(c(0, 1, 0, 0), "lines"),
              panel.background = element_rect(fill = "transparent", color = NA), # Make panel background transparent
              panel.border = if (add_border) element_rect(color = "grey43", fill = NA, size = 0.8) else element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_line(color = "black"), # Set x major gridlines as black
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              axis.line.x = element_blank(), # Hide main x-axis line
              axis.line.y = element_blank(), # Hide main y-axis line
              axis.text.x = element_text(color = "black"),
              axis.text.y = element_text(color = "black"),
              text = element_text(size = text_size))
    
    plot
}

environment(plotClusterExprs3) <- asNamespace('CATALYST')


## ----------------------------------------------------------------------------------------------------------
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(rlang)

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param k PARAM_DESCRIPTION, Default: 'meta20'
#' @param scales PARAM_DESCRIPTION, Default: 'free_y'
#' @param conditions PARAM_DESCRIPTION, Default: NULL
#' @param excluded_clusters PARAM_DESCRIPTION, Default: NULL
#' @param axes PARAM_DESCRIPTION, Default: 'all'
#' @param features PARAM_DESCRIPTION, Default: 'state'
#' @param assay PARAM_DESCRIPTION, Default: 'exprs'
#' @param fun PARAM_DESCRIPTION, Default: c("median", "mean", "sum")
#' @param point_size PARAM_DESCRIPTION, Default: 1
#' @param clusters_order PARAM_DESCRIPTION, Default: NULL
#' @param textsize PARAM_DESCRIPTION, Default: 14
#' @param panel_spacing PARAM_DESCRIPTION, Default: 2
#' @param show_stats PARAM_DESCRIPTION, Default: TRUE
#' @param hide_ns PARAM_DESCRIPTION, Default: TRUE
#' @param facet_by PARAM_DESCRIPTION, Default: c("antigen", "cluster_id")
#' @param color_by PARAM_DESCRIPTION, Default: 'condition'
#' @param mean_or_med PARAM_DESCRIPTION, Default: 'median'
#' @param group_by PARAM_DESCRIPTION, Default: color_by
#' @param shape_by PARAM_DESCRIPTION, Default: NULL
#' @param size_by PARAM_DESCRIPTION, Default: FALSE
#' @param stat_size PARAM_DESCRIPTION, Default: 4
#' @param facet_ratio PARAM_DESCRIPTION, Default: 1.5
#' @param geom PARAM_DESCRIPTION, Default: c("boxes", "bar")
#' @param jitter PARAM_DESCRIPTION, Default: TRUE
#' @param group1keep PARAM_DESCRIPTION, Default: NULL
#' @param nudge PARAM_DESCRIPTION, Default: 0
#' @param step_increase PARAM_DESCRIPTION, Default: 0
#' @param hide_x_labels PARAM_DESCRIPTION, Default: FALSE
#' @param label_parse PARAM_DESCRIPTION, Default: FALSE
#' @param merging_col PARAM_DESCRIPTION, Default: NULL
#' @param swap PARAM_DESCRIPTION, Default: FALSE
#' @param external_stats PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @export 
plotPbExprsDiffCond <- function (
    x, 
    k = "meta20", 
    scales = "free_y", 
    conditions = NULL, 
    excluded_clusters = NULL, 
    axes = "all", 
    features = "state", 
    assay = "exprs", 
    fun = c("median", "mean", "sum"), 
    point_size = 1, 
    clusters_order = NULL, 
    textsize = 14, 
    panel_spacing = 2, 
    show_stats = TRUE, 
    hide_ns = TRUE, 
    facet_by = c("antigen", "cluster_id"), 
    color_by = "condition", 
    mean_or_med = "median", 
    group_by = color_by, 
    shape_by = NULL, 
    size_by = FALSE, 
    stat_size = 4, 
    facet_ratio = 1.5, 
    geom = c("boxes","bar"), 
    jitter = TRUE, 
    group1keep = NULL, 
    nudge = 0, 
    step_increase = 0, 
    hide_x_labels = FALSE, 
    label_parse = FALSE, 
    merging_col = NULL, 
    swap = FALSE,
    external_stats = NULL
) {
    fun <- match.arg(fun)
    facet_by <- match.arg(facet_by)
    stopifnot(is.logical(jitter), length(jitter) == 1)
    
    if (!is.null(merging_col)) {
        cluster_ids <- factor(x[[k]])
    } else {
        .check_sce(x)
        k <- .check_k(x, k)
        cluster_ids <- cluster_ids(x, k)
    }
    
    .check_assay(x, assay)
    .check_cd_factor(x, color_by)
    .check_cd_factor(x, group_by)
    
    x <- x[.get_features(x, features), ]
    x$cluster_id <- cluster_ids
    by <- c("cluster_id", "sample_id")
    ms <- .agg(x, by, fun, assay)
    df <- melt(ms, varnames = c("antigen", by[length(by)]))
    
    if (!is.null(cluster_ids)) {
        df$cluster_id <- df$L1
    }
    
    i <- match(df$sample_id, x$sample_id)
    j <- setdiff(names(colData(x)), c(names(df), "cluster_id"))
    df <- cbind(df, colData(x)[i, j, drop = FALSE])

    DFCOND <<- df
    
    ncs <- table(as.list(colData(x)[by]))
    ncs <- rep(c(t(ncs)), each = nrow(x))
    
    if (size_by) {
        size_by <- "n_cells"
        df$n_cells <- ncs
    } else {
        size_by <- NULL
    }
    
    df <- df[ncs > 0, , drop = FALSE]
    
    if (!is.null(conditions)) {
        df <- df[df[[color_by]] %in% conditions, ]
        df[[color_by]] <- factor(df[[color_by]], levels = conditions, ordered = TRUE)
    }
    
    if (!is.null(excluded_clusters)) {
        df <- df[!df$cluster_id %in% excluded_clusters, ]
    }
    
    if (show_stats) {
        if (!is.null(external_stats)) {
            external_stats <- external_stats %>% filter(antigen %in% features)
            dummy_stat <- df %>% group_by(antigen, cluster_id) %>% wilcox_test(as.formula(paste("value ~", color_by)))
            dummy_stat <- dummy_stat %>% add_y_position(scales = "free", step.increase = step_increase)
            external_stats <- external_stats %>%
                left_join(dummy_stat %>% select(antigen, cluster_id, group1, group2, y.position), 
                          by = c("antigen", "cluster_id", "group1", "group2"))
            stat.test <- external_stats
        } else {
            stat.test <- df %>%
                group_by(antigen, cluster_id) %>%
                wilcox_test(as.formula(paste("value ~", color_by))) %>%
                add_y_position(scales = "free", step.increase = step_increase)
        }
        
        if (!is.null(group1keep)) {
            stat.test <- stat.test[stat.test$group1 %in% group1keep,]
            stat.test <- stat.test %>% add_y_position(scales = "free", step.increase = step_increase)
        }
        
        pbCondStats <<- stat.test
    }
    
    if (mean_or_med == "median") {
        fun = "median"
    } else {
        fun = "mean"
    }
    
    summary_df <- df %>%
        group_by(antigen, cluster_id, !!sym(color_by), condition) %>%
        summarise(
            mean_value = mean(value),
            sem = sd(value) / sqrt(n())
        )
    
    sumout <<- summary_df
    
# Calculate jittered positions
set.seed(123)  # Set seed for reproducibility
summary_df$jittered_x <- jitter(as.numeric(summary_df[[group_by]]), amount = 0.1)

# Create the plot
p <- ggplot(summary_df, aes(x = jittered_x, y = mean_value, color = antigen, group = antigen)) +
    geom_point(size = point_size) +
    geom_errorbar(aes(ymin = mean_value - sem, ymax = mean_value + sem), width = 0.2) +
    geom_line() +
  scale_x_continuous(breaks = unique(as.numeric(summary_df[[group_by]])), labels = unique(summary_df[[group_by]])) +
        labs(y = paste(fun, ifelse(assay == "exprs", "expression", assay))) +
        theme_bw() +
        theme(
            panel.grid = element_blank(),
            text = element_text(size = textsize),
            strip.text = element_text(size = textsize),
            strip.background = element_rect(fill = NA, color = NA),
            strip.placement = "outside",
            legend.text = element_text(size = textsize),
            legend.title = element_text(size = textsize),
            aspect.ratio = facet_ratio,
            axis.ticks = element_line(color = "black"),
            panel.border = element_rect(color = "black", fill = NA, size = 0.5),
            panel.spacing = unit(panel_spacing, "lines"),
            axis.text = element_text(color = "black", size = textsize),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = textsize)
        )
    
    if (hide_x_labels) {
        p <- p + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
    }
    
    if (show_stats) {
        p <- p + stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.025, hide.ns = hide_ns, size = stat_size, bracket.nudge.y = nudge)
    }
    
    if (label_parse) {
        if (!swap) {
            p <- p + facet_grid2(cluster_id ~ condition, scales = scales, axes = axes, labeller = label_parsed)
        } else {
            p <- p + facet_grid2(condition ~ cluster_id, scales = scales, axes = axes, labeller = label_parsed)
        }
    } else {
        if (!swap) {
            p <- p + facet_grid2(cluster_id ~ condition, scales = scales, axes = axes)
        } else {
            p <- p + facet_grid2(condition ~ cluster_id, scales = scales, axes = axes)
        }
    }
    
    return(p)
}

environment(plotPbExprsDiffCond) <- asNamespace("CATALYST")


## ----------------------------------------------------------------------------------------------------------
library(SingleCellExperiment)
library(ggplot2)
library(ggbeeswarm)
library(scales)
library(dplyr)
library(rstatix)
library(rlang)
library(ggh4x)
library(ggrastr)

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param scales PARAM_DESCRIPTION, Default: 'free_y'
#' @param conditions PARAM_DESCRIPTION, Default: NULL
#' @param color_by PARAM_DESCRIPTION, Default: NULL
#' @param excluded_clusters PARAM_DESCRIPTION, Default: NULL
#' @param axes PARAM_DESCRIPTION, Default: 'all'
#' @param features PARAM_DESCRIPTION, Default: 'state'
#' @param assay PARAM_DESCRIPTION, Default: 'exprs'
#' @param fun PARAM_DESCRIPTION, Default: c("median", "mean", "sum")
#' @param point_size PARAM_DESCRIPTION, Default: 1
#' @param clusters_order PARAM_DESCRIPTION, Default: NULL
#' @param textsize PARAM_DESCRIPTION, Default: 14
#' @param panel_spacing PARAM_DESCRIPTION, Default: 2
#' @param show_stats PARAM_DESCRIPTION, Default: TRUE
#' @param line PARAM_DESCRIPTION, Default: FALSE
#' @param hide_ns PARAM_DESCRIPTION, Default: TRUE
#' @param facet_by PARAM_DESCRIPTION, Default: NULL
#' @param x_group PARAM_DESCRIPTION, Default: 'condition'
#' @param mean_or_med PARAM_DESCRIPTION, Default: 'median'
#' @param group_by PARAM_DESCRIPTION, Default: x_group
#' @param shape_by PARAM_DESCRIPTION, Default: NULL
#' @param size_by PARAM_DESCRIPTION, Default: FALSE
#' @param point_alpha PARAM_DESCRIPTION, Default: 1
#' @param stat_size PARAM_DESCRIPTION, Default: 4
#' @param facet_ratio PARAM_DESCRIPTION, Default: 1.5
#' @param geom PARAM_DESCRIPTION, Default: c("boxes", "bar")
#' @param jitter PARAM_DESCRIPTION, Default: TRUE
#' @param geom_type PARAM_DESCRIPTION, Default: c("quasirandom", "point")
#' @param group1keep PARAM_DESCRIPTION, Default: NULL
#' @param nudge PARAM_DESCRIPTION, Default: 0
#' @param step_increase PARAM_DESCRIPTION, Default: 0.1
#' @param hide_x_labels PARAM_DESCRIPTION, Default: FALSE
#' @param label_parse PARAM_DESCRIPTION, Default: FALSE
#' @param merging_col PARAM_DESCRIPTION, Default: NULL
#' @param swap PARAM_DESCRIPTION, Default: FALSE
#' @param external_stats PARAM_DESCRIPTION, Default: NULL
#' @param average_samples PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @export 
plotPbExprsDiffSimp <- function (
    x, 
    scales = "free_y", 
    conditions = NULL, 
    color_by = NULL,
    excluded_clusters = NULL, 
    axes = "all", 
    features = "state", 
    assay = "exprs", 
    fun = c("median", "mean", "sum"), 
    point_size = 1, 
    clusters_order = NULL, 
    textsize = 14, 
    panel_spacing = 2, 
    show_stats = TRUE, 
    line = FALSE,
    hide_ns = TRUE, 
    facet_by = NULL, 
    x_group = "condition", 
    mean_or_med = "median", 
    group_by = x_group, 
    shape_by = NULL, 
    size_by = FALSE, 
    point_alpha = 1,
    stat_size = 4, 
    facet_ratio = 1.5, 
    geom = c("boxes","bar"), 
    jitter = TRUE, 
    geom_type = c("quasirandom", "point"),  # New parameter
    group1keep = NULL, 
    nudge = 0, 
    step_increase = 0.1, 
    hide_x_labels = FALSE, 
    label_parse = FALSE, 
    merging_col = NULL, 
    swap = FALSE,
    external_stats = NULL,  # New parameter for external stats dataframe
    average_samples = FALSE  # New parameter to average across samples
) {
    fun <- match.arg(fun)
    geom <- match.arg(geom)
    geom_type <- match.arg(geom_type)  # Match the new geom_type argument
    stopifnot(is.logical(jitter), length(jitter) == 1)
    
    # Allow plotting from column data
    .check_sce(x)
    .check_assay(x, assay)
    .check_cd_factor(x, x_group)
    .check_cd_factor(x, group_by)
    
    # Retrieve the features
    x <- x[.get_features(x, features), ]
    
    by <- "sample_id"
    ms <- .agg(x, by, fun, assay)
    df <- melt(ms, varnames = c("antigen", by))
    
    DFCOND <<- df
    
    i <- match(df$sample_id, x$sample_id)
    j <- setdiff(names(colData(x)), c(names(df)))
    df <- cbind(df, colData(x)[i, j, drop = FALSE])
    
    if (size_by) {
        size_by <- "n_cells"
        df$n_cells <- table(sample_ids(x))[df$sample_id]
    } else {
        size_by <- NULL
    }
    
    if (!is.null(conditions)) {
        df <- df[df[[x_group]] %in% conditions, ]
        df[[x_group]] <- factor(df[[x_group]], levels = conditions, ordered = TRUE)
    }
    
    # Apply excluded clusters if needed
    if (!is.null(excluded_clusters)) {
        df <- df[!df$cluster_id %in% excluded_clusters, ]
    }
    
    dfout <<- df
    
    # Average across samples if average_samples is TRUE
    if (average_samples) {
        df <- df %>%
            group_by(across(all_of(c(color_by, x_group, "antigen")))) %>%
            summarise(value = mean(value, na.rm = TRUE), .groups = 'drop')
    }
    
    # Check if the data is empty after filtering
    if (nrow(df) == 0) {
        stop("No data available after filtering. Please adjust the parameters.")
    }
    
    dfoutout <<- df
    
    if (show_stats) {
        if (!is.null(external_stats)) {
            # Filter external stats to include only the antigens that are plotted
            external_stats <- external_stats %>% filter(antigen %in% features)
            
            # Use external stats if provided
            dummy_stat <- df %>% group_by(antigen, !!sym(facet_by)) %>% wilcox_test(as.formula(paste("value ~", x_group)))
            dummy_stat <- dummy_stat %>% add_y_position(scales = "free", step.increase = step_increase)
            
            # Merge y.position from dummy_stat to external_stats
            external_stats <- external_stats %>%
                left_join(dummy_stat %>% select(antigen, !!sym(facet_by), group1, group2, y.position), 
                          by = c("antigen", facet_by, "group1", "group2"))
            
            stat.test <- external_stats
        } else {
            # Calculate stats internally using wilcox_test (change to Dunn?)
            stat.test <- df %>%
                group_by(antigen, !!sym(facet_by)) %>%
                wilcox_test(as.formula(paste("value ~", x_group))) %>%
                add_y_position(scales = "free", step.increase = step_increase)
        }
        
        if (!is.null(group1keep)) {
            stat.test <- stat.test[stat.test$group1 %in% group1keep,]
            stat.test <- stat.test %>% add_y_position(scales = "free", step.increase = step_increase)
        }
        
        # Save unfiltered stats globally
        pbCondStats <<- stat.test
    }
    
    if (mean_or_med == "median") {
        fun = "median"
    } else {
        fun = "mean"
    }
    
    # Initialize the ggplot
    p <- ggplot(df, aes_string(x = group_by, y = "value", fill = color_by)) +
        labs(y = paste(fun, ifelse(assay == "exprs", "expression", assay))) + # Plots "median expression"
        theme_bw() +
        theme(
            panel.grid = element_blank(),
            text = element_text(size = textsize),
            strip.text = element_text(size = textsize),
            strip.background = element_rect(fill = NA, color = NA),
            strip.placement = "outside",  # Move facet labels outside
            legend.text = element_text(size = textsize),
            legend.title = element_text(size = textsize),
            aspect.ratio = facet_ratio,
            axis.ticks = element_line(color = "black"),
            panel.border = element_rect(color = "black", fill = NA, size = 0.5),
            panel.spacing = unit(panel_spacing, "lines"),
            axis.text = element_text(color = "black", size = textsize),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = textsize)
        )
    
    # Add geom based on selection
    if (geom == "boxes") {
        p <- p + geom_boxplot(aes_string(group = x_group), fill = "grey84", color = "black", width = 0.75, linewidth = 0.3, alpha = 0.9, outlier.color = NA, show.legend = FALSE)
    } else if (geom == "bar") {
        p <- p + stat_summary(fun = fun, geom = "bar", aes_string(group = x_group), fill = "grey84", color = "black", position = position_dodge(), alpha = 1)
    }
    
    # Define the position adjustment for jitter (unused if geom_type is "point")
    position_adjustment <- position_quasirandom(width = 0.2)
    
    if (line == TRUE) {
        p <- p + geom_line(aes_string(group = color_by), linewidth = 0.4, alpha = 0.6)
    }
    
    # Conditional logic for geom type
    if (geom_type == "quasirandom") {
        if (!is.null(shape_by)) {
            p <- p + geom_quasirandom(aes_string(group = x_group, color = color_by, shape = shape_by), width = 0.4, size = point_size, stroke = 0.5, show.legend = TRUE, alpha = point_alpha)
        } else {
            p <- p + geom_quasirandom(aes_string(group = x_group, fill = color_by), width = 0.4, size = point_size, shape = 21, color = "black", stroke = 0.5, show.legend = TRUE, alpha = point_alpha)
        }
    } else if (geom_type == "point") {
        if (!is.null(shape_by)) {
            p <- p + geom_point(aes_string(group = x_group, color = color_by, shape = shape_by), size = point_size, stroke = 0.5, show.legend = TRUE, alpha = point_alpha)
        } else {
            p <- p + geom_point(aes_string(group = x_group, fill = color_by), size = point_size, shape = 21, color = "black", stroke = 0.5, show.legend = TRUE, alpha = point_alpha)
        }
    }
    
    if (hide_x_labels) {
        p <- p + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
    }
    
    # Add stats
    if (show_stats) {
        p <- p + stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.025, hide.ns = hide_ns, size = stat_size, bracket.nudge.y = nudge)
    }
    
    # Construct the facet formula
    facet_formula <- if (is.null(facet_by)) {
        as.formula("~ antigen")
    } else {
        as.formula(paste(facet_by, "~ antigen"))
    }
    
    p <- p + facet_grid2(facet_formula, scales = scales, axes = axes, independent = "y")
    
    return(p)
}

environment(plotPbExprsDiffSimp) <- asNamespace("CATALYST")


## ----------------------------------------------------------------------------------------------------------
library(shiny)
library(ggplot2)
library(reshape2)
library(ggridges)
library(viridis)
library(CATALYST)

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param sce PARAM_DESCRIPTION
#' @param donor_list PARAM_DESCRIPTION
#' @param subsets_column PARAM_DESCRIPTION, Default: NULL
#' @param condition_column PARAM_DESCRIPTION, Default: NULL
#' @param num_peaks PARAM_DESCRIPTION, Default: 5
#' @param peak_boundary_uncertainty PARAM_DESCRIPTION, Default: 0.2
#' @param moi PARAM_DESCRIPTION, Default: 'CFSE'
#' @param moi2 PARAM_DESCRIPTION, Default: 'Ki67'
#' @param first_valley_init PARAM_DESCRIPTION, Default: 5.8
#' @param peak_width_init PARAM_DESCRIPTION, Default: 0.62
#' @param use_previous_values PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @export 
#' @importFrom reshape2 melt
#' @importFrom scales censor
division_profiler <- function(
  sce, 
  donor_list, 
  subsets_column = NULL, 
  condition_column = NULL, 
  num_peaks = 5,
  peak_boundary_uncertainty = 0.2, 
  moi = "CFSE", 
  moi2 = "Ki67",
  first_valley_init = 5.8, 
  peak_width_init = 0.62,
  use_previous_values = FALSE
) {
  
  # If using previous values, extract them from metadata
  if (use_previous_values && !is.null(metadata(sce)$first_valley_values) && !is.null(metadata(sce)$peak_width_values)) {
    first_valley_vals <- metadata(sce)$first_valley_values
    peak_width_vals <- metadata(sce)$peak_width_values
  } else {
    first_valley_vals <- rep(first_valley_init, length(donor_list))
    peak_width_vals <- rep(peak_width_init, length(donor_list))
  }
  
  # Build 'division_peaks_graphs_list' for the entire donor data 
  division_peaks_graphs_list <- list()
  for (donor in donor_list) {
    sce_temp1 <- filterSCE(sce, patient_id == donor)
    vec_expr <- as.data.frame(t(sce_temp1@assays@data$exprs))[[moi]]
    tmp_melt <- reshape2::melt(vec_expr)
    tmp_melt$culture_condition <- 1
    
    p_base <- ggplot(tmp_melt, aes(x = value)) +
      geom_density(fill = "#00AFBB", alpha = 0.7) +
      xlab(paste(moi, " expression")) +
      ylab("Density") +
      theme_bw(base_size = 14) +
      ggtitle(paste(moi, donor, sep = "_"))
    
    division_peaks_graphs_list[[donor]] <- p_base
  }
  
  # Build dropdown options for subsets and conditions
  subset_options_list <- list()
  condition_options_list <- list()
  
  for (donor in donor_list) {
    sce_subset <- filterSCE(sce, patient_id == donor)
    
    # Subsets
    if (!is.null(subsets_column)) {
      unique_subsets <- sort(as.character(unique(sce_subset[[subsets_column]])))
      subset_options_list[[donor]] <- c("ALL", unique_subsets)
    } else {
      subset_options_list[[donor]] <- "ALL"
    }
    
    # Conditions
    if (!is.null(condition_column)) {
      unique_conds <- sort(as.character(unique(sce_subset[[condition_column]])))
      condition_options_list[[donor]] <- c("ALL", unique_conds)
    } else {
      condition_options_list[[donor]] <- "ALL"
    }
  }
  
  # Define UI
  ui <- fluidPage(
    tags$style(HTML("
      .shiny-output-error { visibility: hidden; }
      .shiny-output-error:before { visibility: hidden; }
      .container-fluid { padding: 0; max-width: 2000px; }  /* Increased max-width */
      .row { margin: 0; }
      .col-sm-6 { padding: 5px; }  /* Wider columns */
      .plot-container { height: 300px; overflow: hidden; }  /* Shorter height */
      .plot-container .plot { height: 100%; width: 100%; }  /* Full width */
    ")),
    titlePanel("Division Profiler"),
    
    fluidRow(
      lapply(seq_along(donor_list), function(i) {
        donor <- donor_list[i]
        column(
          width = 6,  # Wider columns
          selectInput(
            inputId = paste0("subset_select_", i),
            label = paste0("Select subset for ", donor, ":"),
            choices = subset_options_list[[donor]],
            selected = "ALL"
          ),
          selectInput(
            inputId = paste0("condition_select_", i),
            label = paste0("Select condition for ", donor, ":"),
            choices = condition_options_list[[donor]],
            selected = "ALL"
          ),
          numericInput(
            inputId = paste0("first_valley_", i),
            label = paste0("First Valley (Div0-Div1) for ", donor, ":"),
            value = first_valley_vals[i],
            step = 0.05
          ),
          numericInput(
            inputId = paste0("peak_width_", i),
            label = paste0("Peak Width for ", donor, ":"),
            value = peak_width_vals[i],
            step = 0.01
          ),
          div(class = "plot-container",
              plotOutput(outputId = paste0("histogram_", i), height = "300px", width = "100%")  # Histogram
          ),
          div(class = "plot-container",
              plotOutput(outputId = paste0("scatter_", i), height = "300px", width = "100%")  # Scatter plot
          )
        )
      })
    ),
    
    fluidRow(
      column(
        width = 12,
        actionButton("done", "Done")
      )
    )
  )
  
  # Define server logic
  server <- function(input, output, session) {
    
    lapply(seq_along(donor_list), function(i) {
      donor <- donor_list[i]
      
      observeEvent({
        input[[paste0("subset_select_", i)]]
        input[[paste0("condition_select_", i)]]
        input[[paste0("first_valley_", i)]]
        input[[paste0("peak_width_", i)]]
      }, {
        sce_donor <- filterSCE(sce, patient_id == donor)
        
        sub_choice <- input[[paste0("subset_select_", i)]]
        cond_choice <- input[[paste0("condition_select_", i)]]
        
        if (!is.null(subsets_column) && sub_choice != "ALL") {
          sce_donor <- filterSCE(sce_donor, sce_donor[[subsets_column]] == sub_choice)
        }
        if (!is.null(condition_column) && cond_choice != "ALL") {
          sce_donor <- filterSCE(sce_donor, sce_donor[[condition_column]] == cond_choice)
        }
        
        # Check if there are any cells after filtering
        if (ncol(sce_donor) == 0) {
          # Display a message if no data is available
          output[[paste0("histogram_", i)]] <- renderPlot({
            ggplot() +
              annotate("text", x = 4, y = 0.5, label = "No Data for this combination", size = 6) +
              theme_void()
          }, height = 300, width = 600)
          output[[paste0("scatter_", i)]] <- renderPlot({
            ggplot() +
              annotate("text", x = 4, y = 0.5, label = "No Data for this combination", size = 6) +
              theme_void()
          }, height = 300, width = 600)
          return()
        }
        
        expr_vals <- as.vector(sce_donor@assays@data$exprs[moi, ])
        expr_vals_moi2 <- as.vector(sce_donor@assays@data$exprs[moi2, ])
        
        base_plot <- division_peaks_graphs_list[[donor]]
        plot_data <- ggplot_build(base_plot)$data[[1]]
        
        first_valley_val <- input[[paste0("first_valley_", i)]]
        peak_width_val <- input[[paste0("peak_width_", i)]]
        
        pd1 <- subset(plot_data, 
                      x > first_valley_val - peak_boundary_uncertainty & 
                        x < first_valley_val + peak_boundary_uncertainty)
        intersect1 <- if (nrow(pd1) > 0) pd1$x[which.min(pd1$density)] else NA
        
        pd2 <- subset(plot_data, 
                      x > first_valley_val - (peak_width_val + peak_boundary_uncertainty) &
                        x < first_valley_val - (peak_width_val - peak_boundary_uncertainty))
        intersect2 <- if (nrow(pd2) > 0) pd2$x[which.min(pd2$density)] else NA
        
        pd3 <- subset(plot_data, 
                      x > first_valley_val - (2*peak_width_val + peak_boundary_uncertainty) &
                        x < first_valley_val - (2*peak_width_val - peak_boundary_uncertainty))
        intersect3 <- if (nrow(pd3) > 0) pd3$x[which.min(pd3$density)] else NA
        
        ave_ints <- ((intersect1 - intersect2) + (intersect2 - intersect3)) / 2
        line_intervals <- intersect1 + ave_ints - (0:num_peaks)*ave_ints
        line_sorted <- sort(line_intervals)
        
        idx <- findInterval(expr_vals, line_sorted, rightmost.closed = TRUE)
        div_level <- length(line_sorted) - idx
        div_factor <- factor(paste0("Div", div_level))
        
        plot_df <- data.frame(
          value = expr_vals,
          division = div_factor,
          moi2 = expr_vals_moi2
        )
        
        line_blue  <- line_sorted[length(line_sorted)]
        lines_rest <- line_sorted[-length(line_sorted)]
        
        # Histogram plot
        p_hist <- ggplot(plot_df, aes(x = value, fill = division)) +
          geom_histogram(bins = 200, alpha = 1, position = "stack") +
          scale_fill_brewer(palette = "Paired", name = "Division") +
          scale_x_continuous(limits = c(0, 8), oob = scales::censor) +
          geom_vline(xintercept = lines_rest, color = "black", linetype = "dashed") +
          geom_vline(xintercept = line_blue,  color = "blue",  linetype = "dashed") +
          theme_bw(base_size = 14) +
          labs(
            x = paste(moi, " expression"),
            y = "Count",
            title = paste0(
              moi, "_", donor,
              if (sub_choice != "ALL")  paste0("_", sub_choice)  else "",
              if (cond_choice != "ALL") paste0("_", cond_choice) else ""
            )
          )
        
        # Scatter plot with density contours
        p_scatter <- ggplot(plot_df, aes(x = value, y = moi2, color = division)) +
          geom_point(alpha = 0.7, size = 1) +
          #stat_density_2d(aes(fill = ..level..), geom = "polygon", alpha = 0.1, linewidth = 0.4 ,color = "black", linetype = "solid") +
          stat_density_2d(aes(fill = ..level..), geom = "polygon", alpha = 0.01, linewidth = 0.4 ,color = "black", linetype = "solid") +
          scale_color_brewer(palette = "Paired", name = "Division") +
          scale_x_continuous(limits = c(0, 8), oob = scales::censor) +
          theme_bw(base_size = 14) +
          labs(
            x = paste(moi, " expression"),
            y = paste(moi2, " expression"),
            title = paste0(
              moi, "_", donor,
              if (sub_choice != "ALL")  paste0("_", sub_choice)  else "",
              if (cond_choice != "ALL") paste0("_", cond_choice) else ""
            )
          )
        
        output[[paste0("histogram_", i)]] <- renderPlot({ p_hist }, height = 300, width = 600)  # Adjusted dimensions
        output[[paste0("scatter_", i)]] <- renderPlot({ p_scatter }, height = 300, width = 600)  # Adjusted dimensions
      }, ignoreNULL = FALSE)
    })
    
    observeEvent(input$done, {
      first_valley_vals <- sapply(seq_along(donor_list), function(i) input[[paste0("first_valley_", i)]])
      peak_width_vals <- sapply(seq_along(donor_list), function(i) input[[paste0("peak_width_", i)]])
      
      final_values <- data.frame(
        donor                = donor_list,
        first_valley         = first_valley_vals,
        peak_width           = peak_width_vals
      )
      
      # Assign division labels to the entire SCE
      for (i in seq_along(donor_list)) {
        donor <- donor_list[i]
        first_valley    <- first_valley_vals[i]
        peak_width <- peak_width_vals[i]
        
        sce_temp <- sdp(
          sce,
          donor_list = donor_list,
          division_peaks_graphs_list = division_peaks_graphs_list,
          first_valley = first_valley,
          donor = i,
          num_peaks = num_peaks,
          peak_boundary_uncertainty = peak_boundary_uncertainty,
          peak_width = peak_width,
          subsets_column = subsets_column,
          condition_column = condition_column,
          moi = moi,
          moi2 = moi2
        )
        
        # Merge division labels back into the main SCE
        sce$div[sce$patient_id == donor] <- sce_temp$div
      }
      
      # Store final values in metadata
      metadata(sce)$first_valley_values <- first_valley_vals
      metadata(sce)$peak_width_values <- peak_width_vals
      
      stopApp(sce)
    })
  }
  
  # Run the Shiny app
  sce <- runApp(shinyApp(ui = ui, server = server))
  return(sce)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param sce PARAM_DESCRIPTION
#' @param donor_list PARAM_DESCRIPTION
#' @param division_peaks_graphs_list PARAM_DESCRIPTION
#' @param first_valley PARAM_DESCRIPTION
#' @param donor PARAM_DESCRIPTION
#' @param num_peaks PARAM_DESCRIPTION, Default: 5
#' @param peak_boundary_uncertainty PARAM_DESCRIPTION, Default: 0.2
#' @param peak_width PARAM_DESCRIPTION, Default: 1
#' @param subsets_column PARAM_DESCRIPTION, Default: NULL
#' @param condition_column PARAM_DESCRIPTION, Default: NULL
#' @param moi PARAM_DESCRIPTION, Default: 'CFSE'
#' @param moi2 PARAM_DESCRIPTION, Default: 'Ki67'
#' @return OUTPUT_DESCRIPTION
#' @export 
sdp <- function(sce, donor_list, division_peaks_graphs_list, first_valley, donor, 
                num_peaks = 5, peak_boundary_uncertainty = 0.2, 
                peak_width = 1, subsets_column = NULL, 
                condition_column = NULL, 
                moi = "CFSE", moi2 = "Ki67") {
  
  donor <- donor_list[donor]
  sce_temp1 <- filterSCE(sce, patient_id == donor)
  sce_temp1$condition2 <- sce_temp1[[condition_column]] # Condition2 so it doesn't overwrite
  
  plot.data <- ggplot_build(division_peaks_graphs_list[[donor]])
  plot.data <- as.data.frame(plot.data$data)
  
  plot.data2 <- plot.data[plot.data$x > first_valley - peak_boundary_uncertainty & 
                            plot.data$x < first_valley + peak_boundary_uncertainty, ]
  intersect1 <- plot.data2$x[plot.data2$density == min(plot.data2$density)]
  
  plot.data2 <- plot.data[plot.data$x > first_valley - (peak_width + peak_boundary_uncertainty) & 
                            plot.data$x < first_valley - (peak_width - peak_boundary_uncertainty), ]
  intersect2 <- plot.data2$x[plot.data2$density == min(plot.data2$density)]
  
  plot.data2 <- plot.data[plot.data$x > first_valley - (2 * peak_width + peak_boundary_uncertainty) & 
                            plot.data$x < first_valley - (2 * peak_width - peak_boundary_uncertainty), ]
  intersect3 <- plot.data2$x[plot.data2$density == min(plot.data2$density)]
  
  ave_ints <- ((intersect1 - intersect2) + (intersect2 - intersect3)) / 2
  differences <- intersect1 + ave_ints - (0:(num_peaks - 1)) * ave_ints
  differences[1] <- differences[1] + 2
  sorted_differences <- sort(differences)
  
  div_levels <- findInterval(sce_temp1@assays@data$exprs[moi, ], sorted_differences, rightmost.closed = TRUE)
  max_div_level <- length(sorted_differences)
  converted_div_levels <- max_div_level - div_levels
  adjusted_div_levels <- converted_div_levels - 1
  sce_temp1$div <- paste0("Div", adjusted_div_levels)
  
  return(sce_temp1)
}



## ----------------------------------------------------------------------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param div_col PARAM_DESCRIPTION, Default: 'div'
#' @param parse_div PARAM_DESCRIPTION, Default: TRUE
#' @param metric PARAM_DESCRIPTION, Default: c("mean_div", "median_div", "percent_proliferated")
#' @param scales PARAM_DESCRIPTION, Default: 'free_y'
#' @param conditions PARAM_DESCRIPTION, Default: NULL
#' @param facet_by PARAM_DESCRIPTION, Default: NULL
#' @param facet PARAM_DESCRIPTION, Default: TRUE
#' @param excluded_clusters PARAM_DESCRIPTION, Default: NULL
#' @param axes PARAM_DESCRIPTION, Default: 'all'
#' @param line PARAM_DESCRIPTION, Default: FALSE
#' @param fun PARAM_DESCRIPTION, Default: c("median", "mean", "sum")
#' @param point_size PARAM_DESCRIPTION, Default: 1
#' @param textsize PARAM_DESCRIPTION, Default: 14
#' @param panel_spacing PARAM_DESCRIPTION, Default: 2
#' @param show_stats PARAM_DESCRIPTION, Default: TRUE
#' @param hide_ns PARAM_DESCRIPTION, Default: TRUE
#' @param fill_by PARAM_DESCRIPTION, Default: NULL
#' @param x_group PARAM_DESCRIPTION, Default: 'condition'
#' @param group_by PARAM_DESCRIPTION, Default: x_group
#' @param shape_by PARAM_DESCRIPTION, Default: NULL
#' @param stat_size PARAM_DESCRIPTION, Default: 4
#' @param facet_ratio PARAM_DESCRIPTION, Default: 1.5
#' @param geom PARAM_DESCRIPTION, Default: c("boxes", "bar")
#' @param jitter PARAM_DESCRIPTION, Default: TRUE
#' @param geom_type PARAM_DESCRIPTION, Default: c("quasirandom", "point")
#' @param group1keep PARAM_DESCRIPTION, Default: NULL
#' @param nudge PARAM_DESCRIPTION, Default: 0
#' @param step_increase PARAM_DESCRIPTION, Default: 0.1
#' @param hide_x_labels PARAM_DESCRIPTION, Default: FALSE
#' @param label_parse PARAM_DESCRIPTION, Default: FALSE
#' @param merging_col PARAM_DESCRIPTION, Default: NULL
#' @param swap PARAM_DESCRIPTION, Default: FALSE
#' @param external_stats PARAM_DESCRIPTION, Default: NULL
#' @param average_samples PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @export 
plotDivSimple <- function (
    x, 
    div_col = "div",  # Name of the Div column
    parse_div = TRUE,  # Whether to parse Div column (e.g., Div0 → 0)
    metric = c("mean_div", "median_div", "percent_proliferated"),  # Updated metric options
    scales = "free_y", 
    conditions = NULL, 
    facet_by = NULL,
    facet = TRUE,
    excluded_clusters = NULL, 
    axes = "all", 
    line = FALSE,
    fun = c("median", "mean", "sum"), 
    point_size = 1, 
    textsize = 14, 
    panel_spacing = 2, 
    show_stats = TRUE, 
    hide_ns = TRUE, 
    fill_by = NULL,
    x_group = "condition", 
    group_by = x_group, 
    shape_by = NULL, 
    stat_size = 4, 
    facet_ratio = 1.5, 
    geom = c("boxes","bar"), 
    jitter = TRUE,
    geom_type = c("quasirandom", "point"),  # New parameter
    group1keep = NULL, 
    nudge = 0, 
    step_increase = 0.1, 
    hide_x_labels = FALSE, 
    label_parse = FALSE, 
    merging_col = NULL, 
    swap = FALSE,
    external_stats = NULL,  # New parameter for external stats dataframe
    average_samples = FALSE  # New parameter to average across samples
) {
    geom <- match.arg(geom)
    geom_type <- match.arg(geom_type)  # Match the new geom_type argument
    metric <- match.arg(metric)  # Match the metric argument
    stopifnot(is.logical(jitter), length(jitter) == 1)
    
    # Retrieve the div column
    if (!div_col %in% colnames(colData(x))) {
        stop("The specified 'div_col' (", div_col, ") does not exist in the colData of the input object.")
    }
    
    # Extract the div column and convert to numeric if parse_div is TRUE
    div_values <- colData(x)[[div_col]]
    if (parse_div) {
        div_values <- as.numeric(gsub("Div", "", div_values))  # Extract numerical part (e.g., Div0 → 0)
    }
    
    # Create the data frame
    df <- data.frame(div = div_values, colData(x))
    
    # Debugging: Check for non-finite values
    if (any(is.na(df$div)) || any(is.nan(df$div)) || any(is.infinite(df$div))) {
        warning("Non-finite values found in 'div' column. These will be removed.")
    }
    
    if (!is.null(conditions)) {
        df <- df[df[[x_group]] %in% conditions, ]
        df[[x_group]] <- factor(df[[x_group]], levels = conditions, ordered = TRUE)
    }
    
    # Apply excluded clusters if needed
    if (!is.null(excluded_clusters)) {
        df <- df[!df$cluster_id %in% excluded_clusters, ]
    }
    
    # Calculate the metric to plot
    if (metric == "percent_proliferated") {
        # Calculate percentage proliferated (cells not in division 0)
        df_agg <- df %>%
            group_by(sample_id, !!sym(x_group), !!sym(facet_by), patient_id) %>%
            summarise(
                div = sum(div != 0) / n() * 100,  # Percentage of cells not in division 0
                .groups = 'drop'
            )
    } else {
        # Determine the aggregation function based on the metric
        agg_fun <- ifelse(metric == "mean_div", mean, median)
        
        # Aggregate div values by sample_id and condition
        df_agg <- df %>%
            group_by(sample_id, !!sym(x_group), !!sym(facet_by), patient_id) %>%
            summarise(div = agg_fun(div, na.rm = TRUE), .groups = 'drop')
    }
    
    # Debugging: Check the aggregated data
    print(head(df_agg))
    
    # Average across samples if average_samples is TRUE
    if (average_samples) {
        df_agg <- df_agg %>%
            group_by(across(all_of(c(facet_by, x_group)))) %>%
            summarise(div = mean(div, na.rm = TRUE), .groups = 'drop')
    }
    
    # Check if the data is empty after filtering
    if (nrow(df_agg) == 0) {
        stop("No data available after filtering. Please adjust the parameters.")
    }
    
    # Rest of the function remains the same...
    if (show_stats) {
        if (!is.null(external_stats)) {
            # Filter external stats to include only the antigens that are plotted
            external_stats <- external_stats %>% filter(antigen %in% features)
            
            # Use external stats if provided
            dummy_stat <- df_agg %>% group_by(!!sym(x_group)) %>% wilcox_test(as.formula(paste("div ~", x_group)))
            dummy_stat <- dummy_stat %>% add_y_position(scales = "free", step.increase = step_increase)
            
            # Merge y.position from dummy_stat to external_stats
            external_stats <- external_stats %>%
                left_join(dummy_stat %>% select(!!sym(x_group), group1, group2, y.position), 
                          by = c(x_group, "group1", "group2"))
            
            stat.test <- external_stats
        } else {
            # Calculate stats internally using wilcox_test (change to Dunn?)
            stat.test <- df_agg %>%
                group_by(!!sym(x_group)) %>%
                wilcox_test(as.formula(paste("div ~", x_group))) %>%
                add_y_position(scales = "free", step.increase = step_increase)
        }
        
        if (!is.null(group1keep)) {
            stat.test <- stat.test[stat.test$group1 %in% group1keep,]
            stat.test <- stat.test %>% add_y_position(scales = "free", step.increase = step_increase)
        }
        
        # Save unfiltered stats globally
        pbCondStats <<- stat.test
    }

    dfout <<- df_agg
    
    # Initialize the ggplot
    p <- ggplot(df_agg, aes_string(x = group_by, y = "div", fill = fill_by)) +
        labs(y = ifelse(metric == "percent_proliferated", "Percentage Proliferated", 
                        ifelse(metric == "mean_div", "Mean Division", "Median Division"))) +
        theme_bw() +
        theme(
            panel.grid = element_blank(),
            text = element_text(size = textsize),
            strip.text = element_text(size = textsize),
            strip.background = element_rect(fill = NA, color = NA),
            strip.placement = "outside",  # Move facet labels outside
            legend.text = element_text(size = textsize),
            legend.title = element_text(size = textsize),
            aspect.ratio = facet_ratio,
            axis.ticks = element_line(color = "black"),
            panel.border = element_rect(color = "black", fill = NA, size = 0.5),
            panel.spacing = unit(panel_spacing, "lines"),
            axis.text = element_text(color = "black", size = textsize),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = textsize)
        )
    
    # Add geom based on selection
    if (geom == "boxes") {
        p <- p + geom_boxplot(aes_string(group = x_group), fill = "grey84", color = "black", width = 0.75, linewidth = 0.3, alpha = 0.9, outlier.color = NA, show.legend = FALSE)
    } else if (geom == "bar") {
        p <- p + stat_summary(fun = match.fun(fun), geom = "bar", aes_string(group = x_group), fill = "grey84", color = "black", position = position_dodge(), alpha = 1)
    }
    
    # Define the position adjustment for jitter
    position_adjustment <- position_quasirandom(width = 0.2)
    
    if (line == TRUE) {
        p <- p + geom_line(aes_string(group = fill_by), width = 0.4, alpha = 0.4)
    }
    
    # Conditional logic for geom type
    if (geom_type == "quasirandom") {
        if (!is.null(shape_by)) {
            p <- p + geom_quasirandom(aes_string(group = x_group, fill = fill_by, shape = shape_by), width = 0.4, size = point_size, color = "black", stroke = 0.5, show.legend = TRUE, alpha = 1)
        } else {
            p <- p + geom_quasirandom(aes_string(group = x_group, fill = fill_by), width = 0.4, shape = 21, size = point_size, color = "black", stroke = 0.5, show.legend = TRUE, alpha = 1)
        }
    } else if (geom_type == "point") {
        if (!is.null(shape_by)) {
            p <- p + geom_point(aes_string(group = x_group, fill = fill_by, shape = shape_by), size = point_size, color = "black", stroke = 0.5, show.legend = TRUE, alpha = 1)
        } else {
            p <- p + geom_point(aes_string(group = x_group, fill = fill_by), shape = 21, size = point_size, color = "black", stroke = 0.5, show.legend = TRUE, alpha = 1)
        }
    }
    
    if (hide_x_labels) {
        p <- p + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
    }
    
    # Add stats
    if (show_stats) {
        p <- p + stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.025, hide.ns = hide_ns, size = stat_size, bracket.nudge.y = nudge)
    }
    
    if (facet) {
        # Construct the facet formula to include only facet_by
        facet_formula <- as.formula(paste("~", facet_by))
        p <- p + facet_grid(facet_formula, scales = scales)
    }
    
    # Fix y-axis limits for percent_proliferated
    if (metric == "percent_proliferated") {
        p <- p + scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20))
    }
    
    return(p)
}


## ----------------------------------------------------------------------------------------------------------
library(SingleCellExperiment)
library(ggplot2)
library(ggbeeswarm)
library(scales)
library(dplyr)
library(rstatix)
library(rlang)
library(ggh4x)
library(ggrastr)

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param k PARAM_DESCRIPTION, Default: 'meta20'
#' @param scales PARAM_DESCRIPTION, Default: 'free_y'
#' @param conditions PARAM_DESCRIPTION, Default: NULL
#' @param color_by PARAM_DESCRIPTION, Default: NULL
#' @param excluded_clusters PARAM_DESCRIPTION, Default: NULL
#' @param axes PARAM_DESCRIPTION, Default: 'all'
#' @param fun PARAM_DESCRIPTION, Default: c("median", "mean", "sum")
#' @param point_size PARAM_DESCRIPTION, Default: 1
#' @param clusters_order PARAM_DESCRIPTION, Default: NULL
#' @param textsize PARAM_DESCRIPTION, Default: 14
#' @param panel_spacing PARAM_DESCRIPTION, Default: 2
#' @param show_stats PARAM_DESCRIPTION, Default: TRUE
#' @param hide_ns PARAM_DESCRIPTION, Default: TRUE
#' @param facet_by PARAM_DESCRIPTION, Default: 'cluster_id'
#' @param x_group PARAM_DESCRIPTION, Default: 'condition'
#' @param group_by PARAM_DESCRIPTION, Default: x_group
#' @param shape_by PARAM_DESCRIPTION, Default: NULL
#' @param stat_size PARAM_DESCRIPTION, Default: 4
#' @param facet_ratio PARAM_DESCRIPTION, Default: 1.5
#' @param geom PARAM_DESCRIPTION, Default: c("boxes", "bar")
#' @param jitter PARAM_DESCRIPTION, Default: TRUE
#' @param group1keep PARAM_DESCRIPTION, Default: NULL
#' @param nudge PARAM_DESCRIPTION, Default: 0
#' @param step_increase PARAM_DESCRIPTION, Default: 0.1
#' @param hide_x_labels PARAM_DESCRIPTION, Default: FALSE
#' @param label_parse PARAM_DESCRIPTION, Default: FALSE
#' @param merging_col PARAM_DESCRIPTION, Default: NULL
#' @param swap PARAM_DESCRIPTION, Default: FALSE
#' @param external_stats PARAM_DESCRIPTION, Default: NULL
#' @param average_samples PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @export 
plotDivDiff <- function (
    x, 
    k = "meta20", 
    scales = "free_y", 
    conditions = NULL, 
    color_by = NULL,
    excluded_clusters = NULL, 
    axes = "all", 
    fun = c("median", "mean", "sum"), 
    point_size = 1, 
    clusters_order = NULL, 
    textsize = 14, 
    panel_spacing = 2, 
    show_stats = TRUE, 
    hide_ns = TRUE, 
    facet_by = "cluster_id", 
    x_group = "condition", 
    group_by = x_group, 
    shape_by = NULL, 
    stat_size = 4, 
    facet_ratio = 1.5, 
    geom = c("boxes","bar"), 
    jitter = TRUE, 
    group1keep = NULL, 
    nudge = 0, 
    step_increase = 0.1, 
    hide_x_labels = FALSE, 
    label_parse = FALSE, 
    merging_col = NULL, 
    swap = FALSE,
    external_stats = NULL,  # New parameter for external stats dataframe
    average_samples = FALSE  # New parameter to average across samples
) {
    fun <- match.arg(fun)
    geom <- match.arg(geom)
    stopifnot(is.logical(jitter), length(jitter) == 1)
    
    # Use the merging column from colData if provided
    if (!is.null(merging_col)) {
        cluster_ids <- factor(x[[k]])
    } else {
        # If no merging_col, proceed with normal clustering logic
        .check_sce(x)
        k <- .check_k(x, k)
        cluster_ids <- cluster_ids(x, k)
    }
    
    # Retrieve the div column and convert to numeric
    df <- data.frame(div = as.numeric(gsub("div", "", x$div)), colData(x))
    df$cluster_id <- cluster_ids
    
    if (!is.null(conditions)) {
        df <- df[df[[x_group]] %in% conditions, ]
        df[[x_group]] <- factor(df[[x_group]], levels = conditions, ordered = TRUE)
    }
    
    # Apply excluded clusters if needed
    if (!is.null(excluded_clusters)) {
        df <- df[!df$cluster_id %in% excluded_clusters, ]
    }
    
    # Aggregate div values by cluster_id and sample_id
    df_agg <- df %>%
        group_by(cluster_id, sample_id, !!sym(x_group), !!sym(color_by), patient_id) %>%
        summarise(div = match.fun(fun)(div, na.rm = TRUE), .groups = 'drop')
    
    # Average across samples if average_samples is TRUE
    if (average_samples) {
        df_agg <- df_agg %>%
            group_by(across(all_of(c(color_by, x_group, "cluster_id")))) %>%
            summarise(div = mean(div, na.rm = TRUE), .groups = 'drop')
    }
    
    # Check if the data is empty after filtering
    if (nrow(df_agg) == 0) {
        stop("No data available after filtering. Please adjust the parameters.")
    }
    
    if (show_stats) {
        if (!is.null(external_stats)) {
            # Filter external stats to include only the antigens that are plotted
            external_stats <- external_stats %>% filter(antigen %in% features)
            
            # Use external stats if provided
            dummy_stat <- df_agg %>% group_by(!!sym(facet_by)) %>% wilcox_test(as.formula(paste("div ~", x_group)))
            dummy_stat <- dummy_stat %>% add_y_position(scales = "free", step.increase = step_increase)
            
            # Merge y.position from dummy_stat to external_stats
            external_stats <- external_stats %>%
                left_join(dummy_stat %>% select(!!sym(facet_by), group1, group2, y.position), 
                          by = c(facet_by, "group1", "group2"))
            
            stat.test <- external_stats
        } else {
            # Calculate stats internally using wilcox_test (change to Dunn?)
            stat.test <- df_agg %>%
                group_by(!!sym(facet_by)) %>%
                wilcox_test(as.formula(paste("div ~", x_group))) %>%
                add_y_position(scales = "free", step.increase = step_increase)
        }
        
        if (!is.null(group1keep)) {
            stat.test <- stat.test[stat.test$group1 %in% group1keep,]
            stat.test <- stat.test %>% add_y_position(scales = "free", step.increase = step_increase)
        }
        
        # Save unfiltered stats globally
        pbCondStats <<- stat.test
    }
    
    dfout <<- df_agg
    
    # Initialize the ggplot
    p <- ggplot(df_agg, aes_string(x = group_by, y = "div", fill = color_by)) +
        labs(y = paste(fun, "div")) + # Plots "median div"
        theme_bw() +
        theme(
            panel.grid = element_blank(),
            text = element_text(size = textsize),
            strip.text = element_text(size = textsize),
            strip.background = element_rect(fill = NA, color = NA),
            strip.placement = "outside",  # Move facet labels outside
            legend.text = element_text(size = textsize),
            legend.title = element_text(size = textsize),
            aspect.ratio = facet_ratio,
            axis.ticks = element_line(color = "black"),
            panel.border = element_rect(color = "black", fill = NA, size = 0.5),
            panel.spacing = unit(panel_spacing, "lines"),
            axis.text = element_text(color = "black", size = textsize),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = textsize)
        )
    
    # Add geom based on selection
    if (geom == "boxes") {
        p <- p + geom_boxplot(aes_string(group = x_group), fill = "grey84", color = "black", width = 0.75, linewidth = 0.3, alpha = 0.9, outlier.color = NA, show.legend = F)
    } else if (geom == "bar") {
        p <- p + stat_summary(fun = match.fun(fun), geom = "bar", aes_string(group = x_group), fill = "grey84", color = "black", position = position_dodge(), alpha = 1)
    }
    
    # Define the position adjustment for jitter
    position_adjustment <- position_quasirandom(width = 0.2)
    
    p <- p + geom_quasirandom(aes_string(group = x_group, fill = color_by, shape = shape_by), width = 0.4, size = point_size, color = "black", stroke = 0.5, show.legend = TRUE, alpha = 1)
    
    if (hide_x_labels) {
        p <- p + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
    }
    
    # Add stats
    if (show_stats) {
        p <- p + stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.025, hide.ns = hide_ns, size = stat_size, bracket.nudge.y = nudge)
    }
    
    # Construct the facet formula to include both facet_by and color_by
    facet_formula <- as.formula(paste(facet_by, "~", color_by))
    
    if (swap == TRUE) {facet_formula <- as.formula(paste(color_by, "~", facet_by)) }
    
    p <- p + facet_grid(facet_formula, scales = scales)
    
    return(p)
}

environment(plotDivDiff) <- asNamespace("CATALYST")


## ----------------------------------------------------------------------------------------------------------
library(cowplot)
library(grid)
library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param sce PARAM_DESCRIPTION
#' @param master_cluster_column PARAM_DESCRIPTION
#' @param daughter_cluster_column PARAM_DESCRIPTION
#' @param meta PARAM_DESCRIPTION
#' @param sample_id_column PARAM_DESCRIPTION
#' @param cols PARAM_DESCRIPTION, Default: c(brewer.pal(9, "Paired"), "#9b59b6")
#' @param daughter_order PARAM_DESCRIPTION, Default: NULL
#' @param row_order PARAM_DESCRIPTION, Default: NULL
#' @param threshold PARAM_DESCRIPTION, Default: 10
#' @return OUTPUT_DESCRIPTION
#' @export 
plotDonutProportionsN <- function(sce, master_cluster_column, daughter_cluster_column, meta, sample_id_column, cols =  c(brewer.pal(9, "Paired"), "#9b59b6"), daughter_order = NULL, row_order = NULL, threshold = 10) {
  # Convert to data frame
  df <- as.data.frame(colData(sce))

  df[[master_cluster_column]] <- factor(df[[master_cluster_column]], levels = unique(df[[master_cluster_column]]))
  df[[daughter_cluster_column]] <- factor(df[[daughter_cluster_column]], levels = unique(df[[daughter_cluster_column]])) # Ensure consistent factor levels
  df[[sample_id_column]] <- factor(df[[sample_id_column]], levels = unique(df[[sample_id_column]]))

  # Convert the daughter columns to factors with optional custom ordering
  if (!is.null(daughter_order)) {
    df[[daughter_cluster_column]] <- factor(df[[daughter_cluster_column]], levels = daughter_order)
  } else {
    df[[daughter_cluster_column]] <- factor(df[[daughter_cluster_column]], levels = unique(df[[daughter_cluster_column]]))
  }

  # Determine color palette based on unique levels of daughter clustering
  color_palette <- brewer.pal(min(length(levels(df[[daughter_cluster_column]])), 12), "Paired")

  # Calculate counts using dplyr to ensure only existing combinations are included
  count_df <- df %>%
    group_by_at(vars(sample_id_column, master_cluster_column, daughter_cluster_column)) %>%
    summarise(count = n(), .groups = "drop")

  # Ensure all combinations are present with zero counts where necessary
  complete_df <- count_df %>%
    complete(nesting(!!sym(sample_id_column), !!sym(master_cluster_column)), !!sym(daughter_cluster_column), fill = list(count = 0))

  # Calculate total counts for each master cluster within each sample
  master_cluster_counts <- complete_df %>%
    group_by_at(vars(sample_id_column, master_cluster_column)) %>%
    summarise(total_count = sum(count), .groups = "drop")

  # Label clusters that meet the threshold
  labeled_master_clusters <- master_cluster_counts %>%
    mutate(include_in_avg = total_count > threshold)

  # Merge the labels back to the complete_df
  labeled_df <- complete_df %>%
    left_join(labeled_master_clusters, by = c(sample_id_column, master_cluster_column))
  
  labout <<- labeled_df

  # Calculate proportions within each sample_id and master_cluster_column
  proportion_df <- labeled_df %>%
    group_by_at(vars(sample_id_column, master_cluster_column)) %>%
    mutate(proportion = count / sum(count) * 100) %>%
    ungroup()

  # Merge the cell counts with the proportion_df
  summary_df <- proportion_df %>%
    rename(cell_count = count)

  # Output the dataframe to the global environment
  summary_out <<- summary_df
  
  # Calculate the average proportions only for clusters that meet the threshold
  averaged_df <- summary_df %>%
    filter(include_in_avg) %>%
    group_by_at(vars(master_cluster_column, daughter_cluster_column)) %>%
    summarise(mean_proportion = mean(proportion), .groups = "drop") %>%
    ungroup()
  
  average_out <<- averaged_df
  
  # Normalize the averaged proportions to ensure they sum to 100 within each master cluster
  normalized_df <- averaged_df %>%
    group_by_at(vars(master_cluster_column)) %>%
    mutate(normalized_proportion = mean_proportion / sum(mean_proportion) * 100) %>%
    ungroup()

  # Update the factor levels of master_cluster_column based on row_order
  if (!is.null(row_order)) {
    normalized_df[[master_cluster_column]] <- factor(normalized_df[[master_cluster_column]], levels = row_order)
  }

  # Create the plot
  p <- ggplot(normalized_df, aes(x = 2, y = normalized_proportion, fill = get(daughter_cluster_column))) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start = 0) +
    xlim(1, 2.5) +
    scale_fill_manual(values = cols, name = "Isotype") +
    theme_void() +
    facet_grid(rows = vars(get(master_cluster_column))) +
    theme(
      strip.text = element_text(size = 20), # Reduce strip text size
      legend.title = element_text(size = 20), # Reduce legend title size
      legend.text = element_text(size = 20), # Reduce legend text size
      panel.spacing = unit(0.1, "lines"), # Reduce space between panels
      plot.margin = margin(10, 10, 10, 10) # Reduce plot margins
    )

  # Return the final plot
  return(p)
}


## ----------------------------------------------------------------------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param sce PARAM_DESCRIPTION
#' @param master_cluster_column PARAM_DESCRIPTION
#' @param daughter_cluster_column PARAM_DESCRIPTION
#' @param sample_id_column PARAM_DESCRIPTION
#' @param threshold_min PARAM_DESCRIPTION, Default: 1
#' @param threshold_max PARAM_DESCRIPTION, Default: 50
#' @param threshold_step PARAM_DESCRIPTION, Default: 1
#' @param daughter_order PARAM_DESCRIPTION, Default: NULL
#' @param window_size PARAM_DESCRIPTION, Default: 5
#' @return OUTPUT_DESCRIPTION
#' @export 
selectThresholdVarianceFloor <- function(sce, master_cluster_column, daughter_cluster_column, sample_id_column,
                                           threshold_min = 1, threshold_max = 50, threshold_step = 1,
                                           daughter_order = NULL, window_size = 5) {
  library(dplyr)
  library(tidyr)
  
  # Extract data from sce and define factors
  df <- as.data.frame(colData(sce))
  df[[master_cluster_column]] <- factor(df[[master_cluster_column]], levels = unique(df[[master_cluster_column]]))
  
  if (!is.null(daughter_order)) {
    df[[daughter_cluster_column]] <- factor(df[[daughter_cluster_column]], levels = daughter_order)
  } else {
    df[[daughter_cluster_column]] <- factor(df[[daughter_cluster_column]], levels = unique(df[[daughter_cluster_column]]))
  }
  df[[sample_id_column]] <- factor(df[[sample_id_column]], levels = unique(df[[sample_id_column]]))
  
  # Count cells per (sample, master, daughter) combination
  count_df <- df %>%
    group_by_at(vars(sample_id_column, master_cluster_column, daughter_cluster_column)) %>%
    summarise(count = n(), .groups = "drop")
  
  # Ensure all combinations are present (missing ones get count=0)
  complete_df <- count_df %>%
    complete(nesting(!!sym(sample_id_column), !!sym(master_cluster_column)),
             !!sym(daughter_cluster_column), fill = list(count = 0))
  
  # Define thresholds to evaluate
  thresh_vals <- seq(threshold_min, threshold_max, by = threshold_step)
  results_list <- list()
  
  # For each threshold, compute global daughter proportions from valid masters
  for (th in thresh_vals) {
    master_counts <- complete_df %>%
      group_by_at(vars(sample_id_column, master_cluster_column)) %>%
      summarise(total_count = sum(count), .groups = "drop")
    
    valid_masters <- master_counts %>% filter(total_count > th)
    
    filtered_df <- complete_df %>%
      inner_join(valid_masters, by = c(sample_id_column, master_cluster_column))
    
    daughter_totals <- filtered_df %>%
      group_by_at(vars(daughter_cluster_column)) %>%
      summarise(sum_count = sum(count), .groups = "drop")
    
    overall_total <- sum(daughter_totals$sum_count)
    if (overall_total == 0) {
      global_props <- daughter_totals %>%
        mutate(perc = 0) %>%
        select(!!sym(daughter_cluster_column), perc)
    } else {
      global_props <- daughter_totals %>%
        mutate(perc = (sum_count/overall_total)*100) %>%
        select(!!sym(daughter_cluster_column), perc)
    }
    
    # Ensure all daughter levels are present.
    daughter_levels <- levels(df[[daughter_cluster_column]])
    global_props <- global_props %>%
      complete(!!sym(daughter_cluster_column) := daughter_levels, fill = list(perc = 0)) %>%
      arrange(!!sym(daughter_cluster_column))
    
    results_list[[as.character(th)]] <- list(threshold = th,
                                             global_props = global_props)
  }
  
  # Build a matrix where each row corresponds to the global proportions at a threshold.
  # Assumes order of daughter clusters is consistent.
  prop_matrix <- do.call(rbind, lapply(results_list, function(x) x$global_props$perc))
  # Rows correspond to thresh_vals
  
  # For each threshold i, compute a sliding window average of the global proportions.
  # Then calculate the difference metric as sum( abs( proportions[i] - window_average ) ).
  diff_metric <- numeric(nrow(prop_matrix))
  for (i in seq_len(nrow(prop_matrix))) {
    window_end <- min(nrow(prop_matrix), i + window_size - 1)
    window_avg <- colMeans(prop_matrix[i:window_end, , drop = FALSE])
    diff_metric[i] <- sum(abs(prop_matrix[i, ] - window_avg))
    results_list[[as.character(thresh_vals[i])]]$slide_diff_metric <- diff_metric[i]
  }
  
  # Create a summary dataframe.
  summary_df <- data.frame(threshold = thresh_vals, diff_metric = diff_metric)
  
  # Identify the "sweet spot" as the threshold where the diff_metric is minimized.
  sweet_spot <- summary_df$threshold[which.min(summary_df$diff_metric)]
  
  return(list(summary = summary_df, details = results_list, sweet_spot = sweet_spot))
}


## ----------------------------------------------------------------------------------------------------------
library(cowplot)
library(grid)
library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param sce PARAM_DESCRIPTION
#' @param master_cluster_column PARAM_DESCRIPTION
#' @param daughter_cluster_column PARAM_DESCRIPTION
#' @param meta PARAM_DESCRIPTION
#' @param sample_id_column PARAM_DESCRIPTION
#' @param daughter_order PARAM_DESCRIPTION, Default: NULL
#' @param row_order PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @export 
plotDonutProportionsNa <- function(sce, master_cluster_column, daughter_cluster_column, meta, sample_id_column, daughter_order = NULL, row_order = NULL) {
  # Convert to data frame
  df <- as.data.frame(colData(sce))

  df[[master_cluster_column]] <- factor(df[[master_cluster_column]], levels = unique(df[[master_cluster_column]]))
  df[[daughter_cluster_column]] <- factor(df[[daughter_cluster_column]], levels = unique(df[[daughter_cluster_column]])) # Ensure consistent factor levels
  df[[sample_id_column]] <- factor(df[[sample_id_column]], levels = unique(df[[sample_id_column]]))

  # Convert the daughter columns to factors with optional custom ordering
  if (!is.null(daughter_order)) {
    df[[daughter_cluster_column]] <- factor(df[[daughter_cluster_column]], levels = daughter_order)
  } else {
    df[[daughter_cluster_column]] <- factor(df[[daughter_cluster_column]], levels = unique(df[[daughter_cluster_column]]))
  }

  # Determine color palette based on unique levels of daughter clustering
  color_palette <- brewer.pal(min(length(levels(df[[daughter_cluster_column]])), 12), "Paired")

  # Calculate counts using dplyr to ensure only existing combinations are included
  count_df <- df %>%
    group_by_at(vars(sample_id_column, master_cluster_column, daughter_cluster_column)) %>%
    summarise(count = n(), .groups = "drop")

  # Ensure all combinations are present with zero counts where necessary
  complete_df <- count_df %>%
    complete(nesting(!!sym(sample_id_column), !!sym(master_cluster_column)), !!sym(daughter_cluster_column), fill = list(count = 0))

  completeout<<- complete_df
  
  # Calculate proportions within each sample_id and master_cluster_column
  proportion_df <- complete_df %>%
    group_by_at(vars(sample_id_column, master_cluster_column)) %>%
    mutate(proportion = count / sum(count) * 100) %>%
    ungroup()

  # Merge the cell counts with the proportion_df
  summary_df <- proportion_df %>%
    rename(cell_count = count)

  # Output the dataframe to the global environment
  summary_out <<- summary_df
  
  averaged_df <- summary_df %>%
    group_by_at(vars(master_cluster_column, daughter_cluster_column)) %>%
    summarise(mean_proportion = mean(proportion), .groups = "drop") %>%
    ungroup()
  
  average_out <<- averaged_df
  
  # Normalize the averaged proportions to ensure they sum to 100 within each master cluster
  normalized_df <- averaged_df %>%
    group_by_at(vars(master_cluster_column)) %>%
    mutate(normalized_proportion = mean_proportion / sum(mean_proportion) * 100) %>%
    ungroup()

  # Update the factor levels of master_cluster_column based on row_order
  if (!is.null(row_order)) {
    normalized_df[[master_cluster_column]] <- factor(normalized_df[[master_cluster_column]], levels = row_order)
  }

  # Create the plot
  p <- ggplot(normalized_df, aes(x = 2, y = normalized_proportion, fill = get(daughter_cluster_column))) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start = 0) +
    xlim(1, 2.5) +
    scale_fill_manual(values = color_palette, name = "Isotype") +
    theme_void() +
    facet_grid(rows = vars(get(master_cluster_column))) +
    theme(
      strip.text = element_text(size = 20), # Reduce strip text size
      legend.title = element_text(size = 20), # Reduce legend title size
      legend.text = element_text(size = 20), # Reduce legend text size
      panel.spacing = unit(0.1, "lines"), # Reduce space between panels
      plot.margin = margin(10, 10, 10, 10) # Reduce plot margins
    )

  # Return the final plot
  return(p)
}


## ----------------------------------------------------------------------------------------------------------
library(cowplot)
library(grid)
library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param sce PARAM_DESCRIPTION
#' @param master_cluster_column PARAM_DESCRIPTION
#' @param daughter_cluster_column PARAM_DESCRIPTION
#' @param meta PARAM_DESCRIPTION
#' @param facet_column PARAM_DESCRIPTION
#' @param sample_id_column PARAM_DESCRIPTION
#' @param plot_order_df PARAM_DESCRIPTION, Default: NULL
#' @param daughter_order PARAM_DESCRIPTION, Default: NULL
#' @param facet_order_df PARAM_DESCRIPTION, Default: NULL
#' @param cell_threshold PARAM_DESCRIPTION, Default: 100
#' @return OUTPUT_DESCRIPTION
#' @export 
plotDonutProportions <- function(sce, master_cluster_column, daughter_cluster_column, meta, facet_column, sample_id_column, plot_order_df = NULL, daughter_order = NULL, facet_order_df = NULL, cell_threshold = 100) {
  # Convert to data frame
  df <- as.data.frame(colData(sce))

  # Convert the master, daughter, facet, and sample_id columns to factors
  df[[master_cluster_column]] <- factor(df[[master_cluster_column]], levels = unique(df[[master_cluster_column]]))
  df[[daughter_cluster_column]] <- factor(df[[daughter_cluster_column]], levels = unique(df[[daughter_cluster_column]])) # Ensure consistent factor levels
  df[[sample_id_column]] <- factor(df[[sample_id_column]], levels = unique(df[[sample_id_column]]))

  # Convert the facet column to a factor with optional custom ordering
  if (!is.null(facet_order_df)) {
    df[[facet_column]] <- factor(df[[facet_column]], levels = as.character(unlist(facet_order_df)))
  } else {
    df[[facet_column]] <- factor(df[[facet_column]], levels = unique(df[[facet_column]]))
  }

  # Convert the daughter columns to factors with optional custom ordering
  if (!is.null(daughter_order)) {
    df[[daughter_cluster_column]] <- factor(df[[daughter_cluster_column]], levels = daughter_order)
  } else {
    df[[daughter_cluster_column]] <- factor(df[[daughter_cluster_column]], levels = unique(df[[daughter_cluster_column]]))
  }

  # Determine color palette based on unique levels of daughter clustering
  color_palette <- brewer.pal(min(length(levels(df[[daughter_cluster_column]])), 12), "Paired")

  # Calculate counts using dplyr to ensure only existing combinations are included
  count_df <- df %>%
    group_by_at(vars(sample_id_column, master_cluster_column, daughter_cluster_column, facet_column)) %>%
    summarise(count = n(), .groups = "drop")

  # Ensure all combinations are present with zero counts where necessary
  complete_df <- count_df %>%
    complete(nesting(!!sym(sample_id_column), !!sym(master_cluster_column), !!sym(facet_column)), !!sym(daughter_cluster_column), fill = list(count = 0))

  # Calculate proportions within each sample_id and master_cluster_column
  proportion_df <- complete_df %>%
    group_by_at(vars(sample_id_column, master_cluster_column, facet_column)) %>%
    mutate(proportion = count / sum(count) * 100) %>%
    ungroup()

  # Merge the cell counts with the proportion_df
  summary_df <- proportion_df %>%
    rename(cell_count = count)

  # Output the dataframe to the global environment
  summary_out <<- summary_df
  
  averaged_df <- summary_df %>%
    group_by_at(vars(master_cluster_column, daughter_cluster_column, facet_column)) %>%
    summarise(mean_proportion = mean(proportion), total_cells = sum(cell_count), .groups = "drop") %>%
    ungroup()
  
  averaged_out <<- averaged_df
  
  # Filter out combinations of master_cluster_column and facet_column that have a sum of cells less than the threshold
  filtered_averaged_df <- averaged_df %>%
    group_by_at(vars(master_cluster_column, facet_column)) %>%
    filter(sum(total_cells) >= cell_threshold) %>%
    ungroup()
  
  filtout <<- filtered_averaged_df
  
  # Normalize the averaged proportions to ensure they sum to 100 within each master cluster and facet
  normalized_df <- filtered_averaged_df %>%
    group_by_at(vars(master_cluster_column, facet_column)) %>%
    mutate(normalized_proportion = mean_proportion / sum(mean_proportion) * 100) %>%
    ungroup()

  # Determine the order of master clusters for plotting
  plot_order <- if (!is.null(plot_order_df)) {
    as.character(unlist(plot_order_df))
  } else {
    levels(normalized_df[[master_cluster_column]])
  }

  # Create the plot
  p <- ggplot(normalized_df, aes(x = 2, y = normalized_proportion, fill = get(daughter_cluster_column))) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start = 0) +
    xlim(2, 2.5) +
    scale_fill_manual(values = color_palette, name = "Isotype") +
    theme_void() +
    facet_grid(rows = vars(get(master_cluster_column)), cols = vars(get(facet_column))) +
    theme(
      strip.text.x = element_text(size = 20, angle = 90), # Rotate top strip text 90 degrees
      strip.text.y = element_text(size = 20), # Keep side strip text as is
      legend.title = element_text(size = 20), # Reduce legend title size
      legend.text = element_text(size = 20), # Reduce legend text size
      panel.spacing = unit(0.1, "lines"), # Reduce space between panels
      plot.margin = margin(10, 10, 10, 10) # Reduce plot margins
    )

  # Return the final plot
  return(p)
}


## ----------------------------------------------------------------------------------------------------------
# New function to plot scatter plot of cluster abundances between two technical replicates
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param k_abundances PARAM_DESCRIPTION, Default: 'merging1'
#' @param meta PARAM_DESCRIPTION, Default: c("sample_id", "patient_id", "condition", "exp")
#' @param facet_by PARAM_DESCRIPTION, Default: NULL
#' @param text_size PARAM_DESCRIPTION, Default: 16
#' @param k_pal PARAM_DESCRIPTION, Default: CATALYST:::.cluster_cols
#' @return OUTPUT_DESCRIPTION
#' @export 
#' @importFrom CATALYST .cluster_cols
plotAbundanceScatter <- function(x, k_abundances = "merging1", meta = c("sample_id", "patient_id", "condition", "exp"), facet_by = NULL, text_size = 16, k_pal = CATALYST:::.cluster_cols) {
  cluster_ids_abundances <- colData(x)[[k_abundances]]
  
  ns <- table(cluster_id = cluster_ids_abundances, sample_id = colData(x)$sample_id)
  fq <- prop.table(ns, 2) * 100
  df <- as.data.frame(fq)
  m <- match(df$sample_id, colData(x)$sample_id)
  for (i in meta) df[[i]] <- colData(x)[[i]][m]
  
  # Filter to include only the two replicates
  replicates <- unique(df$exp)
  if (length(replicates) != 2) {
    stop("The data must contain exactly two replicates in the 'exp' column.")
  }
  
  df <- df %>% filter(exp %in% replicates)
  
  dfout1 <<- df
  
  # Split the data into two data frames for each replicate
  df1 <- df %>% filter(exp == replicates[1]) %>% rename(Freq1 = Freq)
  df2 <- df %>% filter(exp == replicates[2]) %>% rename(Freq2 = Freq)
  
  # Join the data frames based on patient_id, cluster_id
  df_wide <- left_join(df1, df2, by = c("patient_id", "cluster_id","day"))
  
  dfout2 <<- df_wide
  
  # Plot scatter plot
  p <- ggplot(df_wide, aes(x = Freq1, y = Freq2, color = cluster_id)) +
    geom_point(aes(shape = day), size = 3) +
    labs(x = paste("Replicate", replicates[1]), y = paste("Replicate", replicates[2]), color = "Cluster ID") +
    theme_bw() +
    theme(panel.grid = element_blank(), 
          strip.text = element_text(face = "bold"), 
          strip.background = element_rect(fill = "grey88", color = NA), 
          aspect.ratio = 1,
          axis.text = element_text(color = "black"), 
          text = element_text(size = text_size),
          panel.spacing = unit(1, "lines"),
          legend.key.height = unit(0.8, "lines")) +
    scale_color_manual(values = k_pal)
  
  p <- p + xlim(0,100) + ylim(0,100)
  
  if (!is.null(facet_by)) {
    p <- p + facet_wrap(vars(!!sym(facet_by)), scales = "free_x")
  }
  
  return(p)
}




## ----------------------------------------------------------------------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param div_col PARAM_DESCRIPTION, Default: 'div'
#' @param parse_div PARAM_DESCRIPTION, Default: TRUE
#' @param metric PARAM_DESCRIPTION, Default: c("mean_div", "median_div", "percent_proliferated")
#' @param scales PARAM_DESCRIPTION, Default: 'free_y'
#' @param conditions PARAM_DESCRIPTION, Default: NULL
#' @param condition_pairs PARAM_DESCRIPTION, Default: NULL
#' @param facet_by PARAM_DESCRIPTION, Default: NULL
#' @param facet PARAM_DESCRIPTION, Default: TRUE
#' @param excluded_clusters PARAM_DESCRIPTION, Default: NULL
#' @param axes PARAM_DESCRIPTION, Default: 'all'
#' @param fun PARAM_DESCRIPTION, Default: c("median", "mean", "sum")
#' @param point_size PARAM_DESCRIPTION, Default: 1
#' @param textsize PARAM_DESCRIPTION, Default: 14
#' @param panel_spacing PARAM_DESCRIPTION, Default: 2
#' @param show_stats PARAM_DESCRIPTION, Default: TRUE
#' @param hide_ns PARAM_DESCRIPTION, Default: TRUE
#' @param fill_by PARAM_DESCRIPTION, Default: NULL
#' @param x_group PARAM_DESCRIPTION, Default: 'day'
#' @param shape_by PARAM_DESCRIPTION, Default: NULL
#' @param stat_size PARAM_DESCRIPTION, Default: 4
#' @param facet_ratio PARAM_DESCRIPTION, Default: 1.5
#' @param geom PARAM_DESCRIPTION, Default: c("boxes", "bar")
#' @param jitter PARAM_DESCRIPTION, Default: TRUE
#' @param group1keep PARAM_DESCRIPTION, Default: NULL
#' @param nudge PARAM_DESCRIPTION, Default: 0
#' @param step_increase PARAM_DESCRIPTION, Default: 0.1
#' @param hide_x_labels PARAM_DESCRIPTION, Default: FALSE
#' @param label_parse PARAM_DESCRIPTION, Default: FALSE
#' @param merging_col PARAM_DESCRIPTION, Default: NULL
#' @param swap PARAM_DESCRIPTION, Default: FALSE
#' @param external_stats PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @export 
#' @importFrom purrr imap_dfr
#' @importFrom dplyr filter mutate group_by summarise n case_when
#' @importFrom ggbeeswarm geom_quasirandom
plotDivPairs <- function (
    x, 
    div_col = "div",  
    parse_div = TRUE,  
    metric = c("mean_div", "median_div", "percent_proliferated"),  
    scales = "free_y", 
    conditions = NULL, 
    condition_pairs = NULL, 
    facet_by = NULL,
    facet = TRUE,
    excluded_clusters = NULL, 
    axes = "all", 
    fun = c("median", "mean", "sum"), 
    point_size = 1, 
    textsize = 14, 
    panel_spacing = 2, 
    show_stats = TRUE, 
    hide_ns = TRUE, 
    fill_by = NULL,
    x_group = "day",  
    shape_by = NULL, 
    stat_size = 4, 
    facet_ratio = 1.5, 
    geom = c("boxes","bar"), 
    jitter = TRUE, 
    group1keep = NULL, 
    nudge = 0, 
    step_increase = 0.1, 
    hide_x_labels = FALSE, 
    label_parse = FALSE, 
    merging_col = NULL, 
    swap = FALSE,
    external_stats = NULL
) {
    fun <- match.arg(fun)
    geom <- match.arg(geom)
    metric <- match.arg(metric)
  
    # Verify required column patient_id exists
    if (!"patient_id" %in% colnames(colData(x))) {
        stop("Column 'patient_id' is not found in colData(x).")
    }
  
    # Pull out the div column
    if (!div_col %in% colnames(colData(x))) {
        stop("Specified 'div_col' does not exist.")
    }
    div_values <- colData(x)[[div_col]]
    if (parse_div) {
        div_values <- as.numeric(gsub("Div", "", div_values))
    }
    df <- data.frame(div = div_values, colData(x))
  
    if (!is.null(conditions)) {
        df <- df[df[[x_group]] %in% conditions, ]
        df[[x_group]] <- factor(df[[x_group]], levels = conditions, ordered = TRUE)
    }
    if (!is.null(excluded_clusters)) {
        df <- df[!df$cluster_id %in% excluded_clusters, ]
    }
  
    # If condition_pairs specified, build df with a 'condition_pair' column
    if (!is.null(condition_pairs) && length(condition_pairs) > 0) {
      df <- purrr::imap_dfr(condition_pairs, function(pair, pname) {
        df %>%
          dplyr::filter(!!sym(fill_by) %in% pair) %>%
          dplyr::mutate(condition_pair = pname)
      })
      facet_by <- "condition_pair" 
    }
  
    # Determine summary function
    agg_fun <- if (metric == "mean_div") mean else if (metric == "median_div") median
  
    # Summaries
    if (metric == "percent_proliferated") {
        df_agg <- df %>%
            dplyr::group_by(sample_id, !!sym(x_group), !!sym(facet_by), patient_id, !!sym(fill_by)) %>%
            dplyr::summarise(div = sum(div != 0) / dplyr::n() * 100, .groups = 'drop')
    } else {
        df_agg <- df %>%
            dplyr::group_by(sample_id, !!sym(x_group), !!sym(facet_by), patient_id, !!sym(fill_by)) %>%
            dplyr::summarise(div = agg_fun(div, na.rm = TRUE), .groups = 'drop')
    }
    if (nrow(df_agg) == 0) stop("No data available after filtering.")
  
    # Create grouping variable to connect points from the same patient on the same day
    #df_agg$grouper <- paste(df_agg[[x_group]], df_agg$patient_id, sep = "_")
  
    # Basic plot
    p <- ggplot(df_agg, aes(x = .data[[x_group]], y = div)) +
      labs(y = dplyr::case_when(
        metric == "percent_proliferated" ~ "Percentage Proliferated",
        metric == "mean_div" ~ "Mean Division",
        TRUE ~ "Median Division"
      )) +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        text = element_text(size = textsize),
        strip.text = element_text(size = textsize),
        strip.background = element_rect(fill = NA, color = NA),
        strip.placement = "outside",
        legend.text = element_text(size = textsize),
        legend.title = element_text(size = textsize),
        aspect.ratio = facet_ratio,
        axis.ticks = element_line(color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.spacing = unit(panel_spacing, "lines"),
        axis.text = element_text(color = "black", size = textsize),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_text(size = textsize)
      )
  
    # Box or bar plot
    if (geom == "boxes") {
      p <- p + geom_boxplot(
        aes(
          fill = .data[[fill_by]],
          group = interaction(.data[[x_group]], .data[[fill_by]])
        ),
        position = position_dodge(width = 0.7),
        color = "black",
        width = 0.6,
        alpha = 0.8,
        outlier.color = NA,
        size = 0.5  # Set border size
      )
    } else {
      p <- p + stat_summary(
        fun = match.fun(fun),
        geom = "bar",
        aes(
          fill = .data[[fill_by]],
          group = interaction(.data[[x_group]], .data[[fill_by]])
        ),
        position = position_dodge(width = 0.7),
        width = 0.6,
        color = "black",
        alpha = 0.8,
        size = 0.5  # Set border size
      )
    }
    
    # Points: map shape by 'shape_by' if provided.
    # Default: use custom manual shapes for up to 4 patients: 21, 22, 24, 25.
    if (!is.null(shape_by)) {
      p <- p + ggbeeswarm::geom_quasirandom(
        aes(
          fill = .data[[fill_by]],
          shape = .data[[shape_by]],
          group = interaction(patient_id, .data[[fill_by]])
        ),
        color = "black",
        stroke = 0.5,
        dodge.width = 0.7,
        size = point_size,
        alpha = 0.9
      ) +
        scale_shape_manual(values = c(21, 22, 24, 25))
    } else {
      p <- p + geom_point(
        aes(
          fill = .data[[fill_by]],
          group = interaction(patient_id, .data[[fill_by]])
        ),
        shape = 21,
        color = "black",
        stroke = 0.5,
        position = position_dodge(width = 0.7),
        size = point_size,
        alpha = 0.9
      )
    }
    
   # This just plots vertical lines.  Segments are not plotted on x properly. 
  #  p <- p + geom_line(
    #  aes(
     #   group = interaction(patient_id, .data[[facet_by]], .data[[x_group]])
    #  ),
      #position = position_dodge(width = 0.7),
     # color = "grey60",
    #  alpha = 0.7
   # )
    
    if (hide_x_labels) {
      p <- p + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
    }
    
    if (facet) {
      facet_formula <- as.formula(paste("~", facet_by))
      p <- p + facet_grid(facet_formula, scales = scales)
    }
    if (metric == "percent_proliferated") {
      p <- p + scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20))
    }
  
    return(p)
}


## ----------------------------------------------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(scales)

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param div_col PARAM_DESCRIPTION, Default: 'div'
#' @param parse_div PARAM_DESCRIPTION, Default: TRUE
#' @param x_group PARAM_DESCRIPTION, Default: 'day'
#' @param facet_by PARAM_DESCRIPTION, Default: 'condition'
#' @param colors PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @export 
plotDivisionStackedBar <- function(
    x, 
    div_col = "div",  # Name of the Div column
    parse_div = TRUE,  # Whether to parse Div column (e.g., Div0 → 0)
    x_group = "day",  # Grouping variable for the x-axis (e.g., time points)
    facet_by = "condition",  # Facet variable (e.g., conditions)
    colors = NULL  # Optional: Custom colors for divisions
) {
    # Retrieve the div column
    if (!div_col %in% colnames(colData(x))) {
        stop("The specified 'div_col' (", div_col, ") does not exist in the colData of the input object.")
    }
    
    # Extract the div column and convert to numeric if parse_div is TRUE
    div_values <- colData(x)[[div_col]]
    if (parse_div) {
        div_values <- as.numeric(gsub("Div", "", div_values))  # Extract numerical part (e.g., Div0 → 0)
    }
    
    # Create the data frame
    df <- data.frame(div = div_values, colData(x))
    
    # Calculate proportions of cells in each division for each sample, condition, and time point
    df_prop <- df %>%
        group_by(sample_id, !!sym(x_group), !!sym(facet_by)) %>%
        count(div) %>%
        mutate(prop = n / sum(n)) %>%
        ungroup()
    
    # Aggregate proportions across samples (if needed)
    df_agg <- df_prop %>%
        group_by(!!sym(x_group), !!sym(facet_by), div) %>%
        summarise(prop = mean(prop, na.rm = TRUE), .groups = 'drop')
    
    # Plot
    p <- ggplot(df_agg, aes(x = !!sym(x_group), y = prop, fill = factor(div))) +
        geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.2) +  # Add black border
        labs(
            x = x_group,
            y = "Proportion of Cells",
            fill = "Division"
        ) +
        theme_bw() +
        theme(
            panel.grid = element_blank(),
            strip.background = element_rect(fill = NA, color = NA),
            strip.text = element_text(size = 12),
            panel.border = element_rect(color = "black", fill = NA, size = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            legend.position = "right"
        ) +
        facet_wrap(as.formula(paste("~", facet_by)), scales = "free_x")
    
    # Use the same color palette as division_profiler
    if (is.null(colors)) {
        p <- p + scale_fill_brewer(palette = "Paired")  # Use "Paired" palette
    } else {
        p <- p + scale_fill_manual(values = colors)  # Use custom colors if provided
    }
    
    return(p)
}


## ----------------------------------------------------------------------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param object PARAM_DESCRIPTION
#' @param features PARAM_DESCRIPTION
#' @param cols PARAM_DESCRIPTION, Default: NULL
#' @param pt.size PARAM_DESCRIPTION, Default: NULL
#' @param alpha PARAM_DESCRIPTION, Default: 1
#' @param idents PARAM_DESCRIPTION, Default: NULL
#' @param sort PARAM_DESCRIPTION, Default: FALSE
#' @param assay PARAM_DESCRIPTION, Default: NULL
#' @param group.by PARAM_DESCRIPTION, Default: NULL
#' @param split.by PARAM_DESCRIPTION, Default: NULL
#' @param adjust PARAM_DESCRIPTION, Default: 1
#' @param y.max PARAM_DESCRIPTION, Default: NULL
#' @param same.y.lims PARAM_DESCRIPTION, Default: FALSE
#' @param log PARAM_DESCRIPTION, Default: FALSE
#' @param ncol PARAM_DESCRIPTION, Default: NULL
#' @param slot PARAM_DESCRIPTION, Default: deprecated()
#' @param layer PARAM_DESCRIPTION, Default: NULL
#' @param split.plot PARAM_DESCRIPTION, Default: FALSE
#' @param stack PARAM_DESCRIPTION, Default: FALSE
#' @param combine PARAM_DESCRIPTION, Default: TRUE
#' @param fill.by PARAM_DESCRIPTION, Default: 'feature'
#' @param flip PARAM_DESCRIPTION, Default: FALSE
#' @param add.noise PARAM_DESCRIPTION, Default: TRUE
#' @param raster PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @export 
VlnPlot2 <- function (object, features, cols = NULL, pt.size = NULL, alpha = 1, 
    idents = NULL, sort = FALSE, assay = NULL, group.by = NULL, 
    split.by = NULL, adjust = 1, y.max = NULL, same.y.lims = FALSE, 
    log = FALSE, ncol = NULL, slot = deprecated(), layer = NULL, 
    split.plot = FALSE, stack = FALSE, combine = TRUE, fill.by = "feature", 
    flip = FALSE, add.noise = TRUE, raster = NULL) 
{
    if (is_present(arg = slot)) {
        deprecate_soft(when = "5.0.0", what = "VlnPlot(slot = )", 
            with = "VlnPlot(layer = )")
        layer <- slot %||% layer
    }
    layer.set <- suppressWarnings(Layers(object = object, search = layer %||% 
        "data"))
    if (is.null(layer) && length(layer.set) == 1 && layer.set == 
        "scale.data") {
        warning("Default search for \"data\" layer yielded no results; utilizing \"scale.data\" layer instead.")
    }
    assay.name <- DefaultAssay(object)
    if (is.null(layer.set) & is.null(layer)) {
        warning("Default search for \"data\" layer in \"", assay.name, 
            "\" assay yielded no results; utilizing \"counts\" layer instead.", 
            call. = FALSE, immediate. = TRUE)
        layer.set <- Layers(object = object, search = "counts")
    }
    if (is.null(layer.set)) {
        stop("layer \"", layer, "\" is not found in assay: \"", 
            assay.name, "\"")
    }
    else {
        layer <- layer.set
    }
    if (!is.null(x = split.by) & getOption(x = "Seurat.warn.vlnplot.split", 
        default = TRUE)) {
        message("The default behaviour of split.by has changed.\n", 
            "Separate violin plots are now plotted side-by-side.\n", 
            "To restore the old behaviour of a single split violin,\n", 
            "set split.plot = TRUE.\n      \nThis message will be shown once per session.")
        options(Seurat.warn.vlnplot.split = FALSE)
    }
    return(ExIPlot2(object = object, type = ifelse(test = split.plot, 
        yes = "splitViolin", no = "violin"), features = features, 
        idents = idents, ncol = ncol, sort = sort, assay = assay, 
        y.max = y.max, same.y.lims = same.y.lims, adjust = adjust, 
        pt.size = pt.size, alpha = alpha, cols = cols, group.by = group.by, 
        split.by = split.by, log = log, layer = layer, stack = stack, 
        combine = combine, fill.by = fill.by, flip = flip, add.noise = add.noise, 
        raster = raster))
}

environment(VlnPlot2) <- asNamespace('Seurat')



## ----------------------------------------------------------------------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param object PARAM_DESCRIPTION
#' @param features PARAM_DESCRIPTION
#' @param type PARAM_DESCRIPTION, Default: 'violin'
#' @param idents PARAM_DESCRIPTION, Default: NULL
#' @param ncol PARAM_DESCRIPTION, Default: NULL
#' @param sort PARAM_DESCRIPTION, Default: FALSE
#' @param assay PARAM_DESCRIPTION, Default: NULL
#' @param y.max PARAM_DESCRIPTION, Default: NULL
#' @param same.y.lims PARAM_DESCRIPTION, Default: FALSE
#' @param adjust PARAM_DESCRIPTION, Default: 1
#' @param cols PARAM_DESCRIPTION, Default: NULL
#' @param pt.size PARAM_DESCRIPTION, Default: 0
#' @param alpha PARAM_DESCRIPTION, Default: 1
#' @param group.by PARAM_DESCRIPTION, Default: NULL
#' @param split.by PARAM_DESCRIPTION, Default: NULL
#' @param log PARAM_DESCRIPTION, Default: FALSE
#' @param slot PARAM_DESCRIPTION, Default: deprecated()
#' @param layer PARAM_DESCRIPTION, Default: 'data'
#' @param stack PARAM_DESCRIPTION, Default: FALSE
#' @param combine PARAM_DESCRIPTION, Default: TRUE
#' @param fill.by PARAM_DESCRIPTION, Default: NULL
#' @param flip PARAM_DESCRIPTION, Default: FALSE
#' @param add.noise PARAM_DESCRIPTION, Default: TRUE
#' @param raster PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @export 
ExIPlot2 <- function (object, features, type = "violin", idents = NULL, ncol = NULL, 
    sort = FALSE, assay = NULL, y.max = NULL, same.y.lims = FALSE, 
    adjust = 1, cols = NULL, pt.size = 0, alpha = 1, group.by = NULL, 
    split.by = NULL, log = FALSE, slot = deprecated(), layer = "data", 
    stack = FALSE, combine = TRUE, fill.by = NULL, flip = FALSE, 
    add.noise = TRUE, raster = NULL) 
{
    if (is_present(arg = slot)) {
        layer <- layer %||% slot
    }
    assay <- assay %||% DefaultAssay(object = object)
    DefaultAssay(object = object) <- assay
    cells <- Cells(x = object, assay = NULL)
    if (isTRUE(x = stack)) {
        if (!is.null(x = ncol)) {
            warning("'ncol' is ignored with 'stack' is TRUE", 
                call. = FALSE, immediate. = TRUE)
        }
        if (!is.null(x = y.max)) {
            warning("'y.max' is ignored when 'stack' is TRUE", 
                call. = FALSE, immediate. = TRUE)
        }
    }
    else {
        ncol <- ncol %||% ifelse(test = length(x = features) > 
            9, yes = 4, no = min(length(x = features), 3))
    }
    if (!is.null(x = idents)) {
        cells <- intersect(x = names(x = Idents(object = object)[Idents(object = object) %in% 
            idents]), y = cells)
    }
    data <- FetchData(object = object, vars = features, slot = layer, 
        cells = cells)
    pt.size <- pt.size %||% AutoPointSize(data = object)
    features <- colnames(x = data)
    data <- data[cells, , drop = FALSE]
    idents <- if (is.null(x = group.by)) {
        Idents(object = object)[cells]
    }
    else {
        object[[group.by, drop = TRUE]][cells]
    }
    if (!is.factor(x = idents)) {
        idents <- factor(x = idents)
    }
    if (is.null(x = split.by)) {
        split <- NULL
    }
    else {
        split <- FetchData(object, split.by)[cells, split.by]
        if (!is.factor(x = split)) {
            split <- factor(x = split)
        }
        if (is.null(x = cols)) {
            cols <- hue_pal()(length(x = levels(x = idents)))
            cols <- Interleave(cols, InvertHex(hexadecimal = cols))
        }
        else if (length(x = cols) == 1 && cols == "interaction") {
            split <- interaction(idents, split)
            cols <- hue_pal()(length(x = levels(x = idents)))
        }
        else {
            cols <- Col2Hex(cols)
        }
        if (length(x = cols) < length(x = levels(x = split))) {
            cols <- Interleave(cols, InvertHex(hexadecimal = cols))
        }
        cols <- rep_len(x = cols, length.out = length(x = levels(x = split)))
        names(x = cols) <- levels(x = split)
        if ((length(x = cols) > 2) & (type == "splitViolin")) {
            warning("Split violin is only supported for <3 groups, using multi-violin.")
            type <- "violin"
        }
    }
    if (same.y.lims && is.null(x = y.max)) {
        y.max <- max(data)
    }
    if (isTRUE(x = stack)) {
        return(MultiExIPlot(type = type, data = data, idents = idents, 
            split = split, sort = sort, same.y.lims = same.y.lims, 
            adjust = adjust, cols = cols, pt.size = pt.size, 
            log = log, fill.by = fill.by, add.noise = add.noise, 
            flip = flip))
    }
    plots <- lapply(X = features, FUN = function(x) {
        return(SingleExIPlot2(type = type, data = data[, x, drop = FALSE], 
            idents = idents, split = split, sort = sort, y.max = y.max, 
            adjust = adjust, cols = cols, pt.size = pt.size, 
            alpha = alpha, log = log, add.noise = add.noise, 
            raster = raster))
    })
    label.fxn <- switch(EXPR = type, violin = if (stack) {
        xlab
    } else {
        ylab
    }, splitViolin = if (stack) {
        xlab
    } else {
        ylab
    }, ridge = xlab, stop("Unknown ExIPlot type ", type, call. = FALSE))
    for (i in 1:length(x = plots)) {
        key <- paste0(unlist(x = strsplit(x = features[i], split = "_"))[1], 
            "_")
        obj <- names(x = which(x = Key(object = object) == key))
        if (length(x = obj) == 1) {
            if (inherits(x = object[[obj]], what = "DimReduc")) {
                plots[[i]] <- plots[[i]] + label.fxn(label = "Embeddings Value")
            }
            else if (inherits(x = object[[obj]], what = "Assay") || 
                inherits(x = object[[obj]], what = "Assay5")) {
                next
            }
            else {
                warning("Unknown object type ", class(x = object), 
                  immediate. = TRUE, call. = FALSE)
                plots[[i]] <- plots[[i]] + label.fxn(label = NULL)
            }
        }
        else if (!features[i] %in% rownames(x = object)) {
            plots[[i]] <- plots[[i]] + label.fxn(label = NULL)
        }
    }
    if (combine) {
        plots <- wrap_plots(plots, ncol = ncol)
        if (length(x = features) > 1) {
            plots <- plots & NoLegend()
        }
    }
    return(plots)
}

environment(ExIPlot2) <- asNamespace('Seurat')



## ----------------------------------------------------------------------------------------------------------
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param data PARAM_DESCRIPTION
#' @param idents PARAM_DESCRIPTION
#' @param split PARAM_DESCRIPTION, Default: NULL
#' @param type PARAM_DESCRIPTION, Default: 'violin'
#' @param sort PARAM_DESCRIPTION, Default: FALSE
#' @param y.max PARAM_DESCRIPTION, Default: NULL
#' @param adjust PARAM_DESCRIPTION, Default: 1
#' @param pt.size PARAM_DESCRIPTION, Default: 0
#' @param alpha PARAM_DESCRIPTION, Default: 1
#' @param cols PARAM_DESCRIPTION, Default: NULL
#' @param seed.use PARAM_DESCRIPTION, Default: 42
#' @param log PARAM_DESCRIPTION, Default: FALSE
#' @param add.noise PARAM_DESCRIPTION, Default: TRUE
#' @param raster PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @export 
#' @importFrom ggrastr geom_jitter_rast
SingleExIPlot2 <- function (data, idents, split = NULL, type = "violin", sort = FALSE, 
    y.max = NULL, adjust = 1, pt.size = 0, alpha = 1, cols = NULL, 
    seed.use = 42, log = FALSE, add.noise = TRUE, raster = NULL) 
{
    if (!is.null(x = raster) && isTRUE(x = raster)) {
        if (!PackageCheck("ggrastr", error = FALSE)) {
            stop("Please install ggrastr from CRAN to enable rasterization.")
        }
    }
    if (PackageCheck("ggrastr", error = FALSE)) {
        if ((nrow(x = data) > 1e+05) & is.null(x = raster)) {
            message("Rasterizing points since number of points exceeds 100,000.", 
                "\nTo disable this behavior set `raster=FALSE`")
            raster <- TRUE
        }
    }
    if (!is.null(x = seed.use)) {
        set.seed(seed = seed.use)
    }
    if (!is.data.frame(x = data) || ncol(x = data) != 1) {
        stop("'SingleExIPlot requires a data frame with 1 column")
    }
    feature <- colnames(x = data)
    data$ident <- idents
    if ((is.character(x = sort) && nchar(x = sort) > 0) || sort) {
        data$ident <- factor(x = data$ident, levels = names(x = rev(x = sort(x = tapply(X = data[, 
            feature], INDEX = data$ident, FUN = mean), decreasing = grepl(pattern = paste0("^", 
            tolower(x = sort)), x = "decreasing")))))
    }
    if (log) {
        noise <- rnorm(n = length(x = data[, feature]))/200
        data[, feature] <- data[, feature] + 1
    }
    else {
        noise <- rnorm(n = length(x = data[, feature]))/1e+05
    }
    if (!add.noise) {
        noise <- noise * 0
    }
    if (all(data[, feature] == data[, feature][1])) {
        warning(paste0("All cells have the same value of ", feature, 
            "."))
    }
    else {
        data[, feature] <- data[, feature] + noise
    }
    axis.label <- "Expression Level"
    y.max <- y.max %||% max(data[, feature][is.finite(x = data[, 
        feature])])
    if (type == "violin" && !is.null(x = split)) {
        data$split <- split
        vln.geom <- geom_violin
        fill <- "split"
    }
    else if (type == "splitViolin" && !is.null(x = split)) {
        data$split <- split
        vln.geom <- geom_split_violin
        fill <- "split"
        type <- "violin"
    }
    else {
        vln.geom <- geom_violin
        fill <- "ident"
    }
    switch(EXPR = type, violin = {
        x <- "ident"
        y <- paste0("`", feature, "`")
        xlab <- "Identity"
        ylab <- axis.label
        geom <- list(vln.geom(scale = "width", adjust = adjust, 
            trim = TRUE), theme(axis.text.x = element_text(angle = 45, 
            hjust = 1)))
        if (is.null(x = split)) {
            jitter <- if (isTRUE(x = raster)) {
                ggrastr::geom_jitter_rast(height = 0, size = pt.size,shape = 16, alpha = alpha, show.legend = FALSE, raster.dpi = 600)
            } else {
                geom_jitter(height = 0, size = pt.size, alpha = alpha, show.legend = FALSE)
            }
        } else {
            jitter <- if (isTRUE(x = raster)) {
                ggrastr::geom_jitter_rast(position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.9), shape = 16, raster.dpi = 600, size = pt.size, alpha = alpha, show.legend = FALSE)
            } else {
                geom_jitter(position = position_jitterdodge(jitter.width = 0.4, 
                  dodge.width = 0.9), size = pt.size, alpha = alpha, show.legend = FALSE)
            }
        }
        log.scale <- scale_y_log10()
        axis.scale <- ylim
    }, ridge = {
        x <- paste0("`", feature, "`")
        y <- "ident"
        xlab <- axis.label
        ylab <- "Identity"
        geom <- list(geom_density_ridges(scale = 4), theme_ridges(), 
            scale_y_discrete(expand = c(0.01, 0)), scale_x_continuous(expand = c(0, 
                0)))
        jitter <- geom_jitter(width = 0, size = pt.size, alpha = alpha, 
            show.legend = FALSE)
        log.scale <- scale_x_log10()
        axis.scale <- function(...) {
            invisible(x = NULL)
        }
    }, stop("Unknown plot type: ", type))
    plot <- ggplot(data = data, mapping = aes_string(x = x, y = y, 
        fill = fill)[c(2, 3, 1)]) + labs(x = xlab, y = ylab, 
        title = feature, fill = NULL) + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
    plot <- do.call(what = "+", args = list(plot, geom))
    plot <- plot + if (log) {
        log.scale
    }
    else {
        axis.scale(min(data[, feature]), y.max)
    }
    if (pt.size > 0) {
        plot <- plot + jitter
    }
    if (!is.null(x = cols)) {
        if (!is.null(x = split)) {
            idents <- unique(x = as.vector(x = data$ident))
            splits <- unique(x = as.vector(x = data$split))
            labels <- if (length(x = splits) == 2) {
                splits
            }
            else {
                unlist(x = lapply(X = idents, FUN = function(pattern, 
                  x) {
                  x.mod <- gsub(pattern = paste0(pattern, "."), 
                    replacement = paste0(pattern, ": "), x = x, 
                    fixed = TRUE)
                  x.keep <- grep(pattern = ": ", x = x.mod, fixed = TRUE)
                  x.return <- x.mod[x.keep]
                  names(x = x.return) <- x[x.keep]
                  return(x.return)
                }, x = unique(x = as.vector(x = data$split))))
            }
            if (is.null(x = names(x = labels))) {
                names(x = labels) <- labels
            }
        }
        else {
            labels <- levels(x = droplevels(data$ident))
        }
        plot <- plot + scale_fill_manual(values = cols, labels = labels)
    }
    return(plot)
}

environment(SingleExIPlot2) <- asNamespace('Seurat')

