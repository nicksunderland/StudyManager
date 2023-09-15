#' @title Histogram with outlier bins
#' @description
#' Create a histogram with outlier bins
#' Credit: https://rdrr.io/cran/kim/src/R/histogram_w_outlier_bins.R
#' @param data a data.table or data.frame
#' @param col_name a string col name to extract the numeric data from
#' @param col_group a string group name col
#' @param bin_cutoffs cutoff points for bins
#' @param outlier_bin_left logical. Should the leftmost bin treated as an outlier bin? (default =
#'   TRUE)
#' @param outlier_bin_right logical. Should the rightmost bin treated as an outlier bin? (default =
#'   TRUE)
#' @param x_tick_marks a vector of values at which to place tick marks on the x axis. Note that the
#'   first bar spans from 0.5 to 1.5, second bar from 1.5 to 2.5, ... nth bar from n - 0.5 to n +
#'   0.5. See the example. By default, tick marks will be placed at every cutoff point for bins
#' @param x_tick_mark_labels a character vector to label tick marks. By default, the vector of
#'   cutoff points for bins will also be used as labels.
#' @param y_tick_marks a vector of values at which to place tick marks on the y axis (e.g., setting
#'   \code{y_tick_marks = seq(0, 10, 5)} will put tick marks at 0, 5, and 10.)
#' @param outlier_bin_fill_color color to fill inside of the outlier bins (default = "coral")
#' @param non_outlier_bin_fill_color color to fill inside of the non-outlier bins (default =
#'   "cyan4")
#' @param border_color color for borders of the bins (default = "black")
#' @param y_axis_title_vjust position of the y axis title (default = 0.85).
#' @param x_axis_title title for x axis (default = "Value"). If \code{x_axis_title = FALSE}, x axis
#'   title will be removed from the plot.
#' @param y_axis_title title for y axis. By default, it will be either "Proportion" or "Count".
#' @param notify_na_count if \code{TRUE}, notify how many observations were removed due to missing
#'   values. By default, NA count will be printed only if there are any NA values.
#' @param plot_proportion logical. Should proportions be plotted, as opposed to frequencies?
#'   (default = TRUE)
#' @param plot_frequency logical. Should frequencies be plotted, as opposed to proportions? (default
#'   = FALSE). If \code{plot_frequency = TRUE}, \code{plot_proportion} will switch to be FALSE.
#' @return a ggplot object
#'
#' @importFrom ggplot2 ggplot aes geom_bar theme xlab ylab element_blank
#' @importFrom utils head tail
#' @importFrom rlang .data :=
#' @importFrom dplyr filter select mutate pull all_of
#' @importFrom tibble as_tibble
#' @importFrom forcats fct_relevel
#' @export
#'
histogram_w_outlier_bins <- function(
  data,
  col_name,
  col_group = NULL,
  bin_cutoffs = NULL,
  outlier_bin_left = TRUE,
  outlier_bin_right = TRUE,
  x_tick_marks = NULL,
  x_tick_mark_labels = NULL,
  y_tick_marks = NULL,
  outlier_bin_fill_color = "coral",
  non_outlier_bin_fill_color = "cyan4",
  border_color = "black",
  y_axis_title_vjust = 0.85,
  x_axis_title = NULL,
  y_axis_title = NULL,
  notify_na_count = NULL,
  plot_proportion = TRUE,
  plot_frequency = FALSE) {

  # if no grouping then add single group
  if(is.null(col_group)) {
    col_group <- "group"
    data[ , col_group] <- 1
  }

  # extract the data
  data <- as_tibble(data)
  dat <- data |> select(all_of(c(col_name, col_group)))

  # make sure grouping a factor
  if(!is.factor(dat[[col_group]])) {
    dat <- dat |> mutate(!!as.symbol(col_group) := factor(!!as.symbol(col_group)))
  }
  group_levels <- levels(dat[[col_group]])

  # deal with NA values
  na_flag <- is.na(dat[[col_name]])
  na_count <- sum(na_flag)
  dat <- dat |> filter(!na_flag)

  # by default, notify only if NA values are present
  if (is.null(notify_na_count)) {
    notify_na_count <- ifelse(na_count > 0, TRUE, FALSE)
  }
  if (notify_na_count == TRUE) {
    message(paste0(
      "\n", na_count,
      " observation(s) were removed due to missing values.\n"
    ))
  }

  # check if bin_cutoffs argument is null
  if (is.null(bin_cutoffs)) {
    stop(paste0(
      "Please set cutoff points for bins by entering a numeric vector ",
      "for bin_cutoffs"))
  }

  # bin_cutoffs per group
  groups <- unique(dat[[col_group]])
  dt <- data.frame(
    "bin_number" = numeric(),
    "bin_start" = numeric(),
    "bin_end" = numeric(),
    "count" = numeric(),
    "proportion" = numeric(),
    "group" = factor()
  )
  for(group in groups) {
    cutoffs <- bin_cutoffs

    # data for this group
    d <- dat |> filter(!!as.symbol(col_group) == group) |> select(col_name) |> pull()

    # do bin_cutoffs include min and max values?
    min_bin <- min(bin_cutoffs, na.rm = TRUE)
    max_bin <- max(bin_cutoffs, na.rm = TRUE)
    min_dat <- min(d)
    max_dat <- max(d)

    # if so, add
    if (max_dat > max_bin) {
      cutoffs <- c(cutoffs, max_dat)
    }
    if (min_dat < min_bin) {
      cutoffs <- c(min_dat, cutoffs)
    }

    # characteristics of the histogram
    bin_number <- utils::head(seq_along(cutoffs), - 1)
    bin_start <- utils::head(cutoffs, -1)
    bin_end <- utils::tail(cutoffs, -1)
    n_bins <- max(bin_number)

    # get count of each bin
    count <- vapply(bin_number, function(i) {
      if (i < n_bins) {
        sum(d >= bin_start[i] & d < bin_end[i])
      } else {
        sum(d >= bin_start[i] & d <= bin_end[i])
      }
    }, FUN.VALUE = numeric(1L))

    # get proportion of each bin
    proportion <- count / sum(count)

    # fill colors for bins
    fill_colors <- rep(non_outlier_bin_fill_color, length(bin_number))

    if (outlier_bin_left == TRUE) {
      fill_colors[1] <- outlier_bin_fill_color
    }
    if (outlier_bin_right == TRUE) {
      fill_colors[length(bin_number)] <- outlier_bin_fill_color
    }

    # create a data table
    dt <- rbind(dt,
                data.frame(bin_number,
                           bin_start,
                           bin_end,
                           count,
                           proportion,
                           fill_colors,
                           group ))
  }

  # rename
  names(dt) <- c("bin_number", "bin_start", "bin_end", "count", "proportion", "fill_colors", col_group)

  # ensure levels correct order
  dt[[col_group]] <- forcats::fct_relevel(dt[[col_group]], group_levels)

  # plot frequency or proportion? set a default
  if (plot_frequency == TRUE) {
    plot_proportion <- FALSE
    # message(paste0(
    #   "Plotting frequencies instead of proportions because ",
    #   "plot_frequency = TRUE"))
  }

  # y label
  y <- if(plot_proportion == TRUE){
    "proportion"
  } else if( plot_frequency == TRUE) {
    "count"
  }


  # plot
  g1 <- ggplot2::ggplot(
    data = dt,
    mapping = ggplot2::aes(x = bin_number,
                           y = .data[[y]])) +
    ggplot2::geom_bar(
      stat = "identity",
      color = border_color,
      fill  = dt$fill_colors,
      width = 1)

  # label axes
  if (!is.null(x_axis_title)) {
    if (x_axis_title == FALSE) {
      g1 <- g1 + ggplot2::theme(axis.title.x = ggplot2::element_blank())
    } else {
      g1 <- g1 + ggplot2::xlab(x_axis_title)
    }
  } else {
    g1 <- g1 + ggplot2::xlab("Value")
  }
  if (!is.null(y_axis_title)) {
    g1 <- g1 + ggplot2::ylab(y_axis_title)
  } else {
    g1 <- g1 + ggplot2::ylab(y)   #kim::capitalize(y))
  }

  # adjust x axis tick marks
  if (!is.null(x_tick_marks) & is.null(x_tick_mark_labels)) {
    message("Setting x_tick_mark_labels = x_tick_marks...")
    x_tick_mark_labels <- x_tick_marks
  }
  if (is.null(x_tick_marks)) {
    x_tick_marks <- seq(0.5, n_bins + 0.5, 1)
  }
  if (is.null(x_tick_mark_labels)) {
    x_tick_mark_labels <- bin_cutoffs
  }

  which_bin <- function(x_tick_mark, dt) {
    dt[nrow(dt),"bin_end"] <- dt[nrow(dt),"bin_end"] + 0.00001
    bin_idx <- which( x_tick_mark >= dt$bin_start & x_tick_mark < dt$bin_end )[1]
    return(bin_idx)
  }

  # sort out the x_breaks and labels
  x_breaks <- sapply(x_tick_marks, which_bin, dt=dt)
  if (outlier_bin_left == TRUE) {
    x_breaks[1] <- dt$bin_number[1]
    x_tick_mark_labels[1] <- round(dt$bin_start[1])
  }
  if (outlier_bin_right == TRUE) {
    x_breaks[length(x_breaks)] <- dt$bin_number[nrow(dt)]
    x_tick_mark_labels[length(x_tick_mark_labels)] <- round(dt$bin_end[nrow(dt)])
  }

  g1 <- g1 + ggplot2::scale_x_continuous(
    breaks = x_breaks,
    labels = x_tick_mark_labels
  )

  # update y tick marks
  if (!is.null(y_tick_marks)) {
    g1 <- g1 + ggplot2::scale_y_continuous(
      limits = c(
        min(y_tick_marks, na.rm = TRUE),
        max(y_tick_marks, na.rm = TRUE)),
      breaks = y_tick_marks)
  }

  return(g1)
}
