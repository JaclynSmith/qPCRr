sum_tbl_calc <- function(table, color, house_gene, bars, replicates) {
  #'A function to calculate summary statistics
  #'@param table a tibble
  #'@param color a character or character list, a list of columns to use for color
  #'@param house_gene a character, gene to use as the housekeeping gene
  #'@param bars a character: "std_dev" or "std_error" type of error bars you will get on your graph
  #'@param replicates a character, "no replicates" or the list of columns to use for the replicates
  #'@return summary statistics table

  #group the summary_table by the colors, if given------------------------------------------------
   if ("none" %notin% color) {
     summary_table <- table %>%
       dplyr::group_by_at(dplyr::vars(c("target_name", color)))
   } else {
     summary_table <- table %>%
       dplyr::group_by_at(dplyr::vars("target_name"))
   }

  #calculate summary statistics for data given replicates expts, if provided
   if (replicates != "no_replicates") {
     summary_table <- summary_table %>%
       dplyr::summarise(avg = mean(fold_change), std_dev = stats::sd(fold_change),
                 std_error = sqrt(sum(dd_error^2)))
   } else {
     summary_table <- summary_table %>%
       dplyr::summarise(avg = fold_change, std_dev = 0, std_error = dd_error)
   }

  #label housekeeping gene and calculate ymin and ymax
   summary_table <- summary_table %>%
      dplyr::mutate(control = dplyr::case_when(
                    target_name == house_gene ~ "housekeeping",
                                         TRUE ~ "experimental")) %>%
      dplyr::mutate(ymin := avg - !!rlang::sym(bars),
                    ymax := avg + !!rlang::sym(bars))
  return(summary_table)
}


qPCR_graph <- function(table,
                       x = "target_name",
                       y = "fold_change",
                       house_gene,
                       color = "",
                       xlab = "",
                       ylab = "Fold Change\n",
                       colorlab = "Expt",
                       color_choices =  c("#336699", "#9999FF", "black", "red"),
                       ycoord_max = 5,
                       replicates = "no_replicates",
                       bars = "std_error",
                       hline = "yes",
                       coord_cart = "yes",
                       theme = "jaclyn") {
  #'A function to graph fold change qPCR data
  #'@param table a tibble, columns are expt conditions, fold change, error, target_name
  #'@param x a character, the column title for the variable to be plotted on the x axis, target_name
  #'@param y a character, the column title for the variable to be plotted on the y axis, fold_change
  #'@param house_gene a character, gene to use as the housekeeping gene
  #'@param color a character/character list, the columns to group into color
  #'@param xlab a character, the label for the x axis, default is ""
  #'@param ylab a character, the label for the y axis, default is "Fold Change\\n"
  #'@param colorlab a character, the label for the legend based on color, default is "Expt"
  #'@param color_choices, a character list, the colors you want to use for the graph
  #'@param ycoord_max a numeric, the max value for the fold change displayed on the graph
  #'@param replicates a character, "no replicates" or the list of columns to use for the replicates
  #'@param bars a character: "std_dev" or "std_error" type of error bars you will get on your graph
  #'@param hline a character: "yes" or "no", if yes, a horizontal line for fold change is 1 is added
  #'@param coord_cart a character: "yes" or "no", if no graph scaled to extent of data
  #'@param theme a character: "jaclyn" or "default", if default, the original theme applied and you can add your own theme later
  #'@return summary statistics table is printed and graph object is returned
  #'@export
  #This function calculates summary statistics and graphs data----------------------------------------
  checkmate::assert_tibble(table, types = c("factor", "character", "double"),
                           any.missing = FALSE, min.cols = 4)
  checkmate::assert_subset(c("target_name", "fold_change", "dd_error"),
                           names(table))
  checkmate::assert_subset(c(color), c(names(table), "none"))

  checkmate::assert(
    checkmate::check_character(replicates),
    checkmate::check_subset(replicates, c("no_replicates", names(table))),
    combine = "and"
  )

  checkmate::assert(
    checkmate::check_character(hline),
    checkmate::check_subset(hline, c("yes", "no")),
    combine = "and"
  )

  checkmate::assert(
    checkmate::check_character(bars),
    checkmate::check_subset(bars, c("std_error", "std_dev")),
    combine = "and"
  )

  checkmate::assert(
    checkmate::check_character(hline),
    checkmate::check_subset(hline, c("yes", "no")),
    combine = "and"
  )

  checkmate::assert(
    checkmate::check_character(coord_cart),
    checkmate::check_subset(coord_cart, c("yes", "no")),
    combine = "and"
  )

  checkmate::assert(
    checkmate::check_character(theme),
    checkmate::check_subset(theme, c("jaclyn", "default")),
    combine = "and"
  )

  checkmate::assert_character(x)
  checkmate::assert_character(y)
  checkmate::assert_character(house_gene)
  checkmate::assert_character(color)
  checkmate::assert_character(xlab)
  checkmate::assert_character(ylab)
  checkmate::assert_character(colorlab)
  checkmate::assert_character(color_choices)
  checkmate::assert_numeric(ycoord_max)

  #create summary statistics
  summary_table <- sum_tbl_calc(table, color, house_gene, bars, replicates)
  print("Summary table:")
  print(summary_table)

  replicates <- dplyr::setdiff(names(table),
                               c(color, "target_name", "fold_change", "dd_error"))
  print(paste0(paste0(replicates, collapse =  ", "),
               " - used for the replicates."))

  #graph the plot
  if ("none" %notin% color) {
    total_plot <- ggplot2::ggplot() +
      ggplot2::geom_point(data = table,
                 ggplot2::aes_string(x = x, y = y,
                                     color = paste0("interaction(", paste0(color, collapse =  ", "), ")")),
                 stat = "identity",
                 position = ggplot2::position_dodge(0.5)) +
      ggplot2::labs(x = xlab, y = ylab, color = colorlab) +
      ggplot2::geom_crossbar(data = summary_table,
                    ggplot2::aes_string(x = x,
                               y = "avg",
                               ymin = "avg",
                               ymax = "avg",
                               color = paste0("interaction(", paste0(color, collapse =  ", "), ")")),
                    width = 0.5,
                    position = ggplot2::position_dodge(0.5),
                    fatten = 2.5) +
      ggplot2::geom_errorbar(data = summary_table,
                             ggplot2::aes_string(x = x,
                                                 ymin = "ymin",
                                                 ymax = "ymax",
                                                 color = paste0("interaction(", paste0(color, collapse =  ", "),
                                                                ")")),
                             width = 0.6,
                             position = ggplot2::position_dodge(0.5))
  } else {
    total_plot <- ggplot2::ggplot() +
      ggplot2::geom_point(data = table,
                          ggplot2::aes_string(x = x, y = y),
                          stat = "identity") +
      ggplot2::labs(x = xlab, y = ylab) +
      ggplot2::geom_crossbar(data = summary_table,
                             ggplot2::aes_string(x = "target_name",
                                                 y = "avg",
                                                 ymin = "avg",
                                                 ymax = "avg"),
                             width = 0.5,
                             fatten = 2.5) +
      ggplot2::geom_errorbar(data = summary_table,
                             ggplot2::aes_string(x = "target_name",
                                                 ymin = "ymin",
                                                 ymax = "ymax"),
                             width = 0.6)
  }
  total_plot <- total_plot +
    ggplot2::scale_color_manual(values = color_choices)

  if (theme != "default") {
    total_plot <- total_plot +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 20,
                                                         angle = 90,
                                                         hjust = 1,
                                                         vjust = 0.35),
            axis.text.y = ggplot2::element_text(size = 20),
            axis.title.y = ggplot2::element_text(size = 25),
            plot.background = ggplot2::element_blank(),
            panel.background = ggplot2::element_blank(),
            panel.border = ggplot2::element_blank(),
            axis.line = ggplot2::element_line(color = "black", size = 2),
            axis.ticks.length = ggplot2::unit(0.25, "cm"),
            legend.position = "right",
            strip.text.x = ggplot2::element_text(size = 20),
            legend.title = ggplot2::element_text(size = 15,
                                                 face = "bold",
                                                 hjust = 0.5),
            legend.key = ggplot2::element_rect(fill = "white"),
            legend.text = ggplot2::element_text(size = 15))
  }

  if (hline == "yes") {
    total_plot <- total_plot +
      ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "gray")
  }

  if (coord_cart == "yes") {
    total_plot <- total_plot +
      ggplot2::coord_cartesian(ylim = c(0, ycoord_max))
  }

  return(total_plot)
}
