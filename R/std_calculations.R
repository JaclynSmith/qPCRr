linear_reg <- function(table) {
  #'Summarize data and generate linear regression to fit standards, report linear reg variables
  #'and the plot of the summarized data
  #'@param table a tibble with quantities and corresponding ct values
  #'@return list: (yint, slope, graph of standards)
  #'@export

  #check arguments---------------------------------------------------------------------
  checkmate::assert_tibble(table, types = "double",
                           any.missing = FALSE, ncols = 2)
  checkmate::assert_set_equal(names(table), c("quantity",
                                              "ct"))
  #get summary table--------------------------------------------------------------------
  sum_table <- table %>%
    dplyr::group_by(quantity) %>%
    dplyr::summarise(avg = mean(ct), stdev = sd(ct),
                     std_err = sd(ct) / sqrt(n())) %>%
    dplyr::mutate(log_conc = log(quantity))

  #generate linear regression--------------------------------------------------------
  linear_mod <- stats::lm(formula = avg ~ log_conc, data = sum_table)
  yint <- stats::coef(linear_mod)[[1]]
  slope <- stats::coef(linear_mod)[[2]]
  print(summary(linear_mod))

  #graph summary data and regression line----------------------------------------------------
  standards_plot <- ggplot2::ggplot() +
    ggplot2::geom_point(data = sum_table, aes(x = log_conc, y = avg)) +
    ggplot2::geom_abline(aes(intercept = yint, slope = slope))

  return(list(yint, slope, standards_plot))
}

conc_calc <- function(data_table, yint, slope, avg_length, dilution_fact) {
  #'Summarize data and generate linear regression to fit standards, report linear reg variables
  #'and the plot of the summarized data
  #'@param data_table a tibble with unknown samples, 3 columns: sample_name, target_name, and ct
  #'@param yint a numeric, the y intercept of the regression line
  #'@param slope a numeric, the slope of the regression line
  #'@param avg_length a numeric, the avg length of the molecules in the library
  #'@param dilution_fact a numeric, the dilution factor used for quantification ie for 1:1000, you put 1000
  #'@return a tibble with original library concentrations
  #'@export

  #check arguments---------------------------------------------------------------------
  checkmate::assert_tibble(data_table, types = c("character", "double"),
                           any.missing = FALSE, ncols = 3)
  checkmate::assert_subset(c("sample_name", "target_name", "ct"),
                           names(data_table))
  checkmate::assert_numeric(yint)
  checkmate::assert_numeric(slope)
  checkmate::assert_numeric(avg_length, lower = 0)
  checkmate::assert_numeric(dilution_fact, lower = 0)

  len_kapa_kit_DNA <- 452
  pg_to_ng <- 0.001

  #get summary table--------------------------------------------------------------------
  sum_table <- data_table %>%
    dplyr::group_by(sample_name, target_name) %>%
    dplyr::summarise(avg = mean(ct), stdev = sd(ct), std_err = sd(ct) / n()) %>%
    dplyr::mutate(dil_conc = (10 ^ ((avg - yint) / slope)) *
                    (len_kapa_kit_DNA / avg_length)) %>%
    dplyr::mutate(original_conc = dil_conc * dilution_fact * pg_to_ng)

  return(sum_table)
}
