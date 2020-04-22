check_sep <- function(table, expt_conditions = c(""), sep = "") {
  #'checkmate function to check for the separater
  #'@param table a tibble
  #'@param expt_conditions a character list, a list of columns to combine
  #'@param sep a character, the separator used to combine all the columns
  #'@return true if ok, an error message if false
  #a checkmate function
  #test separator found in each entry in the sample_name column----------------------------------------------
  #is sep in the column sample_name?
  sep_test <- stringr::str_detect(dplyr::pull(table, sample_name), sep)
  #is how many sep are in each entry?
  sep_count <- stringr::str_count(dplyr::pull(table, sample_name), sep)
  #how many times should sep be in each entry?
  len <- length(expt_conditions) - 1
  #are the number of seps in each entry the same or different
  uni_sep <- length(unique(sep_count))

  #combine all the possible separator mistakes and return appropriate errors------------------------------------------
  if (FALSE %in% sep_test) {
    return(paste0("Separater ", sep,
                  " not in every entry in the sample_name column"))
  } else if (sep != "" & uni_sep == 1 & TRUE %in% (sep_count != len)) {
    return(paste0("The number of columns supplied by expt_conditions = ",
                  paste(expt_conditions, collapse = ", "),
                  " does not match number of separators, ",
                  sep))
  } else if (sep != "" & TRUE %in% (sep_count != len)) {
    return(paste0("Separater ", sep,
                  " is not uniformly in all entries in the sample_name column"))
  } else {
    return(TRUE)
  }
}

check_house_gene <- function(table, house_gene) {
  #'checkmate function to check the housekeeping gene
  #'@param table a tibble
  #'@param house_gene a character, the house_keeping gene
  #'@return true if ok, an error message if false
  #a checkmate function
  if (house_gene %notin% dplyr::pull(table, target_name)) {
    return("house_gene: ", house_gene,
           " is not in target_name column. Check you spelled it correctly")
  } else {
    return(TRUE)
  }
}

check_sample_control <- function(table, control, sep) {
  #'checkmate function to check for control condition in the table
  #'@param table a tibble
  #'@param control a character or character list, factors which make up the control condition
  #'@param sep a character, the separator used to combine the control conditions
  #'@return true if ok, an error message if false
  #a checkmate function
  #construct the current sample_name for the control based on user input-----------------------------------------
  sample_control <- paste0(control, collapse = sep)

  #test sample_control is accurately in sample_name--------------------------------------------------------------
  if (sample_control %notin% dplyr::pull(table, sample_name)) {
    return("Control condition is ", sample_control,
           " and is not found in sample_name")
  } else {
    return(TRUE)
  }
}

check_ct_numeric <- function(table, ct) {
  #'checkmate function to check the ct column has proper numeric format
  #'@param table a tibble
  #'@param ct a character, the column title ct
  #'@return true if ok, an error message if false
  #a checkmate function
  ct_numeric_test <- stringr::str_detect(dplyr::pull(table, ct),
                                         "^(([0-9]+(\\.[0-9]*)?)|(\\.[0-9]+))$")
  if (FALSE %in% ct_numeric_test) {
    return("One or more entry in ct is not a number. All values must be real, positive numbers.")
  } else {
    return(TRUE)
  }
}

#checkmate check - tests that a column name is not already in the table---------------------------------------------
check_not_col_name <- function(table, column) {
  #'checkmate function to check that a column name is not already in the column
  #'@param table a tibble
  #'@param column a character, the new column title
  #'@return true if ok, an error message if false
  #a checkmate function
  if (column %in% names(table)) {
    return(paste0(column, " must not be a new column name, not a previous column name"))
  } else {
    return(TRUE)
  }
}

check_replicates <- function(table, replicates) {
  #'checkmate function to check for the replicates
  #'@param table a tibble
  #'@param replicates a character, tests for the replicates being in table from multiple expts
  #'@return true if ok, an error message if false
  #a checkmate function
  #check that no duplicates is actually no duplicates
  #returns columns not used in any way and also says you can say yes and duplicates will be averaged
  #check that there are duplicates?
  #must specify duplicates as yes, or by a specific column,
  if (replicates == "no_duplicates") {

    return(TRUE)
  } else if (replicates %notin% names(table)) {
    return("")
  }
}

are_colors <- function(x) {
  #'function to check that colors are viable
  #'@param x list of colors
  #'@return list of T/F, TRUE if a real color and False if not
  #joshObrien
  #https://stackoverflow.com/questions/13289009/check-if-character-string-is-a-valid-color-representation
  sapply(x, function(X) {
    tryCatch(is.matrix(grDevices::col2rgb(X)),
             error = function(e) FALSE)
  })
}

check_color <- function(table, color, color_choices) {
  #'checkmate function to check for proper colors
  #'@param table a tibble
  #'@param color a character list, a list of columns which to assign the colors by
  #'@param color_choices a character list, list of the colors
  #'@return true if ok, an error message if false
  if ("none" %notin% color) {
    color_test <- table %>%
      dplyr::group_by_at(color)
    colors_needed <- dplyr::n_groups(color_test)
  }
  if (length(color_choices) != colors_needed) {
    return("Number of colors supplied doesn't match the number of groups created. Please provide more colors")
  } else if (FALSE %in% are_colors(color_choices)) {
    return("At least one of the colors supplied is not a valid color. Please supply valid colors")
  } else {
    return(TRUE)
  }
}
