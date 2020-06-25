std_curve_file_clean <- function(file_location, file_name, undetermined = "set40", std_undetermined = "remove") {
  #'Read in qPCR files, clean them and report on controls
  #'@param file_location a character, path to the file, name not included, must have trailing /
  #'@param file_name a character, file name with extension
  #'@param undetermined a character, "set40" or "remove"; set40 takes undetermined values and sets them to 40, otherwise they are removed. default is set40
  #'@return a list; (tibble with cleaned data, tibble with standards),plus printed output on controls
  #'@export
  # this function reads in qPCR files, cleans them and reports on controls.
  # All samples must be properly assigned on the StepOnePlus machine in order
  # for this to work. Controls are water based No Template Controls (NTC) and
  # DNase treated RNA. Standards are also included in this version of the function.
  # output cleaned qPCR file with unknowns only, easy to use names at the top,
  # and will analyze and report on issues in controls.

  checkmate::assert_character(file_location)
  checkmate::assert_character(file_name)
  checkmate::assert_subset(undetermined, c("set40", "remove"))
  checkmate::assert_subset(std_undetermined, c("set0", "remove"))
  checkmate::assert_directory_exists(file_location, access = "rw")

  # read in file-----------------------------------------------------------------------------
  path_to_file <- paste0(file_location, file_name)
  checkmate::assert_file_exists(path_to_file, access = "rw")

  ct_values <- readxl::read_xls(path_to_file, "Results")

  # change titles, remove extraneous columns--------------------------------------------------
  colnames(ct_values) <-
    ct_values[stringr::str_detect(string = tidyr::replace_na(dplyr::pull(ct_values, `Block Type`), ""),
                                  pattern = "Well"), ]
  ct_values <- ct_values %>%
    dplyr::select(Well:Task, "ct" = dplyr::matches("^C.$"), Quantity) %>%
    dplyr::rename(sample_name = 'Sample Name') %>%
    dplyr::rename(target_name = 'Target Name') %>%
    dplyr::rename(quantity = Quantity)

  # create final cleaned file------------------------------------------------------------------
  ct_values_return <- ct_values %>%
    dplyr::filter(Task == "UNKNOWN") %>%
    dplyr::select(sample_name, target_name, ct)
  ct_values_standards <- ct_values %>%
    dplyr::filter(Task == "STANDARD") %>%
    dplyr::select(quantity, ct)

  # this part deals with Undetermined values in the UNKNOWN data------------------------------
  # input default is set40, this sets the undetermined values to 40--------------------------
  if (undetermined == "set40"
      & "Undetermined" %in% dplyr::pull(ct_values_return, ct)) {
    ct_values_return <- ct_values_return %>%
      dplyr::mutate(ct = dplyr::case_when(
        ct == "Undetermined" ~ "40",
                        TRUE ~ ct)) %>%
      dplyr::mutate_at(dplyr::vars(ct), as.numeric)
    print("Warning: Undetermined set to 40")
    # if input is remove instead, it will eliminate undetermined values
  } else if (undetermined == "remove"
             & "Undetermined" %in% dplyr::pull(ct_values_return, ct)) {
    ct_values_return <- ct_values_return %>%
      dplyr::filter(ct != "Undetermined") %>%
      dplyr::mutate_at(dplyr::vars(ct), as.numeric)
    print("Warning: Undetermined removed")
    # if no values are undetermined, this converts them to a numeric
  } else{
    print("No 'Undetermined' detected in experimental samples.")
    ct_values_return <- ct_values_return %>%
      dplyr::mutate_at(dplyr::vars(ct), as.numeric)
  }

  # this part deals with Undetermined values in the STANDARD data------------------------------
  # input default is remove, this removes undetermined entries--------------------------
  if (std_undetermined == "set0"
      & "Undetermined" %in% dplyr::pull(ct_values_standards, ct)) {
    ct_values_standards <- ct_values_standards %>%
      dplyr::mutate(ct = dplyr::case_when(
        ct == "Undetermined" ~ "0",
                        TRUE ~ ct)) %>%
      dplyr::mutate_at(dplyr::vars(c(quantity, ct)), as.numeric)
    print("Warning: Undetermined Standard set to 0")
    # if input is remove instead, it will eliminate undetermined values
  } else if (std_undetermined == "remove"
             & "Undetermined" %in% dplyr::pull(ct_values_standards, ct)) {
    ct_values_standards <- ct_values_standards %>%
      dplyr::filter(ct != "Undetermined") %>%
      dplyr::mutate_at(dplyr::vars(c(quantity, ct)), as.numeric)
    print("Warning: Undetermined Standard removed")
    # if no values are undetermined, this converts them to a numeric
  } else{
    print("No 'Undetermined' detected in standards.")
    ct_values_standards <- ct_values_standards %>%
      dplyr::mutate_at(dplyr::vars(c(quantity,ct)), as.numeric)
  }

  # Test the remaining controls for amplification--------------------------------------------------
  # isloate the NTCs with amplification, replace NA values with n
  controls <- ct_values %>%
    dplyr::filter(Task == "NTC") %>%
    dplyr::filter(ct != "Undetermined") %>%
    replace(is.na(.), "n")

  # The following code deals with messages to the user to indicate the status of-------------------
  # controls---------------------------------------------------------------------------------------
  # no controls have a Ct value
  if (nrow(controls) == 0) {
    print("Congrats! No amplification in your NTC & your DNase treatment looks good.")
  } else {
    # NTCs with water have amplified expression
    if (nrow(controls) > 0 & "n" %in% dplyr::pull(controls, sample_name)) {
      # controls has more than 1 row and no assignment for sampleName (ie water)
      summary <- controls %>%
        dplyr::mutate_at(dplyr::vars(ct), as.numeric) %>%
        dplyr::group_by(target_name) %>%
        dplyr::filter(sample_name == "n") %>%
        dplyr::summarise(avg_ct = mean(ct), count = dplyr::n())
      print(paste0("Warning! NTC has amplification: ",
                   paste(dplyr::pull(summary, target_name),
                         collapse = ", ")))
      print(summary)
    }
    # DNase treated NTCs have DNA still
    if (nrow(controls) > 0 & !("n" %in% dplyr::pull(controls, sample_name))) {
      # controls has more than 1 row and some assignment for sampleName (not water)
      summary <- controls %>%
        dplyr::mutate_at(dplyr::vars(ct), as.numeric) %>%
        dplyr::group_by(sample_name, target_name) %>%
        dplyr::filter(sample_name != "n") %>%
        dplyr::summarise(avg_ct = mean(ct), count = dplyr::n())

      print(paste0("Warning: DNase treatment unsuccessful: ",
                   paste(dplyr::pull(summary, sample_name),
                         collapse = ", ")))
      print(summary)
    }
  }
  # finally, return cleaned data set with no "Undetermined" values -------------------------------
  return(list(ct_values_return, ct_values_standards))
}

remove_outliers <- function(table) {
  #'Remove outliers greater than 0.2 ct values from the triplicate
  #'@param table a tibble with quantities and corresponding ct values
  #'@return a tibble; if range is less than 0.2 return all, drop outliers
  #'@export
  # This function removes outliers from data.

  #check arguments---------------------------------------------------------------------
  checkmate::assert_tibble(table, types = "double",
                           any.missing = FALSE, ncols = 2)
  checkmate::assert_set_equal(names(table), c("quantity",
                                              "ct"))
  #calculate spread--------------------------------------------------------------------
  range_table <- table %>%
    dplyr::group_by(quantity) %>%
    dplyr::summarise(max_val = max(ct),
              min_val = min(ct)) %>%
    dplyr::mutate(data_spread = max_val - min_val)

  #calculate summary statistics--------------------------------------------------------
  summary_table <- table %>%
    dplyr::group_by(quantity) %>%
    dplyr::mutate(avg = mean(ct)) %>%
    dplyr::mutate(dif = abs(avg - ct))

  summary_diff <- summary_table %>%
    dplyr::summarise(max_dif = max(dif))

  #combine tibbles together for easy summarization--------------------------------------
  combined_table <- dplyr::full_join(summary_table,
                                     range_table, by = "quantity")
  combined_table <- dplyr::full_join(combined_table,
                                     summary_diff, by = "quantity")

  outlier <- 0.2

  trimmed_table <- combined_table %>%
    dplyr::group_by(quantity) %>%
    dplyr::filter(data_spread <= outlier | 
                  data_spread > outlier & dif != max_dif) %>%
    dplyr::select(quantity, ct)

  return(trimmed_table)
}
