clean_qPCR_files_check_NTC <- function(file_location, file_name, undetermined = "set40") {
  #'Read in qPCR files, clean them and report on controls
  #'@param file_location a character, path to the file, name not included, must have trailing /
  #'@param file_name a character, file name with extension
  #'@param undetermined a character, "set40" or "remove"; set40 takes undetermined values and sets them to 40, otherwise they are removed. default is set40
  #'@return a tibble with cleaned data, plus printed output on controls
  #'@export
  # this function reads in qPCR files, cleans them and reports on controls.
  # All samples must be properly assigned on the StepOnePlus machine in order
  # for this to work. Controls are water based No Template Controls (NTC) and
  # DNase treated RNA.
  # usage cleanqPCRFiles_checkNTC(file location, fileName)
  # fileLocation = path to file, must have trailing /
  # fileName = the name of your file with the extension
  # undetermined = deault: "set40"; other option: "remove"
  # "set40" - sets the qPCR Value to 40
  # "remove" - removes sample from dataSet
  # output cleaned qPCR file with unknowns only, easy to use names at the top,
  # and will analyze and report on issues in controls.

  checkmate::assert_character(file_location)
  checkmate::assert_character(file_name)
  checkmate::assert_subset(undetermined, c("set40", "remove"))
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
    dplyr::select(Well:Task, "ct" = dplyr::matches("^C.$")) %>%
    dplyr::rename(sample_name = 'Sample Name') %>%
    dplyr::rename(target_name = 'Target Name')

  # create final cleaned file------------------------------------------------------------------
  ct_values_return <- ct_values %>%
    dplyr::filter(Task == "UNKNOWN") %>%
    dplyr::select(sample_name, target_name, ct)

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
  return(ct_values_return)
}

calc_fold_change <- function(table, house_gene, control,
                             expt_conditions = c(""), sep = "",
                             full_table = "no") {
  #'Calculates the fold change using deltadeltaCt
  #'@param table a tibble with columns sample_name, target_name, and ct
  #'@param house_gene a character, the target_name which is the housekeeping gene
  #'@param control a character or character list, conditions which represent the control condition
  #'@param expt_conditions a character list, conditions which make up sample_name, will become titles of new columns
  #'@param sep a character, the separation character used in sample_name to separate conditions
  #'@param full_table a character, "yes" or "no", if yes, will return all intermediate statistics
  #'@return a tibble with fold changes calculated based on house_gene and control conditions
  #'@export
  #This function takes a clean data table with columns sample_name, target_name, and ct.
  #It calculates fold change from these values based on a housekeeping
  #gene and control conditions set by input parameters. After calculations, if specified,
  #the function will split the sample_name column into various experimental conditions.

  #table: the cleaned qPCR table, in tibble format,
  #with columns sample_name, target_name, and ct
  #house_gene: string, the housekeeping gene
  #control: string, control condition found in sample_name OR list of control conditions
  #expt_condition: empty string list or list of strings, this will be used to separate
  #the sample_name column into experimental conditions based on a separator
  #sep: the separator used to name your conditions in your sample_name
  #full_table: if no, will return a compact table, if yes will return all intermediate
  #qPCR calculations

  #test table in proper format-----------------------------------------------------------------------
  checkmate::assert_tibble(table, types = c("character", "double"),
                           any.missing = FALSE, ncols = 3)
  checkmate::assert_set_equal(names(table), c("sample_name",
                                              "target_name",
                                              "ct"))

  #test full_table is either yes or no------------------------------------------------------------------------
  checkmate::assert_subset(full_table, c("yes", "no"))

  #test that sep is type character and is within sample_name the proper number of times------------------------------
  checkmate::assert(
    checkmate::check_character(sep),
    check_sep(table, expt_conditions, sep),
    combine = "and"
  )

  #assert experimental conditions is a list of type character--------------------------------------------------------
  checkmate::assert_character(expt_conditions)

  #assert that house_gene is within column target_name and is a character------------------------------------------
  checkmate::assert(
    checkmate::check_character(house_gene),
    check_house_gene(table, house_gene),
    combine = "and"
  )

  #assert that control is in the sample_name column
  checkmate::assert(
    checkmate::check_character(control),
    check_sample_control(table, control, sep),
    combine = "and"
  )

  #calculate fold change-----------------------------------------------------------------------------------------
  dd_ct_table <- table %>%
    dplyr::group_by(sample_name, target_name) %>%
    dplyr::summarise(avg = mean(ct),
                     std_error = stats::sd(ct) / sqrt(dplyr::n())) %>%
    dplyr::mutate(d_ct = avg - avg[target_name == house_gene],
                  d_error = sqrt(std_error^2 + std_error[target_name == house_gene]^2)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(target_name) %>%
    dplyr::mutate(ddct = d_ct - d_ct[sample_name == sample_control],
                  dd_error = sqrt(d_error^2 + d_error[sample_name == sample_control]^2)) %>%
    dplyr::mutate(fold_change = 2^-ddct)

  #cleanup table if full_table = "no"-------------------------------------------------------------------
  if (full_table == "no") {
    dd_ct_table <- dd_ct_table %>%
      dplyr::select(sample_name, target_name, fold_change, dd_error)
  }

  #break sample_name into multiple columns with experimental conditions-----------------------------------
  if (expt_conditions[1] != "") {
    dd_ct_table <- dd_ct_table %>%
      tidyr::separate(sample_name, expt_conditions, sep = sep)
  }

  #return table-------------------------------------------------------------------------------------------
  return(dd_ct_table)
}

reformat_table <- function(table, sample, target, ct, sep = "_") {
  #'reformat other data tables into one which can be used by calc_fold_change
  #'@param table a data frame or tibble
  #'@param sample a character, the column, or list of columns which make up sample_name
  #'@param target a character, the column which makes up the targeted genes
  #'@param ct a character, the column which makes up the ct values of entries
  #'@param sep a character, the character used to combine columns together
  #'@return a tibble in the proper format for calc_fold_change
  #'@export
  #tests for proper, usable inputs-----------------------------------------------------------------------------
  checkmate::assert_data_frame(table, types = c("factor", "character", "double"),
                               any.missing = FALSE, min.cols = 3)
  checkmate::assert_subset(c(sample, target, ct), names(table))
  checkmate::assert(
    check_ct_numeric(table, ct)
  )
  df <- table
  #tests if data frame is a tibble and returns a tibble if not----------------------------------------------------------------
  if (tibble::is_tibble(table)) {
  } else if (is.data.frame(table)) {
    df <- dplyr::as_tibble(table)
  }

  #Reformat the table into the standard form needed to calculate fold change-----------------------------------------------
  reformat_table <- df %>%
    dplyr::ungroup() %>%
    dplyr::mutate_if(is.factor, as.character) %>%
    dplyr::mutate_at(ct, as.numeric) %>%
    tidyr::unite(sample_name, sample, sep = "_") %>%
    tidyr::unite(target_name, target, sep = "_") %>%
    dplyr::rename(ct = ct) %>%
    dplyr::select(sample_name, target_name, ct)

  #summarise the table to determine number of samples and targets-------------------------------------------------------------
  tbl_summary <- reformat_table %>%
    dplyr::summarise(sample_name = dplyr::n_distinct(sample_name),
                     target_name = dplyr::n_distinct(target_name))

  #report results of sample and target summaries -----------------------------------------------------------------------------------
  print(paste0("There are ", dplyr::pull(tbl_summary, sample_name),
               " samples in final data frame: ",
               paste(unique(dplyr::pull(reformat_table, sample_name)),
                     collapse = ", "), "."))
  print(paste0("There are ", dplyr::pull(tbl_summary, target_name),
               " targets in final data frame.",
               paste(unique(dplyr::pull(reformat_table, target_name)),
                     collapse = ", "), "."))
  print("Remember: all samples from the same source must have the same name.")

  return(reformat_table)
}

