addl_expt_attribute <- function(table, column, contents) {
  #'Create a new column with additional experimental information.
  #'@param table A tibble
  #'@param column A character, the new column name
  #'@param contents A character, the new column contents
  #'@return A tibble with the appended column
  #'@examples
  #'mod_iris <- addl_expt_attribute(dplyr::as_tibble(iris), "flower_type", "iris")
  #'@export
  #tests for proper tibble--------------------------------------------------------------------------------------------
  checkmate::assert_tibble(table)

  #tests that title is a string and not already in the table--------------------------------------------------------------
  checkmate::assert(
    checkmate::check_character(column),
    check_not_col_name(table, column),
    combine = "and"
  )

  #This will add an additional attribute to the table-----------------------------------------------------------------
  table <- table %>%
    dplyr::mutate(column = contents) %>%
    dplyr::rename(!!column := column)

  return(table)
}

reorder_samples <- function(table, column, list) {
  #'Order the entries in a column so that they are displayed in this order for graphs, rather than alphabetically.
  #'@param table A tibble
  #'@param column A character, the title of a column of table to be ordered
  #'@param list A character list, a list of column entries in the correct order for graphing
  #'@return Table is returned with a reordered column stored as a type, factor.
  #'@note list is a set of characters in the order they should appear in graphs. All entires in the list must be unique and all values must appear in the list.
  #'@export
  #tests for proper tibble and column name found within tibble-------------------------------------------------------------------
  checkmate::assert_tibble(table, any.missing = FALSE)
  checkmate::assert_subset(column, names(table))

  #test that the list is equal to unique entries in the column----------------------------------
  checkmate::assert(
    checkmate::check_character(list, unique = TRUE,
                               len = length(unique(dplyr::pull(table, column)))),
    checkmate::check_subset(list, dplyr::pull(table, column)),
    combine = "and"
  )

  #change column specified from character into a factor type, ordered by list-------------------------------------------
  reordered_table <- table %>%
    dplyr::ungroup() %>%
    dplyr::mutate_at(column, as.factor) %>%
    dplyr::mutate(!!column := forcats::fct_relevel(!!rlang::sym(column), list))
  return(reordered_table)
}
