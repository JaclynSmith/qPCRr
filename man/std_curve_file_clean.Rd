% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrangling_fun_std_curve.R
\name{std_curve_file_clean}
\alias{std_curve_file_clean}
\title{Read in qPCR files, clean them and report on controls}
\usage{
std_curve_file_clean(
  file_location,
  file_name,
  undetermined = "set40",
  std_undetermined = "remove"
)
}
\arguments{
\item{file_location}{a character, path to the file, name not included, must have trailing /}

\item{file_name}{a character, file name with extension}

\item{undetermined}{a character, "set40" or "remove"; set40 takes undetermined values and sets them to 40, otherwise they are removed. default is set40}
}
\value{
a list; (tibble with cleaned data, tibble with standards),plus printed output on controls
}
\description{
Read in qPCR files, clean them and report on controls
}
