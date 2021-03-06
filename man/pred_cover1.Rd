% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{pred_cover1}
\alias{pred_cover1}
\title{plot predicted cover prime}
\usage{
pred_cover1(
  poly,
  in_dir = NULL,
  out_dir = NULL,
  rgb_dir = NULL,
  cir_dir = NULL,
  buffer = sf_buff
)
}
\arguments{
\item{poly}{is a spatial polygons object (taxlot or sampling box)}

\item{in_path}{is a character string indicating the path to a directory to rasters}

\item{out_path}{is a character string indicating the path to the output directory}

\item{rgb_path}{is the path to 3-band rgb data}

\item{cir_path}{is the path to 3-band cir data}

\item{buff}{is the 50-ft riparian buffer polygon for `poly`}
}
\value{
a raster plot of predicted cover masked by `poly`
}
\description{
Predicts cover extent within a polygon and returns the results as a raster plot
}
