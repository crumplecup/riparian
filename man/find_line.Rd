% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{find_line}
\alias{find_line}
\title{find line}
\usage{
find_line(pt, lines = prc_strms)
}
\arguments{
\item{pt}{is a coordinate pair (x,y)}

\item{lines}{is a SpatialLinesDataFrame}
}
\value{
coordinate matrix of nearest line to `pt` in `lines`
}
\description{
finds the nearest line in a sp object, returns coordinate matrix
}
