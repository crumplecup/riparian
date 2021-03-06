% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{centerline}
\alias{centerline}
\title{centerline}
\usage{
centerline(find, span, inc, strms, thresh)
}
\arguments{
\item{strms}{is the stream layer `pts` lie along}

\item{pts}{is a spatial points object (stream samples)}
}
\value{
a list of matrices with coordinates for 25 segments along stream path
}
\description{
Given sampling points, find the centerline points of the sampling box.
}
