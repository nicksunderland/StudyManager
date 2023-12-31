% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_col_map.R
\name{set_active}
\alias{set_active}
\alias{set_active,list,character-method}
\alias{set_active,ColMap,character-method}
\title{set_active}
\usage{
set_active(x, col_names, rest.off = TRUE)

\S4method{set_active}{list,character}(x, col_names, rest.off = TRUE)

\S4method{set_active}{ColMap,character}(x, col_names, rest.off = TRUE)
}
\arguments{
\item{x}{a Colmap object}

\item{col_names}{the names to activate}

\item{rest.off}{whether to turn the rest of the columns off}
}
\value{
a Colmap object
}
\description{
set_active
}
