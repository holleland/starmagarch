% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation.R
\name{simSTARMAGARCH}
\alias{simSTARMAGARCH}
\title{Simulation of STARMAGARCH}
\usage{
simSTARMAGARCH(
  parameters,
  n = 100,
  m = c(5, 5),
  W = NULL,
  burnin = 0,
  type = "rook",
  torus = TRUE
)
}
\arguments{
\item{parameters}{List of parameters.}

\item{n}{Temporal dimension (integer).}

\item{m}{Spatial dimension vector of length 1 or 2 (integer).}

\item{W}{Neighbourhood array (optional)}

\item{burnin}{Temporal burnin.}

\item{type}{neighbourhood type (rook or queen)}

\item{torus}{Logical, indicator for circular space.}
}
\description{
Simulation of STARMAGARCH
}
