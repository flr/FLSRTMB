% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{from_logitd}
\alias{from_logitd}
\title{from_logitd()}
\usage{
from_logitd(logit_d, ll = 0.5, ul = 3)
}
\arguments{
\item{logit_d}{logit(depensation)}

\item{ll}{defines lower prior bound of d}

\item{ul}{defines upper prior bound of d}
}
\value{
depensation d
}
\description{
convert depensation d from logit
}
