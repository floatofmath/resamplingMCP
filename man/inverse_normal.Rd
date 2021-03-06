% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/utilities.R
\name{inverse_normal}
\alias{inverse_normal}
\title{Inverse normal combination function}
\usage{
inverse_normal(p1, p2, w1, w2)
}
\arguments{
\item{p1}{First stage p-value}

\item{p2}{Second stage p-value}

\item{w1}{First stage weight}

\item{w2}{Second stage weight}
}
\value{
p-value corresponding to the inverse normal combination z-score
}
\description{
Computes the inverse normal combination (sqrt(w1)*qnorm(1-p1) + sqrt(w2)*qnorm(1-p2)) of two (independent) p-values
}
\author{
Florian Klinglmueller
}

