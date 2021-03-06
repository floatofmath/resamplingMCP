% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/simulation.R
\name{adaptive_npcombtest_os}
\alias{adaptive_npcombtest_os}
\title{adptive NPC test}
\usage{
adaptive_npcombtest_os(x, n1, n, ne, test_statistic,
  combination_function = inverse_normal, perms = 50000, alpha = 0.025)
}
\arguments{
\item{x}{Observations}

\item{n1}{First stage sample size}

\item{n}{Pre-planned total sample size}

\item{ne}{Extended sample size}

\item{test_statistic}{Function that computes the test test statistic}

\item{combination_function}{Function to combine stage-wise (permutation) p-values}

\item{perms}{Maximum number of permutations to use when computing permutation p-values and conditional error rates}

\item{alpha}{}
}
\value{
Logical whether the test rejects
}
\description{
Non-parametric combination of stage-wise test statistics. Combines stage-wise permutation p-values using some combination function; performs the test using the joint conditional permutation distribution of stage wise permutation p-values.
}
\author{
Florian Klinglmueller
}

