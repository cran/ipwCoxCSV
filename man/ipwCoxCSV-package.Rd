\name{ipwCoxCSV-package}
\alias{ipwCoxCSV-package}
\docType{package}
\title{
Inference of Marginal Hazard Ratios (HR) in Inverse Probability Weighted (IPW) Cox Model Using Corrected Sandwich Variance (CSV)
}
\description{
This package is an implementation of the corrected sandwich variance (CSV) estimation method for making inference of marginal hazard ratios (HR) in inverse probability weighted (IPW) Cox model without and with clustered data, under both the conventional inverse probability weights and the stabilized weights, proposed by Shu, Young, Toh, and Wang (2019). Logistic regression model is assumed for propensity score model.
}
\details{
The \code{ipwCoxCSV} package implements the corrected sandwich variance estimation method for making inference of marginal hazard ratios in inverse probability weighted Cox model without and with clustered data, under both the conventional inverse probability weights and the stabilized weights, proposed by Shu, Young, Toh, and Wang (2019). This sandwich type variance estimation is referred to as corrected sandwich variance estimation, because it takes into account the uncertainty in weight estimation. The function \code{\link[ipwCoxCSV]{ipwCoxInd}} implements the corrected sandwich variance estimation method without clustered data (i.e., assuming independence among observations). The function \code{\link[ipwCoxCSV]{ipwCoxCluster}} implements the corrected sandwich variance estimation method with clustered data (i.e., allowing for within-cluster correlation). 
}
\author{
Di Shu and Rui Wang

Maintainer: Di Shu <shudi1991@gmail.com>
}
\references{
Shu D, Young JG, Toh S, Wang R (2019). Variance estimation in inverse probability weighted Cox model.  \emph{Biometrics,} under revision.
}
\keyword{ package }
