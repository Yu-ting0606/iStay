\name{iStay-package}
\alias{iStay-package}
\alias{iStay}
\docType{package}
\title{
  Measuring ecological stability and synchrony.
}
\description{
  \code{iStay} (information-based measures of stability and related synchrony) is an R package for featuring a class of measures of community/ecosystem stability and related synchrony; the corresponding diversity-stability is also assessed if biodiversity data are available. The information-based measures of stability were developed based on the generalized entropy parametrized by an order q > 0; see Chao et al. (2024) for the relevant theory and methodological background. All measures are illustrated by using biomass temporal data from the Jena experiment (Wagg et al. 2022).  \cr

Specifically, iStay features the following measures separately for three cases:

(1) For a single community, it computes stability measures of order q > 0 and provides the stability profile which depicts the stability with respect to order q. When biodiversity data are available, the diversity-stability relationship is also assessed.

(2) For multiple communities, it computes four measures (gamma, alpha and beta stability and synchrony) and provides the corresponding profiles. When biodiversity data are available, the diversity-stability and diversity-synchrony relationships are also assessed.

(3) For hierarchies, it computes three measures (gamma, alpha, and normalized beta stability for data in each hierarchical level) of order q > 0 and provides the stability profile which depicts the stability with respect to order q.\cr


This package contains five main functions: \cr

  1. \code{Stay_Single} calculates stability of the time series data (like biomass, productivity, etc.) for single assemblage. \cr

  2. \code{Stay_Multiple} calculates Gamma, Alpha and Beta stability, and synchrony of the time series data (like biomass, productivity, etc.) for multiple assemblages. \cr

  3. \code{Stay_Hier} calculates stability of the time series data (like biomass, productivity, etc.) for hierarchical structure. \cr

  4. \code{ggStay_qprofile} provides to graph the q-profile of stability (and synchrony if is multiple assemblages) based on the output obtained from the function \code{Stay_Single}, \code{Stay_Multiple} or \code{Stay_Hier}. \cr

  5. \code{ggStay_analysis} provides to graph the relationships between stability (and synchrony if is multiple assemblages) and an additional diversity (or other) variable based on the output obtained from the function \code{Stay_Single} or \code{Stay_Multiple}. \cr

}

\keyword{ package }
