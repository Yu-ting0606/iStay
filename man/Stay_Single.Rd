% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Main_Function.R
\name{Stay_Single}
\alias{Stay_Single}
\title{Calculate stability of the time series data for single assemblage.}
\usage{
Stay_Single(
  data,
  order.q = c(1, 2),
  Alltime = TRUE,
  start_T = NULL,
  end_T = NULL
)
}
\arguments{
\item{data}{can be input as a \code{vector} of time series data, or \code{data.frame} (assemblages by times).}

\item{order.q}{a numerical vector specifying the orders of stability. Default is c(1,2).}

\item{Alltime}{\code{TRUE} or \code{FALSE}, to decide whether to use all the times in the data.}

\item{start_T}{(argument only for \code{Alltime = FALSE}) a positive integer specifying the start column of time in the data.}

\item{end_T}{(argument only for \code{Alltime = FALSE}) a positive integer specifying the end column of time in the data.}
}
\value{
a dataframe with columns "Plot/Community", "Order_q" and "Stability".
}
\description{
\code{Stay_Single} is a function that calculate stability of the time series data (like biomass, productivity, etc.) for single assemblage.
}
\examples{
# Stability of each single plot
data("Jena_plot_biomass_data")
single_plot <- do.call(rbind, Jena_plot_biomass_data)
output_single_plot <- Stay_Single(data=single_plot, order.q=c(1,2), Alltime=TRUE)
output_single_plot

# Stability of each single species in each plot
data("Jena_species_biomass_data")
single_species <- do.call(rbind, Jena_species_biomass_data)
output_single_species <- Stay_Single(data=single_species, order.q=c(1,2), Alltime=TRUE)
output_single_species

}
