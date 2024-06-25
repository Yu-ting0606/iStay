\name{Jena_species_biomass_data}
\alias{Jena_species_biomass_data}
\docType{data}
\title{
16-years biomass data of each speices in the plot in Jena experiment (except 2004)
}
\description{
  In Jena experiment, there are 4 blocks (B1, B2, B3, B4) with a total of 76 plots. For each plot 1, 2, 4, 8 or 16 plants are sown. There are approximately 3 or 4 plots in each block for each number of plant species planted. This data includes 76 lists which are corresponding to 76 plots. In each list, it is a species by time dataframe that record the biomass of each species sown in the plot. However, because Jena experiment only records the total biomass of all species in the plot, biomass of each species are unknown in 2004, the time length of data is about 16-years except 2004. \cr
}

\usage{data(Jena_species_biomass_data)}
\format{
  \code{Jena_species_biomass_data} contains 76 lists. Each list is a dataframe that the rows are species and the columns are the time. For each row, it record the biomass of each species during the time.
 }

\references{
Wagg, C., Roscher, C., Weigelt, A. et al. (2022) Biodiversity–stability relationships strengthen over time in a long-term grassland experiment. Nature Communications 13, 7752. \cr\cr
}

\keyword{datasets}
