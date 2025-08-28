\name{PigDisturbance}
\alias{PigDisturbance}
\docType{data}
\title{
Non-observable information about the disturbances occuring in the \code{PigFarming} dataset
}
\description{
Information about the disturbances relative to the dataset \code{PigFarming}: start, and intensity. An intensity of 0 stands for no disturbance for the considered hierarchical level. Not observable in practice. 
}
\usage{data("PigDisturbance")}
\format{
  A dataframe with 6000 records on the following 12 variables:
  \describe{
    \item{\code{id}}{the identifier of the animal}
    \item{\code{batch}}{the numero of the batch of the animal}
    \item{\code{pen}}{the numero of the batch of the animal}
    \item{\code{int_batch}}{the intensity of the disturbance at the batch level}
    \item{\code{start_batch}}{the starting time of the disturbance at the batch level}
    \item{\code{end_batch}}{the endding time of the disturbance at the batch level}
    \item{\code{int_pen}}{the intensity of the disturbance at the pen level}
    \item{\code{start_pen}}{the starting time of the disturbance at the pen level}
    \item{\code{end_pen}}{the endding time of the disturbance at the pen level}
    \item{\code{int_ind}}{the intensity of the disturbance  at the individual level}
    \item{\code{start_ind}}{the starting time of the disturbance at the individual level}
    \item{\code{end_ind}}{the endding time of the disturbance at the individual level}
  }
}

\references{
Le, Vincent, Tom Rohmer, and Ingrid David. 2022. “Impact of Environmental Disturbances on Estimated Genetic Parameters and Breeding Values for Growth Traits in Pigs.” Animal 16 (4): 100496.
https://doi.org/10.1016/j.animal.2022.100496

}

\examples{
str(PigDisturbance)
PigDisturbance[c(1,6,5405),]
}

\keyword{simulated datasets}

