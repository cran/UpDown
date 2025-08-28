\name{PigFarming}
\alias{PigFarming}
\docType{data}
\title{
Simulated longitudinal phenotypes that mimics a pig-farming dataset subsected to disturbances
}
\description{
Example of a dataset on which UpDown can be applied. It consists in simulated hierarchical data mimicking a pig farming-system dataset subsected to disturbances. The animals (id) were raised in 40 batches and in 15 pens within each batch leading to 15 animals per pen. Hence three hierarchical level are considered: id, pen and batch levels.  Data were simulated following Le et al. 2022. <doi.org/10.1016/j.animal.2022.100496>
}
\usage{data("PigFarming")}
\format{
  A data frame with 578847 individual observations on the following 6 variables:
  \describe{
    \item{\code{id}}{the identifier of the animal}
    \item{\code{batch}}{the numero of the batch}
    \item{\code{pen}}{the numero of the pen}
    \item{\code{age}}{the age (in day) of the animal}
    \item{\code{time}}{the observation times}
    \item{\code{weight}}{the weight (in kg) of the animal}
  }
}


\references{
Le Vincent, Tom Rohmer, and Ingrid David. 2022. “Impact of Environmental Disturbances on Estimated Genetic Parameters and Breeding Values for Growth Traits in Pigs.” Animal 16 (4): 100496.
https://doi.org/10.1016/j.animal.2022.100496

}

\examples{
str(PigFarming)
plot(subset(PigFarming,id==6)$weight)
}
\keyword{simulated datasets}
