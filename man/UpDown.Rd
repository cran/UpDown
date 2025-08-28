\name{UpDown}
\alias{UpDown}
\title{Unsupervised hierarchical classification for elastic or plastic response to a disturbance using mixture of Gaussian  distributions}

\description{Detection and characterisation of disturbances from longitudinal data, organized in hierarchical groups
  }

\usage{
  UpDown(data,levels,obs, vtime, h.int=NULL,
  mixplot=FALSE, correction=NULL, 
  kappa=NULL,thr_va=0.5,
  options=list())
}

\arguments{
  \item{data}{a dataframe containing at minima observations,  a time variable and hierarchical levels (one column per level). The dataframe can also contain other variables. One row per unit and observed time is needed.
  }
  \item{levels}{a vector of character strings specifing the column names corresponding to the considered hierarchical levels appearing in the dataframe \code{data}, ordered from the highest  hierarchical level to the lowest hierarchical level (i.e., the unit level)
  }
  \item{obs}{character string specifing the column names of the considered numeric observations appearing in the data frame \code{data}.
  }
  \item{vtime}{character string specifing the column names of the considered time variable appearing in the data frame \code{data}.
  }
  \item{h.int}{a real parameter specifying the smoothing bandwidth in Nadaraya-Watson's smoothing curves \emph{see ?ksmooth}. The default value is \eqn{\sqrt {n}}{sqrt(n)}, where \eqn{n} is the largest length of observations per unit.
  }
  \item{mixplot}{logical. If TRUE, the mixture curves for each hierarchical levels are plotted. (default value FALSE)
  }
  \item{correction}{an optional character string specifing the column name of a considered time-dependent discrete variable appearing in the dataframe \code{data} for which a correction of the observation \code{obs} is needed. The correction consists in subtracting the median relative to the modalities of the variable. 
  %When \code{correction} has a length greater than 1, corrections are succesively done starting by the first component of the vector \code{correction} up to the last component. 
  %If \code{correction} is \code{null} (default value), no correction is done.
  }
  \item{kappa}{an integer in [0,1) used to eliminate a redundancy disturbance between two distinct levels, based on the estimated starting and end points. It evaluates the overlapping between two considered disturbances. The disturbance is removed in the lowest of the two hierarchical levels. When kappa is not specified (null default value) no redounding disturbances are eliminated. When kappa is equal to 1, only the disturbances with the same rounded starting and ending points are removed. An excessively small values for kappa can lead to wrongly remove disturbances. The suggested value for kappa is 0.75.
  }
  \item{thr_va}{an integer in [0,1) used to validate in the down step the group disturbances. \code{thr_va} control the percentage of starting-times of the underlying elements in a close interval to validate the disturbance. A value of 0 indicate that all the disturbances are validated that is to say no validation is done. The default value is 0.5, i.e., at least 50\% of the elements constituting the group should have a disturbance in a close time to validate the group disturbance.
  }
  \item{options}{A list of options. See the documentation \emph{see ?normalmixEM} for possible options of the mixture model
  }
  
}


\details{Note that unique identifiers are mandatory for the hierarchical levels. Moreover note that UpDown considers that all perturbations have a negative effect on the longitudinal observations. For a positive effect, consider the opposite sign of the observations before using UpDown. Units with less than 20 observations are removed. That can be modified using the option \code{minobs} in \code{options}.
}

\value{
  A list containing the following components:
    \item{data}{the initial dataframe with the supplementary columns (if \code{correction} is not null); \emph{medobs}, the median of observations grouped into the modalities of the qualitative variable used in \code{correction}; and \emph{cobs} the corrected observations. Units with less than \code{minobs} observations are removed.
    }
  \item{levels}{the specified hierarchical levels.
  }
  \item{med_lev}{a list of matrices. The last matrix contains the medians of the unit observations per observed time. The previous one contains the medians of these medians per observed time and so on up to the highest hierarchical level. If there are less than 50\% observations per considered time, the median will not be evaluated (NA).
  }

  \item{mixmdl_lev}{a list of outputs of the mixture models of the hierarchical levels. The first output concerns the lowest hierarchical levels (i.e., unit) and the last output concerns the highest hierarchical levels.
  }
  \item{names_lev}{a list of names of each elements of the hierarchical levels. The first matrix concerns the lowest hierarchical levels (i.e., unit) and the last matrix concerns the highest hierarchical levels.
  }
  \item{sc.x_lev}{a list of matrices giving the time points considered by the smoothing, per identifier and for each hierarchical level. The first matrix concerns the lowest hierarchical levels (i.e., unit) and the last matrix concerns the highest hierarchical levels. 
  }
  \item{sc.y_lev}{a list of matrices of fitted values corresponding to \code{sc.x_lev}. The first matrix concerns the lowest hierarchical levels (i.e., unit) and the last matrix concerns the highest hierarchical levels.
  }
  \item{sc.dx_lev}{a list of matrices giving the time points considered by the derivative smoothing curve, per identifier and for each hierarchical level.  The first matrix concerns the lowest hierarchical levels (i.e., unit) and the last matrix concerns the highest hierarchical levels.
  }
  \item{sc.dy_lev}{a list of matrices of fitted values corresponding to \code{sc.dx_lev}. The first matrix concerns the lowest hierarchical levels (i.e., unit) and the last matrix concerns the highest hierarchical levels.
  }
  \item{Up}{a dataframe that describes for each ids, the type of detected disturbance at the end of the Up-step. '0' stands for no disturbance.
  }
  \item{Down}{a list of matrices that gives the detected disturbances for each hierarchical level and their characteristics.
  }
}

\author{Tom Rohmer, Vincent Le, Ingrid David
}

\references{

Benaglia, T., Chauveau, D., Hunter, D. R., and Young, D. mixtools: An R package for analyzing finite mixture models. Journal of Statistical Software, 32(6):1-29, 2009.

Le, V. 2022. “Nouvelle mesure de la robustesse des animaux d’élevage par utilisation des données de phénotypage haut-débit.” Thesis, INPT Toulouse. https://hal.inrae.fr/tel-03967884.

Nadaraya, E. A. On estimating regression. Theory of Probability & Its Applications, 9(1):141–142, 1964

Scrucca L., Fop M., Murphy T. B. and Raftery A. E. (2016) mclust 5: clustering, classification and density estimation using Gaussian finite mixture models, The R Journal, 8/1, pp. 289-317.

Watson, G. S. Smooth regression analysis. Sankhya: The Indian Journal of Statistics, Series A, pages
359–372, 1964

}


\examples{
\donttest{
# optional arguments
options<-list(maxit=100) 

# considered hierarchical levels
levels=c("batch","pen","id")

UpDown.out<- UpDown(PigFarming, levels=levels, vtime="time", obs="weight",
kappa=0.75, thr_va=0.5, h.int=10, mixplot=FALSE, correction="age", options=options)

UpDown.out$Down$batch
}
}



\seealso{
\link[mixtools]{normalmixEM}, \link[mclust]{Mclust}, \link[stats]{ksmooth}
}


