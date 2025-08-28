\name{UpDownApp}
\alias{UpDownApp}
\title{Shiny app to vizualise the longitudinal intra-group observations and the detected disturbances}

\description{
The function starts a shiny app which visualizes the data organized by the hierarchical levels, and the estimated start and end point of the detected disturbances
  }

\usage{
  UpDownApp(updown.out,obs=NULL,width=1000,height=1000)
}

\arguments{
  \item{updown.out}{Global output of the \code{UpDown} function }
  \item{obs}{(optional) vector of character string specifing the names of the considered longitudinal observations appearing in the dataframe \code{data}.
  }
  \item{width,height}{(optional) interger specifing the width and the height of the plot in the Rshiny.
  }
  
}

\value{No return value, run the application}

\seealso{
\link[shiny]{runApp}, \link[ggplot2]{ggplot}, \link[=UpDown]{UpDown}
}


\examples{
\donttest{
levels=c("batch","pen","id")
UpDown.out<- UpDown(PigFarming,levels=levels, vtime="time", obs="weight",
kappa=0.75, thr_va=0.5, correction="age")
}

if(interactive()){UpDownApp(UpDown.out)}

}

