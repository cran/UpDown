UpDown
==========

In many scientific contexts, longitudinal data are observed and organized in
hierarchical groups. This is the case for example for breeding data from high throughput
sensors, where measurements are automatically performed on animals distributed in pens
belonging to different batches. In these hierarchical systems, disturbances can affect an
individual or a group of individuals by changing their expected trajectories. Being able to
identify a disturbance affecting an individual or a whole group of individuals (which may have
different responses to the disturbance due to their own ‘resilience’) is often a major issue. The
UpDown package consists in using the information of all the groups or subgroups to facilitate
the detection and the characterization of disturbances (beginning, duration, intensity) affecting
the elements constituting the groups. The R package UpDown has been developed to consider as many
hierarchical levels as necessary, thus allowing to adapt to different scientific problems.


Installation
----------
You can install the latest released version from CRAN with:

```
install.packages("UpDown")

```
or the latest development version:

```
#install.packages("remotes")
remotes::install_github("TomRohmer/UpDown")
```

Usage
----------
An example of a dataset on which UpDown could be used is available on the library Updown. It consists of simulated hierarchical data that mimics a pig farming system dataset based subject to disturbances with elastic response.


```
library(UpDown)
data<-get(data(PigFarming))
head(data)

```

Note that unique identifiers are mandatory for the hierarchical levels.  Units with less than 20 observations are removed. That can be modified using the option `minobs` in options as well as the global options of the EM algorithms:

```
options<-list(minobs=25, maxit = 500)

```

Below is a short example of how to execute UpDown:


```
levels=c("batch","pen","id")
UpDown.out<- UpDown(data2,levels=levels, vtime="time", obs="weight",
kappa=0.75, thr_va=0.5, mixplot=FALSE, correction="age",options=options)
```

Then the detected disturbance as well as the estimation of the starting, end points ans the intensity of the disturbance can be visualized using

```
UpDown.out$Down
UpDown.out$Down$batch

```

An Rshiny application
----------

An Rshiny application makes possible to visualize the data organized by the hierarchical levels, and the estimated start and end point of the detected disturbances as well as the 'median observations' and the corresponding smoothing curves.


```
UpDownApp(updown.out)

```


Acknowledgement
----------

We thank Alliance R&D (Axiom, Choice Genetics, Nucleus and IFIP) for the financial support of this research and for testing the \textbf{UpDown} package on real data from pig farms.


References
----------

Le, Vincent. 2022. 'Nouvelle mesure de la robustesse des animaux d’élevage par utilisation des données de phénotypage haut-débit.' 
Thesis, INPT Toulouse. 
https://hal.inrae.fr/tel-03967884.

Le, Vincent, Tom Rohmer, and Ingrid David. 2022. 
'Impact of Environmental Disturbances on Estimated Genetic Parameters and Breeding Values for Growth Traits in Pigs.'
 Animal 16 (4): 100496.
