# Interrogation of Mammalian Protein Complex Structure, Function, and Membership Using Genome-Scale Fitness Screens

## Pan, Meyers, *et al*. Cell Systems. 2018.

Code repository for reproducing analyses and figures in this report.

## Setup

### Clone this repository

From the command line:

```
$ git clone https://github.com/robinmeyers/pan-meyers-et-al
```

### Open an R session from this directory

This can be done a variety of ways.

```
$ cd pan-meyers-et-al
$ R
```

or if you use RStudio

```
$ cd pan-meyers-et-al
$ open Pan_Meyers_et_al.Rproj
```

From the R console, run the SetupProject.R script. This will install any missing packages and download data files from the [figshare record](https://doi.org/10.6084/m9.figshare.6005297). Follow the instructions when the program asks for figshare authentication.

```
> source("./R/SetupProject.R")
```

This will take a while to complete. Once it does, we recommend restarting your R session before continuing.




