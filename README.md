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

### Configure CPU cores

Many of the scripts in this project take advantage of parallelization using the R package `doMC`. The number of threads used is customizable in the `lib/globals.R` file. Edit this file and set to an appopriate number. We recommend allowing 4GB of RAM per core (e.g. set 4 cores on a machine with 16GB of RAM).

### Execute the SetupProject script

From the R console, run the SetupProject.R script. This will install any missing packages and download data files from the [figshare record](https://doi.org/10.6084/m9.figshare.6005297), and do quite a bit of data munging. Follow the instructions when the program asks for figshare authentication.

```
> source("./R/SetupProject.R")
```

This will take a while to complete. Once it does, we recommend restarting your R session before continuing.

*Note: downloading the data files from figshare for the first time requires a web browser to authenticate.* If running on a remote machine without this capability, execute the following commands in an R console locally.

```
> library(rfigshare)
> figshare_article <- fs_details(6005297)
```

This will open a web browser session, log in through figshare, and generate a `.httr-oauth` file in your working directory. Copy this file to the project directory on the remote machine. Running the `SetupProject.R` script as above should no longer require the authentication in a web browser.

## Running analyses

This R project uses a directory structure and functions from the R package, [ProjectTemplate](http://projecttemplate.net/). `run.project()` runs each `*.R` file in the `src/` directory. 

From a fresh R session, run the following commands.

```
library(ProjectTemplate)
run.project()
```

This will take many hours to complete, depending on the number of CPU cores used. The bulk of this time is spent running the `CORUM_Shuffle.R` and `huMAP_Shuffle.R` scripts, which run the network permutation tests for all protein complexes in each dataset. Many of the analyses depend on output files from these scripts. However, any script in the `src/` directory beginning with a lower number prefix (e.g. `00-`, `01-`, `02-`, etc.) will not depend on these output files. Run these scripts individually with `source()` if you're short on time.

