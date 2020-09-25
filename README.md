# sQTLviztools

An R package and associated scripts for visualisation splicing quantitative trait loci (sQTLs).

## Installation

1. clone the sqtlviztools repo

```
git clone https://github.com/jackhump/sqtlviztools.git 
```

2. install the R package and any dependencies

in R:

```{r}
install.packages("remotes")
remotes::install_github("jackhump/sQTLviztools")
```

# Running the example data

I picked a SNP / intron cluster pair at random from the ROSMAP shiny app - https://rajlab.shinyapps.io/sQTLviz_ROSMAP/

In the examples folder there is a toy version of each file you'll need to visualise sQTLs.

1. Create the example Rdata file

on the command line:

```
cd sqtlviztools/
Rscript prepare_example.R
```

2. Run the example Shiny app

```
cd shiny/
Rscript shiny/run_shiny.R
```

Please raise an issue on this github if you can't get the example app working.


