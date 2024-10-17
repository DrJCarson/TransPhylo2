
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TransPhylo2

![GitHub](https://img.shields.io/github/license/DrJCarson/TransPhylo2)

# Introduction

TransPhylo2 is an extension of
[TransPhylo](https://github.com/xavierdidelot/TransPhylo) and
[TransPhyloMulti](https://github.com/DrJCarson/TransPhyloMulti) that
(optionally) incorporates discrete epidemiological data (demes).
Epidemiological parameters (reproduction number and sampling proportion)
may also differ across demes.

For a more formal description of TransPhylo2, see the following paper:

Carson et al. (2024). Incorporating epidemiological data into the
genomic analysis of partially sampled infectious disease outbreaks. In
preparation.

## Installation

You can install TransPhylo2 from github with:

``` r
devtools::install_github("DrJCarson/TransPhylo2")
```

The package can then be loaded using:

``` r
library(TransPhylo2)
```

## Basic usage

You will need as input a dated phylogeny built for example using
[BEAST](https://www.beast2.org/),
[BactDating](https://github.com/xavierdidelot/BactDating) or
[treedater](https://github.com/emvolz/treedater). This dated phylogeny
should be stored in the object `dated_tree` of class phylo from the ape
package. Since this only includes relative rather than absolute dating,
you also need to know the date when the last sample was taken, let’s say
on the 1st July 2022. You will also need a named vector `hosts`
indicating the host from which each of the leaves was sampled. Finally,
if deme data are to be incorpated, they must be added as a vector to the
phylogenetic tree. You can load and plot this data to make sure it looks
correct:

``` r
pt=ptreeFromPhylo(dated_tree,lubridate::decimal_date(as.Date('2022/7/1')),hosts)
pt$demes=demes
plot(pt)
```

You can then infer the transmission tree and associated parameters using
for example:

``` r
res=inferTTree(pt,w.shape=2,w.scale=2,obs.end=lubridate::decimal_date(as.Date('2023/1/1')))
plot(res)
```

if all demes use the same model parameters, or

``` r
res=inferTTreeMutli(pt,w.shape=2,w.scale=2,obs.end=lubridate::decimal_date(as.Date('2023/1/1')))
plot(res)
```

otherwise. Note that the parameters `w.shape` and `w.scale` specify a
gamma distribution that needs to be appropriate to represent the
generation time distribution of the pathogen under investigation. It is
also important to input the date when sampling of cases ended, for
example the 1st January 2023 in the code above.

## More information and getting help

For more detailed examples of how to use TransPhylo2, see the vignettes
[here](https://github.com/DrJCarson/TransPhyloMulti/tree/master/vignettes).
We also include the code used to generate the results in the manuscript
[here](https://github.com/DrJCarson/TransPhyloMulti/tree/master/run).
See also the help included in the package using the R command
`help(package='TransPhylo2')`.

If you have any problem or question please create an issue
[here](https://github.com/DrJCarson/TransPhylo2/issues)
