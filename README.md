#Modern Measures of Differentiation

`mmod` is an R package for calculating modern population divergence statistics. 

##Quickstart

###Install

mmod isn't on CRAN, so you can install the latest stable version using 
`install.packages("mmod")`. This github repository may be running ahead of 
the version on CRAN, if you really want the latest version you can download
 `mmod_*.tar.gz` and install via the command line (or with the GUI if you'd prefer):

        $R CMD INSTALL mmod_*.tar.gz 

###Usage
Once it's up an running all you need is genepop file with your data

        >library(mmod)
        >my_data <- read.genepop("my_file.gen")
        >diff_stats(my_data)
   
##Overview

Population geneticists have traditionally used Nei's Gst (often confusingly called 
Fst...) to measure divergence between populations. It turns out, Gst doesn't really
measure divergence so, [a set of new measures have been developed]
(http://www.molecularecologist.com/2011/03/should-i-use-fst-gst-or-d-2/)
  
mmod is a package that brings two of these measures, Hedricks (2008) G'st 
and Jost's (2008) D to R, along with an implementation of Nei's Gst that
uses nearly unbiased estimators for Hs and Ht, the two key parameters from
which all these stats are calculated. All these functions work on `genind`
objects from the library `adegenet` so data can be read in from standard
`genepop`files. An overview of a typical usage is provided in a vignette
called "demo", acessable from `vignette("demo", package="mmod")`

Briefly yhere are functions for each of these measures which give values for 
each locus in a `genind` object and a global estimate:`D_Jost()`, 
`Gst_Hedrick()`, `Gst_Nei()`. Because most of the heavy-lifting in calculating 
all these stats is finding Hs and Ht, a function, `diff_stats()` is 
provided to calculate each at once.

Each of the stats can be calculated for each pairwise comparison of populations
 in a dataset: `pairwise_D()`, `pairwise_Gst_Hedrick()`, `pairswise_Gst_Nei()`.

Finally, the function `jacknife_population()` can, as the name suggests,
calculate these statistics in a sample jacknifed across populations

mmod is still very much in development, so I'm happy to receive and
suggestions contributions or bugs you might find

