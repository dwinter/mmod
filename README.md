#Modern Measures of Differentiation

`mmod` is an R package for calculating modern population divergence statistics. 

##Quickstart

###Install

mmod is on CRAN, so you can install the latest stable version using 
`install.packages("mmod")`. This github repository may be running ahead of 
the version on CRAN, if you really want the latest version you can download
 `mmod_*.tar.gz` and install via the command line (or with the GUI if you'd 
 prefer):

        $R CMD INSTALL mmod_*.tar.gz 

###Usage
Once it's up an running all you need is genepop (or fstat) file with your data

        >library(mmod)
        >my_data <- read.genepop("my_file.gen")
        >diff_stats(my_data)
   
##Overview

Population geneticists have traditionally used Nei's Gst (often confusingly 
called Fst...) to measure divergence between populations. It turns out, Gst 
doesn't really measure divergence so, [a set of new measures have been developed]
(http://www.molecularecologist.com/2011/03/should-i-use-fst-gst-or-d-2/)
  
mmod is a package that brings two of these measures; Hedricks (2005, 2011) G''st 
 Jost's (2008) D and Meirman's (2005) φ'st to R, along with a function that 
 calculates Nei's Gst using nearly unbiased estimators for Hs and Ht  (the two 
 key parameters from which most of these stats are calculated). All these 
 functions work on `genind` objects from the library `adegenet` so data can be 
 read in from standard `genepop` for `fstat` files. An overview of a typical 
 usage is provided in a vignette called "demo", acessable from 
 `vignette("demo", package="mmod")`, I suggest new users read this before that 
 start. 

##Help
All functions are documented, there are [pretty versions of the help pages here]
(http://dwinter.github.com/mmod/).



