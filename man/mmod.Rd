\docType{package}
\name{mmod}
\alias{mmod}
\alias{"mmod-package"}
\title{Modern Measures of Divergence}
\description{
  Population geneticists have traditionally used Nei's Gst
  (often confusingly called Fst...) to measure divergence
  between populations. It turns out, Gst doesn't really
  measure divergence so, a set of new measures have been
  developed.
}
\details{
  mmod is a package that brings two of these mesures,
  Hedricks (2008) G'st and Jost's (2008) D to R, along wiht
  an implementation of Nei's Gst that uses nearly unbiased
  estimators for Hs and Ht, the two key paramaters from
  which all these stats are calculated. All these functions
  work on genind objects from the libary adegenet so data
  can be read in from standard genepop files. See
  vignette("demo", package="mmod") for an example usage.
}

