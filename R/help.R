#' Modern Measures of Divergence
#'
#'
#'Population geneticists have traditionally used Nei's Gst (often confusingly called 
#'Fst...) to measure divergence between populations. It turns out, Gst doesn't really
#'measure divergence so, a set of new measures have been developed.
#'
#' 
#'mmod is a package that brings two of these mesures, Hedricks (2008) G'st 
#'and Jost's (2008) D to R, along wiht an implementation of Nei's Gst that
#'uses nearly unbiased estimators for Hs and Ht, the two key paramaters from
#'which all these stats are calculated. All these functions work on genind
#'objects from the libary adegenet so data can be read in from standard
#'genepop files.
#'
#' @docType package
#' @name mmod
#' @aliases mmod mmod-package
NULL
