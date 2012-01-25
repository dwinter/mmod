pkgname <- "mmod"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('mmod')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("D_Jost")
### * D_Jost

flush(stderr()); flush(stdout())

### Name: D_Jost
### Title: Calculate Jost's D
### Aliases: D_Jost

### ** Examples

data(nancycats)
D_Jost(nancycats)



cleanEx()
nameEx("Gst_Hedrick")
### * Gst_Hedrick

flush(stderr()); flush(stdout())

### Name: Gst_Hedrick
### Title: Calculate Nei's Gst using estimators for Hs and Ht
### Aliases: Gst_Hedrick

### ** Examples

data(nancycats)
Gst_Hedrick(nancycats)



cleanEx()
nameEx("Gst_Nei")
### * Gst_Nei

flush(stderr()); flush(stdout())

### Name: Gst_Nei
### Title: Calculate Nei's Gst using estimators for Hs and Ht
### Aliases: Gst_Nei

### ** Examples

data(nancycats)
Gst_Nei(nancycats)



cleanEx()
nameEx("chao_bootstrap")
### * chao_bootstrap

flush(stderr()); flush(stdout())

### Name: chao_bootstrap
### Title: Produce bootstrap samples from each subpopulation of a genind
###   object
### Aliases: chao_bootstrap

### ** Examples

## Not run: 
##D data(nancycats)
##D bs <- chao_bootstrap(nancycats)
##D bs_D <- sapply(bs, D_Jost)
##D hist(unlist(bs_D[2,]))
##D quantile( unlist(bs_D[2,]), c(0.025, 0.5, 0.975) )
## End(Not run)



cleanEx()
nameEx("diff_stats")
### * diff_stats

flush(stderr()); flush(stdout())

### Name: diff_stats
### Title: Calculate differentiation statistics for a genind objects
### Aliases: diff_stats

### ** Examples

data(nancycats)
diff_stats(nancycats)



cleanEx()
nameEx("diff_test")
### * diff_test

flush(stderr()); flush(stdout())

### Name: diff_test
### Title: An exact test of population differentiation for Genind objects
### Aliases: diff_test

### ** Examples

data(nancycats)
diff_test(seploc(nancycats)[[2]], nreps=100)



cleanEx()
nameEx("harmonic_mean")
### * harmonic_mean

flush(stderr()); flush(stdout())

### Name: harmonic_mean
### Title: Harmonic mean
### Aliases: harmonic_mean

### ** Examples

data(nancycats)
pop.sizes <- table(pop(nancycats))
harmonic_mean(pop.sizes)



cleanEx()
nameEx("jacknife_populations")
### * jacknife_populations

flush(stderr()); flush(stdout())

### Name: jacknife_populations
### Title: Calculate differentiation stats for a jacknife sample of a
###   Genind opject
### Aliases: jacknife_populations

### ** Examples

## Not run: 
##D data(nancycats)
##D obs <- diff_stats(nancycats)
##D jn <- jacknife_populations(nancycats)
##D D_sampled <-jn[5,]
##D hist(D_sampled)
##D abline(h=obs$global)
## End(Not run)



cleanEx()
nameEx("pairwise_D")
### * pairwise_D

flush(stderr()); flush(stdout())

### Name: pairwise_D
### Title: Calculates pairwise values of Jost's D
### Aliases: pairwise_D

### ** Examples

data(nancycats)
pairwise_D(nancycats[1:26,])



cleanEx()
nameEx("pairwise_Gst_Hedrick")
### * pairwise_Gst_Hedrick

flush(stderr()); flush(stdout())

### Name: pairwise_Gst_Hedrick
### Title: Calculates pairwise values of Hedrick's G'st
### Aliases: pairwise_Gst_Hedrick

### ** Examples

data(nancycats)
pairwise_Gst_Hedrick(nancycats[1:26,])



cleanEx()
nameEx("pairwise_Gst_Nei")
### * pairwise_Gst_Nei

flush(stderr()); flush(stdout())

### Name: pairwise_Gst_Nei
### Title: Calculates pairwise values of Nei's Gst
### Aliases: pairwise_Gst_Nei

### ** Examples

data(nancycats)
pairwise_Gst_Nei(nancycats[1:26,])



cleanEx()
nameEx("summarise_bootsrap")
### * summarise_bootsrap

flush(stderr()); flush(stdout())

### Name: summarise_bootsrap
### Title: Apply a differentiation statistic to a bootstrap sample sample
### Aliases: summarise_bootsrap

### ** Examples

## Not run: 
##D data(nancycats)
##D bs <- chao_bootstrap(nancycats)
##D summarise_bootstrap(bs, D_Jost)
## End(Not run)



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
