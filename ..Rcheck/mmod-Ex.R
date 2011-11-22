pkgname <- "mmod"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('mmod')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("D.Jost")
### * D.Jost

flush(stderr()); flush(stdout())

### Name: D.Jost
### Title: Calculate Jost's D using
### Aliases: D.Jost

### ** Examples

data(nancycats)
D.Jost(nancycats)



cleanEx()
nameEx("Gst.Hedrick")
### * Gst.Hedrick

flush(stderr()); flush(stdout())

### Name: Gst.Hedrick
### Title: Calculate Nei's Gst using estimators for Hs and Ht
### Aliases: Gst.Hedrick

### ** Examples

data(nancycats)
Gst.Hedrick(nancycats)



cleanEx()
nameEx("Gst.Nei")
### * Gst.Nei

flush(stderr()); flush(stdout())

### Name: Gst.Nei
### Title: Calculate Nei's Gst using estimators for Hs and Ht
### Aliases: Gst.Nei

### ** Examples

data(nancycats)
Gst.Nei(nancycats)



cleanEx()
nameEx("diff.stats")
### * diff.stats

flush(stderr()); flush(stdout())

### Name: diff.stats
### Title: Calculate differentiation statistics for a genetic dataset
### Aliases: diff.stats

### ** Examples

data(nancycats)
diff.stats(nancycats)



cleanEx()
nameEx("harmonic.mean")
### * harmonic.mean

flush(stderr()); flush(stdout())

### Name: harmonic.mean
### Title: Harmonic mean
### Aliases: harmonic.mean

### ** Examples

data(nancycats)
pop.sizes <- table(pop(nancycats))
harmonic.mean(pop.sizes)



cleanEx()
nameEx("jacknife.populations")
### * jacknife.populations

flush(stderr()); flush(stdout())

### Name: jacknife.populations
### Title: Calculate differentiation stats for a jacknife sample of a
###   Genind opject
### Aliases: jacknife.populations

### ** Examples

## Not run: 
##D data(nancycats)
##D jacknife.populations(nancycats)
## End(Not run)



cleanEx()
nameEx("pairwise.D")
### * pairwise.D

flush(stderr()); flush(stdout())

### Name: pairwise.D
### Title: Calculates pairwise values of Jost's D
### Aliases: pairwise.D

### ** Examples

data(nancycats)
pairwise.D(nancycats[1:26,])



cleanEx()
nameEx("pairwise.Gst.Hedrick")
### * pairwise.Gst.Hedrick

flush(stderr()); flush(stdout())

### Name: pairwise.Gst.Hedrick
### Title: Calculates pairwise values of Hedrick's G'st
### Aliases: pairwise.Gst.Hedrick

### ** Examples

data(nancycats)
pairwise.Gst.Hedrick(nancycats[1:26,])



cleanEx()
nameEx("pairwise.Gst.Nei")
### * pairwise.Gst.Nei

flush(stderr()); flush(stdout())

### Name: pairwise.Gst.Nei
### Title: Calculates pairwise values of Nei's Gst
### Aliases: pairwise.Gst.Nei

### ** Examples

data(nancycats)
pairwise.Gst.Nei(nancycats[1:26,])



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
