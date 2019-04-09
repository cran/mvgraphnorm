### R code from vignette source 'mvgraphnorm-intro.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: style-Sweave
###################################################
library(tools)
#BiocStyle::latex()


###################################################
### code chunk number 2: preliminaries
###################################################
library(mvgraphnorm)
options(SweaveHooks=list(twofig=function() {par(mfrow=c(1,2))},
                         twofig2=function() {par(mfrow=c(2,1))},
                         onefig=function() {par(mfrow=c(1,1))}))


###################################################
### code chunk number 3: plotunif
###################################################
############ Figure A #############
library("mvgraphnorm")
## input network ##
g <- barabasi.game(100, directed=FALSE)
zz1 <- rmvggm(net.str=g, method="htf")
zz2 <- rmvggm(net.str=g, method="ipf")
zz3 <- rmvggm(net.str=g, method="kim")
names(zz1)
head(zz1[[2]][1:10,1:10])


###################################################
### code chunk number 4: plotunif
###################################################
############ Figure A #############
library("mvgraphnorm")
## input network ##
g <- barabasi.game(50, directed=F)

## Generate a random correlation matrix ##
m <- matrix(0,50,50)
m[upper.tri(m)] <- runif(min=-.8, max=.8, 25*49)
m <- m+t(m)
diag(m) <- 1

zz <- rmvggm(20, net.str=g, method="htf"
, cor=m)
head(zz[[2]][1:10,1:10])



###################################################
### code chunk number 5: plotunif
###################################################
############ Figure A #############
library("mvgraphnorm")
## input network ##
g <- barabasi.game(50, directed=F)

## Generate a random correlation matrix ##
m <- matrix(0,50,50)
m[upper.tri(m)] <- runif(min=-.8, max=.8, 25*49)
m <- m+t(m)
diag(m) <- 1
zz <- rmvggm(20, net.str=g, method="htf"
, cor=m)

p <- NULL
p <- viz.rmvggm(zz, net=FALSE)
dev.off()
print(p$covp)
print(p$covsmp)




###################################################
### code chunk number 6: plotunif
###################################################
############ Figure A #############
library("mvgraphnorm")
## input network ##
g <- barabasi.game(50, directed=F)

## Generate a random correlation matrix ##
m <- matrix(0,50,50)
m[upper.tri(m)] <- runif(min=-.8, max=.8, 25*49)
m <- m+t(m)
diag(m) <- 1
zz <- rmvggm(20, net.str=g, method="htf", cor=m)

p <- NULL
p <- viz.rmvggm(zz, net=TRUE, test="cor", undirected = TRUE)
print(p$netp)




###################################################
### code chunk number 7: mvgraphnorm-intro.Rnw:190-191
###################################################
sessionInfo()


