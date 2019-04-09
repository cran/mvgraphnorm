\name{rmvggm}
\alias{rmvggm}
\title{
	Generate samples from a multivariate gaussian graphical distribution
}
\description{
	Generate samples from a multivariate gaussian graphical distribution using two algorithms  (1) Castelo & Roverato, 2006, 2009 (2) Kim & van de Wiel, 2008.
}
\usage{
rmvggm(n.samples=100, 
  net.str=NULL, method=c("ipf", "htf", "kim" ), 
  cor=0.0, mean.vector=NULL, m.var=100)
}

\arguments{
  \item{n.samples}{
   is a numeric variable which represents total 
   numbers of samples to be generated. 
}
  \item{net.str}{
  is a undirected graph of igraph-class object or can be a
  binary-matrix representation of a graph and should be symmetric.
}
  \item{method}{
    Three methods are used. "ipf", "kim"  and "htf". Default method is "ipf". 
}
  \item{cor}{
     is a numeric variable or an input symmetric correlation matrix,  
     which represents the average correlation of edge-components of 
     the graph given as an input. This is required for "IPF and "HTF"
     method. If the input matrix is not positive definite matrix
     we use  Higham (2002) algorithm to convert the matrix into the
     nearest positive definite matrix.
}
  \item{mean.vector}{
 is a vector which represents the mean of variables of 
 gaussian graphical model. This option is used with "castelo" method.
}
  \item{m.var}{
	number of additional auxilary variables to be used to generate 
  a covariance matrix using "kim" method, deafult value of 'm.var' is 50. 
  This option is not used with 'castelo' method. 
  The minimum value of 'm.var' option should be 0.  
}
}
\details{
 random samples are geneerated from a gaussian graphical 
 model using covariance matrix generated from two algorithms:
  (1) Castelo & Roverato, 2006, 2009 (2) Kim & van de Wiel, 2008.  

}
\value{
 returns a list object of "rmvggm" class. 
 The first component of the list object is a matrix 
 of the sampled data. The rows are samples and columns 
 represent the total number of valriables. 
 Second component is a covariance matrix which 
 is used to generate samples from the multivariate 
 normal distribution. 
}
\references{
 (1) Kim, K. I. & van de Wiel, M. (2008). Effects of dependence 
 in high-dimensional multiple testing problems. BMC Bioinformatics, 9 (1), 114.

 (2) Castelo, R. & Roverato, A. (2006). A robust procedure for 
 Gaussian Graphical Model search from microarray 
 data with p larger than n. J. Mach. Learn. Res., 7, 2621–2650.

 (3) Higham, Nick (2002) Computing the nearest correlation matrix - 
 a problem from finance; IMA Journal of Numerical Analysis 22, 329–343.
}
\author{
Shailesh Tripathi, Frank Emmert-Streib
}


\seealso{
 'pseudoinverse', 'qpG2Sigma', 'rmvnorm' 
}
\examples{
 	library(mvgraphnorm)
 	g <- barabasi.game(100, directed=FALSE)
 	kk <- rmvggm(net.str=g)
	kk <- rmvggm(net.str=g, n.samples=1000, 
    cor=0, method="ipf")
	summary(kk)
	kk1 <- rmvggm(net.str=g, n.samples=1000, 
    cor=0, method="kim")
	summary(kk1)
###kk$dat is the data where columns 
###represent variables and rows represent samples. 
###Total 1000 samples are generated for each variable
###kk$sigma is a covariance matrix used to generate samples from gaussian graphical model.
###kk$net.str is a binary matrix of graph object given as an input.

}
