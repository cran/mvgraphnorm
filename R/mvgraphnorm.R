#library(mvtnorm)
#library(igraph)
#library(qpgraph)
#library(corpcor)
rmvggm <- function(n.samples=100, 
net.str=NULL, method=c("ipf", "htf", "kim" ), 
cor=0.0, mean.vector=NULL, m.var=100){

	if (is.null(net.str)) {
           		stop("Please input an undirected graph object of igraph or a symmetric binary matrix")

	}
	if ((!is.igraph(net.str)) && (!is.matrix(net.str)) &&  (class(net.str)!="dgCMatrix")) {
			stop("Not a graph object or a binary matrix")
			
	}
	
	if(is.igraph(net.str)){
		if(is.directed(net.str)){
			stop("Not an undirected graph object")
		}
		net.str <- get.adjacency(net.str)
		if(!isSymmetric(as.matrix(net.str))) {
			stop("Please input an undirected graph object of igraph or a symmetric binary matrix")
		}
	}
	
	smp <- .gensample(n=n.samples+nrow(net.str), 
		bin.mat=net.str, method=method[1], 
		cor=cor, mean.vector=mean.vector, m.var=m.var)
	smp
}

.gensample <- function(n=NULL, bin.mat=NULL, 
method=NULL, cor=0.0, mean.vector=NULL, m.var=100){
	
	print(method)
	if(method[1]=="ipf"){
		sigval <- .GRCM_n(diag_=bin.mat, rho=cor, method="IPF")
	}
	if(method[1]=="htf"){
		sigval <- .GRCM_n(diag_=bin.mat, rho=cor, method="HTF")
	}
        if(method[1]=="kim"){
	           sigval <- .GRCM3(diag_=bin.mat, m=m.var)
	}
	if(is.null(mean.vector)||(length(mean.vector)!=ncol(bin.mat))){
		mean.vector=rep(0, nrow(bin.mat))
	}
	if(n<2 || is.null(n)){
		d <- NULL
	}
	else{
		d <- rmvnorm(n=n, mean=mean.vector, sigma=sigval[[1]], method="svd")
	}
        gendata <- list(dat=d, sigma=sigval[[1]], net.str=bin.mat)
        class(gendata) <- "rmvggm"
        gendata  	
}
viz.rmvggm <- function(x, col=c("red", "blue"), net=FALSE, ...){
 m <- g <- p3 <- NULL
 if(net){
 	pc <- as.data.frame(x$dat)
 	zz <- gs(pc,...)
 	ed <- zz$arcs
 	g <- simplify(graph.edgelist(ed, directed=FALSE))
 	p3 <- plot(g, vertex.size=5)
 	p3 <- recordPlot()
 	
 }
 x1 <- pseudoinverse(x$sigma)
 x2 <- x$net.str
 x1 <- x1[upper.tri(x1)]
 x2 <- x2[upper.tri(x2)]
 
 y1 <- pseudoinverse(cov(x$dat))
 y1 <- y1[upper.tri(y1)]
 
 dist1 <- list("1"=x1[which(x2==1)], "0"=x1[which(x2==0)])
 dist2 <- list("1"= y1[which(x2==1)], "0"= y1[which(x2==0)] )
 
 boxplot(dist1, xlab="distribution of edge and non-edge components \n of inverse of  covariance matrix from \n generated covariance matrix using algorithm",
 	cex.lab=.65, col=c("red", "blue"))
 p1 <- recordPlot()
 boxplot(dist2, xlab="distribution of edge and non-edge components \n of inverse of covariance matrix of sampled-data", 
 col=col, cex.lab=.7)
 p2 <- recordPlot()
 invisible(list("covp"=p1,"covsmp"=p2,
 	"netp"=p3, "infnet" = g))
}
