
.GRCM<-function(diag_, m = 50)
{
        mat<-c()
	N <- nrow(diag_) + m
        temp_mat<-matrix(rnorm(N*ncol(diag_)) , nrow=N, ncol=ncol(diag_) )
        for(i in 2:nrow(diag_)){
                for(j in 1:(i-1) ){
                        if(diag_[j,i]==0)
                        {
                                mat<-append(mat,j)

                        }

                }
                if(length(mat)!=0)
                {
                        temp<-temp_mat[,(mat)]
                        tm<-.calc_P(temp)
                        temp_mat[,i]<-temp_mat[,i]-tm%*%temp_mat[,i]

                }
                mat<-c()

        }
        w<-pseudoinverse(t(temp_mat)%*%temp_mat)
        d_d<-diag(diag(w)^(-1/2),nrow(diag_))
        dg<-d_d%*%w%*%d_d
        retval <- list(sigcov=w, sigcor=dg)
        return (retval)

}
###############################################################
.GRCM2 <- function(diag_, m=50)
{
        mat<-c()
        N <- nrow(diag_) + m
        Z <- matrix(rnorm(N*ncol(diag_)) , nrow=N, ncol=ncol(diag_) )
	#print(nrow(Z))
        for(i in 2:nrow(diag_)){
          mat <- which(diag_[1:(i-1), i]==0)
                if(length(mat) != 0)
                {
                  z.I.tilde <- Z[,(mat)]
                  aux <- t(z.I.tilde) %*% z.I.tilde
                  aux2 <- t(z.I.tilde) %*% Z[,i]
                  aux3 <- pseudoinverse(aux) %*% aux2
                  aux4 <- z.I.tilde %*% aux3
                  Z[,i] <- Z[,i] - aux4
                 
                }
                mat<-c()
        }
        aux <- t(Z)%*%Z
        w<-pseudoinverse( aux )
        d_d<-diag(diag(w)^(-1/2),nrow(diag_))
        dg<-d_d%*%w%*%d_d
        retval <- list(sigcov=w, sigcor=dg)
        return (retval)
}

###############################################################
.GRCM3 <- function(diag_, m=100)
{
        mat<-c()
        N <- nrow(diag_) + m
        Z <- matrix(rnorm(N*ncol(diag_)) , nrow=N, ncol=ncol(diag_) )
	#print(nrow(Z))
        for(i in 2:nrow(diag_)){
          mat <- which(diag_[1:(i-1), i]==0)
                if(length(mat) != 0)
                {
                  z.I.tilde <- Z[,(mat)]
                  aux <- t(z.I.tilde) %*% z.I.tilde
                  aux2 <- t(z.I.tilde) %*% Z[,i]
                  aux3 <- solve(aux) %*% aux2
                  aux4 <- z.I.tilde %*% aux3
                  Z[,i] <- Z[,i] - aux4
                 
                }
                mat<-c()
        }
        aux <- t(Z)%*%Z
        w<-solve( aux )
        d_d<-diag(diag(w)^(-1/2),nrow(diag_))
        dg<-d_d%*%w%*%d_d
        retval <- list(sigcov=w, sigcor=dg)
        return (retval)
}

#########################This is called internally by GRCM function #####################################################
.calc_P<-function(mat)
{
        temp_mat<-(mat%*%(pseudoinverse(t(mat)%*%mat)))%*%t(mat)
        return (temp_mat)
}
###################################################################################################
.GRCM_n <- function(diag_, rho = 0, method= "IPF"){

    if(class(rho)=="matrix"){
        diag(rho) <- 1
        if(!is.positive.definite(rho)){
            rho <- as.matrix(nearPD(rho, corr=T)$mat)
        }
    }

	#print(method)
	w <- as.matrix(qpG2Sigma(diag_, rho=rho, matrix.completion=method))
	d_d<-diag(diag(w)^(-1/2),nrow(diag_))
        dg<-d_d%*%w%*%d_d
        retval <- list(sigcov=w, sigcor=dg)
        return (retval)
	
}
###############################################################################










