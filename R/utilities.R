##' Make contrast matrix for all (pairwise) many-to-one comparisons
##'
##' @param g Vector of treatment assignments
##' @title Contrasts for many-to-one comparisons
##' @export
##' @author Florian Klinglmueller
make_pw_contrasts <- function(g,control=0){
    levels <- unique(g)
    sapply(levels[levels!=control],function(l) (g == l) - (g == control))
}

##'  Simulate normal heteroscedastic data with scaling of the variance drawn from a chi-square distribution with \
##'
##' @title Generate heteroskedastic data
##' @param n  number of observations 
##' @param df degrees of freedom for the variance distribution
##' @param delta non-centrality parameter
##' @export
rnorm_hetero <- function(n,df=1,ncp=0){
  rnorm(n)/sqrt(rchisq(n,df)/df) + ncp
}


##' Simulate from a tail contaminated normal distribution
##'
##' @title Generate contaminated normally distributed data
##' @param n number of observations
##' @param mean mean value
##' @param sd standard deviation 
##' @param cprop proportion of contaminated samples
##' @param csd standard deviation of contaminated samples
##' @return vector with \code{n} pseudo random numbers
##' @author Florian Klinglmueller
##' @export
rnorm_tcont <- function(n,mean=0,sd=1,cprop=.1,csd=3){
    if(n == 1)
        return(rnorm(1,mean=mean,sd=sample(c(sd,csd),1)))
    out <- sample(c(rnorm(ceiling(n*(1-cprop)),mean=mean,sd=sd),rnorm(floor(n*cprop),mean=mean,sd=csd)))
    out
}

##' Simulate from a location shift contaminated normal distribution
##'
##' @title Generate contaminated normally distributed data
##' @param n number of observations
##' @param mean mean value
##' @param sd standard deviation 
##' @param cprop proportion of contaminated samples
##' @param cmean standard deviation of contaminated samples
##' @return vector with \code{n} pseudo random numbers
##' @author Florian Klinglmueller
##' @export
rnorm_scont <- function(n,mean=0,sd=1,cprop=.1,cmean=3,csd=sd){
    if(n == 1)
        return(rnorm(1,mean=mean,sd=sample(c(sd,csd),1)))
    m <- ceiling(n*(1-cprop))
    k <- n-m
    out <- sample(c(rnorm(m,mean=mean,sd=sd),
                    sample(c(-1,1),k,replace=TRUE)*rnorm(k,mean=mean+cmean,sd=csd)))
    out
}

##' Template function to compute the Type 1 error of a trial design
##'
##' @title Type 1 error 
##' @param X Random number generater for first stage observations
##' @param X2 Random number generator for second stage observations (if NULL first stage generator will be used)
##' @param n pre-planned total sample size
##' @param n1 first stagen sample size
##' @param B number of permutations to be used to compute the permutation null distribution
##' @param MCMC number of simulation runs for Type I error estimation
##' @param one_sample should a one-sample test be simulated
##' @param ... additional functions to conditional error rate function
##' @return matrix of p-values from different test procedures
##' @author Florian Klinglmueller
##'
##' @export
type1 <- function(X=rnorm,X2=NULL,n,n1,stat=meandiff,B,MCMC,one_sample=FALSE,...){
  x <- X(n)
  pm <- dplyr::bind_rows(lapply(1:MCMC,function(i) {
    g <- sample(rep(c(-1,1),each=n1/2))
    y <- if(is.null(X2)) x else c(x[1:n1],X2(n-n1))
        list(preplanned=permutation_CER(x[1:n1],g,x[(n1+1):n],B=B,stat=stat,one_sample=one_sample),
             influenced=permutation_CER(y[1:n1],g,y[(n1+1):n],B=B,stat=stat,one_sample=one_sample),
             normal=normal_CER(x[1:n1],g,n,one_sample=one_sample,...),
             ttest=t_CER(x[1:n1],g,n,one_sample=one_sample,...))
        }))
    pm
}


##' Template function to compute the Type 1 error of a trial design
##'
##' @title Type 1 error 
##' @param X Random number generater for first stage observations
##' @param X2 Random number generator for second stage observations (if NULL first stage generator will be used)
##' @param ncp Non-centrality (shift) paramater 
##' @param n pre-planned total sample size
##' @param n1 first stagen sample size
##' @param stat 
##' @param B number of permutations to be used to compute the permutation null distribution
##' @param MCMC number of simulation runs for Type I error estimation
##' @param one_sample should a one-sample test be simulated
##' @param ... additional functions to conditional error rate function
##' @return matrix of p-values from different test procedures
##' @author Florian Klinglmueller
##'
##' @export
compare_power <- function(X=rnorm,X2=NULL,ncp,n,n1,stat=meandiff,B,MCMC,one_sample=FALSE,...){
  x <- X(n)
  pm <- dplyr::bind_rows(lapply(1:MCMC,function(i) {
    g <- sample(rep(c(-1,1),each=n1/2))
    y <- if(is.null(X2)) x else c(x[1:n1],X2(n-n1))
        list(preplanned=permutation_CER(x[1:n1],g,x[(n1+1):n],B=B,stat=stat,one_sample=one_sample),
             influenced=permutation_CER(y[1:n1],g,y[(n1+1):n],B=B,stat=stat,one_sample=one_sample),
             normal=normal_CER(x[1:n1],g,n,one_sample=one_sample,...),
             ttest=t_CER(x[1:n1],g,n,one_sample=one_sample,...))
        }))
    pm
}

##' Simulate trial data
##'
##' @title Simulate trials
##' @param MCMC Number of trials to simulate
##' @param n1 First stage sample size
##' @param n Total preplanned sample size
##' @param rfs First stage random number generator
##' @param rss Second stage random number generator
##' @param ncp Non-centrality parameter
##' @param n3 Third stage sample size (maybe a function of \code{x},\code{g2})
##' @param cond Condition on first and second stage observations ("both"), only on first stage observations ("first"), or unconditional of observations
##' @return list object with simulated data
##' @export
##' @author Florian Klinglmueller
simulate_trials <- function(MCMC,n1,n,n3=0,r=1/2,rfs=rnorm,rss=rnorm,ncp=0,cond=c("both","first","none"),restricted=TRUE){
    cond <- match.arg(cond)[1]
    n2 <- n-n1
    ns <- c(n1,n2)
    if(!isTRUE(all.equal(ns*r,as.integer(ns*r),check.attributes=FALSE))) { warning(paste("Relative group size",r,"do not permit whole numbered group sizes:",paste(ns*r,collapse=' ')," will be rounded!")) }
    ks <- floor(ns*r)
    if(ncp != 0 && cond != "none") { stop("Can not condition on observations under the alternative") }
    if(cond == "both"){
        x1 <- rfs(n1)
        x2 <- rss(n2)
        ans <- list(x=c(x1,x2),g=strat_reassignments(ns,ks,restricted=restricted,B=MCMC))
    } else if(cond == "first") {
        g1 <- strat_reassignments(ns[1],ks[1],restricted=restricted,B=MCMC)
        g2 <- if(restricted) matrix(c(rep(0L,ns[2]-ks[2]),rep(1L,ks[2])),ncol=ncol(g1),nrow=ns[2]) else matrix(sample(0L:1L,ns[2]*ncol(g1),replace=TRUE,prob=c(1-r,r)),nrow=ns[2])
        x1 <- rfs(n1)
        x2 <- matrix(rss(n2*ncol(g1)),nrow=n2)
        ans <- list(x=rbind(x1,x2),
                    g=rbind(g1,g2))
    } else if(cond == "none") {
        g1 <- if(restricted) matrix(c(rep(0L,ns[1]-ks[1]),rep(1L,ks[2])),nrow=ns[1],ncol=MCMC) else matrix(sample(0L:1L,ns[1]*MCMC,replace=TRUE,prob=c(1-r,r)),nrow=ns[1])
        g2 <- if(restricted) matrix(c(rep(0L,ns[2]-ks[2]),rep(1L,ks[2])),nrow=ns[2],ncol=MCMC) else matrix(sample(0L:1L,ns[2]*MCMC,replace=TRUE,prob=c(1-r,r)),nrow=ns[2])
        x1 <- matrix(rfs(n1*MCMC),nrow=n1)
        x2 <- matrix(rss(n2*MCMC),nrow=n2)
        x1 <- x1 + ncp*g1
        x2 <- x2 + ncp*g2
        ans <- list(x=rbind(x1,x2),
                    g=rbind(g1,g2))
        
    }
    if(is.function(n3)){
        n3 <- n3(x1,ans$g[,1:n1])
        k3 <- floor(n3*r)
        ans  <- within(ans,
                       {
                           if(restricted){
                               g3 <- lapply(1:ncol(g),function(i) c(rep(0L,n3[i]-k3[i])),rep(1L,k3))
                           } else {
                               g3 <- lapply(1:ncol(g),function(i) sample(0L:1L,n3[i],replace=TRUE,prob=c(1-r,r)))
                           }
                           x3 <- lapply(1:ncol(g),rss(n3[i])+ncp*g[i])
                       })
    } else if(n3>0) {
        k3 <- floor(n3*r)
        ans <- within(ans,
                      {
                          x3 <- matrix(rss(n3*ncol(g)),nrow=n3)+ncp*g3
                          g3 <- if(restricted) matrix(c(rep(0L,n3-k3),rep(1L,k3)),nrow=n3,ncol=ncol(ans$g)) else matrix(sample(0L:1L,n3*ncol(ans$g),replace=TRUE,prob=c(1-r,r)),nrow=n3)
                      })
    }
    return(ans)
}

## n1=3
## n=5
## B=10
## expect_warning(s1 <- simulate_trials(B,n1,n))
## expect_equal(length(s1$x),n)
## expect_equal(nrow(s1$g),n)
## expect_equal(ncol(s1$g),pmin(prod(choose(c(n1,n-n1),floor(1/2*c(n1,n-n1)))),B))
## expect_equal(unique(colSums(s1$g)),sum(floor(1/2*c(n1,n-n1))))
## expect_equal(unique(colSums(s1$g[1:n1,])),floor(1/2*n1))
## expect_equal(unique(colSums(s1$g[(n1+1):n,])),floor(1/2*(n-n1)))
