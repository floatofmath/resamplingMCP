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
rnorm_hetero <- function(n,df=1,ncp=0){
  rnorm(n)/sqrt(rchisq(n,df)/variance_df) + ncp
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
type1 <- function(X=rnorm,X2=NULL,n,n1,B,MCMC,one_sample=FALSE,...){
  x <- X(n)
  pm <- rbindlist(lapply(1:MCMC,function(i) {
    g <- sample(rep(c(-1,1),each=n1/2))
    y <- if(is.null(X2)) x else c(x[1:n1],X2(n-n1))
        list(preplanned=permutation_CER(x[1:n1],g,x[(n1+1):n],B=B,one_sample=one_sample),
             influenced=permutation_CER(y[1:n1],g,y[(n1+1):n],B=B,one_sample=one_sample),
             normal=normal_CER(x[1:n1],g,n,one_sample=one_sample,...),
             ttest=t_CER(x[1:n1],g,n,one_sample=one_sample,...))
        }))
    pm
}
