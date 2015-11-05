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
  rnorm(n)/sqrt(rchisq(n,df)/df) + ncp
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

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param MCMC Number of trials to simulate
##' @param n1 First stage sample size
##' @param n Total preplanned sample size
##' @param rfs First stage random number generator
##' @param rss Second stage random number generator
##' @param ncp Non-centrality parameter
##' @param n3 Third stage sample size (maybe a function of \code{x},\code{g2})
##' @param cond Condition on first and second stage observations ("both"), only on first stage observations ("first"), or unconditional of observations
##' @return 
##' @author Florian Klinglmueller
simulate_trials <- function(MCMC,n1,n,n3=0,r=1/2,rfs=rnorm,rss=rnorm,ncp=0,cond=c("both","first","none")){
    cond <- match.arg(cond)[1]
    restricted = TRUE
    n2 <- n-n1
    ns <- c(n1,n2)
    if(!isTRUE(all.equal(ns,as.integer(ns),check.attributes=FALSE))) { warning(cat("Relative group size",r,"do not permit whole numbered group sizes:",n*r," will be rounded!")) }
    ks <- floor(ns*r)
    if(ncp != 0 && cond != "none") { stop("Can not condition on observations under the alternative") }
    if(cond == "both"){
        x1 <- rfs(n1)
        x2 <- rss(n2)
        ans <- list(x1=x1,x2=x2,x3=numeric(0),g=strat_reassignments(ns,ks,restricted=restricted,B=MCMC))
        if(is.function(n3)){
            n3 <- n3(x1,ans$g[,1:n1])
            ans  <- within(ans,
                           {
                               x3 <- lapply(n3,rss)
                               if(restricted){
                                   g <- lapply(1:MCMC,function(i) c(g[,i],sample(c(rep(0L,n3[i]-ceiling((1-r)*n3[i])),rep(1L,n3[i]-floor(n3[i]*r))))))
                               } else {
                                   g <- lapply(1:MCMC,function(i) c(g[,i],sample(0L:1L,n3[i],replace=TRUE,prob=c(1-r,r))))
                               }
                           })
        } else if(n3>0) {
            ans <- within(ans,
                          {
                              x3 <- rss(n3)
                              g = rbind(g,strat_reassignments(n3,floor(n3*r),B=MCMC))
                          })
        }
    } else if(cond == "first") {
        x1 <- rfs(n1)
        x2 <- matrix(rss(n2*MCMC),ncol=MCMC,nrow=n2)
        ans <- list(x1=x1,x2=x2,x3=numeric(0),g=rbind(strat_reassignments(ns[1],ks[1],restricted=restricted,B=MCMC),matrix(rep(0L,ns[2]-ks[2]),rep(1L,ks[2]),nrow=MCMC,ncol=n2,byrow=T)))
        if(is.function(n3)){
            n3 <- n3(x1,ans$g[,1:n1])
            ans  <- within(ans,
                           {
                               x3 <- lapply(n3,rss)
                               if(restricted){
                                   g <- lapply(1:MCMC,function(i) c(g[i,],sample(c(rep(0L,n3[i]-ceiling((1-r)*n3[i])),rep(1L,n3[i]-floor(n3[i]*r))))))
                               } else {
                                   g <- lapply(1:MCMC,function(i) c(g[i,],sample(0L:1L,n3[i],replace=TRUE,prob=c(1-r,r))))
                               }
                           })
        } else if(n3>0) {
            ans <- within(ans,
                          {
                              x3 <- rss(n3)
                              g = cbind(g,strat_reassignments(n3,floor(n3*r),B=MCMC))
                          })
        }

    }
    return(ans)
}
