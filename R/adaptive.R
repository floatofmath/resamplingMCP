##' Computes the conditional type I error rate of a pre-planned z-test in a two-stage parallel group adaptive design. We condition on the observed first stage data and treatment assignment as well as the pre-planned second stage sample size.
##' 
##' @title z - Test Conditional Error Rate
##' @param x1 First stage observations
##' @param g1 First stage treatment assignments (1 for control, -1 for treatment)
##' @param n Preplanned overal sample size
##' @param alpha pre-fixed significance level
##' @return numeric value of the conditional error rate
##' @author Florian Klinglmueller
##'
##' @export
normal_CER <- function(x1,g1,n,alpha=0.025,sigma=1,one_sample=F){
    n1 <- length(x1)
    z  <- zstat(x1,g1,sigma,one_sample=one_sample)
    pnorm((qnorm(alpha,lower=F)-sqrt(n1/n)*z)/sqrt(1-n1/n),lower=F)
}


##' Computes the conditional type I error rate of a pre-planned t-test (using inverse normal combination of stage-wise t-tests) in a two-stage parallel group adaptive design. We condition on the observed first stage data and treatment assignment as well as the pre-planned second stage sample size.
##' 
##' @title t - Test Conditional Error Rate
##' @param x1 First stage observations
##' @param g1 First stage treatment assignments (1 for control, -1 for treatment)
##' @param n Preplanned overal sample size
##' @param alpha pre-fixed significance level
##' @return numeric value of the conditional error rate
##' @author Florian Klinglmueller
##'
##' @export
t_CER <- function(x1,g1,n,alpha=0.025,one_sample=FALSE,...){
    n1 <- length(x1)
    if(one_sample){
        z <- qnorm(t.test(x1,alternative='less')$p.value,lower=FALSE)
    } else {
        z <- qnorm(t.test(x1~g1,alternative='less',var.equal=T)$p.value,lower=FALSE)
    }
    pnorm((qnorm(alpha,lower=F)-sqrt(n1/n)*z)/sqrt(1-n1/n),lower=F)
}


##' Computes the p-value of the t-test (assuming equal variances).
##'
##' @title t-test p-value
##' @param x1 First stage observations
##' @param x2 Second stage observations
##' @param g1 First stage treatment assignments
##' @param g2 Second stage treatment assignments
##' @param x3 Third stage observations
##' @param g3 Third stage treatment assignments
##' @return p-value of the t-test
##' @author Florian Klinglmueller
##'
##' @export
t_test <- function(x1,x2,g1,g2,x3=NULL,g3=NULL){
    x <- c(x1,x2,x3)
    y <- c(g1,g2,g3)
    t.test(x~y,alternative='less',var.equal=TRUE)$p.value
}


##' Computes the p-value of the z-test (assuming known variances).
##'
##' @title z-test p-value
##' @param x1 First stage observations
##' @param x2 Second stage observations
##' @param g1 First stage treatment assignments
##' @param g2 Second stage treatment assignments
##' @param sigma Common standard deviation
##' @param x3 Third stage observations
##' @param g3 Third stage treatment assignments
##' @return p-value of the z-test
##' @author Florian Klinglmueller
##'
##' @export
z_test <- function(x1,x2,g1,g2,sigma=1,x3=NULL,g3=NULL){
    pnorm(zstat(c(x1,x2,x3),c(g1,g2,g3)),lower=F)
}

##' Computes the conditional type I error rate of a pre-planned permutation test in a two-stage adaptive design. We condition on the observed first stage data and treatment assignment as well as the observed second stage data - which we assume are obtained when the experiment reaches its preplanned sample size.
##'
##' Based on the first stage data and treatment assignments one may perform sample size reassassment - and possibly other trial modifications - as long as the (preplanned) second stage sample size is not reduced.
##'
##' \code{stat} needs to be a function of the form \code{function(x,g,...)} returning a numeric of length one. Possible options are \code{\link{sumdiff}}, \code{\link{meandiff}}, \code{\link{zstat}}
##'
##' For the moment, we assume that observations are randomized using random allocation blocked by stages, (i.e. we resample using \code{sample (g1)}). \code{g2} does not have to be the actual second stage treatment assignments but just one possible example randomization, that fixes the treatment group sizes. 
##'
##' @title Permutation Conditional Error Rate
##' @param x1 vector of preplanned first stage observations
##' @param g1 vector of first stage treatment assignments
##' @param x2 vector of preplanned second stage observations
##' @param stat function computing the test statistic (see Details)
##' @param permutations number of permutations (rerandomizations) used to compute unconditional and conditional permutation distributions
##' @param alpha pre-fixed significance level
##' @param g2 template vector for second stage treatment assignments
##' @param each 
##' @param ... additional options to \code{stat}
##' @return numeric value of the conditional error rate
##' 
##' @author Florian Klinglmueller
##' 
##' @export
permutation_CER <- function(x1,g1,x2,stat=sumdiff,
                            permutations = 1000,alpha=.025,
                            g2=rep(c(-1,1),each=length(x2)/2),one_sample=FALSE,
                            restricted=TRUE,...){
  n1 <- length(g1)
  n <- length(c(x1,x2))
  n2 <- n-n1
  ## right parameter setting for one_sample
  if(one_sample&restricted) warning("Only un-restricted sampling is allowed for one-sample setting!")
  if(one_sample) {
    restricted <- FALSE 
    g1 <- sign(x1)
    x1 <- abs(x1)
    g2 <- sign(x2)
    x2 <- abs(x2)
  }
  ## balanced second stage!!
  dist <- perm_dist(x1,x2,g1,g2,stat,permutations,restricted=restricted,...)
  cdist <- cond_dist(x1,x2,g1,g2,stat,permutations,restricted=restricted,...)
  m <- length(dist)
  talpha <- min(dist[rank(dist,ties.method='min')>=ceiling((1-alpha)*m)])
  cer1 <- mean(cdist>talpha)
  ##    pvals <- unlist(lapply(cdist,function(x) sum(dist>=x)/B))
  ##    cer2 <- sum(pvals<=alpha)/B
  c(cer1)
}


##' Computes the conditional type I error rate of a pre-planned permutation test in a two-stage adaptive design. We condition on the observed first stage data and treatment assignment as well as the observed second stage data as well as the observations from the extension of the trial.
##'
##' Based on the first stage data and treatment assignments one may perform sample size reassassment - and possibly other trial modifications - as long as the (preplanned) second stage sample size is not reduced.
##'
##' \code{stat} needs to be a function of the form \code{function(x,g,...)} returning a numeric of length one. Possible options are \code{\link{sumdiff}}, \code{\link{meandiff}}, \code{\link{zstat}}
##'
##' For the moment, we assume that observations are randomized using random allocation blocked by stages, (i.e. we resample using \code{sample (g1)}). \code{g2} does not have to be the actual second stage treatment assignments but just one possible example randomization, that fixes the treatment group sizes. 
##'
##' @title Permutation Conditional Error Rate
##' @param x1 vector of preplanned first stage observations
##' @param g1 vector of first stage treatment assignments
##' @param x2 vector of preplanned second stage observations
##' @param x3 vector of observations from the extended trial
##' @param stat function computing the test statistic (see Details)
##' @param permutations number of permutations (rerandomizations) used to compute unconditional and conditional permutation distributions
##' @param subsamples number of (second stage) subsamples used for conditional error rate estimation
##' @param alpha pre-fixed significance level
##' @param g2 template vector for second stage treatment assignments
##' @param each 
##' @param one_sample 
##' @param restricted 
##' @param ... additional options to \code{stat}
##' @return numeric value of the conditional error rate
##' 
##' @author Florian Klinglmueller
##' 
##' @export
permutation_CER2 <- function(x1,g1,x2,x3,stat=sumdiff,
                             permutations = 1000, subsamples = 1000, alpha=.025,
                             g2=rep(c(-1,1),each=length(x2)/2),one_sample=FALSE,
                             restricted=TRUE,...){
    n1 <- length(g1)
    n <- length(c(x1,x2))
    n2 <- n-n1
    nt <- length(c(x1,x2,x3))
    ne <- nt-n
    if(ne == 0) stop("Conditional should not be computed if sample size is not increased")
    ## right parameter setting for one_sample
    if(one_sample&restricted) warning("Only un-restricted sampling is allowed for one-sample setting!")
    if(one_sample) {
        restricted <- FALSE 
        g1 <- sign(x1)
        x1 <- abs(x1)
        g2 <- sign(x2)
        x2 <- abs(x2)
    }
    ## balanced second stage!!
    
    if(subsamples < choose(n2+ne,n2)){
        x2rs <- random_samples_cpp(c(x2,x3),n2,subsamples)
    } else {
        x2rs <- subsamples_cpp(c(x2,x3),n2)
    }
    cers <- apply(x2rs,2,cer,x1=x1,g1=g1,g2=g2,stat=stat,permutations=permutations,restricted=restricted,alpha=alpha,...)
    ##    pvals <- unlist(lapply(cdist,function(x) sum(dist>=x)/B))
    ##    cer2 <- sum(pvals<=alpha)/B
    mean(cers)
}




cer <- function(x1,x2,g1,g2,stat,permutations,restricted,alpha,...){
    dist <- perm_dist(x1,x2,g1,g2,stat,permutations,restricted=restricted,...)
    m <- length(dist)
    talpha <- min(dist[rank(dist,ties.method='min')>=ceiling((1-alpha)*m)])
    cdist <- cond_dist(x1,x2,g1,g2,stat,permutations,restricted=restricted,...)
    mean(cdist>talpha)
}
