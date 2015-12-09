##' Adaptive permutation test one-sample problems
##'
##' @title One-sample adaptive permutation test
##' @template onesample_sims
##' @param combination_function Function to combine stage-wise (permutation) p-values
##' @param perms Maximum number of permutations to use when computing permutation p-values and conditional error rates
##' @author Florian Klinglmueller
##' @export
adaptive_permtest_os <- function(x,n1,n,ne,test_statistic,perms=50000,alpha=0.025){
    if(ne>n){
        xs <- split(x,rep(1:3,c(n1,n-n1,ne-n)))
        gs <- split(sign(x)>0,rep(1:3,c(n1,n-n1,ne-n)))
        A <- permutation_CER(xs[[1]],gs[[1]],xs[[2]],test_statistic,one_sample=TRUE,restricted=FALSE,
                             permutations=perms,alpha=alpha)
        q <- perm_test(xs[[2]],xs[[3]],gs[[2]],gs[[3]],test_statistic,restricted=FALSE,B=perms)
        return(A>=q)
    } else {
        xs <- split(x,rep(1:2,c(n1,ne-n1)))
        gs <- split(sign(x)>0,rep(1:2,c(n1,ne-n1)))
        return(alpha>=perm_test(xs[[1]],xs[[2]],gs[[1]],gs[[2]],test_statistic,restricted=FALSE,B=perms))
    }
}


##' Adaptive t-test as described in Timmesfeld et al. (2007)
##'
##' Warning the current implementation seems to be numerically instable
##' @title Adaptive t-test
##' @param mpfr Whether to use high precision numbers from \code{Rmpfr}
##' @template onesample_sims
##' @author Florian Klinglmueller
##' @export
adaptive_ttest_os <- function(x,n1,n,ne,alpha=0.025,mpfr=FALSE) {
    if(n == ne){
        return(0.025 >= t.test(x,alternative='greater')$p.value)
    }
    xs <- split(x,rep(1:2,c(n1,ne-n1)))
    V1 <- sum(xs[[1]])
    U <- sum(xs[[1]]^2)
    tU <- sum(xs[[2]]^2)
    A <- clev(tU,U,V1,ne-n1,n,n1,alpha=alpha,mpfr=mpfr)
    A >= t.test(xs[[2]],alternative='greater')$p.value
}

##' Adaptive combination test of stage-wise t-tests using the inverse normal combination function. 
##'
##' @title Inverse normal adaptive t-test
##' @template onesample_sims
##' @author Florian Klinglmueller
##' @export
adaptive_invnormtest_os <- function(x,n1,n,ne,alpha=0.025){
    xs <- split(x,rep(1:2,c(n1,ne-n1)))
    p1 <- t.test(xs[[1]],alternative='greater')$p.value
    p2 <- t.test(xs[[2]],alternative='greater')$p.value
    pnorm(sum(sqrt(c(n1,n-n1)/n) * qnorm(c(p1,p2),lower=F)),lower=FALSE) <= alpha 
}

##' Non-parametric combination of stage-wise test statistics. Combines stage-wise permutation p-values using some combination function; performs the test using the joint conditional permutation distribution of stage wise permutation p-values.
##'
##' @title adptive NPC test
##' @template onesample_sims
##' @param test_statistic Function that computes the test test statistic
##' @param combination_function Function to combine stage-wise (permutation) p-values
##' @param perms Maximum number of permutations to use when computing permutation p-values and conditional error rates
##' @author Florian Klinglmueller
##' @export
adaptive_npcombtest_os <- function(x,n1,n,ne,test_statistic,combination_function=inverse_normal,perms=50000,alpha=0.025) {
    xs <- split(x,rep(1:2,c(n1,ne-n1)))
    gs <- split(sign(x)>0,rep(1:2,c(n1,ne-n1)))
    G <- omega(gs[[1]],gs[[2]],restricted=FALSE,B=perms)
    rB <- ncol(G)
    p1 <- 1-(rank(test_statistic(xs[[1]],G[1:n1,]))/(rB+1))
    p2 <- 1-(rank(test_statistic(xs[[2]],G[(n1+1):ne,]))/(rB+1))
    ct <- combination_function(p1,p2,n1/n,(n-n1)/n)
    sum(ct[1]>=ct)/length(ct) <= alpha
}



adaptive_invnormtest_2s <- function(x,y,n1,n,ne,m1,m,me,alpha=0.025){
    xs <- split(x,rep(1:2,c(n1,ne-n1)))
    p1 <- t.test(xs[[1]],alternative='greater')$p.value
    p2 <- t.test(xs[[2]],alternative='greater')$p.value
    alpha >= pnorm(sum(sqrt(c(n1,n-n1)/n) * qnorm(c(p1,p2),lower=F)),lower=FALSE) 
}



##' Compare adaptive test procedures for general two-sample cases
##'
##' \code{rdist} needs to take the number of samples to return as its first argument.
##' 
##' @title Compare two-sample adaptive tests
##' @param n1 first stage sample size (control group)
##' @param n preplanned total sample size (control group)
adaptive_permtest_2s <- function(x,y,n1,n,ne,m1,m,me,test_statistic,perms=50000,alpha=0.025){
    if(ne>n){
        xs <- split(x,rep(1:3,c(n1,n-n1,ne-n)))
    } else {
        if(me>m) stop('Stages with controls only not supported')
        xs <- split(x,rep(1:2,c(n1,ne-n1)))
    }
    if(me>m){
        ys <- split(y,rep(1:3,c(m1,m-m1,me-m)))
    } else {
        if(ne>n) stop('Stages with treatments only not supported')
        ys <- split(y,rep(1:2,c(m1,me-m1)))
    }
    gs <- lapply(1:length(xs),function(i) rep(0:1,c(length(xs[[i]],ys[[i]]))))
    xs <- lapply(1:length(xs),function(i) c(xs[[i]],ys[[i]]))
    A <- permutation_CER(xs[[1]],gs[[1]],xs[[2]],test_statistic,permutations=perms,alpha=alpha)
    q <- perm_test(xs[[2]],xs[[3]],gs[[2]],gs[[3]],test_statistic,B=perms)
    A>=q
}


adaptive_invnormtest_2s <- function(x,y,n1,n,ne,m1,m,me,alpha=0.025){
    xs <- split(x,rep(1:2,c(n1,ne-n1)))
    p1 <- t.test(xs[[1]],alternative='greater')$p.value
    p2 <- t.test(xs[[2]],alternative='greater')$p.value
    alpha >= pnorm(sum(sqrt(c(n1,n-n1)/n) * qnorm(c(p1,p2),lower=F)),lower=FALSE) 
}



