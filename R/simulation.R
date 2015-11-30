adaptive_permtest_os <- function(x,n1,n,ne,test_statistic,perms=50000,alpha=0.025){
    if(ne>n){
        xs <- split(x,rep(1:3,c(n1,n-n1,ne-n)))
        gs <- split(sign(x)>0,rep(1:3,c(n1,n-n1,ne-n)))
    } else {
        xs <- split(x,rep(1:2,c(n1,ne-n1)))
        gs <- split(sign(x)>0,rep(1:2,c(n1,ne-n1)))
        gs[[3]] <- xs[[3]] <- numeric(0)
    }
    A <- permutation_CER(xs[[1]],gs[[1]],xs[[2]],test_statistic,one_sample=TRUE,restricted=FALSE,B=perms,alpha=alpha)
    q <- perm_test(xs[[2]],xs[[3]],gs[[2]],gs[[3]],test_statistic,restricted=FALSE,B=perms)
    A>=q
}

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
    A <- permutation_CER(xs[[1]],gs[[1]],xs[[2]],test_statistic,B=perms,alpha=alpha)
    q <- perm_test(xs[[2]],xs[[3]],gs[[2]],gs[[3]],test_statistic,B=perms)
    A>=q
}



adaptive_ttest_os <- function(x,n1,n,ne,alpha=0.025) {
    if(n == ne){
        return(0.025 >= t.test(x,alternative='greater')$p.value)
    }
    xs <- split(x,rep(1:2,c(n1,ne-n1)))
    V1 <- sum(xs[[1]])
    U <- sum(xs[[1]]^2)
    tU <- sum(xs[[2]]^2)
    A <- clev(tU,U,V1,ne-n1,n,n1,alpha=alpha)
    A >= t.test(xs[[2]],alternative='greater')$p.value
}

adaptive_invnormtest_os <- function(x,n1,n,ne,alpha=0.025){
    xs <- split(x,rep(1:2,c(n1,ne-n1)))
    p1 <- t.test(xs[[1]],alternative='greater')$p.value
    p2 <- t.test(xs[[2]],alternative='greater')$p.value
    alpha >= {sqrt(c(n1,n-n1)/n) * qnorm(c(p1,p2),lower=F)} %>% sum() %>% pnorm(lower=FALSE) 
}

adaptive_invnormtest_2s <- function(x,y,n1,n,ne,m1,m,me,alpha=0.025){
    xs <- split(x,rep(1:2,c(n1,ne-n1)))
    p1 <- t.test(xs[[1]],alternative='greater')$p.value
    p2 <- t.test(xs[[2]],alternative='greater')$p.value
    alpha >= {sqrt(c(n1,n-n1)/n) * qnorm(c(p1,p2),lower=F)} %>% sum() %>% pnorm(lower=FALSE) 
}

inverse_normal <- function(p1,p2,w1,w2){
    (sqrt(w1) * qnorm(p1,lower=F) + sqrt(w2) * qnorm(p2,lower=F)) %>% pnorm(lower=FALSE)
}

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

compare_adaptive_tests <- function(n1,n,rule,rdist,test_statistic,...){
    x <- rdist(n,...)
    ne <- rule(x[1:n1])
    if(ne>n){
        x <- c(x,rdist(ne-n,...))
    } else {
        ne <- n
    }
    list(permtest = adaptive_permtest_os(x,n1,n,ne,test_statistic),
         ttest = adaptive_ttest_os(x,n1,n,ne),
         invnorm = adaptive_invnormtest_os(x,n1,n,ne),
         npcomb = adaptive_npcombtest_os(x,n1,n,ne,test_statistic))
}
    

##' Compare adaptive test procedures for general two-sample cases
##'
##' \code{rdist} needs to take the number of samples to return as its first argument.
##' 
##' @title Compare two-sample adaptive tests
##' @param n1 first stage sample size (control group)
##' @param n preplanned total sample size (control group)
##' @param rule adaptive sample size rule
##' @param rdist random number generator for the data
##' @param test_statistic function that computes the test statistic
##' @param control_opts list of options past to \code{rdist} for the control group
##' @param treatment_opts list of options past to \code{rdist} for the treatment group
##' @param m1 first stage sample size (control group)
##' @param m preplanned total sample size (control group)
##' @param ... 
##' @return 
##' @author Florian Klinglmueller
compare_adaptive_tests_2s <- function(n1,n,rule,rdist,
                                      test_statistic,
                                      control_opts=NULL,
                                      treatment_opts=control_opts,
                                      m1=n1,m=n,...){
    x <- do.call(match.fun(rdist),c(n=n,control_opts))
    y <- do.call(match.fun(rdist),c(n=n,treatment_opts))
    nes <- rule(x[1:n1],y[1:m1])
    if(nes[1]>n){
        ne <- nes[1]
        x <- c(x,rdist(ne[1]-n,...))
    } else {
        ne <- n
    }
    if(ne[2]>m){
        me <- nes[2]
        y <- c(x,rdist(ne[2]-m,...))
    } else {
        me <- m
    }
    list(permtest = adaptive_permtest_2s(x,y,n1,n,ne,m1,m,me,test_statistic),
#         ttest = adaptive_ttest_os(x,n1,n,ne),
         invnorm = adaptive_invnormtest_2s(x,y,n1,n,ne,m1,m,me),
         npcomb = adaptive_npcombtest_2s(x,y,n1,n,ne,m1,m,me,test_statistic))
}


