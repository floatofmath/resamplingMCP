##' Compute the conditional permutation distribution given first stage group assignments
##'
##' Second and third stage treatment assignments are only passed to define the second stage sample size and the sizes of the corresponding (treatment) groups. If the number of requested permutations \code{B} is larger than the number of all possible permutations, only the latter will be used.
##' 
##' @title Conditional permutations distribution
##' @param x1 First stage data
##' @param x2 Second stage data
##' @param g1 First stage treatment assignments
##' @param g2 (Dummy) second stage treatment assignments (see Details)
##' @param stat Function that computes the test statistic (signature \code{function(x,g)})
##' @param B Number of permutations 
##' @param x3 Third stage data (e.g. sample size increase)
##' @param g3 (Dummy) third stage treatment assignments (see Details)
##' @return numeric vector 
##' @author Florian Klinglmueller
cond_dist <- function(x1,x2,g1,g2,stat,B,x3=NULL,g3=NULL){
    if(B >= choose(length(g2),sum(g2>0))*choose(length(g3),sum(g3>0))){
        omega <- t(apply(omega(g2,g3),1,function(x) c(g1,x)))
        unlist(lapply(1:nrow(omega),function(i) stat(c(x1,x2,x3),omega[i,])))
    } else {
        unlist(lapply(1:B,function(n) stat(c(x1,x2,x3),c(g1,sample(g2),sample(g3)))))
    }
}



##' Compute the conditional permutation distribution given first stage group assignments
##'
##' First, second and third stage treatment assignments are only passed to define the second stage sample size and the sizes of the corresponding (treatment) groups. If the number of requested permutations \code{B} is larger than the number of all possible permutations, only the latter will be used.
##' 
##' @title Conditional permutations distribution
##' @param x1 First stage data
##' @param x2 Second stage data
##' @param g1 (Dummy) first stage treatment assignments
##' @param g2 (Dummy) second stage treatment assignments (see Details)
##' @param stat Function that computes the test statistic (signature \code{function(x,g)})
##' @param B Number of permutations 
##' @param x3 (Dummy) third stage data (e.g. sample size increase)
##' @param g3 (Dummy) third stage treatment assignments (see Details)
##' @param restricted Should group sizes be considered fixed
##' @return numeric vector
##' @author Florian Klinglmueller
perm_dist <- function(x1,x2,g1,g2,stat,B,x3=NULL,g3=NULL,restricted=TRUE){
  if(restricted){  
    if(B >= choose(length(g1),sum(g1>0))*choose(length(g2),sum(g2>0))*choose(length(g3),sum(g3>0))){
        omega <- omega(g1,g2,g3,B)
        unlist(lapply(1:nrow(omega),function(i) stat(c(x1,x2,x3),omega[i,])))
    } else {
        unlist(lapply(1:B,function(n) stat(c(x1,x2,x3),c(sample(g1),sample(g2),sample(g3)))))
    }
  } else {
      return(NULL)
  }
}

##' Returns a logical matrix with all possible assignments (e.g. of patients) of \code{n} objects to two groups (e.g. treatment and control) given that \code{k} objects are in one group
##'
##' @title All combinations of two group assignments
##' @param n number of observations
##' @param k number of observations in (e.g. treatment) group
##' @return logical matrix
##' @author Florian Klinglmueller
all_reassignments <- function(n,k){
    N <- choose(n,k)
    M <- matrix(FALSE,N,n)
    M[cbind(rep(1:N,k),as.integer(gtools::combinations(n,k)))] <- TRUE
    M
}

##' Returns a logical matrix with all possible assignments (e.g. of patients) of \code{n} objects to two groups (e.g. treatment and control)
##'
##' @title All assignments
##' @param n number of observations
##' @return logical matrix
##' @author Florian Klinglmueller
all_assignments <- function(n){
    e1071::bincombinations(n)>0
}

##' Returns a logical matrix with B random assignments (e.g. of patients) of \code{n} objects to two groups (e.g. treatment and control)
##'
##' @title Random combinations of two group assignments
##' @param n number of observations
##' @param B number of random combinations
##' @return logical matrix
##' @author Florian Klinglmueller
random_assignments <- function(n,B){
    matrix(sample(c(T,F),n*B,rep=T),B,n)
}

##' Returns a logical matrix with B random assignments (e.g. of patients) of \code{n} objects to two groups (e.g. treatment and control) given that \code{k} objects are in one group
##'
##' @title Random reassignments to two groups
##' @param n number of observations
##' @param k number of observations in (e.g. treatment) group
##' @param B number of random combinations
##' @return logical matrix
##' @author Florian Klinglmueller
random_reassignments <- function(n,k,B){
    t(replicate(B,sample(c(rep(TRUE,k),rep(FALSE,n-k)))))
}
##library(microbenchmark)

## only faster for small c1
## join_to_row1 <- function(c1,c2){
##     do.call('cbind',lapply(1:ncol(c1),function(i) rbind(matrix(c1[,i],nc=ncol(c2),nr=nrow(c1)),c2)))
## }


## slightly slower maybe due to magrittr
## join_to_row2 <- function(c1,c2){
##     list(x=c1,y=c2) %>%
##         lapply(FUN=function(combs) { combs %>% t  %>% as.data.frame %>% mutate(idx=1)}) %$%
##             full_join(x,y,by='idx') %>% mutate(idx=NULL)
## }

##' Joins two matrices repeats each row of the first matrix by the
##' number of rows of the second matrix and column binds them together to a 
##' 
##' @title Join two matrices by row
##' @param c1 a n1 \times m1 matrix
##' @param c2 a n2 \times m2 matrix
##' @return n1 + n2 \times m1*m2 matrix
##' @author Florian Klinglmueller
join_to_row <- function(c1,c2){
    out <- dplyr::mutate(full_join(dplyr::mutate(as.data.frame(c1),idx=1),
                                   dplyr::mutate(as.data.frame(c2),idx=1),by='idx'),idx=NULL)
    dimnames(out) <- list(1:nrow(out),1:ncol(out))
    out
}

## c1 <- all_reassignments(10,7)
## c2 <- all_reassignments(10,5)
## microbenchmark(join_to_row1(c1,c2),join_to_row2(c1,c2),join_to_row3(c1,c2))

##' Returns all possible assignments to two groups stratified by stages.
##'
##' Number of observations \code{ns} and observations from one group \code{ks} will be recycled to match the longer argument. 
##'
##' @title Stagewise two-group assignments
##' @param ns vector of number of observations per stage
##' @param ks vector of number of observations in group (e.g. treatment) per stage
##' @param random whether all or a limited number of random assignments are desired
##' @param B number of assignments
##' @return logical matrix of dimension \code{sum(ns)} \times \code{prod(choose(ns,ks))}
##' @author Florian Klinglmueller
strat_reassignments <- function(ns,ks,restricted=FALSE,random=FALSE,B=NULL){
    ns <- bt88.03.704::recycle(ns,ks)
    ks <- bt88.03.704::recycle(ks,ns)
    combinations <- if(random) {
        if(restricted) {
            function(n,k) random_reassignments(n,k,B)
        } else {
            function(n,k) random_assignments(n,B)
        }
    } else {
        if(restricted) {
            all_reassignments
        } else {
            function(n,k) all_assignments(n)
        }
    }
    cs <- lapply(which(ns>0),function(i) combinations(ns[i],ks[i]))
    Reduce(join_to_row,cs)
}

##' Reference/permutation space of a two-stage adaptive permutation test
##'
##' @title Reference space
##' @param g1 First stage treatment assignments
##' @param g2 Second stage treatment assignments
##' @param g3 Third stage treatment assignments
##' @param restricted Whether group sizes are considered fixed
##' @param B Number of permutations to be used (if smaller than all permutations)
##' @return logical matrix 
##' @author Florian Klinglmueller
##'
##' @export
omega <- function(g1,g2=NULL,g3=NULL,restricted = TRUE,B=1000){
    ns <- sapply(list(g1,g2,g3),length)
    ks  <- sapply(list(g1,g2,g3),function(g) sum(g>0))
    n_combs <- prod(choose(ns,ks))
    random <- (n_combs > B)
    ps <- strat_reassignments(ns,ks,restricted,random,B)
    if(random) rbind(c(g1,g2,g3),ps) else ps
}

##' Perform a permutation test (stratified by stages)
##'
##' @title perform permutation test
##' @param x1 First stage observations
##' @param x2 Second stage observations
##' @param g1 First stage group assignments
##' @param g2 Second stage group assignments
##' @param stat Function that computes test statistic 
##' @param B Number of permutations to use
##' @param x3 Third stage observations
##' @param g3 Third stage group assignments
##' @return p-value of the permutation test
##' @author Florian Klinglmueller
##'
##' @export
perm_test <- function(x1,x2,g1,g2,stat,B,x3=NULL,g3=NULL){
    cdist <- perm_dist(x1,x2,g1,g2,stat,B,x3,g3)
    mean(cdist>=stat(c(x1,x2,x3),c(g1,g2,g3)))
}
