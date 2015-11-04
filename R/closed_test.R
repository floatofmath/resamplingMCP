##' make a binary vector out of a number
##'
##' @title Convert to binary
##' @param int integer
##' @param n length of binary vector
##' @return binary vector
##' @examples
##' to.binom(3,4)
##'
##' @author Florian Klinglmueller
to.binom <- function (int, n = floor(log2(int)) + 1) 
{
    if (n + 2 <= floor(log2(int))) {
        stop("Vector length to small to hold binary number")
    }
    ((int)%/%2^((n:1) - 1))%%2
}

##' turn a binary vector into a decimal number
##'
##' @title Binary vector to decimal number
##' @param binom 
##' @return Integer
##' 
##' @author Florian Klinglmueller
from.binom <- function (binom){
    sum(2L^(which(rev(binom) > 0)-1))
}

##' Computes and returns all numbers that as binary numbers of length \code{m} that have a 1 at the \code{i}-th position. Can be used to represent intersection hypotheses for efficient implementation of closed testing procedures. 
##'
##' @title Contains i
##' @param i position 
##' @param m number of digits
##' @return integer vector
##' @author Florian Klinglmueller
##'
##' @examples
##'
##' ## which binary numbers (of length 3) have a 1 at the 3rd position
##' out <- contains(3,3)
##' out
##' ## check if true
##' bins <- lapply(out,to.binom,3)
##' bins
##' all(lapply(bins,`[`,3) == 1)
contains <- function(i,m){
    ## computes binary numbers with a 1 at the i'th place
    i <- m-i+1
    n <- 2L^m-1
    if(i == 1){
        return((1:ceiling(n/2)*2)-1)
    }
    first <- 2^(i-1)
    breakpoints <- (1:ceiling(n/first))*first
    breaks <- ceiling(length(breakpoints)/2)
    unlist(lapply(1:breaks,function(i) breakpoints[i*2-1]:(breakpoints[i*2]-1)))
}





##' Performs a closed test procedure using a binary representation of intersection hypotheses 
##'
##' Intersection hypothesis - of a subset of m overall elementary hypotheses - can be represented by integer numbers of length m. \code{closure} takes a matrix of local significance levels ordered by the value of the binary representation of the corresponding intersection hypothesis.
##' 
##' @title Perform a closed test procedure
##' @param p vector of unadjusted p-values
##' @param A matrix of local critical levels sorted by binary representation of intersection hypotheses (see details)
##' @return logical vector wether elementary hypotheses can be rejected
##' @author Florian Klinglmueller
closure <- function(p,A){
    ## this assumes that each line corresponds
    ## to the intersection hypotheses given
    ## the line numbers binary representation
    m <- length(p)
    n <- nrow(A)
    if(n != 2L^m-1) stop("Too few local levels for a closed procedure")
    d <- sweep(A,2,p,">=")
    D <- rep(NA,m)
    for(i in 1:m){
        index <- contains(i,m)
        D[m-i+1] <- all(apply(d[index,],1,any))
    }
    return(D)
}


##' compute ajusted p-value for an elementary hypothesis based on a closed test procedure 
##'
##' @title adjusted closed test pvalue
##' @param i elementary hypothesis to test
##' @param p list of local p-values (see details)
##' @param m overall number of elementary hypotheses
##' @return adjusted p-value
##' 
##' @author Florian Klinglmueller
##'
##' @export
closed_pval <- function(i,p,m=log2(length(p)+1)){
    max(sapply(contains(i,m),function(J) p[J]))
}

##' compute test decision for an elementary hypothesis based on a closed test procedure
##'
##' @title closed test decision
##' @param i elementary hypothesis to test
##' @param p list of local p-values (see details)
##' @param alpha significance level
##' @return logical whether hypothesis can be rejected at FWER level \code{alpha}
##'
##' @author Florian Klinglmueller
##' @export
closed_test  <- function(i,p,m=log2(length(p)+1),alpha=.025){
    all(lapply(contains(i,m),function(J) p[J] <= alpha))
}


setGeneric('closed_pval')
setMethod(closed_pval,signature(p="flip.object"),
          definition = function(i,p,m){
              p <- p@res$`p-value`
              if(m != log2(length(p)+1)) stop("Need p-values for all 2^m-1 intersection hypotheses")
              closed_pval(i,p,m)
          })
