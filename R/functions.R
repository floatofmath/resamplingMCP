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



##' Matrix computation of difference of means between two groups
##'
##' @template matrix_stats_details
##' 
##' @title Difference of means
##'
##' @export
##' 
##' @author Florian Klinglmueller
meandiff <- function(x,g){
    if(is.matrix(g)){
        if(is.matrix(x)){
            stop("Only one of g or x may be passed as a matrix")
        }
        (x %*% {g>0})/colSums({g>0}) - (x%*%({g<=0}))/colSums({g<=0})
    } else if(is.matrix(x)){
        colMeans(x[g>0,]) - colMeans(x[g<=0,])
    } else {
        sum(x*{g>0})/sum({g>0}) - sum(x*({g<=0}))/sum({g<=0})
    }
}




##' Matrix computation of pooled variances for two groups
##'
##' @template matrix_stats_details
##' 
##' @title Pooled variances
##' @export
##' @author Florian Klinglmueller
pooled_variance <- function(x,g){
    require(matrixStats)
    if(is.matrix(g)){
        stop("Pooled variances not yet implemented for matrix treatment assignments")
    } else if(is.matrix(x)){
        n1 <- sum(g>0)
        n2 <- sum(g<=0)
        ((n1-1)*colVars(x[g>0,])+(n1-1)*colVars(x[g<=0,]))/(n1+n2-2)
    } else {
        n1 <- sum(g>0)
        n2 <- sum(g<=0)
        ((n1-1)*var(x[g>0])+(n2-1)*var(x[g<=0]))/(n1+n2-2)
    }
}


##' Matrix computation of z-scores
##'
##' @template matrix_stats_details
##' @param one_sample Whether a one or two sample test should be performed
##' @title z-scores
##' @export
##' @author Florian Klinglmueller
zstat <- function(x,g,sigma=1,one_sample=FALSE){
    if(one_sample){
        if(is.matrix(x)){
            colMeans(x)/(sigma*sqrt(nrow(x)))
        } else {
            mean(x)/sigma*sqrt(nrow(x))
        }
    } else {
        meandiff(x,g)/(sigma*sqrt(1/sum(g>0)+1/sum(g<=0)))
    }
}

##' Matrix computation of t-statistics
##'
##' @template matrix_stats_details
##' @param one_sample Whether a one or two sample test should be performed
##' @title t-statistics
##' @export
##' @author Florian Klinglmueller
tstat <- function(x,g,one_sample=FALSE){
    if(one_sample){
        if(is.matrix(x)){
            sigma  <- colSds(x)
        } else {
            sigma  <- sd(x)
        }
    } else {
        sigma <- sqrt(pooled_variance(x,g))
    }
    zstat(x,g,sigma,one_sample)
}

##' Matrix computation of all pairwise statistics between (treatment) groups and one control group.
##'
##' @template matrix_stats_details
##' @param stat function that computes the test statistic
##' @param control Label that defines the control group
##' @title Pairwise (many-to-one) statistics
##' @export
##' @author Florian Klinglmueller
pw_stat <- function(x,g,stat,control=0,...){
    levels = unique(g[g!=control])
    if(is.matrix(g)){
        if(is.matrix(x)){
            stop('Vectorized rerandomization not implemented for matrix outcomes')
        }
        if(!all(all(sweep(apply(g,2,table),1,table(g[,1]),`==`)))){
            stop('Group sizes need to be constant across columns of g')
        }
#        stop("Vectorized resampling not yet implemented for pairwise statistics")
        nc <- ncol(g)
        sapply(levels,function(l){
            idx <- apply(g,2,`%in%`,c(control,l))
            gl <- matrix(g[idx],ncol=nc)>0
            xl <- apply(g,2,function(i) x[i])
            stat(t(xl),gl,...)
        })
    } else {
        sapply(levels,function(l){
            idx <- g %in% c(control,l)
            stat(subset(x,idx),g[idx]>0,...)
        })
    }
}

##' Matrix computation of all pairwise mean differences between (treatment) groups and one control group.
##'
##' @template matrix_stats_details
##' @param control Label that defines the control group
##' @title Pairwise (many-to-one) mean differences
##' @export
##' @author Florian Klinglmueller
pw_meandiff <- function(x,g,control=0){
    pw_stat(x,g,meandiff,control)
}



##' Matrix computation of all pairwise z-scores between (treatment) groups and one control group.
##'
##' @template matrix_stats_details
##' @param control Label that defines the control group
##' @title Pairwise (many-to-one) z-scores
##' @export
##' @author Florian Klinglmueller
pw_zstat <- function(x,g,control=0,sigma=1){
    pw_stat(x,g,zstat,control,sigma=sigma)
}

##' Matrix computation of all pairwise t-statistics between (treatment) groups and one control group.
##'
##' @template matrix_stats_details
##' @param control Label that defines the control group
##' @title Pairwise (many-to-one) t-statistics
##' @export
##' @author Florian Klinglmueller
pw_tstat <- function(x,g,control=0){
    pw_stat(x,g,tstat,control)
}


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



##' Perform a graph-based closed test procedure using permutation tests
##'
##' @title graph based multiple comparsions procedure using permutation tests
##' @param permTP flip.object
##' @param G graphMCP
##' @param ... additional arguments to \code{\link{flip::npc}}
##' @return vector with adjustde p-values
##'
##' @importFrom flip npc
##' @importFrom gMCP generateWeights
##' @export
##' @author Florian Klinglmueller
gMCfliP <- function(permTP,G,...){
    weights <- gMCP:::generateWeights(G@m,G@weights)
    n2 <- ncol(weights)
    weight_list <- alply(weights,1,tail,n=n2/2)
    weight_list <- lapply(weight_list,`names<-`,names(permTP))
    subset_list <- alply(weights,1,function(w) names(permTP)[head(w,n=n2/2)>0])
    weight_list <- lapply(1:length(subset_list),function(i) weight_list[[i]][subset_list[[i]]])
    names(weight_list) <- names(subset_list)
    p <- npc(permTP,subsets=subset_list,weights=weight_list,...)
    adjP <- sapply(1:length(permTP),closed_pval,p,m=length(permTP))
    names(adjP) <- names(permTP)
    return(adjP)
}
