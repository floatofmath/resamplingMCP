##' @param tu2 Sum of squares of extended second stage 
##' @param u1 First stage sum of squares
##' @param v1 First stage sum of observations
##' @param tn2 Extended second stage sample size
##' @param n Preplanned total sample size
##' @param n1 First stage sample size
##' @param alpha Significance level
##' @examples
##'
##' ## From Timmesfeld et al. (2007) - should be 0.042 
##' tu2 = 49968.5
##' u1 = 7989
##' v1 = 82.5
##' tn2 = 104
##' n1 = 15
##' n = 65
##' clev(tu2,u1,v1,tn2,n,n1)
##' tccv(tu2,u1,v1,tn2,n,n1)
