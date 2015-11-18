v <- function(x){
    sum(x)
}

u <- function(x){
    sum(x*x)
}

cv <- function(u,n,alpha){
    cv <- qt(alpha,n-1,lower.tail=FALSE)
    cv*sqrt(n*u/(n-1+cv^2))
}

ccv <- function(u2,u1,v1,n,alpha){
    cv(u1+u2,n,alpha)-v1
}

ccer <- function(u2,u1,v1,n,n1,alpha){
    n2 <- n-n1
    .ccv <- ccv(u2,u1,v1,n,alpha)
    ifelse(.ccv^2 >= n2*u2,.ccv>0,pt(sqrt(n2-1)*.ccv/sqrt(abs((n2*u2) - (.ccv^2))),n2-1,lower.tail=FALSE))
}
    
dutu <- function(u,tu,n,tn){
    gamma(tn/2)/(gamma(n/2)*gamma((tn-n)/2)) * (u^(n/2-1)*(tu-u)^((tn-n)/2-1))/tu^(tn/2-1)
}


##' Computes the conditional level of the adaptive one-sided one-sample t-test for the null hypothesis mu=0 against the alternative mu > 0, given first stage sum of observations, sum of squares, and adapted secon stage sum of squares
##' @title Condition level function adaptive t-test
##' @template timmesfeld
##' @return Conditional level for the adapted test
##' @export
clev <- function(tu2,u1,v1,tn2,n,n1,alpha=0.025){
    integrate(function(u2) ccer(u2,u1,v1,n,n1,alpha) * dutu(u2,tu2,n,tn2),lower=0,upper=tu2)$value
}


##' Computes the conditional critical value of the adaptive one-sided one-sample t-test for the null hypothesis mu=0 against the alternative mu > 0, given first stage sum of observations, sum of squares, and adapted secon stage sum of squares 
##' @title Conditional critical value adaptive t-test
##' @template timmesfeld
##' @return Conditional level for the adapted test
##' @author Florian Klinglmueller
##' @export
tccv <- function(tu2,u1,v1,tn2,n,n1,alpha=0.025){
    qt(clev(tu2,u1,v1,tn2,n,n1,alpha=0.025),tn2-1,lower.tail=FALSE)
}
