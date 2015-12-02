compute_fraction <- function(x,p,y,q){
    if(q>p) (x/y)^p /y^(q-p) else (x/y)^q *x^(p-q)
}
        

v <- function(x,g=rep.int(1,length(x))){
    sum(x[g>0])
}

uu <- function(x){
    c(sum(x),sum(x*x))
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

cv_twosample <- function(uu_1,uu_2,m,n,alpha){
    cv <- qt(alpha,n+m-2,lower.tail=FALSE)
    (cv* sqrt(m*n*((m+n)*uu_2-uu_1^2)/(m+n-2+cv^2))+m*uu_1)/(m+n)
}

ccv_twosample <- function(uu1_1,uu1_2,uu2_1,uu2_2,v1,m,n,alpha){
    cv_twosample(uu1_1+uu2_1,uu1_2+uu2_2,m,n,alpha)-v1
}

ccer <- function(u2,u1,v1,n,n1,alpha){
    n2 <- n-n1
    .ccv <- ccv(u2,u1,v1,n,alpha)
    ifelse(.ccv^2 >= n2*u2,.ccv>0,pt(sqrt(n2-1)*.ccv/sqrt(abs((n2*u2) - (.ccv^2))),n2-1,lower.tail=FALSE))
}



kmn <- function(x,uu_1,uu_2,m,n){
    ## ! not clear whether they mean (x*(m+n)/(m*n)) or (m+n)/(m*n*x)
    (sqrt(m+n-2) * sqrt(m*n/(m+n)) * ((m+n)/(m*n)*x - uu_1/n))/
        sqrt(abs(uu_2-uu_1^2/(m+n) - (m*n)/(m+n) * ((m+n)/(m*n)*x - uu_1/n)^2))
}

ccer_twosample <- function(uu2_1,uu2_2,uu1_1,uu1_2,v1,m,n,m1,n1,alpha){
    n2 <- n-n1
    m2 <- m-m1
    .ccv <- ccv_twosample(uu1_1,uu1_2,uu2_1,uu2_2,v1,m,n,alpha)
    ifelse((.ccv-m2/(m2+n2)*uu2_1)^2 > abs((uu2_2-uu2_1^2/(m2+n2))*m2*n2/(m2+n2)),
           (.ccv-m2/(m2+n2)*uu2_1) > 0,
           pt(kmn(.ccv,uu2_1,uu2_2,m2,n2),m2+n2-2,lower.tail=FALSE))
}

## slightly faster and numerically ok
dutu_cf <- function(u,tu,n,tn){
    gamma(tn/2)/(gamma(n/2)*gamma((tn-n)/2)) * compute_fraction(u,n/2-1,tu,tn/2-1)*(tu-u)^((tn-n)/2-1)
                                        #(u^(n/2-1)*(tu-u)^((tn-n)/2-1))/tu^(tn/2-1)
}

dutu_const <- function(n,tn){
    gamma(tn/2)/(gamma(n/2)*gamma((tn-n)/2))
}
dutu_var <- function(u,tu,n,tn,mpfr){
    if(mpfr){
        u <- mpfr(u,300)
        tu <- mpfr(tu,300)
        (u^(n/2-1)*(tu-u)^((tn-n)/2-1))/tu^(tn/2-1)
    } else {
        return(compute_fraction(u,n/2-1,tu,tn/2-1)*(tu-u)^((tn-n)/2-1))
    }
}

dutu_mpfr <- function(u,tu,n,tn){
    u <- mpfr(u,500)
    tu <- mpfr(tu,500)
    gamma(tn/2)/(gamma(n/2)*gamma((tn-n)/2)) * (u^(n/2-1)*(tu-u)^((tn-n)/2-1))/tu^(tn/2-1)
}


## numerically unstable
dutu_chisq <- function(u,tu,n,tn){
    dchisq(tu-u,tn-n)*dchisq(u,n)/dchisq(tu,tn)
}

dutu  <- function(u,tu,n,tn,mpfr){
    if(mpfr){
        return(dutu_mpfr(u,tu,n,tn))
    } else {
        return(dutu_cf(u,tu,n,tn))
    }
}
duutuu <- function(uu_1,uu_2,tuu_1,tuu_2,m,n,tm,tn){
    uu_1 <- mpfr(uu_1,10000)
    uu_2 <- mpfr(uu_2,10000)
    tuu_1 <- mpfr(tuu_1,10000)
    tuu_2 <- mpfr(tuu_2,10000)
    g <- m+n
    tg <- tm+tn
    ## (sqrt(tg)*gamma((tg-1)/2)*(uu_2-uu_1^2/g)^((g-3)/2)*(tuu_2-uu_2-(tuu_1-tuu_1)^2/(tg-g))^((tg-g-3)/2))/
    ##     (sqrt(pi)*sqrt(g)*gamma((tg-1)/2)*sqrt(tg-g)*gamma((tg-g-1)/2)*(tuu_2-tuu_1^2/tg)^((tg-3)/2))
    as.numeric({ sqrt(tg)*gamma((tg-1)/2) / (sqrt(pi*g * (tg - g))* gamma((g-1)/2) *gamma((tg-g-1)/2)) } *
        compute_fraction((tuu_2-uu_2-(tuu_1-uu_1)^2/(tg-g)) * (uu_2-uu_1^2/g),(tg-g-3)/2,tuu_2-tuu_1^2/tg,(tg-3)/4) *
                                             compute_fraction((uu_2-uu_1^2/g),(2*g-tg)/2,tuu_2-tuu_1^2/tg,(tg-3)/4))
}

## original <- function(){
##     { sqrt(tg) * gamma((tg-1)/2) *
##          (uu_2-uu_1^2/g)^((g-3)/2) * (tuu_2 - uu_2 - (tuu_1-uu_1)^2/(tg-g))^((tg-g-3)/2) } /
##              ## note that instead of repeating gamma((tg-1)/2) as in the appendix of Timmesfeld et al. (2007) we
##              ## use gamma((g-1)/2) in the denominator
##              { sqrt(pi*g*(tg-g)) * gamma((g-1)/2) * gamma((tg-g-1)/2) *
##                   (tuu_2-tuu_2^2/tg)^((tg-3)/2) }
## }


## new <- function(){
##     { sqrt(tg) * gamma((tg-1)/2) / (sqrt(pi*g*(tg-g)) * gamma((g-1)/2) * gamma((tg-g-1)/2)) } *
##          {(uu_2-uu_1^2/g)^((g-3)/2) * (tuu_2 - uu_2 - (tuu_1-uu_1)^2/(tg-g))^((tg-g-3)/2) } /
##              ## note that instead of repeating gamma((tg-1)/2) as in the appendix of Timmesfeld et al. (2007) we
##              ## use gamma((g-1)/2) in the denominator
##              { (tuu_2-tuu_1^2/tg)^((tg-3)/2) }
## }    

## prec <- function(){
##     uu_1 <- mpfr(uu_1,2000)
##     uu_2 <- mpfr(uu_2,2000)
##     tuu_1 <- mpfr(tuu_1,2000)
##     tuu_2 <- mpfr(tuu_2,2000)
##     { sqrt(tg) * gamma((tg-1)/2) *
##          (uu_2-uu_1^2/g)^((g-3)/2) * (tuu_2 - uu_2 - (tuu_1-uu_1)^2/(tg-g))^((tg-g-3)/2) } /
##              ## note that instead of repeating gamma((tg-1)/2) as in the appendix of Timmesfeld et al. (2007) we
##              ## use gamma((g-1)/2) in the denominator
##              { sqrt(pi*g*(tg-g)) * gamma((g-1)/2) * gamma((tg-g-1)/2) *
##                   (tuu_2-tuu_2^2/tg)^((tg-3)/2) }
## }





##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Conditional level of extended two sample t-test
##' @template timmesfeld
##' @param tuu2 Two dimensional vector containing the sum and the sum of squares of the extended second stage observations
##' @param uu1 Two dimensional vector containing the  sum and the sum of squares of the first stage observations
##' @param v1 Sum of first stage treatment group observations 
##' @param tm2 Extended second stage sample size treatment group
##' @param tn2 Extended second stage sample size control group
##' @param m Preplanned total sample size treatment group
##' @param n Preplanned total sample size control group
##' @param m1 First stage sample size treatment group
##' @param n1 First stage sample size control group
##' @param alpha Significance level
##' @author Florian Klinglmueller
##' @return Conditional level of an extended two sample t-test
clev_twosample <- function(tuu2,uu1,v1,tm2,tn2,m,n,m1,n1,alpha=0.025){
# ##' @import pracma
    m2 <- m-m1
    n2 <- n-n1
    f <- function(x,y) {
        ccer_twosample(y,x,uu1[1],uu1[2],v1,m,n,m1,n1,alpha) *  duutuu(y,x,tuu2[1],tuu2[2],m2,n2,tm2,tn2)
    }
    integral2(f,xmin=0,xmax=tuu2[2],
              ymin=function(x) -sqrt((m2+n2)*x),
              ymax = function(x) sqrt((m2+n2)*x))
}




## m <- n <- 75
## m1 <- n1 <- 20
## v1 <- 112
## uu1 <- c(102,8955.72)
## tm2 <- tn2 <- 76
## tuu2 <- c(421.8,31107)
## clev_twosample(tuu2,uu1,v1,tm2,tn2,m,n,m1,n1)



##' Computes the conditional level of the adaptive one-sided one-sample t-test for the null hypothesis mu=0 against the alternative mu > 0, given first stage sum of observations, sum of squares, and adapted secon stage sum of squares
##' @title Condition level function adaptive t-test
##' @template timmesfeld
##' @return Conditional level for the adapted test
##' @export
clev <- function(tu2,u1,v1,tn2,n,n1,alpha=0.025,mpfr=F){
# ##' @import Rmpfr
    f <- function(u2) ccer(u2,u1,v1,n,n1,alpha) * dutu_var(u2,tu2,n-n1,tn2,mpfr) 
                                        #    tol = ifelse(mpfr,.5*.Machine$double.eps^0.25,.Machine$double.eps^0.25)
    if(mpfr){
        expectation <- try(integrateR(f,lower=0,upper=tu2,rel.tol=1e-20)$value)
    } else {
        expectation <- try(integrate(f,lower=0,upper=tu2)$value)
    }
    if(class(expectation) == "try-error"){
        browser()
    }  
    dutu_const(n-n1,tn2) * expectation
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
