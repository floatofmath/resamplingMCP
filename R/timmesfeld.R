compute_fraction <- function(x,p,y,q){
    ifelse(q>p,(x/y)^p /y^(q-p),(x/y)^q *x^(p-q))
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
    ifelse((.ccv-m2/(m2+n2)*uu2_1)^2 > (m2*n/(m2+n2)*(uu2_2-uu2_1^2/(m2+n2))),(.ccv-m2/(m2+n2)*uu2_1)>0,
           pt(kmn(.ccv,uu2_1,uu2_2,m2,n2),m2+n2-2,lower.tail=FALSE))
}

dutu <- function(u,tu,n,tn){
    gamma(tn/2)/(gamma(n/2)*gamma((tn-n)/2)) * compute_fraction(u,n/2-1,tu,tn/2-1)*(tu-u)^((tn-n)/2-1)
                                        #(u^(n/2-1)*(tu-u)^((tn-n)/2-1))/tu^(tn/2-1)
}

duutuu <- function(uu_1,uu_2,tuu_1,tuu_2,m,n,tm,tn){
    g <- m+n
    tg <- tm+tn
    ## (sqrt(tg)*gamma((tg-1)/2)*(uu_2-uu_1^2/g)^((g-3)/2)*(tuu_2-uu_2-(tuu_1-tuu_1)^2/(tg-g))^((tg-g-3)/2))/
    ##     (sqrt(pi)*sqrt(g)*gamma((tg-1)/2)*sqrt(tg-g)*gamma((tg-g-1)/2)*(tuu_2-tuu_1^2/tg)^((tg-3)/2))
    sqrt(tg)*gamma((tg-1)/2)*(tuu_2-uu_2-(tuu_1-tuu_1)^2/(tg-g))^((tg-g-3)/2)*
        compute_fraction(uu_2-uu_1^2/g,(g-3)/2,tuu_2-tuu_1^2/tg,(tg-3)/2)/
            (sqrt(pi)*sqrt(g)*gamma((tg-1)/2)*sqrt(tg-g)*gamma((tg-g-1)/2))
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param tuu2_1 
##' @param tuu2_2 
##' @param uu1_1 
##' @param uu1_2 
##' @param v1 
##' @param tm2 
##' @param tn2 
##' @param m 
##' @param n 
##' @param m1 
##' @param n1 
##' @param alpha 
##' @return
##' @author Florian Klinglmueller
##' @export
##' @import pracma
clev_twosample <- function(tuu2,uu1,v1,tm2,tn2,m,n,m1,n1,alpha=0.025){
    m2 <- m-m1
    n2 <- n-n1
    f <- function(x,y) {
        ccer_twosample(y,x,uu1[1],uu1[2],v1,m,n,m1,n1,alpha) *  duutuu(y,x,tuu2[1],tuu2[2],m2,n2,tm2,tn2)
    }
    integral2(f,xmin=0,xmax=tuu2[2],
              ymin=function(x) -sqrt((m2+n2)*x),
              ymax = function(x) sqrt((m2+n2)*x))
}




m <- n <- 75
m1 <- n1 <- 20
v1 <- 112
uu1 <- c(102,8955.72)
tm2 <- tn2 <- 76
tuu2 <- c(421.8,31107)
clev_twosample(tuu2,uu1,v1,tm2,tn2,m,n,m1,n1)



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
