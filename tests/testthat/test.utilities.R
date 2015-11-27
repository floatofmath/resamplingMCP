context("Utilities")
expect_aboutequal <- function(object,expected,digits=4,...){
    ro <- round(object,digits)
    re <- round(expected,digits)
    expect_equivalent(ro,re,...)
}

## test mean diff
n <- 30
m <- 100
x <- rnorm(2*n)
X <- matrix(rnorm(2*n*m),nc=m)
largeX <- matrix(rnorm(2*n*10^5),nc=10^5)
g <- rep(c(0,1),each=n)
g2 <- rep(0:2,each=2*n/3)

test_that("Test of meandiff function",
          {
              expect_equivalent(meandiff(x,g),mean(x[(n+1):(2*n)])-mean(x[1:n]))
              expect_equivalent(meandiff(X,g),colMeans(X[(n+1):(2*n),])-colMeans(X[1:n,]))
          })

test_that("Test of pooled variances function",
          {
              n1 <- sum(g>0)
              n2 <- sum(g<=0)
              pv <- function(x,g)(var(x[g>0])*(n1-1)+var(x[g<=0])*(n2-1))/(n1+n2-2)
              expect_equivalent(pooled_variance(x,g),pv(x,g))
              expect_equivalent(pooled_variance(X,g),apply(X,2,function(x) pv(x,g)))
              x <- rnorm(4)
              G <- omega(rep(0:1,each=2))
              X <- apply(G,2,function(i) c(x[i>0],x[i<=0]))
              expect_equivalent(pooled_variance(X,G[,1]),pooled_variance(x,G))
              expect_equivalent(sumdiff(X,G[,1]),sumdiff(x,G))
              expect_equivalent(meandiff(X,G[,1]),meandiff(x,G))
              expect_equivalent(zstat(X,G[,1]),zstat(x,G))
              expect_equivalent(tstat(X,G[,1]),tstat(x,G))
          })

test_that("Test of tstat function",
          {
              expect_equivalent(tstat(x,g),-t.test(x~g,var.equal=T)$statistic)
              expect_equivalent(tstat(X,g),apply(X,2,function(x) -t.test(x~g,var.equal=T)$statistic))
          })

test_that("Test of zstat function",
          {
              expect_equivalent(zstat(x,g),sqrt(n/2)*(mean(x[(n+1):(2*n)])-mean(x[1:n])))
              expect_equivalent(zstat(X,g),sqrt(n/2)*(colMeans(X[(n+1):(2*n),])-colMeans(X[1:n,])))
              expect_aboutequal(sd(zstat(largeX,g)),1,2)
              expect_aboutequal(mean(zstat(largeX,g)),0,2)
          })


test_that("Test pairwise stuff",
          {
              expect_equivalent(pw_meandiff(x,g2),c(mean(x[g2==1])-mean(x[g2==0]),mean(x[g2==2])-mean(x[g2==0])))
              expect_equivalent(pw_meandiff(X,g2),cbind(colMeans(X[g2==1,])-colMeans(X[g2==0,]),colMeans(X[g2==2,])-colMeans(X[g2==0,])))
              expect_equivalent(pw_zstat(x,g2),c(sqrt(n/3)*(mean(x[g2==1])-mean(x[g2==0])),sqrt(n/3)*(mean(x[g2==2])-mean(x[g2==0]))))
              expect_equivalent(pw_zstat(X,g2),sqrt(n/3)*cbind(colMeans(X[g2==1,])-colMeans(X[g2==0,]),colMeans(X[g2==2,])-colMeans(X[g2==0,])))
              expect_aboutequal(pw_meandiff(x,g2)*n*2/3,flip(x,make_pw_contrasts(g2),statTest='sum')@res$Stat)
              expect_aboutequal(as.numeric(t(pw_meandiff(X,g2)*n*2/3)),flip(X,make_pw_contrasts(g2),statTest='sum')@res$Stat)
          })
          
