context("Permutation Spaces")
## test omega
n1 <- 9
g1 <- sample(c(1,-1),n1,rep=T)



test_that("Omega one stage, restricted",
          {
              o1 <- omega(g1)
              expect_equal(sum(g1>0),unique(rowSums(o1 > 0)))
              expect_equal(nrow(o1),choose(n1,sum(g1>0)))
              expect_equal(ncol(o1),length(g1))
          })


test_that("Omega one stage unrestricted",
          {
              o2 <- omega(g1,restricted=FALSE)
              expect_equal(nrow(o2),2^length(g1))
              expect_equal(ncol(o2),length(g1))
              o2 <- omega(g1,restricted=FALSE,B=10^6)
              expect_equal(dim(o2),c(2^length(g1),length(g1)))
              o2 <- omega(g1,restricted=FALSE,B=3)
              expect_equal(dim(o2),c(4,length(g1)))
          })


n2 <- 7
g2 <- sample(c(1,-1),n2,rep=T)

test_that("Omega two stages restricted",
          {
              o3 <- omega(g1,g2,B = 10^6)
              expect_equal(sum(g1>0) + sum(g2>0),unique(rowSums(o3 >0)))
              expect_equal(dim(o3),c(choose(n1,sum(g1>0))*choose(n2,sum(g2>0)),length(g2) + length(g1)))
          })
## test conddist:
test_that("Conditional permutation distribution",
          {
              expect_true(all(cond_dist(rnorm(4),rnorm(4),rep(c(-1,1),each=2),rep(c(-1,1),each=2),function(y,x) paste(which(x>0),collapse=''),100000) %in% apply(expand.grid('34',apply(4+gtools::combinations(4,2),1,paste,collapse='')),1,paste,collapse='')))
          })

test_that("Permutation distribution",
          {
              expect_true(all(perm_dist(rnorm(4),rnorm(4),rep(c(-1,1),each=2),rep(c(-1,1),each=2),function(y,x) paste(which(x>0),collapse=''),100000) %in% apply(expand.grid(apply(gtools::combinations(4,2),1,paste,collapse=''),apply(4+gtools::combinations(4,2),1,paste,collapse='')),1,paste,collapse='')))
          }
          )

n <- 8
n1 <- 4
x1 <- rep(c(1,-1),n1/2)
g1 <- rep(c(-1,1),each=n1/2)
x2 <- rep(c(1,-1),n1/2)
g2 <- rep(c(-1,1),each=n1/2)
test_that("Permutation test",
          {
              expect_equivalent(perm_test(x1,x2,g1,g2,sumdiff,100000),.75)
          })
          
## test t_test
n <- 8
n1 <- 4
x1 <- rep(c(1,-1),n1/2)
g1 <- rep(c(-1,1),each=n1/2)
x2 <- rep(c(1,-1),n1/2)
g2 <- rep(c(-1,1),each=n1/2)
test_that("t-test",
          {
              expect_equivalent(t_test(x1,x2,g1,g2),.5)
          }
          )

## test permutation cer
n <- 12
n1 <- 6
x1 <- rep(c(-1,1),each=n1/2)
g1 <- rep(c(-1,1),each=n1/2)
x2 <- rep(c(-1,1),each=n1/2)
g2 <- rep(c(-1,1),each=n1/2)
test_that("permutation CER",
          {
              expect_equivalent(permutation_CER(x1,g1,x2,B=100000),0.5)
          })
