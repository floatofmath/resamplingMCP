context("Adaptive procedures")

## Parameter values from Timmesfeld et al.
test_that("Reproduce results from Timmesfeld et al. (2007)",
          {
              tu2 = 49968.5
              u1 = 7989
              v1 = 82.5
              tn2 = 104
              n1 = 15
              n = 65
              expect_equivalent(round(clev(tu2,u1,v1,tn2,n,n1),3),0.042)
              expect_equivalent(round(tccv(tu2,u1,v1,tn2,n,n1),2),1.74)
          })

library(Rmpfr)
n <- 100
tn <- 120
u <- rchisq(10,n)
tu <- u+rchisq(10,tn-n)
pu <- mpfr(u,200)
ptu <- mpfr(tu,200)
uu_1 <- sum(x <- rnorm(n))
puu_1 <- mpfr(uu_1,200)
uu_2 <- sum(x^2)
puu_2 <- mpfr(uu_2,200)
tuu_1 <- sum(c(x,y <- rnorm(tn-n)))
ptuu_1 <- mpfr(tuu_1,200)
tuu_2 <- sum(c(x^2,y^2))
ptuu_2 <- mpfr(tuu_2,200)
g <- 2*n
tg <- 2*tn
test_that("Compute fractions of large numbers",{
    expect_equal(compute_fraction(u,n/2-1,tu,tn/2-1) * (tu-u)^((tn-n)/2-1),
                 as.numeric((pu^(n/2-1)*(ptu-pu)^((tn-n)/2-1))/ptu^(tn/2-1)),expected.label='native',label='compute_fraction',tolerance=0.0001)
    expect_equal(sqrt(tg)*gamma((tg-1)/2)*(tuu_2-uu_2-(tuu_1-tuu_1)^2/(tg-g))^((tg-g-3)/2)*
                     compute_fraction(uu_2-uu_1^2/g,(g-3)/2,tuu_2-tuu_1^2/tg,(tg-3)/2)/
                         (sqrt(pi)*sqrt(g)*gamma((tg-1)/2)*sqrt(tg-g)*gamma((tg-g-1)/2)),
                 as.numeric((sqrt(tg)*gamma((tg-1)/2)*(puu_2-puu_1^2/g)^((g-3)/2)*(ptuu_2-puu_2-(ptuu_1-ptuu_1)^2/(tg-g))^((tg-g-3)/2))/
                                (sqrt(pi)*sqrt(g)*gamma((tg-1)/2)*sqrt(tg-g)*gamma((tg-g-1)/2)*(ptuu_2-ptuu_1^2/tg)^((tg-3)/2))),
                 tolerance=0.0001)
})
