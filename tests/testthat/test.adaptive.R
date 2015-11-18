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
