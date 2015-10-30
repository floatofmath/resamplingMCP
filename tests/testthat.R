library(testthat)
library(codepackage)

test_check("codepackage")

expect_aboutequal <- function(object,expected,digits=4,...){
    ro <- round(object,digits)
    re <- round(expected,digits)
    expect_equivalent(ro,re,...)
}
