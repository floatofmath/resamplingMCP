##' @details Group assignments can be passed via \code{g} using either a numeric/integer factor, vector, or matrix. All elements with values below or equal to zero are treated as belonging to the control, all above as belonging to the treatment group. 
##'
##' If either x or g are passed as matrices test statistics each computed for each column of outcomes or treatment assignments.  Passing both as matrices does not work!
##'
##' @param x either a numeric vector or column matrix of outcomes (see Details)
##' @param g group assignments (see Details)
##' @return numeric vector of test statistics
