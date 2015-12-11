##' Mean of positive entries
##' 
##' @title Positive sum
##' @template matrix_stats_details
##' @author Florian Klinglmueller
##' @export
possum <- function(x,g){
    if(is.matrix(g)){
        if(is.matrix(x)){
            stop("Only one of g or x may be passed as a matrix")
        }
        colSums(abs(x) * (g>0))
    } else if(is.matrix(x)){
        colSums(abs(x)[g>0,])
    } else {
        sum(abs(x)[g>0])
    }
}
    

##' Mean of paired differences. Usefull in one-sample tests. 
##'
##' @template matrix_stats_details
##' @title Mean of paired differences
##' @author Florian Klinglmueller
##' @export
diffmean <- function(x,g){
    if(is.matrix(g)){
        if(is.matrix(x)){
            stop("Only one of g or x may be passed as a matrix")
        }
        colMeans(abs(x) * (((g>0)*2)-1))
    } else if(is.matrix(x)){
        colMeans(abs(x) * (((g>0)*2)-1))
    } else {
        mean(abs(x) * (((g>0)*2)-1))
    }
}

##' Sum of signed ranks, equivalent to the wilcoxon signed rank test.
##'
##' @template matrix_stats_details
##' @title Sum of signed ranks
##' @author Florian Klinglmueller
##' @export
signedranks <- function(x,g){
    if(is.matrix(g)){
        if(is.matrix(x)){
            stop("Only one of g or x may be passed as a matrix")
        }
        colSums(rank(abs(x)) * (((g>0)*2)-1))
    } else if(is.matrix(x)){
        colSums(colRanks(abs(x)) * (((g>0)*2)-1))
    } else {
        sum(rank(abs(x)) * (((g>0)*2)-1))
    }
}


##' Median of paired differences. Usefull in one-sample tests. 
##'
##' @template matrix_stats_details
##' @title Median of paired differences
##' @author Florian Klinglmueller
##' @export
diffmedian <- function(x,g){
    if(is.matrix(g)){
        if(is.matrix(x)){
            stop("Only one of g or x may be passed as a matrix")
        }
        colMedians(abs(x) * (((g>0)*2)-1))
    } else if(is.matrix(x)){
        colMedians(abs(x) * (((g>0)*2)-1))
    } else {
        median(abs(x) * (((g>0)*2)-1))
    }
}




##' Matrix computation of difference of means between two groups
##'
##' @template matrix_stats_details
##' 
##' @title Difference of means
##' 
##' @author Florian Klinglmueller
##'
##' @export
meandiff <- function(x,g){
    if(is.matrix(g)){
        if(is.matrix(x)){
            stop("Only one of g or x may be passed as a matrix")
        }
        as.numeric(meandiffC(x,g>0))#as.numeric((x %*% {g>0})/colSums({g>0}) - ((x %*% {g<=0})/colSums({g<=0})))
    } else if(is.matrix(x)){
        colMeans(x[g>0,]) - colMeans(x[g<=0,])
    } else {
        sum(x*{g>0})/sum({g>0}) - sum(x*({g<=0}))/sum({g<=0})
    }
}


##' Matrix computation of difference of sums between two groups
##'
##' @template matrix_stats_details
##'
##' @title Sum of differences
##' @author Florian Klinglmueller
##' @export
sumdiff <- function(x,g,...){
    if(is.matrix(x)){
        if(length(dim(g))>1){
            stop("Only one of g orx may be passed as a matrix")
        }
        colSums(x[g>0,])- colSums(x[g<=0,])
    } else if(length(dim(g))>1){
        as.numeric(sumdiffC(x,g>0))#colSums(x * {g>0}) - colSums(x * {g<=0})
    } else {
        sum(x*(g>0)) - sum(x*(g<=0))
    }
}


##' Matrix computation of pooled variances for two groups
##'
##' @template matrix_stats_details
##' 
##' @title Pooled variances
##' @author Florian Klinglmueller
##'
##' @export
pooled_variance <- function(x,g){
    require(matrixStats)
    if(length(dim(g))>1){
        n1 <- colSums(g>0)
        n2 <- colSums(g<=0)
        x1 <- x * (g>0)
        x1[g<=0] <- NA
        x2 <- x * (g<=0)
        x2[g>0] <- NA
        ((n1-1)*colVars(x1,na.rm=T)+(n1-1)*colVars(x2,na.rm=T))/(n1+n2-2)
    } else if(is.matrix(x)){
        n1 <- sum(g>0)
        n2 <- sum(g<=0)
        ((n1-1)*colVars(x[g>0,])+(n1-1)*colVars(x[g<=0,]))/(n1+n2-2)
    } else {
        n1 <- sum(g>0)
        n2 <- sum(g<=0)
        ((n1-1)*var(x[g>0])+(n2-1)*var(x[g<=0]))/(n1+n2-2)
    }
}



##' Matrix computation of z-scores
##'
##' @template matrix_stats_details
##' @param one_sample Whether a one or two sample test should be performed
##' @title z-scores
##' @author Florian Klinglmueller
##'
##' @export
zstat <- function(x,g,sigma=1,one_sample=FALSE){
    if(one_sample){
        if(is.matrix(x)){
            colMeans(x)/(sigma*sqrt(nrow(x)))
        } else {
            mean(x)/sigma*sqrt(nrow(x))
        }
    } else {
        if(is.matrix(g)){
            meandiff(x,g)/(sigma*sqrt(1/colSums(g>0)+1/colSums(g<=0)))
        } else {
            meandiff(x,g)/(sigma*sqrt(1/sum(g>0)+1/sum(g<=0)))
        }
    }
}

##' Matrix computation of t-statistics
##'
##' @template matrix_stats_details
##' @param one_sample Whether a one or two sample test should be performed
##' @title t-statistics
##' @author Florian Klinglmueller
##'
##' @export
tstat <- function(x,g,one_sample=FALSE){
    if(one_sample){
        if(is.matrix(x)){
            sigma  <- colSds(x)
        } else {
            sigma  <- sd(x)
        }
    } else {
        sigma <- sqrt(pooled_variance(x,g))
    }
    zstat(x,g,sigma,one_sample)
}

##' Matrix computation of all pairwise statistics between (treatment) groups and one control group.
##'
##' @template matrix_stats_details
##' @param stat function that computes the test statistic
##' @param control Label that defines the control group
##' @title Pairwise (many-to-one) statistics
##' @author Florian Klinglmueller
##' @export
pw_stat <- function(x,g,stat,control=0,...){
    levels = unique(g[g!=control])
    if(length(dim(g))>1){
        if(is.matrix(x)){
            stop('Vectorized rerandomization not implemented for matrix outcomes')
        }
        if(!all(all(sweep(apply(g,2,table),1,table(g[,1]),`==`)))){
            stop('Group sizes need to be constant across columns of g')
        }
#        stop("Vectorized resampling not yet implemented for pairwise statistics")
        nc <- ncol(g)
        sapply(levels,function(l){
            idx <- apply(g,2,`%in%`,c(control,l))
            gl <- matrix(g[idx],ncol=nc)>0
            xl <- apply(g,2,function(i) x[i])
            stat(t(xl),gl,...)
        })
    } else {
        sapply(levels,function(l){
            idx <- g %in% c(control,l)
            stat(subset(x,idx),g[idx]>0,...)
        })
    }
}

##' Matrix computation of all pairwise mean differences between (treatment) groups and one control group.
##'
##' @template matrix_stats_details
##' @param control Label that defines the control group
##' @title Pairwise (many-to-one) mean differences
##' @author Florian Klinglmueller
##'
##' @export
pw_meandiff <- function(x,g,control=0){
    pw_stat(x,g,meandiff,control)
}



##' Matrix computation of all pairwise z-scores between (treatment) groups and one control group.
##'
##' @template matrix_stats_details
##' @param control Label that defines the control group
##' @title Pairwise (many-to-one) z-scores
##' @author Florian Klinglmueller
##'
##' @export
pw_zstat <- function(x,g,control=0,sigma=1){
    pw_stat(x,g,zstat,control,sigma=sigma)
}

##' Matrix computation of all pairwise t-statistics between (treatment) groups and one control group.
##'
##' @template matrix_stats_details
##' @param control Label that defines the control group
##' @title Pairwise (many-to-one) t-statistics
##' @author Florian Klinglmueller
##'
##' @export
pw_tstat <- function(x,g,control=0){
    pw_stat(x,g,tstat,control)
}
